#!/bin/bash

#SBATCH --job-name=93c0_coloc_audit
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=93c0_coloc_audit_%j.out
#SBATCH --error=93c0_coloc_audit_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"
OUT="${ROOT}/95_mobilome/93C0_colocalization_audit"

mkdir -p "$OUT"

FILES=(
"${ROOT}/25_abricate_consolidated/plasmid_replicon_loci_90_90.tsv"
"${ROOT}/25_abricate_consolidated/amr_loci_90_90.tsv"
"${ROOT}/25_abricate_consolidated/virulence_loci_90_90.tsv"
"${ROOT}/25_abricate_consolidated/abricate_nonredundant_loci_90_90.tsv"
"${ROOT}/80_resistome_community_consolidated/92"*".tsv"
)

echo "======================================================================"
echo "93C0 - AUDITORIA PARA CO-LOCALIZACION"
echo "Inicio: $(date)"
echo "======================================================================"

python3 - "$ROOT" "$OUT" <<'PY'
import csv
import glob
import os
import sys
from pathlib import Path

root = Path(sys.argv[1])
out = Path(sys.argv[2])

patterns = [
    root / "25_abricate_consolidated" / "plasmid_replicon_loci_90_90.tsv",
    root / "25_abricate_consolidated" / "amr_loci_90_90.tsv",
    root / "25_abricate_consolidated" / "virulence_loci_90_90.tsv",
    root / "25_abricate_consolidated" / "abricate_nonredundant_loci_90_90.tsv",
]

# Search all plausible community-resistome tables rather than assuming
# one filename.
for p in [
    root / "80_resistome_community_consolidated",
    root / "81_resistome_ARG_abundance",
]:
    if p.exists():
        patterns.extend(
            Path(x)
            for x in glob.glob(str(p / "*.tsv"))
        )

files = []
seen = set()

for p in patterns:
    p = Path(p)
    if p.is_file() and p not in seen:
        files.append(p)
        seen.add(p)

if not files:
    raise RuntimeError("No encontré archivos candidatos")


interesting = (
    "mag", "contig", "sequence", "seq", "scaffold",
    "gene", "symbol", "product", "database", "db",
    "start", "end", "strand", "identity", "coverage",
    "type", "element", "method", "protein"
)

rows = []

for path in files:

    with path.open(errors="replace") as fh:

        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []

        examples = []
        n = 0

        for row in reader:
            n += 1

            if len(examples) < 5:

                vals = []

                for field in fields:
                    lf = field.lower()

                    if any(k in lf for k in interesting):

                        value = str(row.get(field, ""))

                        if len(value) > 120:
                            value = value[:117] + "..."

                        vals.append(
                            f"{field}={value}"
                        )

                examples.append(
                    " | ".join(vals)
                )

    rows.append({
        "file": path.name,
        "path": str(path),
        "rows": n,
        "n_columns": len(fields),
        "columns": ";".join(fields),
        "example1": examples[0] if len(examples) > 0 else "",
        "example2": examples[1] if len(examples) > 1 else "",
        "example3": examples[2] if len(examples) > 2 else "",
        "example4": examples[3] if len(examples) > 3 else "",
        "example5": examples[4] if len(examples) > 4 else "",
    })


outfile = out / "93C0_input_schema_audit.tsv"

with outfile.open("w", newline="") as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(rows[0].keys()),
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(rows)

PY


# ======================================================================
# Reconstruir mapa traceable <-> original por orden de contigs
#
# Solo auditoría: verificamos que pueda construirse sin ambigüedad.
# ======================================================================

python3 - "$ROOT" "$OUT" <<'PY'
import sys
from pathlib import Path

root = Path(sys.argv[1])
out = Path(sys.argv[2])

magdir = (
    root
    / "22_final_representative_mags"
    / "fastas"
)

trace = (
    root
    / "22_final_representative_mags"
    / "representative18_traceable.fasta"
)

if not trace.is_file():
    raise RuntimeError(f"Falta {trace}")


def fasta_headers(path):
    vals = []

    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                vals.append(
                    line[1:].strip().split()[0]
                )

    return vals


magfiles = sorted(
    list(magdir.glob("*.fa"))
    + list(magdir.glob("*.fna"))
    + list(magdir.glob("*.fasta"))
)

if len(magfiles) != 18:
    raise RuntimeError(
        f"Esperaba 18 MAGs; encontré {len(magfiles)}"
    )


trace_ids = fasta_headers(trace)

trace_by_mag = {}

for tid in trace_ids:

    if "___contig" not in tid:
        raise RuntimeError(
            f"Formato traceable inesperado: {tid}"
        )

    mag = tid.split(
        "___contig",
        1
    )[0]

    trace_by_mag.setdefault(
        mag,
        []
    ).append(tid)


outfile = (
    out
    / "93C0_traceable_original_mapping.tsv"
)

auditfile = (
    out
    / "93C0_traceable_mapping_audit.tsv"
)


with outfile.open("w") as oh, auditfile.open("w") as ah:

    oh.write(
        "MAG\tcoassembly\ttraceable_contig\t"
        "original_contig\tfull_original_header\n"
    )

    ah.write(
        "MAG\toriginal_contigs\ttraceable_contigs\tstatus\n"
    )


    for fasta in magfiles:

        mag = fasta.stem
        coassembly = mag.split("__", 1)[0]

        original_headers = []
        original_ids = []

        with fasta.open() as fh:

            for line in fh:

                if line.startswith(">"):

                    full = line[1:].strip()
                    oid = full.split()[0]

                    original_headers.append(full)
                    original_ids.append(oid)


        tids = trace_by_mag.get(
            mag,
            []
        )


        if len(tids) != len(original_ids):

            status = "COUNT_MISMATCH"

        else:

            # traceable contig numbers should be sequential
            expected = [
                f"{mag}___contig{i:06d}"
                for i in range(
                    1,
                    len(tids) + 1
                )
            ]

            status = (
                "PASS"
                if tids == expected
                else "ORDER_REVIEW"
            )


        ah.write(
            f"{mag}\t"
            f"{len(original_ids)}\t"
            f"{len(tids)}\t"
            f"{status}\n"
        )


        if status == "PASS":

            for tid, oid, full in zip(
                tids,
                original_ids,
                original_headers,
            ):

                oh.write(
                    f"{mag}\t"
                    f"{coassembly}\t"
                    f"{tid}\t"
                    f"{oid}\t"
                    f"{full}\n"
                )

PY


# ======================================================================
# geNomad candidate names
# ======================================================================

python3 - "$ROOT" "$OUT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
out = Path(sys.argv[2])

catalog = (
    root
    / "95_mobilome"
    / "93B_plasmid_catalog"
    / "93B_plasmid_catalog_all.tsv"
)

if not catalog.is_file():
    raise RuntimeError(f"Falta {catalog}")

rows = []

with catalog.open() as fh:

    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in r:
        rows.append(row)

outfile = (
    out
    / "93C0_genomad_plasmid_id_examples.tsv"
)

with outfile.open("w") as fh:

    fh.write(
        "seq_name\tcoassembly\tlength\t"
        "conjugation_genes\tamr_genes\n"
    )

    for row in rows[:30]:

        fh.write(
            f"{row['seq_name']}\t"
            f"{row['coassembly']}\t"
            f"{row['length']}\t"
            f"{row['conjugation_genes']}\t"
            f"{row['amr_genes']}\n"
        )

PY


# ======================================================================
# Final
# ======================================================================

printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"$SLURM_JOB_ID" \
> "${OUT}/93C0_COMPLETE.ok"


echo
echo "======================================================================"
echo "INPUT SCHEMAS"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93C0_input_schema_audit.tsv"


echo
echo "======================================================================"
echo "TRACEABLE MAPPING AUDIT"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93C0_traceable_mapping_audit.tsv"


echo
echo "======================================================================"
echo "geNomad ID EXAMPLES"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93C0_genomad_plasmid_id_examples.tsv"


echo
echo "======================================================================"
echo "93C0 FINALIZO CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "======================================================================"

