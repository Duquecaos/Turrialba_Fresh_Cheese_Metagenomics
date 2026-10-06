#!/bin/bash

#SBATCH --job-name=95a_16S_shotgun
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=95a_16S_shotgun_%j.out
#SBATCH --error=95a_16S_shotgun_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"
BARRNAP="${BARRNAP:-/opt/ohpc/pub/containers/BIO/barrnap-0.9.sif}"
SAMTOOLS="${SAMTOOLS:-/opt/ohpc/pub/containers/BIO/samtools-1.21.sif}"
BOWTIE2="${BOWTIE2:-/opt/ohpc/pub/containers/BIO/bowtie2-2.5.4.sif}"
BLAST="${BLAST:-/opt/ohpc/pub/containers/BIO/blast-2.16.0.sif}"

ASSEMBLY="${ROOT}/05_assembly"
HOSTFREE="${ROOT}/04_host_removed"
MAPPING="${ROOT}/06_mapping"
FINAL18="${ROOT}/22_final_representative_mags/fastas"

OUT="${ROOT}/97_16S_shotgun/95A_audit"

rm -rf "$OUT"
mkdir -p "$OUT"

export LC_ALL=C
export LANG=C
export LANGUAGE=C

echo "============================================================"
echo "95A - AUDITORIA 16S SHOTGUN"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

# ------------------------------------------------------------
# 1. Herramientas
# ------------------------------------------------------------

{
    echo -e "tool\tstatus\tpath"

    for SPEC in \
        "apptainer:${APPTAINER}" \
        "barrnap:${BARRNAP}" \
        "samtools:${SAMTOOLS}" \
        "bowtie2:${BOWTIE2}" \
        "blast:${BLAST}"
    do
        NAME="${SPEC%%:*}"
        PATHX="${SPEC#*:}"

        if [[ -e "$PATHX" ]]; then
            echo -e "${NAME}\tFOUND\t${PATHX}"
        else
            echo -e "${NAME}\tNOT_FOUND\t${PATHX}"
        fi
    done
} > "${OUT}/95A_tools.tsv"

[[ -x "$APPTAINER" ]] || {
    echo "ERROR: Apptainer no disponible" >&2
    exit 2
}

[[ -s "$BARRNAP" ]] || {
    echo "ERROR: Barrnap no disponible" >&2
    exit 3
}

"$APPTAINER" exec \
    "$BARRNAP" \
    barrnap --version \
    > "${OUT}/95A_barrnap_version.txt" \
    2> "${OUT}/95A_barrnap_version.stderr" \
    || true

# ------------------------------------------------------------
# 2. Coensamblajes
# ------------------------------------------------------------

python3 - \
    "$ASSEMBLY" \
    "$OUT" <<'PY'

import csv
import sys
from pathlib import Path

assembly = Path(sys.argv[1])
out = Path(sys.argv[2])

groups = ["L1", "L2", "L3", "M1", "M2", "M3"]

rows = []

for group in groups:

    fasta = assembly / group / "final.contigs.fa"

    rows.append({
        "coassembly": group,
        "fasta": str(fasta),
        "exists": "YES" if fasta.is_file() else "NO",
        "size_bytes": fasta.stat().st_size if fasta.is_file() else 0,
    })

if sum(x["exists"] == "YES" for x in rows) != 6:
    raise RuntimeError(
        "No están disponibles los seis coensamblajes"
    )

with (out / "95A_coassemblies.tsv").open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=list(rows[0].keys()),
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)

PY

# ------------------------------------------------------------
# 3. Final18 MAGs
# ------------------------------------------------------------

python3 - \
    "$FINAL18" \
    "$OUT" <<'PY'

import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
out = Path(sys.argv[2])

files = []

for pattern in ("*.fa", "*.fna", "*.fasta"):
    files.extend(root.glob(pattern))

files = sorted(set(files))

rows = []

for f in files:

    mag = f.stem

    coassembly = (
        mag.split("__", 1)[0]
        if "__" in mag
        else ""
    )

    rows.append({
        "MAG": mag,
        "coassembly": coassembly,
        "fasta": str(f),
        "size_bytes": f.stat().st_size,
    })

if len(rows) != 18:
    raise RuntimeError(
        f"Esperaba 18 MAGs; observé {len(rows)}"
    )

with (out / "95A_final18.tsv").open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=[
            "MAG",
            "coassembly",
            "fasta",
            "size_bytes",
        ],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)

PY

# ------------------------------------------------------------
# 4. FASTQ host-free
# ------------------------------------------------------------

python3 - \
    "$HOSTFREE" \
    "$OUT" <<'PY'

import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
out = Path(sys.argv[2])

files = sorted(
    f for f in root.rglob("*.fastq.gz")
    if "hostfree" in f.name.lower()
)

rows = []

for f in files:

    rows.append({
        "sample_directory": f.parent.name,
        "filename": f.name,
        "path": str(f),
        "size_bytes": f.stat().st_size,
    })

with (out / "95A_hostfree_fastqs.tsv").open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=[
            "sample_directory",
            "filename",
            "path",
            "size_bytes",
        ],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)

sample_dirs = sorted(
    {
        x["sample_directory"]
        for x in rows
    }
)

with (out / "95A_hostfree_summary.tsv").open("w") as fh:

    fh.write("metric\tvalue\n")
    fh.write(f"hostfree_fastq_files\t{len(rows)}\n")
    fh.write(f"sample_directories\t{len(sample_dirs)}\n")

PY

# ------------------------------------------------------------
# 5. BAM existentes
# ------------------------------------------------------------

python3 - \
    "$MAPPING" \
    "$OUT" <<'PY'

import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
out = Path(sys.argv[2])

bams = sorted(
    root.rglob("*.sorted.bam")
)

rows = []

for bam in bams:

    rel = bam.relative_to(root)
    parts = rel.parts

    coassembly = (
        parts[0]
        if len(parts) >= 1
        else ""
    )

    sample = (
        parts[1]
        if len(parts) >= 2
        else bam.stem
    )

    bai1 = Path(str(bam) + ".bai")
    bai2 = bam.with_suffix(".bai")

    rows.append({
        "coassembly": coassembly,
        "sample": sample,
        "bam": str(bam),
        "index_exists": (
            "YES"
            if bai1.is_file() or bai2.is_file()
            else "NO"
        ),
        "size_bytes": bam.stat().st_size,
    })

with (out / "95A_existing_BAMs.tsv").open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=[
            "coassembly",
            "sample",
            "bam",
            "index_exists",
            "size_bytes",
        ],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)

with (out / "95A_existing_BAM_summary.tsv").open("w") as fh:

    fh.write("metric\tvalue\n")
    fh.write(f"sorted_BAMs\t{len(rows)}\n")
    fh.write(
        "BAMs_with_index\t"
        f"{sum(x['index_exists']=='YES' for x in rows)}\n"
    )

PY

# ------------------------------------------------------------
# 6. Buscar SILVA solamente como referencia taxonómica
# ------------------------------------------------------------

: > "${OUT}/95A_SILVA_candidates.txt"

for D in \
    "$HOME" \
    "/scratch/global/$USER"
do
    [[ -d "$D" ]] || continue

    find "$D" \
        -maxdepth 5 \
        -type f \
        \( \
            -iname '*silva*.fa' \
            -o -iname '*silva*.fasta' \
            -o -iname '*silva*.fna' \
            -o -iname '*silva*.qza' \
            -o -iname '*silva*.rds' \
            -o -iname '*silva*.tsv' \
            -o -iname '*silva*.txt' \
        \) \
        2>/dev/null \
        >> "${OUT}/95A_SILVA_candidates.txt" \
        || true
done

sort -u \
    "${OUT}/95A_SILVA_candidates.txt" \
    -o "${OUT}/95A_SILVA_candidates.txt"

# ------------------------------------------------------------
# 7. Alcance metodológico congelado
# ------------------------------------------------------------

cat > "${OUT}/95A_methodological_scope.tsv" <<'EOF'
field	value
analysis_scope	shotgun_16S_only
amplicon_analysis	NO
amplicon_comparison	NO
primary_16S_detection	Barrnap_on_six_coassemblies
secondary_16S_detection	Barrnap_on_final18_MAGs
coassembly_role	primary_community_16S_locus_discovery
MAG_role	exact_bin_association_when_16S_is_retained_in_MAG
rRNA_not_present_in_MAG_equals_taxon_absence	NO
fragmented_16S_negative_equals_absence	NO
taxonomy_reference	prefer_SILVA_if_available
VFDB_or_GTDB_reference_used_for_16S_taxonomy	NO
read_based_quantification_planned	YES
16S_multimapping_expected	YES
multiple_16S_operons_per_genome_possible	YES
assembled_16S_copy_count_equals_organism_abundance	NO
EOF

# ------------------------------------------------------------
# 8. Validación
# ------------------------------------------------------------

for F in \
    95A_tools.tsv \
    95A_coassemblies.tsv \
    95A_final18.tsv \
    95A_hostfree_fastqs.tsv \
    95A_hostfree_summary.tsv \
    95A_existing_BAMs.tsv \
    95A_existing_BAM_summary.tsv \
    95A_SILVA_candidates.txt \
    95A_methodological_scope.tsv
do
    [[ -e "${OUT}/${F}" ]] || {
        echo "ERROR: falta ${F}" >&2
        exit 4
    }
done

printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"${SLURM_JOB_ID}" \
> "${OUT}/95A_COMPLETE.ok"

echo
echo "============================================================"
echo "TOOLS"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/95A_tools.tsv"

echo
echo "============================================================"
echo "BARRNAP"
echo "============================================================"

cat "${OUT}/95A_barrnap_version.txt" \
    2>/dev/null \
    || true

echo
echo "============================================================"
echo "HOSTFREE"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/95A_hostfree_summary.tsv"

echo
echo "============================================================"
echo "BAMS"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/95A_existing_BAM_summary.tsv"

echo
echo "============================================================"
echo "SILVA"
echo "============================================================"

head -n 50 \
    "${OUT}/95A_SILVA_candidates.txt"

echo
echo "============================================================"
echo "95A FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"

