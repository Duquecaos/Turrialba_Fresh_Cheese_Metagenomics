#!/bin/bash

#SBATCH --job-name=94a_VFDB_audit
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=94a_VFDB_audit_%j.out
#SBATCH --error=94a_VFDB_audit_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"
ABRICATE="${ABRICATE:-/opt/ohpc/pub/containers/BIO/abricate-1.0.1.sif}"

ASSEMBLY="${ROOT}/05_assembly"
REP18="${ROOT}/24_abricate_rep18/abricate_all_hits.tsv"

OUT="${ROOT}/96_virulence/94A_audit"

rm -rf "$OUT"
mkdir -p "$OUT"

export LC_ALL=C
export LANG=C

echo "======================================================================"
echo "94A - AUDITORIA DE VIRULENCIA / VFDB"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "======================================================================"

# ======================================================================
# 1. Validar herramientas
# ======================================================================

[[ -x "$APPTAINER" ]] || {
    echo "ERROR: Apptainer no disponible" >&2
    exit 2
}

[[ -s "$ABRICATE" ]] || {
    echo "ERROR: contenedor ABRicate no disponible: $ABRICATE" >&2
    exit 3
}

"$APPTAINER" exec \
    "$ABRICATE" \
    abricate --version \
    > "${OUT}/94A_abricate_version.txt" \
    2> "${OUT}/94A_abricate_version.stderr"

"$APPTAINER" exec \
    "$ABRICATE" \
    abricate --list \
    > "${OUT}/94A_abricate_databases.txt" \
    2> "${OUT}/94A_abricate_databases.stderr"

# ======================================================================
# 2. Confirmar VFDB
# ======================================================================

awk '
    NR > 1 && tolower($1) == "vfdb" {
        print
    }
' "${OUT}/94A_abricate_databases.txt" \
> "${OUT}/94A_vfdb_database_record.txt"

if [[ ! -s "${OUT}/94A_vfdb_database_record.txt" ]]
then
    echo "ERROR: VFDB no aparece en abricate --list" >&2
    exit 4
fi

# ======================================================================
# 3. Auditar los hits final18 existentes
# ======================================================================

[[ -s "$REP18" ]] || {
    echo "ERROR: falta tabla final18: $REP18" >&2
    exit 5
}

python3 - \
    "$REP18" \
    "$OUT" <<'PY'

import csv
import sys
from pathlib import Path

src = Path(sys.argv[1])
out = Path(sys.argv[2])

with src.open(errors="replace") as fh:
    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )
    rows = list(reader)
    fields = reader.fieldnames or []

if not fields:
    raise RuntimeError(
        "No se pudo recuperar el header de abricate_all_hits.tsv"
    )

db_candidates = [
    x for x in fields
    if x.lower().lstrip("#") == "database"
]

db_col = (
    db_candidates[0]
    if db_candidates
    else None
)

if db_col is not None:
    vfdb = [
        r for r in rows
        if r.get(db_col, "").strip().lower() == "vfdb"
    ]
else:
    # Rescate conservador si la columna DATABASE no existe.
    vfdb = [
        r for r in rows
        if any(
            str(v).strip().lower() == "vfdb"
            for v in r.values()
        )
    ]

with (
    out / "94A_existing_final18_vfdb.tsv"
).open("w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        extrasaction="ignore",
    )
    writer.writeheader()
    writer.writerows(vfdb)

with (
    out / "94A_existing_final18_summary.tsv"
).open("w") as fh:

    fh.write("metric\tvalue\n")
    fh.write(f"existing_total_abricate_hits\t{len(rows)}\n")
    fh.write(f"existing_VFDB_hits\t{len(vfdb)}\n")
    fh.write(f"database_column\t{db_col or 'NOT_DETECTED'}\n")
    fh.write(f"n_columns\t{len(fields)}\n")

with (
    out / "94A_existing_final18_header.txt"
).open("w") as fh:
    fh.write("\t".join(fields) + "\n")

print(f"EXISTING_VFDB_HITS={len(vfdb)}")

PY

# ======================================================================
# 4. Auditar las seis coensamblajes
# ======================================================================

python3 - \
    "$ASSEMBLY" \
    "$OUT" <<'PY'

import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
out = Path(sys.argv[2])

groups = [
    "L1",
    "L2",
    "L3",
    "M1",
    "M2",
    "M3",
]

def fasta_stats(path):

    n = 0
    total = 0
    current = 0
    lengths = []

    with path.open() as fh:
        for line in fh:

            if line.startswith(">"):

                if current:
                    lengths.append(current)
                    total += current

                n += 1
                current = 0

            else:
                current += len(
                    line.strip()
                )

        if current:
            lengths.append(current)
            total += current

    return {
        "contigs": n,
        "total_bp": total,
        "contigs_ge_1000": sum(
            x >= 1000
            for x in lengths
        ),
        "contigs_ge_5000": sum(
            x >= 5000
            for x in lengths
        ),
        "longest_contig_bp": (
            max(lengths)
            if lengths
            else 0
        ),
    }

rows = []

for group in groups:

    fasta = (
        root
        / group
        / "final.contigs.fa"
    )

    if not fasta.is_file():
        raise RuntimeError(
            f"Falta coensamblaje: {fasta}"
        )

    stats = fasta_stats(fasta)

    rows.append({
        "coassembly": group,
        "fasta": str(fasta),
        **stats,
        "file_size_bytes": fasta.stat().st_size,
    })

with (
    out / "94A_coassembly_inputs.tsv"
).open("w", newline="") as fh:

    fields = [
        "coassembly",
        "fasta",
        "contigs",
        "total_bp",
        "contigs_ge_1000",
        "contigs_ge_5000",
        "longest_contig_bp",
        "file_size_bytes",
    ]

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(rows)

with (
    out / "94A_coassembly_global_summary.tsv"
).open("w") as fh:

    fh.write("metric\tvalue\n")
    fh.write(f"coassemblies\t{len(rows)}\n")
    fh.write(
        "total_contigs\t"
        f"{sum(x['contigs'] for x in rows)}\n"
    )
    fh.write(
        "total_bp\t"
        f"{sum(x['total_bp'] for x in rows)}\n"
    )
    fh.write(
        "total_contigs_ge_1000\t"
        f"{sum(x['contigs_ge_1000'] for x in rows)}\n"
    )
    fh.write(
        "total_contigs_ge_5000\t"
        f"{sum(x['contigs_ge_5000'] for x in rows)}\n"
    )

print("COASSEMBLY_AUDIT=PASS")

PY

# ======================================================================
# 5. Test de visibilidad de una coensamblaje
# ======================================================================

TEST_FASTA="${ASSEMBLY}/L1/final.contigs.fa"

"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$ABRICATE" \
    test -r "$TEST_FASTA"

# ======================================================================
# 6. Metodología congelada
# ======================================================================

cat > "${OUT}/94A_methodological_scope.tsv" <<'EOF'
field	value
primary_database	VFDB
community_scope	six_coassemblies
final18_scope	existing_ABRicate_results
community_detection_method	ABRicate_VFDB
gene_homology_equals_pathogenic_phenotype	NO
virulence_gene_detection_equals_pathogenic_organism	NO
negative_result_equals_absence_from_community	NO
MAG_assignment_requires_exact_contig_membership	YES
community_contig_hit_without_MAG_assignment	community_level_only
planned_primary_filter_percent_identity	90
planned_primary_filter_percent_coverage	90
EOF

# ======================================================================
# 7. Validación
# ======================================================================

REQUIRED=(
    "94A_abricate_version.txt"
    "94A_abricate_databases.txt"
    "94A_vfdb_database_record.txt"
    "94A_existing_final18_vfdb.tsv"
    "94A_existing_final18_summary.tsv"
    "94A_coassembly_inputs.tsv"
    "94A_coassembly_global_summary.tsv"
    "94A_methodological_scope.tsv"
)

for F in "${REQUIRED[@]}"
do
    [[ -s "${OUT}/${F}" ]] || {
        echo "ERROR: salida faltante/vacía ${OUT}/${F}" >&2
        exit 6
    }
done

printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"${SLURM_JOB_ID}" \
> "${OUT}/94A_COMPLETE.ok"

echo
echo "======================================================================"
echo "ABRICATE VERSION"
echo "======================================================================"

cat "${OUT}/94A_abricate_version.txt"

echo
echo "======================================================================"
echo "VFDB DATABASE"
echo "======================================================================"

cat "${OUT}/94A_vfdb_database_record.txt"

echo
echo "======================================================================"
echo "EXISTING FINAL18"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/94A_existing_final18_summary.tsv"

echo
echo "======================================================================"
echo "COASSEMBLY SUMMARY"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/94A_coassembly_global_summary.tsv"

echo
echo "======================================================================"
echo "94A FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "======================================================================"

