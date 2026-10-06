#!/bin/bash

#SBATCH --job-name=94b1_VFDB
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --array=0-5%3

#SBATCH --mail-type=END,FAIL

#SBATCH --output=94b1_VFDB_%A_%a.out
#SBATCH --error=94b1_VFDB_%A_%a.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"
ABRICATE="${ABRICATE:-/opt/ohpc/pub/containers/BIO/abricate-1.0.1.sif}"

ASSEMBLY="${ROOT}/05_assembly"
OUT="${ROOT}/96_virulence/94B1_VFDB_community"

mkdir -p "$OUT"

export LC_ALL=C
export LANG=C
export LANGUAGE=C

GROUPS_LIST=(
    "L1"
    "L2"
    "L3"
    "M1"
    "M2"
    "M3"
)

IDX="${SLURM_ARRAY_TASK_ID}"
GROUP="${GROUPS_LIST[$IDX]}"

FASTA="${ASSEMBLY}/${GROUP}/final.contigs.fa"

TASKOUT="${OUT}/${GROUP}"

mkdir -p "$TASKOUT"

RAW="${TASKOUT}/${GROUP}_VFDB_raw.tsv"
SCREEN="${TASKOUT}/${GROUP}_VFDB_screen80_80.tsv"
HIGH="${TASKOUT}/${GROUP}_VFDB_high90_90.tsv"
SUMMARY="${TASKOUT}/${GROUP}_VFDB_summary.tsv"

echo "======================================================================"
echo "94B1 - VFDB COMMUNITY SCREEN"
echo "Grupo: ${GROUP}"
echo "Job: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Inicio: $(date)"
echo "======================================================================"

# ======================================================================
# 1. Validar entradas
# ======================================================================

[[ -x "$APPTAINER" ]] || {
    echo "ERROR: Apptainer no disponible" >&2
    exit 2
}

[[ -s "$ABRICATE" ]] || {
    echo "ERROR: falta contenedor ABRicate" >&2
    exit 3
}

[[ -s "$FASTA" ]] || {
    echo "ERROR: falta FASTA $FASTA" >&2
    exit 4
}

"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$ABRICATE" \
    test -r "$FASTA" || {
        echo "ERROR: el contenedor no puede leer $FASTA" >&2
        exit 5
    }

# ======================================================================
# 2. Ejecutar ABRicate/VFDB
#
# Recuperamos primero el cribado 80/80 para conservar comparabilidad
# con los resultados existentes de final18.
# ======================================================================

echo
echo "Ejecutando ABRicate VFDB sobre ${GROUP}..."

"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$ABRICATE" \
    abricate \
        --db vfdb \
        --minid 80 \
        --mincov 80 \
        "$FASTA" \
    > "$RAW"

[[ -s "$RAW" ]] || {
    echo "ERROR: ABRicate no produjo una tabla legible" >&2
    exit 6
}

# ======================================================================
# 3. Parsear y clasificar
# ======================================================================

python3 - \
    "$RAW" \
    "$SCREEN" \
    "$HIGH" \
    "$SUMMARY" \
    "$GROUP" <<'PY'

import csv
import sys
from pathlib import Path

raw = Path(sys.argv[1])
screen_out = Path(sys.argv[2])
high_out = Path(sys.argv[3])
summary_out = Path(sys.argv[4])
group = sys.argv[5]


with raw.open(errors="replace") as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    rows = list(reader)
    fields = reader.fieldnames or []


if not fields:
    raise RuntimeError(
        "No se recuperó header de ABRicate"
    )


required = {
    "SEQUENCE",
    "START",
    "END",
    "STRAND",
    "GENE",
    "%COVERAGE",
    "%IDENTITY",
    "DATABASE",
    "ACCESSION",
    "PRODUCT",
}

missing = sorted(
    required - set(fields)
)

if missing:
    raise RuntimeError(
        "Columnas faltantes: "
        + ",".join(missing)
    )


def number(value):

    try:
        return float(value)
    except Exception:
        return None


screen = []

high = []


for i, row in enumerate(
    rows,
    start=1
):

    cov = number(
        row.get("%COVERAGE", "")
    )

    ident = number(
        row.get("%IDENTITY", "")
    )

    if cov is None or ident is None:
        raise RuntimeError(
            f"Hit {i}: cobertura/identidad no numérica"
        )

    if cov < 80 or ident < 80:
        raise RuntimeError(
            f"ABRicate devolvió hit por debajo de 80/80: "
            f"cov={cov}, id={ident}"
        )

    rec = {
        "community_vfdb_locus_id":
            f"VFDB_{group}_{i:05d}",
        "coassembly": group,
        "screen_class": (
            "high_90_90"
            if cov >= 90 and ident >= 90
            else "screen_80_80"
        ),
        **row,
    }

    screen.append(rec)

    if cov >= 90 and ident >= 90:
        high.append(rec)


out_fields = [
    "community_vfdb_locus_id",
    "coassembly",
    "screen_class",
] + fields


def write_table(path, data):

    with path.open(
        "w",
        newline=""
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=out_fields,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(data)


write_table(
    screen_out,
    screen
)

write_table(
    high_out,
    high
)


genes_all = {
    x.get("GENE", "").strip()
    for x in screen
    if x.get("GENE", "").strip()
}

genes_high = {
    x.get("GENE", "").strip()
    for x in high
    if x.get("GENE", "").strip()
}

contigs_all = {
    x.get("SEQUENCE", "").strip()
    for x in screen
    if x.get("SEQUENCE", "").strip()
}

contigs_high = {
    x.get("SEQUENCE", "").strip()
    for x in high
    if x.get("SEQUENCE", "").strip()
}


with summary_out.open("w") as fh:

    fh.write("metric\tvalue\n")

    fh.write(
        f"coassembly\t{group}\n"
    )

    fh.write(
        f"VFDB_hits_80_80\t{len(screen)}\n"
    )

    fh.write(
        f"VFDB_hits_high_90_90\t{len(high)}\n"
    )

    fh.write(
        f"unique_gene_labels_80_80\t{len(genes_all)}\n"
    )

    fh.write(
        f"unique_gene_labels_high_90_90\t{len(genes_high)}\n"
    )

    fh.write(
        f"contigs_with_VFDB_80_80\t{len(contigs_all)}\n"
    )

    fh.write(
        f"contigs_with_VFDB_high_90_90\t{len(contigs_high)}\n"
    )


print(
    f"GROUP={group}"
)

print(
    f"VFDB_80_80={len(screen)}"
)

print(
    f"VFDB_90_90={len(high)}"
)

PY

# ======================================================================
# 4. Provenance
# ======================================================================

{
    echo -e "field\tvalue"
    echo -e "tool\tABRicate"
    echo -e "version\t1.0.1"
    echo -e "database\tVFDB"
    echo -e "database_record\tvfdb 2597 nucl 2024-Dec-15"
    echo -e "coassembly\t${GROUP}"
    echo -e "input\t${FASTA}"
    echo -e "primary_screen_min_identity_percent\t80"
    echo -e "primary_screen_min_coverage_percent\t80"
    echo -e "high_similarity_identity_percent\t90"
    echo -e "high_similarity_coverage_percent\t90"
    echo -e "VFDB_hit_equals_pathogenic_phenotype\tNO"
    echo -e "VFDB_hit_equals_pathogenic_organism\tNO"
} > "${TASKOUT}/${GROUP}_VFDB_metadata.tsv"

sha256sum \
    "$RAW" \
    "$SCREEN" \
    "$HIGH" \
> "${TASKOUT}/${GROUP}_VFDB_sha256.tsv"

# ======================================================================
# 5. Marker
# ======================================================================

printf \
"COMPLETED\t%s\tjob=%s_%s\tcoassembly=%s\n" \
"$(date -Iseconds)" \
"${SLURM_ARRAY_JOB_ID}" \
"${SLURM_ARRAY_TASK_ID}" \
"${GROUP}" \
> "${TASKOUT}/94B1_COMPLETE.ok"

echo
echo "======================================================================"
echo "SUMMARY ${GROUP}"
echo "======================================================================"

column -t -s $'\t' \
    "$SUMMARY"

echo
echo "======================================================================"
echo "94B1 ${GROUP} FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "======================================================================"

