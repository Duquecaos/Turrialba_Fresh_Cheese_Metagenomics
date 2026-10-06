#!/bin/bash

#SBATCH --job-name=93d3_isescan
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=23:00:00
#SBATCH --array=0-5%2

#SBATCH --mail-type=END,FAIL

#SBATCH --output=93d3_isescan_%A_%a.out
#SBATCH --error=93d3_isescan_%A_%a.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"
IS_SIF="${ROOT}/software/mobilome/isescan-1.7.3.sif"

OUTROOT="${ROOT}/95_mobilome/93D3_ISEScan"

GROUP_LIST=(
    L1
    L2
    L3
    M1
    M2
    M3
)

GROUP="${GROUP_LIST[$SLURM_ARRAY_TASK_ID]}"

INPUT="${ROOT}/05_assembly/${GROUP}/final.contigs.fa"
OUT="${OUTROOT}/${GROUP}"

mkdir -p "$OUT"

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"


echo "======================================================================"
echo "93D3 - ISEScan"
echo "Grupo: $GROUP"
echo "Inicio: $(date)"
echo "Job: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "======================================================================"


[[ -x "$APPTAINER" ]] || {
    echo "ERROR: Apptainer no disponible" >&2
    exit 2
}

[[ -s "$IS_SIF" ]] || {
    echo "ERROR: falta $IS_SIF" >&2
    exit 3
}

[[ -s "$INPUT" ]] || {
    echo "ERROR: falta $INPUT" >&2
    exit 4
}


# Confirmar visibilidad dentro del contenedor
"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$IS_SIF" \
    test -r "$INPUT" || {
        echo "ERROR: contenedor no puede leer $INPUT" >&2
        exit 5
    }


cd "$OUT"


echo
echo "Ejecutando ISEScan..."

srun --exact \
    -n 1 \
    -c "$SLURM_CPUS_PER_TASK" \
    "$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$IS_SIF" \
    isescan.py \
    --seqfile "$INPUT" \
    --output "$OUT" \
    --nthread "$SLURM_CPUS_PER_TASK" \
    > "${OUT}/93D3_isescan.stdout" \
    2> "${OUT}/93D3_isescan.stderr"


# ======================================================================
# Localizar TSV principal
# ======================================================================

mapfile -t TSVFILES < <(
    find "$OUT" \
        -type f \
        -name 'final.contigs.fa.tsv' \
        -print \
    | sort
)

if [[ "${#TSVFILES[@]}" -ne 1 ]]
then
    echo "ERROR: esperaba 1 TSV principal; encontré ${#TSVFILES[@]}" >&2
    printf '%s\n' "${TSVFILES[@]:-}" >&2
    exit 6
fi

TSV="${TSVFILES[0]}"

echo -e "result_tsv\t${TSV}" \
    > "${OUT}/93D3_result_file.tsv"


# ======================================================================
# Parseo y QC
# ======================================================================

python3 - \
    "$TSV" \
    "$OUT" \
    "$GROUP" <<'PY'

import csv
import sys
from collections import Counter
from pathlib import Path

tsv = Path(sys.argv[1])
out = Path(sys.argv[2])
group = sys.argv[3]

with tsv.open(errors="replace") as fh:
    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )
    rows = list(reader)
    fields = reader.fieldnames or []


required = {
    "seqID",
    "family",
    "cluster",
    "isBegin",
    "isEnd",
    "isLen",
    "orfBegin",
    "orfEnd",
    "strand",
    "type",
}

missing = required - set(fields)

if missing:
    raise RuntimeError(
        f"Faltan columnas: {sorted(missing)}"
    )


type_counts = Counter()
family_counts = Counter()
contigs = set()
complete_contigs = set()
partial_contigs = set()

for row in rows:

    typ = str(
        row.get("type", "")
    ).strip()

    fam = str(
        row.get("family", "")
    ).strip()

    seqid = str(
        row.get("seqID", "")
    ).strip()

    type_counts[typ] += 1

    if fam:
        family_counts[fam] += 1

    if seqid:
        contigs.add(seqid)

        if typ == "c":
            complete_contigs.add(seqid)

        elif typ == "p":
            partial_contigs.add(seqid)


unknown_types = {
    x
    for x in type_counts
    if x not in {"c", "p"}
}

if unknown_types:
    raise RuntimeError(
        f"Tipos ISEScan inesperados: {unknown_types}"
    )


with (
    out / "93D3_group_summary.tsv"
).open("w") as fh:

    fh.write("metric\tvalue\n")
    fh.write(f"coassembly\t{group}\n")
    fh.write(
        f"IS_predictions_total\t{len(rows)}\n"
    )
    fh.write(
        f"complete_IS_c\t{type_counts.get('c', 0)}\n"
    )
    fh.write(
        f"partial_IS_p\t{type_counts.get('p', 0)}\n"
    )
    fh.write(
        f"contigs_with_any_IS\t{len(contigs)}\n"
    )
    fh.write(
        f"contigs_with_complete_IS\t{len(complete_contigs)}\n"
    )
    fh.write(
        f"contigs_with_partial_IS\t{len(partial_contigs)}\n"
    )
    fh.write(
        f"IS_families\t{len(family_counts)}\n"
    )


with (
    out / "93D3_family_counts.tsv"
).open("w") as fh:

    fh.write(
        "family\tcount\n"
    )

    for family, count in sorted(
        family_counts.items(),
        key=lambda x: (-x[1], x[0])
    ):
        fh.write(
            f"{family}\t{count}\n"
        )


with (
    out / "93D3_schema.tsv"
).open("w") as fh:

    fh.write(
        "columns\t"
        + ";".join(fields)
        + "\n"
    )


preview = rows[:10]

with (
    out / "93D3_preview.tsv"
).open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()

    if preview:
        writer.writerows(preview)


print(
    f"{group}\t"
    f"total={len(rows)}\t"
    f"complete={type_counts.get('c',0)}\t"
    f"partial={type_counts.get('p',0)}\t"
    f"families={len(family_counts)}"
)

PY


# ======================================================================
# Provenance
# ======================================================================

{
    echo -e "field\tvalue"
    echo -e "tool\tISEScan"
    echo -e "version\t1.7.3"
    echo -e "input\t${INPUT}"
    echo -e "coassembly\t${GROUP}"
    echo -e "threads\t${SLURM_CPUS_PER_TASK}"
    echo -e "complete_and_partial_IS_retained\tYES"
    echo -e "removeShortIS_used\tNO"
    echo -e "complete_type_code\tc"
    echo -e "partial_type_code\tp"
    echo -e "IS_prediction_equals_demonstrated_transposition\tNO"
    echo -e "cross_coassembly_detection_equals_prevalence\tNO"
} > "${OUT}/93D3_run_metadata.tsv"


sha256sum "$TSV" \
    > "${OUT}/93D3_primary_result_sha256.tsv"


for F in \
    "${OUT}/93D3_group_summary.tsv" \
    "${OUT}/93D3_family_counts.tsv" \
    "${OUT}/93D3_schema.tsv" \
    "${OUT}/93D3_preview.tsv" \
    "${OUT}/93D3_run_metadata.tsv" \
    "${OUT}/93D3_primary_result_sha256.tsv"
do
    [[ -s "$F" ]] || {
        echo "ERROR: salida faltante/vacía $F" >&2
        exit 7
    }
done


printf \
"COMPLETED\t%s\tjob=%s_%s\tgroup=%s\n" \
"$(date -Iseconds)" \
"$SLURM_ARRAY_JOB_ID" \
"$SLURM_ARRAY_TASK_ID" \
"$GROUP" \
> "${OUT}/93D3_COMPLETE.ok"


echo
echo "======================================================================"
echo "RESULTADO $GROUP"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93D3_group_summary.tsv"

echo
echo "======================================================================"
echo "93D3 ${GROUP} FINALIZO CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "======================================================================"

