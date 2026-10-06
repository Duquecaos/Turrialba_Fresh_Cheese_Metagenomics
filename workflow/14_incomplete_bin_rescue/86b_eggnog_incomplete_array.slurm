#!/bin/bash

#SBATCH --partition=shared_64
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=24:00:00

#SBATCH --array=1-32%2

#SBATCH --job-name=86b_eggnog
#SBATCH --output=86b_eggnog_%A_%a.out
#SBATCH --error=86b_eggnog_%A_%a.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

BASE="${ROOT}/71_incomplete_eggnog"
CHUNKS="${BASE}/chunks"
RESULTS="${BASE}/results"
TMPBASE="${BASE}/tmp"

MANIFEST="${BASE}/chunk_manifest.tsv"

DB="${ROOT}/databases/eggnog_5.0.2"

EGGNOG="${EGGNOG:-/opt/ohpc/pub/containers/BIO/eggnog-mapper-2.1.12.sif}"

TASK="${SLURM_ARRAY_TASK_ID}"
IDX=$(printf "%02d" "${TASK}")

PREFIX="chunk_${IDX}"

INPUT="${CHUNKS}/${PREFIX}.faa"

RESULT="${RESULTS}/${PREFIX}"
TMP="${TMPBASE}/${PREFIX}"

ANNOT="${RESULT}/${PREFIX}.emapper.annotations"
HITS="${RESULT}/${PREFIX}.emapper.hits"
SEEDS="${RESULT}/${PREFIX}.emapper.seed_orthologs"


echo "============================================================"
echo "PASO 86b - eggNOG BINS INCOMPLETOS"
echo "============================================================"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task: ${TASK}"
echo "Chunk: ${PREFIX}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memoria solicitada: 128G"
echo "Input: ${INPUT}"
echo "============================================================"
echo


# ============================================================
# VALIDACIONES
# ============================================================

if [[ ! -s "${INPUT}" ]]; then
    echo "ERROR: falta input: ${INPUT}" >&2
    exit 1
fi

if [[ ! -s "${MANIFEST}" ]]; then
    echo "ERROR: falta manifiesto: ${MANIFEST}" >&2
    exit 1
fi

if [[ ! -s "${EGGNOG}" ]]; then
    echo "ERROR: falta contenedor: ${EGGNOG}" >&2
    exit 1
fi

REQUIRED_DB=(
    "${DB}/eggnog.db"
    "${DB}/eggnog_proteins.dmnd"
    "${DB}/eggnog.taxa.db"
    "${DB}/EGGNOG_DB_5.0.2_READY"
)

for F in "${REQUIRED_DB[@]}"; do
    if [[ ! -e "${F}" ]]; then
        echo "ERROR: falta componente eggNOG: ${F}" >&2
        exit 1
    fi
done


# ============================================================
# VALIDAR CONTEO DEL CHUNK
# ============================================================

EXPECTED=$(
    awk -F'\t' \
        -v idx="${TASK}" \
        'NR > 1 && $1 == idx {print $3}' \
        "${MANIFEST}"
)

if [[ -z "${EXPECTED}" ]]; then
    echo "ERROR: no se encontró task ${TASK} en manifiesto." >&2
    exit 1
fi

OBSERVED=$(
    grep -c '^>' "${INPUT}" || true
)

echo "Proteínas esperadas: ${EXPECTED}"
echo "Proteínas observadas: ${OBSERVED}"

if [[ "${EXPECTED}" -ne "${OBSERVED}" ]]; then
    echo "ERROR: conteo del chunk no coincide." >&2
    exit 1
fi


# ============================================================
# SI YA ESTÁ COMPLETO, NO REPETIR
# ============================================================

if [[ -s "${ANNOT}" && -s "${HITS}" && -s "${SEEDS}" ]]; then

    echo
    echo "============================================================"
    echo "${PREFIX} YA ESTÁ COMPLETO"
    echo "No se repetirá eggNOG."
    echo "============================================================"

    NANNOT=$(
        grep -vc '^#' "${ANNOT}" || true
    )

    echo "Registros de anotación existentes: ${NANNOT}"
    echo "PASO 86b ${PREFIX} YA ESTABA COMPLETO; ES SEGURO SALIR."

    exit 0
fi


# ============================================================
# EVITAR SOBREESCRIBIR RESULTADOS INCOMPLETOS
# ============================================================

if [[ -d "${RESULT}" ]]; then

    if find "${RESULT}" \
        -mindepth 1 \
        -maxdepth 1 \
        -type f \
        | grep -q .
    then
        echo "ERROR: existe resultado incompleto previo:" >&2
        echo "${RESULT}" >&2
        echo "No se sobrescribirá automáticamente." >&2
        exit 2
    fi
fi


rm -rf "${TMP}"

mkdir -p \
    "${RESULT}" \
    "${TMP}"


# ============================================================
# PREFLIGHT DEL CONTENEDOR
# ============================================================

echo
echo "===== PREFLIGHT eggNOG ====="

apptainer exec \
    --cleanenv \
    --bind "${DB}:/eggnog_data:ro" \
    "${EGGNOG}" \
    emapper.py \
        --data_dir /eggnog_data \
        --version


# ============================================================
# EGGNOG-MAPPER
# ============================================================

echo
echo "============================================================"
echo "INICIANDO eggNOG-mapper"
echo "Entrada: proteínas previamente predichas"
echo "Tipo: proteins"
echo "Búsqueda: DIAMOND"
echo "============================================================"
echo

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    apptainer exec \
        --cleanenv \
        --bind "${DB}:/eggnog_data:ro" \
        --bind "${CHUNKS}:/input:ro" \
        --bind "${BASE}:/work" \
        "${EGGNOG}" \
        emapper.py \
            --cpu "${SLURM_CPUS_PER_TASK}" \
            -i "/input/${PREFIX}.faa" \
            --itype proteins \
            -m diamond \
            --data_dir /eggnog_data \
            -o "${PREFIX}" \
            --output_dir "/work/results/${PREFIX}" \
            --temp_dir "/work/tmp/${PREFIX}"


# ============================================================
# VALIDAR OUTPUTS
# ============================================================

echo
echo "============================================================"
echo "VALIDANDO OUTPUTS ${PREFIX}"
echo "============================================================"

for F in \
    "${ANNOT}" \
    "${HITS}" \
    "${SEEDS}"
do
    if [[ ! -s "${F}" ]]; then
        echo "ERROR: output faltante/vacío: ${F}" >&2
        exit 1
    fi
done


NANNOT=$(
    grep -vc '^#' "${ANNOT}" || true
)

NHITS=$(
    grep -vc '^#' "${HITS}" || true
)

NSEEDS=$(
    grep -vc '^#' "${SEEDS}" || true
)

echo "Proteínas input:              ${OBSERVED}"
echo "Registros emapper.annotations: ${NANNOT}"
echo "Registros emapper.hits:        ${NHITS}"
echo "Registros seed_orthologs:      ${NSEEDS}"


if [[ "${NANNOT}" -le 0 ]]; then
    echo "ERROR: annotations no contiene registros." >&2
    exit 1
fi


# ============================================================
# MARCADOR DE COMPLETADO
# ============================================================

{
    echo -e "array_task\t${TASK}"
    echo -e "chunk\t${PREFIX}"
    echo -e "input_proteins\t${OBSERVED}"
    echo -e "annotation_records\t${NANNOT}"
    echo -e "hit_records\t${NHITS}"
    echo -e "seed_records\t${NSEEDS}"
    echo -e "completed_at\t$(date --iso-8601=seconds)"
} > "${RESULT}/COMPLETED.tsv"


rm -rf "${TMP}"


echo
echo "============================================================"
echo "PASO 86b ${PREFIX} FINALIZÓ CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "============================================================"
