#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00

#SBATCH --job-name=gunc_post_mdm
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

INPUT="${ROOT}/16_mdmcleaner_kept/fastas"
OUT="${ROOT}/18_gunc_post_mdm"

GUNC="${GUNC:-/opt/ohpc/pub/containers/BIO/gunc-1.0.6.sif}"

echo "============================================================"
echo "GUNC POST-MDMCLEANER"
echo "============================================================"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Input: ${INPUT}"
echo "Output: ${OUT}"
echo "============================================================"
echo


# ============================================================
# VALIDAR CONTENEDOR
# ============================================================

if [[ ! -s "${GUNC}" ]]; then
    echo "ERROR: no existe contenedor GUNC:" >&2
    echo "${GUNC}" >&2
    exit 1
fi


# ============================================================
# VALIDAR 8 MAGs
# ============================================================

N_MAG=$(find "${INPUT}" \
    -maxdepth 1 \
    -type f \
    -name '*.fa' \
    | wc -l)

echo "MAGs detectados: ${N_MAG}"

if [[ "${N_MAG}" -ne 8 ]]; then
    echo "ERROR: se esperaban 8 MAGs." >&2
    exit 1
fi

echo
find "${INPUT}" \
    -maxdepth 1 \
    -type f \
    -name '*.fa' \
    -printf '%f\t%s bytes\n' \
    | sort


# ============================================================
# LOCALIZAR BASE proGenomes 2.1
# ============================================================

DB=""

CANDIDATES=(
    "${ROOT}/databases/gunc/gunc_db_progenomes2.1.dmnd"
    "${ROOT}/databases/gunc_db_progenomes2.1.dmnd"
    "${GUNC_DB:-/home/public/met-workshop/databases/gunc_db_progenomes2.1.dmnd}"
)

for F in "${CANDIDATES[@]}"; do
    if [[ -s "${F}" ]]; then
        DB="${F}"
        break
    fi
done

if [[ -z "${DB}" ]]; then
    echo "ERROR: no se localizó gunc_db_progenomes2.1.dmnd" >&2
    exit 1
fi

echo
echo "===== BASE GUNC ====="
ls -lh "${DB}"


# ============================================================
# EVITAR SOBREESCRITURA
# ============================================================

if [[ -e "${OUT}" ]]; then
    echo "ERROR: ya existe ${OUT}" >&2
    echo "No se sobrescribirá automáticamente." >&2
    exit 1
fi

mkdir -p "${OUT}"


# ============================================================
# VERSIÓN
# ============================================================

echo
echo "===== GUNC ====="

apptainer exec \
    --cleanenv \
    "${GUNC}" \
    gunc --version || true


# ============================================================
# EJECUTAR GUNC
# ============================================================

echo
echo "============================================================"
echo "INICIANDO GUNC"
echo "============================================================"
echo

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    apptainer exec \
    --cleanenv \
    --bind "${INPUT}:/input:ro" \
    --bind "${OUT}:/output" \
    --bind "${DB}:/guncdb/gunc_db_progenomes2.1.dmnd:ro" \
    "${GUNC}" \
    gunc run \
    -d /input \
    -r /guncdb/gunc_db_progenomes2.1.dmnd \
    -t "${SLURM_CPUS_PER_TASK}" \
    -o /output \
    -e .fa


# ============================================================
# VALIDAR RESULTADO
# ============================================================

REPORT="${OUT}/GUNC.progenomes_2.1.maxCSS_level.tsv"

if [[ ! -s "${REPORT}" ]]; then
    echo "ERROR: no se generó el reporte maxCSS esperado:" >&2
    echo "${REPORT}" >&2
    exit 1
fi

N_RESULT=$(awk '
NR > 1 && NF > 0 {
    n++
}
END {
    print n+0
}' "${REPORT}")

echo
echo "MAGs reportados por GUNC: ${N_RESULT}"

if [[ "${N_RESULT}" -ne 8 ]]; then
    echo "ERROR: se esperaban 8 resultados GUNC." >&2
    exit 1
fi


# ============================================================
# MOSTRAR RESULTADO
# ============================================================

echo
echo "===== GUNC POST-MDMCLEANER ====="

column -t -s $'\t' "${REPORT}"


echo
echo "============================================================"
echo "GUNC POST-MDMCLEANER COMPLETADO"
echo "Reporte: ${REPORT}"
echo "Fin: $(date)"
echo "============================================================"
