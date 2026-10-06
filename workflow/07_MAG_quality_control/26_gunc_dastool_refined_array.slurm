#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00

#SBATCH --array=1-6%2
#SBATCH --job-name=gunc_refined

#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

# Locate repository configuration portably.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PROJECT_CONFIG="${PROJECT_CONFIG:-${REPO_ROOT}/config/config.sh}"

if [[ ! -r "$PROJECT_CONFIG" ]]; then
    echo "ERROR: configuration file not found: $PROJECT_CONFIG" >&2
    echo "Create config/config.sh from config/config.example.sh before running this script." >&2
    exit 2
fi

source "$PROJECT_CONFIG"

COASSEMBLY_NAMES=(L1 L2 L3 M1 M2 M3)

INDEX=$((SLURM_ARRAY_TASK_ID - 1))
GROUP="${COASSEMBLY_NAMES[$INDEX]}"

OUT="${SCRATCH_DIR}/14_gunc_refined/${GROUP}"
LIST="${OUT}/${GROUP}_candidates_50_10.list"
REPORT="${OUT}/GUNC.progenomes_2.1.maxCSS_level.tsv"
DONE="${OUT}/done"

TMP_HOST="${SCRATCH_DIR}/tmp/gunc_ref_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
TMP_CONTAINER="/gunc_tmp"

mkdir -p "${OUT}" "${TMP_HOST}"

cleanup() {
    echo
    echo "===== CLEANING TEMPORARY FILES ====="
    rm -rf "${TMP_HOST}"
}

trap cleanup EXIT

echo "============================================================"
echo "GUNC DAS TOOL REFINED"
echo "============================================================"
echo "Group:    ${GROUP}"
echo "Inicio:   $(date)"
echo "Nodo:     $(hostname)"
echo "Lista:    ${LIST}"
echo "Database: ${GUNC_DB}"
echo "============================================================"
echo

if [[ ! -s "${GUNC}" ]]; then
    echo "ERROR: contenedor GUNC no encontrado." >&2
    exit 1
fi

if [[ ! -s "${GUNC_DB}" ]]; then
    echo "ERROR: base GUNC no encontrada." >&2
    exit 1
fi

if [[ ! -s "${LIST}" ]]; then
    echo "ERROR: lista de candidatos vacía." >&2
    exit 1
fi

N_EXPECTED=$(wc -l < "${LIST}")

echo "Candidatos esperados: ${N_EXPECTED}"

while read -r FASTA
do
    if [[ ! -s "${FASTA}" ]]; then
        echo "ERROR: FASTA no encontrado: ${FASTA}" >&2
        exit 1
    fi
done < "${LIST}"

if [[ -f "${DONE}" && -s "${REPORT}" ]]; then

    echo "GUNC ya completado para ${GROUP}."

else

    rm -f "${DONE}"

    echo
    echo "===== RUNNING GUNC ====="
    echo

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        apptainer exec \
        --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
        --bind "${TMP_HOST}:${TMP_CONTAINER}" \
        --env TMPDIR="${TMP_CONTAINER}" \
        "${GUNC}" \
        gunc run \
        -f "${LIST}" \
        -r "${GUNC_DB}" \
        -t "${SLURM_CPUS_PER_TASK}" \
        -o "${OUT}" \
        --temp_dir "${TMP_CONTAINER}" \
        --detailed_output

fi

if [[ ! -s "${REPORT}" ]]; then
    echo "ERROR: reporte GUNC no encontrado." >&2
    exit 1
fi

N_REPORTED=$(( $(wc -l < "${REPORT}") - 1 ))

echo
echo "Esperados: ${N_EXPECTED}"
echo "Reportados: ${N_REPORTED}"

if [[ "${N_EXPECTED}" -ne "${N_REPORTED}" ]]; then
    echo "ERROR: GUNC no reportó todos los candidatos." >&2
    exit 1
fi

echo
echo "===== GUNC COUNTS ====="

awk -F '\t' '
NR==1 {
    for(i=1;i<=NF;i++)
        if($i=="pass.GUNC")
            passcol=i
    next
}

{
    total++

    if($passcol=="True")
        pass++
    else if($passcol=="False")
        fail++
}

END {
    print "Total =", total+0
    print "GUNC pass =", pass+0
    print "GUNC fail =", fail+0
}
' "${REPORT}"

echo
echo "===== FIRST RESULTS ====="

head -n 10 "${REPORT}" | column -s $'\t' -t

touch "${DONE}"

echo
echo "============================================================"
echo "GUNC REFINED COMPLETADO"
echo "Group: ${GROUP}"
echo "Bins evaluados: ${N_REPORTED}"
echo "Fin: $(date)"
echo "============================================================"
