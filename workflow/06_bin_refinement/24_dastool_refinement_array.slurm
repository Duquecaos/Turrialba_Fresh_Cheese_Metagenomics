#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00

#SBATCH --array=1-6%2
#SBATCH --job-name=dastool_refine

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

ASSEMBLY="${SCRATCH_DIR}/05_assembly/${GROUP}/final.contigs.fa"

INPUT_DIR="${SCRATCH_DIR}/12_dastool/${GROUP}/inputs"

META="${INPUT_DIR}/${GROUP}_metabat2_contigs2bin.tsv"
MAXBIN="${INPUT_DIR}/${GROUP}_maxbin2_contigs2bin.tsv"
CONCOCT="${INPUT_DIR}/${GROUP}_concoct_contigs2bin.tsv"

OUTDIR="${SCRATCH_DIR}/12_dastool/${GROUP}/refined"
OUTBASE="${OUTDIR}/${GROUP}_DASTool"

BINDIR="${OUTBASE}_DASTool_bins"
LOG="${OUTBASE}_DASTool.log"
C2B="${OUTBASE}_DASTool_contig2bin.tsv"
EVAL="${OUTBASE}_allBins.eval"

DONE="${OUTDIR}/done"

TMP_HOST="${SCRATCH_DIR}/tmp/dastool_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
TMP_CONTAINER="/dastool_tmp"

mkdir -p "${OUTDIR}"
mkdir -p "${TMP_HOST}"

cleanup() {
    echo
    echo "===== CLEANING TEMPORARY FILES ====="
    rm -rf "${TMP_HOST}"
}

trap cleanup EXIT


echo "============================================================"
echo "DAS TOOL REFINEMENT"
echo "============================================================"
echo "Group:       ${GROUP}"
echo "Inicio:      $(date)"
echo "Nodo:        $(hostname)"
echo "Job ID:      ${SLURM_JOB_ID}"
echo "Array ID:    ${SLURM_ARRAY_TASK_ID}"
echo "CPUs:        ${SLURM_CPUS_PER_TASK}"
echo "Memoria:     32G"
echo
echo "Assembly:"
echo "  ${ASSEMBLY}"
echo
echo "MetaBAT2:"
echo "  ${META}"
echo
echo "MaxBin2:"
echo "  ${MAXBIN}"
echo
echo "CONCOCT:"
echo "  ${CONCOCT}"
echo
echo "Container:"
echo "  ${DASTOOL}"
echo
echo "Output:"
echo "  ${OUTBASE}"
echo "============================================================"
echo


# ============================================================
# VALIDACIONES
# ============================================================

for F in \
    "${ASSEMBLY}" \
    "${META}" \
    "${MAXBIN}" \
    "${CONCOCT}" \
    "${DASTOOL}"
do

    if [[ ! -s "${F}" ]]; then
        echo "ERROR: archivo no encontrado o vacío:" >&2
        echo "${F}" >&2
        exit 1
    fi

done


echo "===== DEPENDENCIES ====="

apptainer exec \
    "${DASTOOL}" \
    bash --noprofile --norc -c '
    for X in DAS_Tool prodigal diamond pullseq ruby
    do
        printf "%-20s " "$X"

        command -v "$X" || {
            echo "NOT_FOUND"
            exit 1
        }
    done
'

echo


# ============================================================
# EJECUCION
# ============================================================

if [[ -f "${DONE}" && -d "${BINDIR}" && -s "${C2B}" ]]; then

    echo "DAS Tool ya fue completado para ${GROUP}."
    echo "Se reutilizarán los resultados."

else

    rm -f "${OUTBASE}"*
    rm -rf "${BINDIR}"
    rm -f "${DONE}"

    echo "===== RUNNING DAS TOOL ====="
    echo

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        apptainer exec \
        --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
        --bind "${TMP_HOST}:${TMP_CONTAINER}" \
        --env TMPDIR="${TMP_CONTAINER}" \
        "${DASTOOL}" \
        DAS_Tool \
        -i "${META},${MAXBIN},${CONCOCT}" \
        -l "metabat2,maxbin2,concoct" \
        -c "${ASSEMBLY}" \
        -o "${OUTBASE}" \
        --search_engine diamond \
        --write_bin_evals \
        --write_bins \
        --threads "${SLURM_CPUS_PER_TASK}"

fi


# ============================================================
# VALIDACION
# ============================================================

echo
echo "===== VALIDATING DAS TOOL OUTPUT ====="

if [[ ! -s "${LOG}" ]]; then
    echo "ERROR: DAS Tool log no generado." >&2
    exit 1
fi

if [[ ! -s "${C2B}" ]]; then
    echo "ERROR: DAS Tool contig2bin no generado." >&2
    exit 1
fi

if [[ ! -d "${BINDIR}" ]]; then
    echo "ERROR: DAS Tool bins directory no generado." >&2
    exit 1
fi

shopt -s nullglob
BINS=( "${BINDIR}"/*.fa )
shopt -u nullglob

N_BINS=${#BINS[@]}

if [[ "${N_BINS}" -eq 0 ]]; then
    echo "ERROR: DAS Tool no seleccionó bins." >&2
    exit 1
fi


# ============================================================
# ORIGEN DE LOS BINS
# ============================================================

META_SELECTED=0
MAX_SELECTED=0
CON_SELECTED=0
SUB_SELECTED=0

for BIN in "${BINS[@]}"
do

    NAME=$(basename "${BIN}" .fa)

    if [[ "${NAME}" == L*_bin.* ]]; then
        META_SELECTED=$((META_SELECTED + 1))

    elif [[ "${NAME}" == M*_bin.* ]]; then
        META_SELECTED=$((META_SELECTED + 1))

    elif [[ "${NAME}" == *maxbin2* ]]; then
        MAX_SELECTED=$((MAX_SELECTED + 1))

    elif [[ "${NAME}" == concoct_* ]]; then
        CON_SELECTED=$((CON_SELECTED + 1))
    fi

    if [[ "${NAME}" == *_sub ]]; then
        SUB_SELECTED=$((SUB_SELECTED + 1))
    fi

done


echo "Bins seleccionados: ${N_BINS}"
echo "MetaBAT2:            ${META_SELECTED}"
echo "MaxBin2:             ${MAX_SELECTED}"
echo "CONCOCT:             ${CON_SELECTED}"
echo "Sub-bins (_sub):     ${SUB_SELECTED}"


# ============================================================
# RESUMEN
# ============================================================

SUMMARY="${OUTDIR}/${GROUP}_dastool_selection_summary.tsv"

printf \
"coassembly\ttotal_selected\tmetabat2\tmaxbin2\tconcoct\tsubbins\n" \
> "${SUMMARY}"

printf \
"%s\t%s\t%s\t%s\t%s\t%s\n" \
"${GROUP}" \
"${N_BINS}" \
"${META_SELECTED}" \
"${MAX_SELECTED}" \
"${CON_SELECTED}" \
"${SUB_SELECTED}" \
>> "${SUMMARY}"


echo
echo "===== FIRST DAS TOOL EVALUATIONS ====="

if [[ -s "${EVAL}" ]]; then
    head -n 12 "${EVAL}" | column -s $'\t' -t
fi


touch "${DONE}"


echo
echo "============================================================"
echo "DASTOOL REFINEMENT COMPLETADO"
echo "Group: ${GROUP}"
echo "Bins seleccionados: ${N_BINS}"
echo "Fin: $(date)"
echo "============================================================"
