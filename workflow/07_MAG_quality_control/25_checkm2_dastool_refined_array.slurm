#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00

#SBATCH --array=1-6%2
#SBATCH --job-name=cm2_das_ref

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

BIN_DIR="${SCRATCH_DIR}/12_dastool/${GROUP}/refined/${GROUP}_DASTool_DASTool_bins"

OUT="${SCRATCH_DIR}/13_checkm2_refined/${GROUP}"

RESULT="${OUT}/checkm2_output"
REPORT="${RESULT}/quality_report.tsv"
SUMMARY="${OUT}/${GROUP}_dastool_checkm2.tsv"
UNASSESSED="${OUT}/${GROUP}_checkm2_unassessed.tsv"
DONE="${OUT}/done"

TMP_HOST="${SCRATCH_DIR}/tmp/cm2_das_ref_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
TMP_CONTAINER="/cm2tmp"

mkdir -p "${OUT}"
mkdir -p "${TMP_HOST}"

cleanup() {
    echo
    echo "===== CLEANING TEMPORARY FILES ====="
    rm -rf "${TMP_HOST}"
}

trap cleanup EXIT


echo "============================================================"
echo "CHECKM2 DAS TOOL REFINED"
echo "============================================================"
echo "Group:       ${GROUP}"
echo "Inicio:      $(date)"
echo "Nodo:        $(hostname)"
echo "Job ID:      ${SLURM_JOB_ID}"
echo "Array ID:    ${SLURM_ARRAY_TASK_ID}"
echo "CPUs:        ${SLURM_CPUS_PER_TASK}"
echo "Memoria:     64G"
echo
echo "Bins:"
echo "  ${BIN_DIR}"
echo
echo "CheckM2:"
echo "  ${CHECKM2}"
echo
echo "Database:"
echo "  ${CHECKM2_DB}"
echo
echo "Output:"
echo "  ${OUT}"
echo "============================================================"
echo


# ============================================================
# VALIDACIONES
# ============================================================

if [[ ! -d "${BIN_DIR}" ]]; then
    echo "ERROR: directorio de bins refinados no encontrado." >&2
    exit 1
fi

if [[ ! -s "${CHECKM2}" ]]; then
    echo "ERROR: contenedor CheckM2 no encontrado." >&2
    exit 1
fi

if [[ ! -s "${CHECKM2_DB}" ]]; then
    echo "ERROR: base CheckM2 no encontrada." >&2
    exit 1
fi


N_BINS=$(
    find "${BIN_DIR}" \
        -maxdepth 1 \
        -type f \
        -name "*.fa" \
        | wc -l
)

echo "Bins detectados: ${N_BINS}"

if [[ "${N_BINS}" -eq 0 ]]; then
    echo "ERROR: no se encontraron bins refinados." >&2
    exit 1
fi


# ============================================================
# CHECKM2
# ============================================================

if [[ -f "${DONE}" && -s "${REPORT}" ]]; then

    echo
    echo "CheckM2 ya fue completado para ${GROUP}."
    echo "Se reutilizará el resultado."

else

    rm -rf "${RESULT}"
    rm -f "${SUMMARY}" "${UNASSESSED}" "${DONE}"

    echo
    echo "===== RUNNING CHECKM2 ====="
    echo

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        apptainer exec \
        --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
        --bind "${TMP_HOST}:${TMP_CONTAINER}" \
        --env CHECKM2DB="${CHECKM2_DB}" \
        --env TMPDIR="${TMP_CONTAINER}" \
        "${CHECKM2}" \
        checkm2 predict \
        --input "${BIN_DIR}" \
        --output-directory "${RESULT}" \
        --threads "${SLURM_CPUS_PER_TASK}" \
        -x fa

fi


# ============================================================
# VALIDAR REPORTE
# ============================================================

if [[ ! -s "${REPORT}" ]]; then
    echo "ERROR: quality_report.tsv no generado." >&2
    exit 1
fi


N_REPORT=$(( $(wc -l < "${REPORT}") - 1 ))

if [[ "${N_REPORT}" -gt "${N_BINS}" ]]; then
    echo "ERROR: CheckM2 reportó más bins que la entrada." >&2
    exit 1
fi

N_SKIPPED=$((N_BINS - N_REPORT))


echo
echo "===== CHECKM2 COUNTS ====="
echo "Bins input:         ${N_BINS}"
echo "Bins evaluados:     ${N_REPORT}"
echo "Bins no evaluables: ${N_SKIPPED}"


# ============================================================
# RESUMEN CHECKM2
# ============================================================

awk -F '\t' \
    -v OFS='\t' \
    -v group="${GROUP}" '

    NR == 1 {

        for (i=1; i<=NF; i++) {

            if ($i=="Name")
                name=i

            if ($i=="Completeness")
                comp=i

            if ($i=="Contamination")
                contam=i
        }

        if (!name || !comp || !contam) {
            print "ERROR: columnas CheckM2 no encontradas" > "/dev/stderr"
            exit 2
        }

        print "coassembly",
              "bin",
              "completeness",
              "contamination"

        next
    }

    {
        print group,
              $name,
              $comp,
              $contam
    }

' "${REPORT}" > "${SUMMARY}"


# ============================================================
# IDENTIFICAR NO EVALUABLES
# ============================================================

INPUT_IDS="${OUT}/input_bins.tmp"
REPORT_IDS="${OUT}/reported_bins.tmp"

find "${BIN_DIR}" \
    -maxdepth 1 \
    -type f \
    -name "*.fa" \
    -printf "%f\n" \
    | sed 's/\.fa$//' \
    | sort \
    > "${INPUT_IDS}"

awk -F '\t' '
    NR == 1 {
        for(i=1;i<=NF;i++)
            if($i=="Name")
                name=i
        next
    }

    {
        print $name
    }
' "${REPORT}" \
    | sort \
    > "${REPORT_IDS}"


printf "coassembly\tbin\treason\n" > "${UNASSESSED}"

comm -23 "${INPUT_IDS}" "${REPORT_IDS}" |
while read -r BIN
do
    printf "%s\t%s\t%s\n" \
        "${GROUP}" \
        "${BIN}" \
        "Not reported by CheckM2; inspect CheckM2 log" \
        >> "${UNASSESSED}"
done

rm -f "${INPUT_IDS}" "${REPORT_IDS}"


# ============================================================
# CALIDAD
# ============================================================

echo
echo "===== QUALITY COUNTS ====="

awk -F '\t' '

NR > 1 {

    total++

    comp=$3+0
    contam=$4+0

    if (comp >= 50 && contam < 10)
        candidate++

    if (comp >= 90 && contam < 5)
        strong++
}

END {

    print "Evaluados =", total+0
    print ">=50/<10  =", candidate+0
    print ">=90/<5   =", strong+0
}

' "${SUMMARY}"


echo
echo "===== FIRST RESULTS ====="

head -n 12 "${SUMMARY}" | column -s $'\t' -t


if [[ "${N_SKIPPED}" -gt 0 ]]; then

    echo
    echo "===== UNASSESSED BINS ====="

    column -s $'\t' -t "${UNASSESSED}"

fi


touch "${DONE}"


echo
echo "============================================================"
echo "CHECKM2 DAS TOOL REFINED COMPLETADO"
echo "Group: ${GROUP}"
echo "Bins input: ${N_BINS}"
echo "Bins evaluados: ${N_REPORT}"
echo "No evaluables: ${N_SKIPPED}"
echo "Fin: $(date)"
echo "============================================================"
