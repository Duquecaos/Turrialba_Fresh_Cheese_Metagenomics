#!/bin/bash

#SBATCH --partition=shared_64
#SBATCH --qos=sm

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=192G
#SBATCH --time=23:50:00

#SBATCH --array=1-6%2

#SBATCH --job-name=megahit_coasm

#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ============================================================
# CONFIGURACION
# ============================================================

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

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

GROUP_FILE="${PROJECT_HOME}/metadata/coassembly_groups.tsv"

ASSEMBLY_ROOT="${SCRATCH_DIR}/05_assembly"

mkdir -p "${ASSEMBLY_ROOT}"


# ============================================================
# DEFINIR COASSEMBLY
#
# Array:
# 1 = L1
# 2 = L2
# 3 = L3
# 4 = M1
# 5 = M2
# 6 = M3
#
# IMPORTANTE:
# No usar la variable GROUPS.
# GROUPS es una variable especial/reservada de Bash.
# ============================================================

COASSEMBLY_NAMES=(L1 L2 L3 M1 M2 M3)

ARRAY_INDEX=$((SLURM_ARRAY_TASK_ID - 1))

GROUP_NAME="${COASSEMBLY_NAMES[${ARRAY_INDEX}]}"

OUT="${ASSEMBLY_ROOT}/${GROUP_NAME}"


# ============================================================
# VALIDAR METADATA
# ============================================================

if [[ ! -s "${GROUP_FILE}" ]]; then
    echo "ERROR: no existe metadata de coassembly:"
    echo "${GROUP_FILE}"
    exit 1
fi


# ============================================================
# OBTENER LAS 3 MUESTRAS DE LA UNIDAD BIOLOGICA
# ============================================================

mapfile -t GROUP_SAMPLES < <(
    awk -F '\t' \
        -v group="${GROUP_NAME}" \
        'NR > 1 && $1 == group {print $2}' \
        "${GROUP_FILE}"
)


# ============================================================
# INFORMACION PRELIMINAR
# ============================================================

echo "============================================================"
echo "MEGAHIT LONGITUDINAL COASSEMBLY"
echo "============================================================"
echo "Array task:    ${SLURM_ARRAY_TASK_ID}"
echo "Array index:   ${ARRAY_INDEX}"
echo "Coassembly:    ${GROUP_NAME}"
echo "Metadata:      ${GROUP_FILE}"
echo "Samples found: ${#GROUP_SAMPLES[@]}"
echo "============================================================"
echo


# ============================================================
# VALIDAR NUMERO DE MUESTRAS
# ============================================================

if [[ "${#GROUP_SAMPLES[@]}" -ne 3 ]]; then

    echo "ERROR: ${GROUP_NAME} no tiene exactamente 3 muestras."
    echo
    echo "Muestras encontradas:"

    printf '  %s\n' "${GROUP_SAMPLES[@]}"

    echo
    echo "Metadata utilizada:"
    cat "${GROUP_FILE}"

    exit 1
fi


# ============================================================
# CONSTRUIR LISTAS R1 / R2
# ============================================================

R1_FILES=()
R2_FILES=()

for SAMPLE in "${GROUP_SAMPLES[@]}"
do

    R1="${SCRATCH_DIR}/04_host_removed/${SAMPLE}/${SAMPLE}_1.hostfree.fastq.gz"
    R2="${SCRATCH_DIR}/04_host_removed/${SAMPLE}/${SAMPLE}_2.hostfree.fastq.gz"

    if [[ ! -s "${R1}" ]]; then
        echo "ERROR: no existe R1 para ${SAMPLE}:"
        echo "${R1}"
        exit 1
    fi

    if [[ ! -s "${R2}" ]]; then
        echo "ERROR: no existe R2 para ${SAMPLE}:"
        echo "${R2}"
        exit 1
    fi

    R1_FILES+=("${R1}")
    R2_FILES+=("${R2}")

done


# ============================================================
# MEGAHIT ACEPTA MULTIPLES LIBRERIAS SEPARADAS POR COMAS
# ============================================================

R1_LIST=$(IFS=, ; echo "${R1_FILES[*]}")
R2_LIST=$(IFS=, ; echo "${R2_FILES[*]}")


# ============================================================
# INFORMACION COMPLETA DEL JOB
# ============================================================

echo "============================================================"
echo "COASSEMBLY CONFIGURATION"
echo "============================================================"
echo "Coassembly:    ${GROUP_NAME}"
echo "Job ID:        ${SLURM_JOB_ID}"
echo "Array ID:      ${SLURM_ARRAY_TASK_ID}"
echo "Nodo:          $(hostname)"
echo "Inicio:        $(date)"
echo "CPUs:          ${SLURM_CPUS_PER_TASK}"
echo "Memoria:       192G"
echo
echo "Muestras:"
printf '  %s\n' "${GROUP_SAMPLES[@]}"
echo
echo "R1:"
printf '  %s\n' "${R1_FILES[@]}"
echo
echo "R2:"
printf '  %s\n' "${R2_FILES[@]}"
echo
echo "Output:"
echo "  ${OUT}"
echo "============================================================"
echo


# ============================================================
# TAMANOS DE INPUT
# ============================================================

echo "===== INPUT FILE SIZES ====="

for FILE in "${R1_FILES[@]}" "${R2_FILES[@]}"
do
    ls -lh "${FILE}"
done

echo


# ============================================================
# SI YA EXISTE UN ASSEMBLY TERMINADO, NO REPETIR
# ============================================================

if [[ -s "${OUT}/final.contigs.fa" ]]; then

    echo "============================================================"
    echo "COASSEMBLY YA COMPLETADO"
    echo "${OUT}/final.contigs.fa"
    echo "Se omite MEGAHIT."
    echo "============================================================"

    exit 0
fi


# ============================================================
# MEGAHIT
#
# - coassembly longitudinal
# - tres tiempos de la misma muestra biologica
# - meta-sensitive
# - contigs finales >= 1000 bp
# ============================================================

echo "===== RUNNING MEGAHIT ====="
echo


# ============================================================
# SI EXISTE OUTPUT PARCIAL DE MEGAHIT, INTENTAR CONTINUAR
# ============================================================

if [[ -d "${OUT}" && ! -s "${OUT}/final.contigs.fa" ]]; then

    echo "Output parcial detectado:"
    echo "${OUT}"
    echo
    echo "Intentando continuar con --continue..."
    echo

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        apptainer exec \
        --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
        "${MEGAHIT}" \
        megahit \
        --continue \
        -o "${OUT}"

else

    echo "Iniciando assembly nuevo..."
    echo

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        apptainer exec \
        --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
        "${MEGAHIT}" \
        megahit \
        -1 "${R1_LIST}" \
        -2 "${R2_LIST}" \
        -o "${OUT}" \
        --presets meta-sensitive \
        --min-contig-len 1000 \
        --num-cpu-threads "${SLURM_CPUS_PER_TASK}" \
        --memory 170000000000

fi

echo


# ============================================================
# VALIDAR ASSEMBLY
# ============================================================

CONTIGS="${OUT}/final.contigs.fa"

if [[ ! -s "${CONTIGS}" ]]; then
    echo "ERROR: MEGAHIT no produjo final.contigs.fa" >&2
    exit 1
fi


# ============================================================
# ESTADISTICAS BASICAS DEL ASSEMBLY
# ============================================================

echo "===== CALCULATING ASSEMBLY STATISTICS ====="

N_CONTIGS=$(grep -c '^>' "${CONTIGS}")

TOTAL_BP=$(
    awk '
        /^>/ {next}
        {
            gsub(/[[:space:]]/, "")
            total += length($0)
        }
        END {
            print total
        }
    ' "${CONTIGS}"
)

MAX_CONTIG=$(
    awk '
        /^>/ {
            if (len > max) max=len
            len=0
            next
        }
        {
            gsub(/[[:space:]]/, "")
            len += length($0)
        }
        END {
            if (len > max) max=len
            print max
        }
    ' "${CONTIGS}"
)

SUMMARY="${OUT}/assembly_summary.tsv"

printf "coassembly\tn_samples\tn_contigs\ttotal_bp\tmax_contig_bp\n" \
    > "${SUMMARY}"

printf "%s\t3\t%s\t%s\t%s\n" \
    "${GROUP_NAME}" \
    "${N_CONTIGS}" \
    "${TOTAL_BP}" \
    "${MAX_CONTIG}" \
    >> "${SUMMARY}"


# ============================================================
# RESULTADOS
# ============================================================

echo
echo "===== ASSEMBLY OUTPUT ====="
echo

ls -lh "${CONTIGS}"

echo

echo "===== ASSEMBLY SUMMARY ====="
echo

column -t "${SUMMARY}"

echo


# ============================================================
# FINAL
# ============================================================

echo "============================================================"
echo "MEGAHIT COASSEMBLY COMPLETADO"
echo "Coassembly: ${GROUP_NAME}"
echo "Fin: $(date)"
echo "============================================================"
