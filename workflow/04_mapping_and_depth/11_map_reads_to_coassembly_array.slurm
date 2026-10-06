#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=06:00:00

#SBATCH --array=1-18%4
#SBATCH --job-name=map_coasm

#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err

#SBATCH --mail-type=FAIL

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

LINE=$((SLURM_ARRAY_TASK_ID + 1))

SAMPLE=$(awk -F '\t' -v line="${LINE}" \
    'NR==line {print $1}' "${SAMPLES}")

if [[ -z "${SAMPLE}" ]]; then
    echo "ERROR: no se pudo obtener sample_id para array ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

# Ejemplos:
# L1_1 -> L1
# L1_2 -> L1
# M3_3 -> M3

GROUP_NAME="${SAMPLE%_*}"

# ============================================================
# INPUTS
# ============================================================

R1="${SCRATCH_DIR}/04_host_removed/${SAMPLE}/${SAMPLE}_1.hostfree.fastq.gz"
R2="${SCRATCH_DIR}/04_host_removed/${SAMPLE}/${SAMPLE}_2.hostfree.fastq.gz"

ASSEMBLY_DIR="${SCRATCH_DIR}/05_assembly/${GROUP_NAME}"

CONTIGS="${ASSEMBLY_DIR}/final.contigs.fa"

INDEX_PREFIX="${ASSEMBLY_DIR}/bowtie2_index/${GROUP_NAME}"

# ============================================================
# OUTPUT
# ============================================================

MAP_ROOT="${SCRATCH_DIR}/06_mapping"
OUT="${MAP_ROOT}/${GROUP_NAME}/${SAMPLE}"

mkdir -p "${OUT}"

BAM="${OUT}/${SAMPLE}.sorted.bam"
BAI="${BAM}.bai"

BT2_LOG="${OUT}/${SAMPLE}_bowtie2.log"
FLAGSTAT="${OUT}/${SAMPLE}_flagstat.txt"
IDXSTATS="${OUT}/${SAMPLE}_idxstats.txt"
SUMMARY="${OUT}/${SAMPLE}_mapping_summary.tsv"

# ============================================================
# INFORMACION
# ============================================================

echo "============================================================"
echo "MAP HOST-FREE READS TO COASSEMBLY"
echo "============================================================"
echo "Sample:       ${SAMPLE}"
echo "Coassembly:   ${GROUP_NAME}"
echo "Array ID:     ${SLURM_ARRAY_TASK_ID}"
echo "Job ID:       ${SLURM_JOB_ID}"
echo "Nodo:         $(hostname)"
echo "Inicio:       $(date)"
echo "CPUs:         ${SLURM_CPUS_PER_TASK}"
echo
echo "R1:           ${R1}"
echo "R2:           ${R2}"
echo "Contigs:      ${CONTIGS}"
echo "Index:        ${INDEX_PREFIX}"
echo "BAM:          ${BAM}"
echo "============================================================"
echo

# ============================================================
# VALIDACIONES
# ============================================================

if [[ ! -s "${R1}" ]]; then
    echo "ERROR: R1 no existe: ${R1}" >&2
    exit 1
fi

if [[ ! -s "${R2}" ]]; then
    echo "ERROR: R2 no existe: ${R2}" >&2
    exit 1
fi

if [[ ! -s "${CONTIGS}" ]]; then
    echo "ERROR: contigs no existen: ${CONTIGS}" >&2
    exit 1
fi

if [[ ! -s "${INDEX_PREFIX}.1.bt2" && \
      ! -s "${INDEX_PREFIX}.1.bt2l" ]]; then
    echo "ERROR: índice Bowtie2 no encontrado para ${GROUP_NAME}" >&2
    exit 1
fi

echo "===== INPUT FILES ====="
ls -lh "${R1}" "${R2}" "${CONTIGS}"
echo

# ============================================================
# SI YA EXISTE BAM COMPLETO, NO REPETIR
# ============================================================

if [[ -s "${BAM}" && -s "${BAI}" && -s "${BT2_LOG}" ]]; then

    echo "BAM e índice ya existen para ${SAMPLE}."
    echo "Se omite el mapeo."

else

    rm -f "${BAM}" "${BAI}"

    echo "===== RUNNING BOWTIE2 + SAMTOOLS SORT ====="
    echo

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        bash -c "
        set -o pipefail

        apptainer exec \
            --bind '${SCRATCH_DIR}:${SCRATCH_DIR}' \
            '${BOWTIE2}' \
            bowtie2 \
            --very-sensitive \
            --threads 6 \
            --no-unal \
            -x '${INDEX_PREFIX}' \
            -1 '${R1}' \
            -2 '${R2}' \
            2> '${BT2_LOG}' \
        |
        apptainer exec \
            --bind '${SCRATCH_DIR}:${SCRATCH_DIR}' \
            '${SAMTOOLS}' \
            samtools sort \
            -@ 2 \
            -m 4G \
            -o '${BAM}' \
            -
        "

    echo

    # ========================================================
    # INDEXAR BAM
    # ========================================================

    echo "===== INDEXING BAM ====="

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        apptainer exec \
        --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
        "${SAMTOOLS}" \
        samtools index \
        -@ "${SLURM_CPUS_PER_TASK}" \
        "${BAM}"

fi

# ============================================================
# VALIDAR BAM
# ============================================================

if [[ ! -s "${BAM}" ]]; then
    echo "ERROR: no se generó BAM para ${SAMPLE}" >&2
    exit 1
fi

if [[ ! -s "${BAI}" ]]; then
    echo "ERROR: no se generó BAI para ${SAMPLE}" >&2
    exit 1
fi

# ============================================================
# BAM QC
# ============================================================

echo
echo "===== FLAGSTAT ====="

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    apptainer exec \
    --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
    "${SAMTOOLS}" \
    samtools flagstat \
    -@ "${SLURM_CPUS_PER_TASK}" \
    "${BAM}" \
    > "${FLAGSTAT}"

cat "${FLAGSTAT}"

echo

echo "===== IDXSTATS ====="

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    apptainer exec \
    --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
    "${SAMTOOLS}" \
    samtools idxstats \
    "${BAM}" \
    > "${IDXSTATS}"

# ============================================================
# EXTRAER ALIGNMENT RATE DE BOWTIE2
# ============================================================

ALIGN_RATE=$(awk '
    /overall alignment rate/ {
        gsub("%","",$1)
        print $1
    }
' "${BT2_LOG}")

if [[ -z "${ALIGN_RATE}" ]]; then
    ALIGN_RATE="NA"
fi

# ============================================================
# SUMMARY
# ============================================================

printf "sample\tcoassembly\tbowtie2_alignment_rate\tbam_size_bytes\n" \
    > "${SUMMARY}"

BAM_SIZE=$(stat -c%s "${BAM}")

printf "%s\t%s\t%s\t%s\n" \
    "${SAMPLE}" \
    "${GROUP_NAME}" \
    "${ALIGN_RATE}" \
    "${BAM_SIZE}" \
    >> "${SUMMARY}"

echo
echo "===== MAPPING SUMMARY ====="

column -t "${SUMMARY}"

echo

echo "===== OUTPUT FILES ====="

ls -lh "${OUT}"

echo

echo "============================================================"
echo "MAPPING COMPLETADO"
echo "Sample:     ${SAMPLE}"
echo "Coassembly: ${GROUP_NAME}"
echo "Fin: $(date)"
echo "============================================================"
