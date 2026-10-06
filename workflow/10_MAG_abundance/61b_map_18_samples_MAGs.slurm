#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=0-17%4

#SBATCH --job-name=map_18_MAGs
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail


# ============================================================
# CONFIGURACIÓN
# ============================================================

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"
BASE="${ROOT}/39_mag_abundance"

MANIFEST="${BASE}/samples_hostfree.tsv"

REF="${BASE}/reference/representative_18.prefixed.fna"
INDEX_PREFIX="${BASE}/reference/bowtie2_index/representative_18"

BOWTIE2_SIF="${BOWTIE2_SIF:-/opt/ohpc/pub/containers/BIO/bowtie2-2.5.4.sif}"

SAMTOOLS_SIF=$(
    find "${CONTAINER_ROOT:-/opt/ohpc/pub/containers/BIO}" \
        -maxdepth 1 \
        -type f \
        -iname 'samtools*.sif' \
        | sort -V \
        | tail -n 1
)

BOWTIE_THREADS=5
SORT_THREADS=2


# ============================================================
# OBTENER MUESTRA DEL ARRAY
# ============================================================

LINE_NUMBER=$((SLURM_ARRAY_TASK_ID + 2))

LINE=$(
    sed -n "${LINE_NUMBER}p" \
        "${MANIFEST}"
)

if [[ -z "${LINE}" ]]; then
    echo "ERROR: no existe línea ${LINE_NUMBER} en ${MANIFEST}" >&2
    exit 1
fi

IFS=$'\t' read -r \
    SAMPLE \
    PRODUCER \
    BIOLOGICAL_UNIT \
    WEEK \
    R1 \
    R2 \
    <<< "${LINE}"


FINAL="${BASE}/mapping/${SAMPLE}"
WORK="${BASE}/mapping/.working_${SAMPLE}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

BAM="${WORK}/${SAMPLE}.MAGs.sorted.bam"
BT2LOG="${WORK}/${SAMPLE}.bowtie2.log"


echo "============================================================"
echo "MAPEO CONTRA 18 MAGs"
echo "Array job: ${SLURM_ARRAY_JOB_ID}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Muestra: ${SAMPLE}"
echo "Productor: ${PRODUCER}"
echo "Unidad biológica: ${BIOLOGICAL_UNIT}"
echo "Semana: ${WEEK}"
echo "Nodo: $(hostname)"
echo "Inicio: $(date)"
echo "============================================================"


# ============================================================
# SI YA ESTÁ COMPLETADO, NO REPETIR
# ============================================================

if [[ -f "${FINAL}/MAPPING_COMPLETED" ]]; then
    echo "Muestra ya completada: ${SAMPLE}"
    echo "Se omite sin repetir el mapeo."
    echo "MAPEO ${SAMPLE} YA ESTABA COMPLETADO"
    exit 0
fi

if [[ -e "${FINAL}" ]]; then
    echo "ERROR: existe ${FINAL} pero no contiene MAPPING_COMPLETED." >&2
    exit 1
fi


# ============================================================
# PREFLIGHT
# ============================================================

test -s "${MANIFEST}" || {
    echo "ERROR: falta manifiesto." >&2
    exit 1
}

test -s "${REF}" || {
    echo "ERROR: falta referencia." >&2
    exit 1
}

test -s "${R1}" || {
    echo "ERROR: falta R1: ${R1}" >&2
    exit 1
}

test -s "${R2}" || {
    echo "ERROR: falta R2: ${R2}" >&2
    exit 1
}

test -s "${BOWTIE2_SIF}" || {
    echo "ERROR: falta contenedor Bowtie2." >&2
    exit 1
}

test -n "${SAMTOOLS_SIF}" || {
    echo "ERROR: no se encontró contenedor Samtools." >&2
    exit 1
}

test -s "${SAMTOOLS_SIF}" || {
    echo "ERROR: contenedor Samtools inválido: ${SAMTOOLS_SIF}" >&2
    exit 1
}

for SUFFIX in \
    1.bt2 \
    2.bt2 \
    3.bt2 \
    4.bt2 \
    rev.1.bt2 \
    rev.2.bt2
do
    test -s "${INDEX_PREFIX}.${SUFFIX}" || {
        echo "ERROR: falta ${INDEX_PREFIX}.${SUFFIX}" >&2
        exit 1
    }
done


# ============================================================
# WORKDIR
# ============================================================

mkdir -p "${BASE}/mapping"

rm -rf "${WORK}"
mkdir -p "${WORK}"


# ============================================================
# MAPEO BOWTIE2 -> BAM ORDENADO
# ============================================================

echo
echo "============================================================"
echo "BOWTIE2 -> SAMTOOLS SORT"
echo "============================================================"

export ROOT
export BOWTIE2_SIF
export SAMTOOLS_SIF
export INDEX_PREFIX
export R1
export R2
export BAM
export BT2LOG
export SAMPLE
export BOWTIE_THREADS
export SORT_THREADS

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    bash -c '
        set -euo pipefail

        apptainer exec \
            --bind "${ROOT}:${ROOT}" \
            "${BOWTIE2_SIF}" \
            bowtie2 \
                --very-sensitive \
                --threads "${BOWTIE_THREADS}" \
                --rg-id "${SAMPLE}" \
                --rg "SM:${SAMPLE}" \
                --no-unal \
                -x "${INDEX_PREFIX}" \
                -1 "${R1}" \
                -2 "${R2}" \
                2> "${BT2LOG}" \
        | \
        apptainer exec \
            --bind "${ROOT}:${ROOT}" \
            "${SAMTOOLS_SIF}" \
            samtools sort \
                -@ "${SORT_THREADS}" \
                -m 2G \
                -o "${BAM}" \
                -
    '


# ============================================================
# VALIDAR BAM
# ============================================================

test -s "${BAM}" || {
    echo "ERROR: BAM vacío para ${SAMPLE}." >&2
    exit 1
}

srun --exact \
    -n 1 \
    -c 1 \
    apptainer exec \
        --bind "${ROOT}:${ROOT}" \
        "${SAMTOOLS_SIF}" \
        samtools quickcheck \
            -v "${BAM}"

echo "quickcheck: OK"


# ============================================================
# INDEXAR
# ============================================================

srun --exact \
    -n 1 \
    -c 3 \
    apptainer exec \
        --bind "${ROOT}:${ROOT}" \
        "${SAMTOOLS_SIF}" \
        samtools index \
            -@ 2 \
            "${BAM}"

test -s "${BAM}.bai" || {
    echo "ERROR: falta índice BAM." >&2
    exit 1
}


# ============================================================
# QC
# ============================================================

srun --exact \
    -n 1 \
    -c 2 \
    apptainer exec \
        --bind "${ROOT}:${ROOT}" \
        "${SAMTOOLS_SIF}" \
        samtools flagstat \
            -@ 1 \
            "${BAM}" \
        > "${WORK}/${SAMPLE}.flagstat.txt"


srun --exact \
    -n 1 \
    -c 1 \
    apptainer exec \
        --bind "${ROOT}:${ROOT}" \
        "${SAMTOOLS_SIF}" \
        samtools idxstats \
            "${BAM}" \
        > "${WORK}/${SAMPLE}.idxstats.tsv"


# ============================================================
# VALIDAR LAS 2463 REFERENCIAS
# ============================================================

N_SQ=$(
    apptainer exec \
        --bind "${ROOT}:${ROOT}" \
        "${SAMTOOLS_SIF}" \
        samtools view -H "${BAM}" \
    | awk '
        $1=="@SQ" {n++}
        END {print n+0}
    '
)

if [[ "${N_SQ}" -ne 2463 ]]; then
    echo "ERROR: ${SAMPLE} tiene ${N_SQ} @SQ; esperábamos 2463." >&2
    exit 1
fi


# ============================================================
# EXTRAER RESUMEN DE BOWTIE2
# ============================================================

TOTAL_PAIRS=$(
    awk '
        /^[0-9]+ reads; of these:/ {
            print $1
            exit
        }
    ' "${BT2LOG}"
)

OVERALL_RATE=$(
    awk '
        /overall alignment rate/ {
            print $1
            exit
        }
    ' "${BT2LOG}"
)

CONCORDANT_1=$(
    awk '
        /aligned concordantly exactly 1 time/ {
            print $1
            exit
        }
    ' "${BT2LOG}"
)

CONCORDANT_MULTI=$(
    awk '
        /aligned concordantly >1 times/ {
            print $1
            exit
        }
    ' "${BT2LOG}"
)


# ============================================================
# NÚMERO DE REGISTROS MAPEADOS EN BAM
# ============================================================

MAPPED_RECORDS=$(
    apptainer exec \
        --bind "${ROOT}:${ROOT}" \
        "${SAMTOOLS_SIF}" \
        samtools view \
            -c \
            "${BAM}"
)


# ============================================================
# RESUMEN TABULAR POR MUESTRA
# ============================================================

SUMMARY="${WORK}/${SAMPLE}.mapping_summary.tsv"

printf \
'sample\tproducer\tbiological_unit\tweek\tread_pairs\tbowtie2_overall_alignment_rate\tconcordant_unique_pairs\tconcordant_multimapping_pairs\tmapped_BAM_records\n' \
> "${SUMMARY}"

printf \
'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${SAMPLE}" \
    "${PRODUCER}" \
    "${BIOLOGICAL_UNIT}" \
    "${WEEK}" \
    "${TOTAL_PAIRS}" \
    "${OVERALL_RATE}" \
    "${CONCORDANT_1}" \
    "${CONCORDANT_MULTI}" \
    "${MAPPED_RECORDS}" \
    >> "${SUMMARY}"


# ============================================================
# PROVENANCE
# ============================================================

{
    echo "sample=${SAMPLE}"
    echo "producer=${PRODUCER}"
    echo "biological_unit=${BIOLOGICAL_UNIT}"
    echo "week=${WEEK}"
    echo "R1=${R1}"
    echo "R2=${R2}"
    echo "reference=${REF}"
    echo "index=${INDEX_PREFIX}"
    echo "bowtie2_container=${BOWTIE2_SIF}"
    echo "samtools_container=${SAMTOOLS_SIF}"
    echo "bowtie2_mode=very-sensitive"
    echo "bowtie2_threads=${BOWTIE_THREADS}"
    echo "samtools_sort_threads=${SORT_THREADS}"
    echo "unmapped_records_written=no"
    echo "job_id=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
    echo "completed=$(date)"
} > "${WORK}/mapping_provenance.txt"

touch "${WORK}/MAPPING_COMPLETED"


# ============================================================
# CERRAR ATÓMICAMENTE
# ============================================================

mv "${WORK}" "${FINAL}"


echo
echo "============================================================"
echo "RESUMEN ${SAMPLE}"
echo "============================================================"

column -t -s $'\t' \
    "${FINAL}/${SAMPLE}.mapping_summary.tsv"

echo
echo "BAM:"
ls -lh \
    "${FINAL}/${SAMPLE}.MAGs.sorted.bam" \
    "${FINAL}/${SAMPLE}.MAGs.sorted.bam.bai"

echo
echo "============================================================"
echo "MAPEO ${SAMPLE} COMPLETADO"
echo "Fin: $(date)"
echo "============================================================"
