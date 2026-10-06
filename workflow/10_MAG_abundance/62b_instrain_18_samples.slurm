#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --array=0-17%4

#SBATCH --job-name=instrain_18
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail


# ============================================================
# CONFIGURACIÓN
# ============================================================

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

MAPBASE="${ROOT}/39_mag_abundance"
OUTBASE="${ROOT}/40_instrain_MAG_profiles"

MANIFEST="${MAPBASE}/samples_hostfree.tsv"

REF="${MAPBASE}/reference/representative_18.prefixed.fna"
STB="${MAPBASE}/reference/representative_18.instrain.stb"

INSTRAIN_SIF="${INSTRAIN_SIF:-/opt/ohpc/pub/containers/BIO/instrain-1.9.0.sif}"


# ============================================================
# IDENTIFICAR MUESTRA
# ============================================================

LINE_NUMBER=$((SLURM_ARRAY_TASK_ID + 2))

LINE=$(
    sed -n "${LINE_NUMBER}p" \
        "${MANIFEST}"
)

if [[ -z "${LINE}" ]]; then
    echo "ERROR: no existe línea ${LINE_NUMBER} del manifiesto." >&2
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

BAM="${MAPBASE}/mapping/${SAMPLE}/${SAMPLE}.MAGs.sorted.bam"
BAI="${BAM}.bai"

OUT="${OUTBASE}/${SAMPLE}.IS"


echo "============================================================"
echo "inStrain - 18 MAGs REPRESENTATIVOS"
echo "Muestra: ${SAMPLE}"
echo "Productor: ${PRODUCER}"
echo "Unidad biológica: ${BIOLOGICAL_UNIT}"
echo "Semana: ${WEEK}"
echo "Array: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Nodo: $(hostname)"
echo "Inicio: $(date)"
echo "============================================================"


# ============================================================
# PREFLIGHT
# ============================================================

for F in \
    "${MANIFEST}" \
    "${REF}" \
    "${STB}" \
    "${BAM}" \
    "${BAI}" \
    "${INSTRAIN_SIF}"
do
    test -s "${F}" || {
        echo "ERROR: falta archivo: ${F}" >&2
        exit 1
    }
done


N_REF=$(
    grep -c '^>' "${REF}"
)

N_STB=$(
    wc -l < "${STB}"
)

N_STB="${N_STB//[[:space:]]/}"

if [[ "${N_REF}" -ne 2463 ]]; then
    echo "ERROR: referencia tiene ${N_REF} contigs; esperábamos 2463." >&2
    exit 1
fi

if [[ "${N_STB}" -ne 2463 ]]; then
    echo "ERROR: STB tiene ${N_STB} entradas; esperábamos 2463." >&2
    exit 1
fi


mkdir -p "${OUTBASE}"


# ============================================================
# RECONOCER PILOTO YA COMPLETADO
# ============================================================

if [[ -d "${OUT}" ]]; then

    EXISTING_GENOME_INFO=$(
        find "${OUT}" \
            -type f \
            -name '*genome_info.tsv' \
            -print \
            -quit
    )

    if [[ -n "${EXISTING_GENOME_INFO}" ]] \
        && [[ -s "${EXISTING_GENOME_INFO}" ]] \
        && [[ -s "${OUT}/Turrialba_profile_provenance.txt" ]]
    then
        touch "${OUT}/PROFILE_COMPLETED"

        echo "Perfil existente validado: ${SAMPLE}"
        echo "No se repetirá."

        echo
        echo "============================================================"
        echo "inStrain ${SAMPLE} YA ESTABA COMPLETADO"
        echo "Fin: $(date)"
        echo "============================================================"

        exit 0
    fi

    echo "ERROR: existe ${OUT}, pero parece incompleto." >&2
    exit 1
fi


# ============================================================
# VERSIÓN
# ============================================================

echo
echo "============================================================"
echo "VERSIÓN"
echo "============================================================"

apptainer exec \
    --bind "${ROOT}:${ROOT}" \
    "${INSTRAIN_SIF}" \
    inStrain profile --version


# ============================================================
# PERFIL
# ============================================================

echo
echo "============================================================"
echo "PERFILANDO ${SAMPLE}"
echo "============================================================"

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    apptainer exec \
        --bind "${ROOT}:${ROOT}" \
        "${INSTRAIN_SIF}" \
        inStrain profile \
            "${BAM}" \
            "${REF}" \
            -o "${OUT}" \
            -p "${SLURM_CPUS_PER_TASK}" \
            -s "${STB}" \
            --min_read_ani 0.95 \
            --min_mapq 2 \
            --pairing_filter paired_only \
            --skip_mm_profiling \
            --skip_plot_generation


# ============================================================
# VALIDAR genome_info
# ============================================================

GENOME_INFO=$(
    find "${OUT}" \
        -type f \
        -name '*genome_info.tsv' \
        -print \
        -quit
)

if [[ -z "${GENOME_INFO}" ]]; then
    echo "ERROR: ${SAMPLE} no generó genome_info.tsv." >&2
    exit 1
fi

test -s "${GENOME_INFO}" || {
    echo "ERROR: genome_info vacío para ${SAMPLE}." >&2
    exit 1
}


N_GENOMES=$(
    awk '
        END {print NR-1}
    ' "${GENOME_INFO}"
)

echo "MAGs reportados por inStrain: ${N_GENOMES}"


# ============================================================
# GUARDAR PROVENANCE
# ============================================================

{
    echo "sample=${SAMPLE}"
    echo "producer=${PRODUCER}"
    echo "biological_unit=${BIOLOGICAL_UNIT}"
    echo "week=${WEEK}"
    echo "bam=${BAM}"
    echo "reference=${REF}"
    echo "stb=${STB}"
    echo "container=${INSTRAIN_SIF}"
    echo "min_read_ani=0.95"
    echo "min_mapq=2"
    echo "pairing_filter=paired_only"
    echo "skip_mm_profiling=yes"
    echo "skip_plot_generation=yes"
    echo "database_mode=no"
    echo "genomes_reported=${N_GENOMES}"
    echo "job_id=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
    echo "completed=$(date)"
} > "${OUT}/Turrialba_profile_provenance.txt"


touch "${OUT}/PROFILE_COMPLETED"


# ============================================================
# RESUMEN
# ============================================================

echo
echo "============================================================"
echo "RESULTADO ${SAMPLE}"
echo "============================================================"

echo "genome_info:"
echo "${GENOME_INFO}"

echo
echo "MAGs reportados: ${N_GENOMES}"

echo
echo "Primeras filas:"
sed -n '1,6p' "${GENOME_INFO}"


echo
echo "============================================================"
echo "PASO 62b COMPLETADO: inStrain ${SAMPLE}"
echo "Fin: $(date)"
echo "============================================================"
