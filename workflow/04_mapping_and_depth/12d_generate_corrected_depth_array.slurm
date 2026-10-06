#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=03:00:00

#SBATCH --array=1-6%3
#SBATCH --job-name=depth_corrected

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

export OMP_NUM_THREADS=3

NEW_METABAT="${SCRATCH_DIR}/containers/metabat2-2.18-r3.sif"

COASSEMBLY_NAMES=(L1 L2 L3 M1 M2 M3)

ARRAY_INDEX=$((SLURM_ARRAY_TASK_ID - 1))

GROUP="${COASSEMBLY_NAMES[${ARRAY_INDEX}]}"


# ============================================================
# INPUT
# ============================================================

ASSEMBLY="${SCRATCH_DIR}/05_assembly/${GROUP}/final.contigs.fa"

MAPROOT="${SCRATCH_DIR}/06_mapping/${GROUP}"

S1="${GROUP}_1"
S2="${GROUP}_2"
S3="${GROUP}_3"

BAM1="${MAPROOT}/${S1}/${S1}.sorted.bam"
BAM2="${MAPROOT}/${S2}/${S2}.sorted.bam"
BAM3="${MAPROOT}/${S3}/${S3}.sorted.bam"


# ============================================================
# OUTPUT
# ============================================================

OUT="${SCRATCH_DIR}/07_depth_corrected/${GROUP}"

mkdir -p "${OUT}"

DEPTH="${OUT}/${GROUP}_depth.txt"


# ============================================================
# INFORMACION
# ============================================================

echo "============================================================"
echo "CORRECTED CONTIG DEPTH"
echo "============================================================"
echo "Group:       ${GROUP}"
echo "Job ID:      ${SLURM_JOB_ID}"
echo "Array ID:    ${SLURM_ARRAY_TASK_ID}"
echo "Nodo:        $(hostname)"
echo "Inicio:      $(date)"
echo
echo "Container:"
echo "  ${NEW_METABAT}"
echo
echo "Assembly:"
echo "  ${ASSEMBLY}"
echo
echo "BAMs:"
echo "  ${BAM1}"
echo "  ${BAM2}"
echo "  ${BAM3}"
echo
echo "Output:"
echo "  ${DEPTH}"
echo "============================================================"
echo


# ============================================================
# VALIDAR CONTENEDOR
# ============================================================

if [[ ! -s "${NEW_METABAT}" ]]; then
    echo "ERROR: no existe el nuevo contenedor MetaBAT2:" >&2
    echo "${NEW_METABAT}" >&2
    exit 1
fi


# ============================================================
# VALIDAR ASSEMBLY
# ============================================================

if [[ ! -s "${ASSEMBLY}" ]]; then
    echo "ERROR: no existe assembly para ${GROUP}" >&2
    exit 1
fi


# ============================================================
# VALIDAR BAMS
# ============================================================

for BAM in "${BAM1}" "${BAM2}" "${BAM3}"
do

    if [[ ! -s "${BAM}" ]]; then
        echo "ERROR: BAM no encontrado:" >&2
        echo "${BAM}" >&2
        exit 1
    fi

    if [[ ! -s "${BAM}.bai" ]]; then
        echo "ERROR: BAM index no encontrado:" >&2
        echo "${BAM}.bai" >&2
        exit 1
    fi

done


# ============================================================
# MOSTRAR VERSION EXACTA
# ============================================================

echo "===== JGI VERSION ====="

apptainer exec \
    --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
    "${NEW_METABAT}" \
    jgi_summarize_bam_contig_depths \
    --help 2>&1 | head -2 || true

echo


# ============================================================
# SI YA EXISTE UN RESULTADO VALIDADO, NO REPETIR
# ============================================================

if [[ -s "${DEPTH}" ]]; then

    echo "Depth existente detectado."
    echo "Será validado antes de decidir si repetir."

    EXISTING_NEG=$(
        awk '
            NR > 1 {
                for (i=3; i<=NF; i++) {
                    if (($i+0) < 0) n++
                }
            }
            END {print n+0}
        ' "${DEPTH}"
    )

    if [[ "${EXISTING_NEG}" -eq 0 ]]; then

        echo "Depth existente sin valores negativos."
        echo "Se reutiliza."
        SKIP_JGI=1

    else

        echo "Depth existente inválido."
        echo "Se recalculará."
        rm -f "${DEPTH}"
        SKIP_JGI=0

    fi

else

    SKIP_JGI=0

fi


# ============================================================
# JGI SUMMARIZE
# ============================================================

if [[ "${SKIP_JGI}" -eq 0 ]]; then

    echo "===== RUNNING JGI ====="
    echo

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        apptainer exec \
        --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
        "${NEW_METABAT}" \
        jgi_summarize_bam_contig_depths \
        --referenceFasta "${ASSEMBLY}" \
        --outputDepth "${DEPTH}" \
        "${BAM1}" \
        "${BAM2}" \
        "${BAM3}"

fi

echo


# ============================================================
# VALIDAR NUMERO DE CONTIGS
# ============================================================

N_ASSEMBLY=$(grep -c '^>' "${ASSEMBLY}")

N_DEPTH=$(( $(wc -l < "${DEPTH}") - 1 ))

echo "===== CONTIG COUNTS ====="
echo "Assembly: ${N_ASSEMBLY}"
echo "Depth:    ${N_DEPTH}"

if [[ "${N_ASSEMBLY}" -ne "${N_DEPTH}" ]]; then

    echo "ERROR: número de contigs no coincide." >&2
    exit 1

fi

echo


# ============================================================
# VALIDAR VALORES NEGATIVOS
# ============================================================

NEG=$(
    awk '
        NR > 1 {
            for (i=3; i<=NF; i++) {
                if (($i+0) < 0) n++
            }
        }
        END {print n+0}
    ' "${DEPTH}"
)

echo "Negative numeric values: ${NEG}"

if [[ "${NEG}" -ne 0 ]]; then

    echo "ERROR: valores negativos detectados." >&2
    exit 1

fi


# ============================================================
# VALIDAR NAN / INF
# ============================================================

NONFINITE=$(
    awk '
        NR > 1 {
            for (i=3; i<=NF; i++) {
                x=tolower($i)
                if (x ~ /nan/ || x ~ /inf/) n++
            }
        }
        END {print n+0}
    ' "${DEPTH}"
)

echo "NaN / Inf values: ${NONFINITE}"

if [[ "${NONFINITE}" -ne 0 ]]; then

    echo "ERROR: valores no finitos detectados." >&2
    exit 1

fi


# ============================================================
# MAXIMOS DE COBERTURA
#
# Columnas:
# 3 = totalAvgDepth
# 4 = sample 1 depth
# 6 = sample 2 depth
# 8 = sample 3 depth
#
# NO usamos columnas 5,7,9 porque son varianzas.
# ============================================================

echo
echo "===== MAX DEPTH VALUES ====="

awk '
    NR > 1 {

        if ($3 > maxt) maxt=$3
        if ($4 > max1) max1=$4
        if ($6 > max2) max2=$6
        if ($8 > max3) max3=$8

    }

    END {

        printf "max totalAvgDepth = %.6f\n", maxt
        printf "max sample_1      = %.6f\n", max1
        printf "max sample_2      = %.6f\n", max2
        printf "max sample_3      = %.6f\n", max3

    }
' "${DEPTH}"

echo


# ============================================================
# PRIMERAS FILAS
# ============================================================

echo "===== FIRST ROWS ====="

head -n 3 "${DEPTH}"

echo


# ============================================================
# OUTPUT
# ============================================================

echo "===== OUTPUT ====="

ls -lh "${DEPTH}"

echo

echo "============================================================"
echo "CORRECTED CONTIG DEPTH COMPLETADO"
echo "Group: ${GROUP}"
echo "Fin: $(date)"
echo "============================================================"
