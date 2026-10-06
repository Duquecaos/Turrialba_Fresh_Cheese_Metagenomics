#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00

#SBATCH --array=1-6%3
#SBATCH --job-name=metabat2_bins

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

NEW_METABAT="${SCRATCH_DIR}/containers/metabat2-2.18-r3.sif"

COASSEMBLY_NAMES=(L1 L2 L3 M1 M2 M3)

ARRAY_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
GROUP="${COASSEMBLY_NAMES[${ARRAY_INDEX}]}"


# ============================================================
# INPUT
# ============================================================

ASSEMBLY="${SCRATCH_DIR}/05_assembly/${GROUP}/final.contigs.fa"

DEPTH="${SCRATCH_DIR}/07_depth_corrected/${GROUP}/${GROUP}_depth.txt"


# ============================================================
# OUTPUT
# ============================================================

BIN_ROOT="${SCRATCH_DIR}/08_binning"

OUT="${BIN_ROOT}/${GROUP}/metabat2"

mkdir -p "${OUT}"

BIN_PREFIX="${OUT}/${GROUP}_bin"

SUMMARY="${OUT}/${GROUP}_metabat2_summary.tsv"


# ============================================================
# INFORMACION
# ============================================================

echo "============================================================"
echo "METABAT2 BINNING"
echo "============================================================"
echo "Group:       ${GROUP}"
echo "Job ID:      ${SLURM_JOB_ID}"
echo "Array ID:    ${SLURM_ARRAY_TASK_ID}"
echo "Nodo:        $(hostname)"
echo "Inicio:      $(date)"
echo "CPUs:        ${SLURM_CPUS_PER_TASK}"
echo
echo "Assembly:"
echo "  ${ASSEMBLY}"
echo
echo "Depth:"
echo "  ${DEPTH}"
echo
echo "Output:"
echo "  ${OUT}"
echo "============================================================"
echo


# ============================================================
# VALIDACIONES
# ============================================================

if [[ ! -s "${NEW_METABAT}" ]]; then
    echo "ERROR: no existe el contenedor MetaBAT2 nuevo." >&2
    exit 1
fi

if [[ ! -s "${ASSEMBLY}" ]]; then
    echo "ERROR: no existe assembly para ${GROUP}" >&2
    exit 1
fi

if [[ ! -s "${DEPTH}" ]]; then
    echo "ERROR: no existe depth corregido para ${GROUP}" >&2
    exit 1
fi


# ============================================================
# VALIDAR DEPTH ANTES DEL BINNING
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

echo "Negative depth values: ${NEG}"

if [[ "${NEG}" -ne 0 ]]; then
    echo "ERROR: depth inválido para ${GROUP}" >&2
    exit 1
fi

echo


# ============================================================
# MOSTRAR VERSION
# ============================================================

echo "===== METABAT2 VERSION ====="

apptainer exec \
    --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
    "${NEW_METABAT}" \
    metabat2 --version || true

echo


# ============================================================
# LIMPIAR SOLO RESULTADOS INCOMPLETOS PREVIOS
# ============================================================

DONE="${OUT}/done"

if [[ -f "${DONE}" ]]; then

    echo "Binning previamente completado para ${GROUP}."
    echo "Se omite."

else

    rm -f "${BIN_PREFIX}".*.fa
    rm -f "${SUMMARY}"

    echo "===== RUNNING METABAT2 ====="
    echo

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        apptainer exec \
        --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
        "${NEW_METABAT}" \
        metabat2 \
        -i "${ASSEMBLY}" \
        -a "${DEPTH}" \
        -o "${BIN_PREFIX}" \
        -m 1500 \
        -t "${SLURM_CPUS_PER_TASK}"

    echo

fi


# ============================================================
# CONTAR BINS
# ============================================================

N_BINS=$(find "${OUT}" \
    -maxdepth 1 \
    -name "${GROUP}_bin.*.fa" \
    -type f | wc -l)

echo "===== BIN COUNT ====="
echo "Bins generated: ${N_BINS}"
echo


# ============================================================
# RESUMEN DE TAMAÑO DE LOS BINS
# ============================================================

printf "coassembly\tbin\tn_contigs\ttotal_bp\tmax_contig_bp\n" \
    > "${SUMMARY}"

for BIN in "${OUT}"/${GROUP}_bin.*.fa
do

    [[ -e "${BIN}" ]] || continue

    BINNAME=$(basename "${BIN}")

    awk -v group="${GROUP}" \
        -v bin="${BINNAME}" '
        /^>/ {
            if (len > 0) {
                total += len
                n++
                if (len > max) max=len
            }
            len=0
            next
        }

        {
            len += length($0)
        }

        END {
            if (len > 0) {
                total += len
                n++
                if (len > max) max=len
            }

            printf "%s\t%s\t%d\t%d\t%d\n",
                   group, bin, n, total, max
        }
    ' "${BIN}" >> "${SUMMARY}"

done


# ============================================================
# VALIDACION FINAL
# ============================================================

if [[ "${N_BINS}" -eq 0 ]]; then

    echo "WARNING: MetaBAT2 no generó bins para ${GROUP}."

else

    touch "${DONE}"

fi


# ============================================================
# OUTPUT
# ============================================================

echo "===== FIRST BIN SUMMARY ROWS ====="

head -n 10 "${SUMMARY}"

echo

echo "===== OUTPUT SIZE ====="

du -sh "${OUT}"

echo

echo "============================================================"
echo "METABAT2 BINNING COMPLETADO"
echo "Group: ${GROUP}"
echo "Bins:  ${N_BINS}"
echo "Fin: $(date)"
echo "============================================================"
