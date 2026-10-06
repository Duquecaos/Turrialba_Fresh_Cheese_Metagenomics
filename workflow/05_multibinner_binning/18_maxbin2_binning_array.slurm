#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00

#SBATCH --array=1-6%2
#SBATCH --job-name=maxbin2_bins

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

export LC_ALL=C
export LANG=C

COASSEMBLY_NAMES=(L1 L2 L3 M1 M2 M3)

ARRAY_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
GROUP="${COASSEMBLY_NAMES[${ARRAY_INDEX}]}"

ASSEMBLY="${SCRATCH_DIR}/05_assembly/${GROUP}/final.contigs.fa"
DEPTH="${SCRATCH_DIR}/07_depth_corrected/${GROUP}/${GROUP}_depth.txt"

OUT="${SCRATCH_DIR}/08_binning/${GROUP}/maxbin2"
DONE="${OUT}/done"

TMP_HOST="${SCRATCH_DIR}/tmp/maxbin2_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
TMP_CONTAINER="/mbtmp"

mkdir -p "${TMP_HOST}"

cleanup() {
    echo
    echo "===== CLEANING TEMPORARY FILES ====="
    rm -rf "${TMP_HOST}"
}

trap cleanup EXIT

echo "============================================================"
echo "MAXBIN2 BINNING"
echo "============================================================"
echo "Group:       ${GROUP}"
echo "Inicio:      $(date)"
echo "Nodo:        $(hostname)"
echo "Job ID:      ${SLURM_JOB_ID}"
echo "Array ID:    ${SLURM_ARRAY_TASK_ID}"
echo "CPUs:        ${SLURM_CPUS_PER_TASK}"
echo "Memoria:     64G"
echo
echo "Assembly:"
echo "  ${ASSEMBLY}"
echo
echo "Depth:"
echo "  ${DEPTH}"
echo
echo "Container:"
echo "  ${MAXBIN2}"
echo
echo "Output:"
echo "  ${OUT}"
echo "============================================================"
echo

if [[ ! -s "${MAXBIN2}" ]]; then
    echo "ERROR: contenedor MaxBin2 no encontrado." >&2
    exit 1
fi

if [[ ! -s "${ASSEMBLY}" ]]; then
    echo "ERROR: assembly no encontrado." >&2
    exit 1
fi

if [[ ! -s "${DEPTH}" ]]; then
    echo "ERROR: tabla de profundidad no encontrada." >&2
    exit 1
fi

if [[ -f "${DONE}" ]]; then
    echo "MaxBin2 ya fue completado previamente para ${GROUP}."
    echo "Se omite ejecución."
else

    rm -rf "${OUT}"
    mkdir -p "${OUT}/abundance"

    ABUND_LIST="${OUT}/abundance_files.txt"

    : > "${ABUND_LIST}"

    echo "===== GENERATING ABUNDANCE FILES ====="

    for SAMPLE in "${GROUP}_1" "${GROUP}_2" "${GROUP}_3"
    do

        ABUND="${OUT}/abundance/${SAMPLE}.abund.tsv"

        awk -F '\t' \
            -v OFS='\t' \
            -v target="${SAMPLE}.sorted" '

            NR == 1 {

                for (i=1; i<=NF; i++) {
                    if ($i == target)
                        col=i
                }

                if (!col) {
                    print "ERROR: column not found: " target > "/dev/stderr"
                    exit 2
                }

                next
            }

            {
                print $1, $col
            }

        ' "${DEPTH}" > "${ABUND}"

        if [[ ! -s "${ABUND}" ]]; then
            echo "ERROR: no se generó abundancia para ${SAMPLE}" >&2
            exit 1
        fi

        N_ASSEMBLY=$(grep -c '^>' "${ASSEMBLY}")
        N_ABUND=$(wc -l < "${ABUND}")

        if [[ "${N_ASSEMBLY}" -ne "${N_ABUND}" ]]; then
            echo "ERROR: contigs y abundancias no coinciden para ${SAMPLE}" >&2
            echo "Assembly=${N_ASSEMBLY} Abundance=${N_ABUND}" >&2
            exit 1
        fi

        echo "${ABUND}" >> "${ABUND_LIST}"

    done

    echo
    echo "Abundance files:"
    cat "${ABUND_LIST}"
    echo

    PREFIX="${OUT}/${GROUP}_maxbin2"

    echo "===== RUNNING MAXBIN2 ====="
    echo

    srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    apptainer exec \
    --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
    --bind "${TMP_HOST}:${TMP_CONTAINER}" \
    --pwd "${TMP_CONTAINER}" \
    --env TMPDIR="${TMP_CONTAINER}" \
    "${MAXBIN2}" \
    run_MaxBin.pl \
        -contig "${ASSEMBLY}" \
        -abund_list "${ABUND_LIST}" \
        -out "${PREFIX}" \
        -min_contig_length 1500 \
        -thread "${SLURM_CPUS_PER_TASK}"

    echo

    shopt -s nullglob
    BIN_FILES=( "${PREFIX}".*.fasta )
    shopt -u nullglob

    N_BINS=${#BIN_FILES[@]}

    if [[ "${N_BINS}" -eq 0 ]]; then
        echo "ERROR: MaxBin2 no generó bins para ${GROUP}." >&2
        exit 1
    fi

    SUMMARY="${OUT}/${GROUP}_maxbin2_bins.tsv"

    printf "coassembly\tbin\tn_contigs\ttotal_bp\tmax_contig_bp\n" > "${SUMMARY}"

    for BIN in "${BIN_FILES[@]}"
    do

        BINNAME=$(basename "${BIN}")

        awk \
            -v OFS='\t' \
            -v group="${GROUP}" \
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

                print group, bin, n, total, max
            }

        ' "${BIN}" >> "${SUMMARY}"

    done

    touch "${DONE}"

fi

shopt -s nullglob
FINAL_BINS=( "${OUT}/${GROUP}_maxbin2".*.fasta )
shopt -u nullglob

N_FINAL=${#FINAL_BINS[@]}

echo "===== MAXBIN2 RESULTS ====="
echo "Bins generated: ${N_FINAL}"
echo

if [[ -s "${OUT}/${GROUP}_maxbin2_bins.tsv" ]]; then
    head -n 10 "${OUT}/${GROUP}_maxbin2_bins.tsv" | column -t
fi

echo
du -sh "${OUT}"

echo
echo "============================================================"
echo "MAXBIN2 COMPLETADO"
echo "Group: ${GROUP}"
echo "Bins:  ${N_FINAL}"
echo "Fin: $(date)"
echo "============================================================"
