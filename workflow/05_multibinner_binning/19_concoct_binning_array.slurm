#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00

#SBATCH --array=1-6%2
#SBATCH --job-name=concoct_bins

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

MAPROOT="${SCRATCH_DIR}/06_mapping/${GROUP}"

BAM1="${MAPROOT}/${GROUP}_1/${GROUP}_1.sorted.bam"
BAM2="${MAPROOT}/${GROUP}_2/${GROUP}_2.sorted.bam"
BAM3="${MAPROOT}/${GROUP}_3/${GROUP}_3.sorted.bam"

OUT="${SCRATCH_DIR}/08_binning/${GROUP}/concoct"

DONE="${OUT}/done"

TMP_HOST="${SCRATCH_DIR}/tmp/concoct_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
TMP_CONTAINER="/concoct_tmp"

mkdir -p "${TMP_HOST}"

cleanup() {
    echo
    echo "===== CLEANING TEMPORARY FILES ====="
    rm -rf "${TMP_HOST}"
}

trap cleanup EXIT

echo "============================================================"
echo "CONCOCT BINNING"
echo "============================================================"
echo "Group:       ${GROUP}"
echo "Inicio:      $(date)"
echo "Nodo:        $(hostname)"
echo "Job ID:      ${SLURM_JOB_ID}"
echo "Array ID:    ${SLURM_ARRAY_TASK_ID}"
echo "CPUs:        ${SLURM_CPUS_PER_TASK}"
echo "Memoria:     96G"
echo
echo "Assembly:"
echo "  ${ASSEMBLY}"
echo
echo "BAMs:"
echo "  ${BAM1}"
echo "  ${BAM2}"
echo "  ${BAM3}"
echo
echo "Container:"
echo "  ${CONCOCT}"
echo
echo "Output:"
echo "  ${OUT}"
echo "============================================================"
echo

if [[ ! -s "${CONCOCT}" ]]; then
    echo "ERROR: contenedor CONCOCT no encontrado." >&2
    exit 1
fi

if [[ ! -s "${ASSEMBLY}" ]]; then
    echo "ERROR: assembly no encontrado." >&2
    exit 1
fi

for BAM in "${BAM1}" "${BAM2}" "${BAM3}"
do

    if [[ ! -s "${BAM}" ]]; then
        echo "ERROR: BAM no encontrado: ${BAM}" >&2
        exit 1
    fi

    if [[ ! -s "${BAM}.bai" ]]; then
        echo "ERROR: índice BAM no encontrado: ${BAM}.bai" >&2
        exit 1
    fi

done

if [[ -f "${DONE}" ]]; then

    echo "CONCOCT ya fue completado previamente para ${GROUP}."
    echo "Se omite ejecución."

else

    rm -rf "${OUT}"
    mkdir -p "${OUT}"

    BED="${OUT}/contigs_10K.bed"
    CUT_FASTA="${OUT}/contigs_10K.fa"
    COVERAGE="${OUT}/coverage_table.tsv"

    CLUSTER_DIR="${OUT}/clustering"
    BINS_DIR="${OUT}/fasta_bins"

    mkdir -p "${CLUSTER_DIR}"
    mkdir -p "${BINS_DIR}"

    echo "===== CUTTING CONTIGS INTO 10 KB PARTS ====="
    echo

    srun --exact \
        -n 1 \
        -c 1 \
        apptainer exec \
        --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
        "${CONCOCT}" \
        cut_up_fasta.py \
        "${ASSEMBLY}" \
        -c 10000 \
        -o 0 \
        --merge_last \
        -b "${BED}" \
        > "${CUT_FASTA}"

    if [[ ! -s "${BED}" || ! -s "${CUT_FASTA}" ]]; then
        echo "ERROR: falló cut_up_fasta.py" >&2
        exit 1
    fi

    echo "===== GENERATING COVERAGE TABLE ====="
    echo

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        apptainer exec \
        --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
        --bind "${TMP_HOST}:${TMP_CONTAINER}" \
        --env TMPDIR="${TMP_CONTAINER}" \
        "${CONCOCT}" \
        concoct_coverage_table.py \
        "${BED}" \
        "${BAM1}" \
        "${BAM2}" \
        "${BAM3}" \
        > "${COVERAGE}"

    if [[ ! -s "${COVERAGE}" ]]; then
        echo "ERROR: no se generó coverage_table.tsv" >&2
        exit 1
    fi

    echo "===== RUNNING CONCOCT ====="
    echo

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        apptainer exec \
        --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
        --bind "${TMP_HOST}:${TMP_CONTAINER}" \
        --env TMPDIR="${TMP_CONTAINER}" \
        "${CONCOCT}" \
        concoct \
        --composition_file "${CUT_FASTA}" \
        --coverage_file "${COVERAGE}" \
        -b "${CLUSTER_DIR}/" \
        -t "${SLURM_CPUS_PER_TASK}" \
        -l 1500 \
        -r 151 \
        -s 1

    echo

    shopt -s nullglob
    CLUSTER_FILES=( "${CLUSTER_DIR}"/clustering_gt*.csv )
    shopt -u nullglob

    if [[ "${#CLUSTER_FILES[@]}" -ne 1 ]]; then
        echo "ERROR: se esperaba un archivo clustering_gt*.csv." >&2
        printf '%s\n' "${CLUSTER_FILES[@]:-}"
        exit 1
    fi

    CLUSTERING="${CLUSTER_FILES[0]}"
    MERGED="${CLUSTER_DIR}/clustering_merged.csv"

    echo "Clustering file:"
    echo "${CLUSTERING}"
    echo

    echo "===== MERGING CUT-UP CLUSTERING ====="
    echo

    srun --exact \
        -n 1 \
        -c 1 \
        apptainer exec \
        --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
        "${CONCOCT}" \
        merge_cutup_clustering.py \
        "${CLUSTERING}" \
        > "${MERGED}"

    if [[ ! -s "${MERGED}" ]]; then
        echo "ERROR: no se generó clustering_merged.csv" >&2
        exit 1
    fi

    echo "===== EXTRACTING FASTA BINS ====="
    echo

    srun --exact \
        -n 1 \
        -c 1 \
        apptainer exec \
        --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
        "${CONCOCT}" \
        extract_fasta_bins.py \
        "${ASSEMBLY}" \
        "${MERGED}" \
        --output_path "${BINS_DIR}"

    shopt -s nullglob
    BIN_FILES=( "${BINS_DIR}"/*.fa )
    shopt -u nullglob

    N_BINS=${#BIN_FILES[@]}

    if [[ "${N_BINS}" -eq 0 ]]; then
        echo "ERROR: CONCOCT no generó bins." >&2
        exit 1
    fi

    SUMMARY="${OUT}/${GROUP}_concoct_bins.tsv"

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

BINS_DIR="${OUT}/fasta_bins"

shopt -s nullglob
FINAL_BINS=( "${BINS_DIR}"/*.fa )
shopt -u nullglob

N_FINAL=${#FINAL_BINS[@]}

echo "===== CONCOCT RESULTS ====="
echo "Bins generated: ${N_FINAL}"
echo

if [[ -s "${OUT}/${GROUP}_concoct_bins.tsv" ]]; then
    head -n 10 "${OUT}/${GROUP}_concoct_bins.tsv" | column -t
fi

echo
du -sh "${OUT}"

echo
echo "============================================================"
echo "CONCOCT COMPLETADO"
echo "Group: ${GROUP}"
echo "Bins:  ${N_FINAL}"
echo "Fin: $(date)"
echo "============================================================"
