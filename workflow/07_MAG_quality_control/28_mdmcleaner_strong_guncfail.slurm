#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=24:00:00

#SBATCH --job-name=mdmcleaner_8MAGs
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

# Resolve repository-relative support files.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"
DASTOOL="${ROOT}/12_dastool"

OUT_PARENT="${ROOT}/15_mdmcleaner_rescue"
OUT_HOST="${OUT_PARENT}/strong_gunc_fail"
STAGE="${OUT_PARENT}/inputs_writable"

DB_HOST="${ROOT}/databases/mdmcleaner"

MDMCLEANER="${ROOT}/containers/mdmcleaner-0.8.7.sif"
SAMTOOLS="${SAMTOOLS:-/opt/ohpc/pub/containers/BIO/samtools-1.21.sif}"

CONFIG_HOST="${MDMCLEANER_CONFIG_DIR:-${REPO_ROOT}/config/mdmcleaner}"
CONFIG_FILE="${MDMCLEANER_CONFIG:-${CONFIG_HOST}/mdmcleaner.config}"

PATCH_PREP="${REPO_ROOT}/patches/mdmcleaner-0.8.7/prepare_mdmcleaner_patch.sh"
PATCH="${MDMCLEANER_PATCH:-${STAGE}/read_gtdb_taxonomy.py}"
PATCH_TARGET="/usr/local/lib/python3.11/site-packages/mdmcleaner/read_gtdb_taxonomy.py"

EXPECTED_PATCH_SHA256="e9b683ce506941e89638691ac57bd66da6620d7e2a0defb8fa528df54cae69f3"

DONE="${OUT_HOST}/MDMCLEANER_8MAGS_COMPLETED"

echo "============================================================"
echo "MDMCLEANER - 8 MAGs FUERTES GUNC-FAIL"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "============================================================"
echo

# ============================================================
# VALIDACIONES
# ============================================================

for F in \
    "${MDMCLEANER}" \
    "${SAMTOOLS}" \
    "${CONFIG_FILE}" \
    "${PATCH_PREP}"
do
    [[ -s "${F}" ]] || {
        echo "ERROR: falta ${F}" >&2
        exit 1
    }
done

if [[ -e "${DONE}" ]]; then
    echo "ERROR: existe marcador de ejecución completa:"
    echo "${DONE}"
    exit 1
fi

mkdir -p \
    "${OUT_HOST}" \
    "${STAGE}/L1" \
    "${STAGE}/L2" \
    "${STAGE}/M3"


# ============================================================
# PREPARE REPRODUCIBLE MDMCLEANER COMPATIBILITY PATCH
# ============================================================

echo
echo "Preparing MDMcleaner compatibility source..."

srun --exact \
    -n 1 \
    -c 1 \
    env \
    MDMCLEANER_PATCH="${PATCH}" \
    "${PATCH_PREP}" \
    "${MDMCLEANER}"

[[ -s "${PATCH}" ]] || {
    echo "ERROR: patched MDMcleaner source was not generated:" >&2
    echo "${PATCH}" >&2
    exit 1
}

ACTUAL_PATCH_SHA256="$(
    sha256sum "${PATCH}" |
    awk '{print $1}'
)"

echo
echo "MDMcleaner patched source SHA-256:"
echo "${ACTUAL_PATCH_SHA256}"

if [[ "${ACTUAL_PATCH_SHA256}" != "${EXPECTED_PATCH_SHA256}" ]]; then
    echo "ERROR: generated MDMcleaner patch does not match" >&2
    echo "the exact source used in the thesis." >&2
    echo >&2
    echo "Expected: ${EXPECTED_PATCH_SHA256}" >&2
    echo "Observed: ${ACTUAL_PATCH_SHA256}" >&2
    exit 1
fi

echo "MDMCLEANER_PATCH_EXACT_HISTORICAL_MATCH=YES"


# ============================================================
# FASTA ORIGINALES REFINADOS
# ============================================================

SRC_L1_A="${DASTOOL}/L1/refined/L1_DASTool_DASTool_bins/L1_maxbin2.051.fa"
SRC_L1_B="${DASTOOL}/L1/refined/L1_DASTool_DASTool_bins/concoct_20.fa"
SRC_L1_C="${DASTOOL}/L1/refined/L1_DASTool_DASTool_bins/concoct_24.fa"

SRC_L2_A="${DASTOOL}/L2/refined/L2_DASTool_DASTool_bins/L2_maxbin2.004_sub.fa"
SRC_L2_B="${DASTOOL}/L2/refined/L2_DASTool_DASTool_bins/concoct_55.fa"
SRC_L2_C="${DASTOOL}/L2/refined/L2_DASTool_DASTool_bins/concoct_72.fa"
SRC_L2_D="${DASTOOL}/L2/refined/L2_DASTool_DASTool_bins/concoct_91.fa"

SRC_M3_A="${DASTOOL}/M3/refined/M3_DASTool_DASTool_bins/concoct_37.fa"

SOURCES=(
    "${SRC_L1_A}"
    "${SRC_L1_B}"
    "${SRC_L1_C}"
    "${SRC_L2_A}"
    "${SRC_L2_B}"
    "${SRC_L2_C}"
    "${SRC_L2_D}"
    "${SRC_M3_A}"
)

for F in "${SOURCES[@]}"; do
    [[ -s "${F}" ]] || {
        echo "ERROR: falta MAG original: ${F}" >&2
        exit 1
    }
done


# ============================================================
# COPIAS ESCRIBIBLES
# ============================================================

cp -f "${SRC_L1_A}" "${STAGE}/L1/L1_maxbin2.051.fa"
cp -f "${SRC_L1_B}" "${STAGE}/L1/concoct_20.fa"
cp -f "${SRC_L1_C}" "${STAGE}/L1/concoct_24.fa"

cp -f "${SRC_L2_A}" "${STAGE}/L2/L2_maxbin2.004_sub.fa"
cp -f "${SRC_L2_B}" "${STAGE}/L2/concoct_55.fa"
cp -f "${SRC_L2_C}" "${STAGE}/L2/concoct_72.fa"
cp -f "${SRC_L2_D}" "${STAGE}/L2/concoct_91.fa"

cp -f "${SRC_M3_A}" "${STAGE}/M3/concoct_37.fa"

INPUTS=(
    "${STAGE}/L1/L1_maxbin2.051.fa"
    "${STAGE}/L1/concoct_20.fa"
    "${STAGE}/L1/concoct_24.fa"
    "${STAGE}/L2/L2_maxbin2.004_sub.fa"
    "${STAGE}/L2/concoct_55.fa"
    "${STAGE}/L2/concoct_72.fa"
    "${STAGE}/L2/concoct_91.fa"
    "${STAGE}/M3/concoct_37.fa"
)


# ============================================================
# TRAZABILIDAD
# ============================================================

echo "===== SHA256 ORIGINALES ====="

sha256sum "${SOURCES[@]}" \
    > "${OUT_PARENT}/strong_gunc_fail_originals.sha256"

cat "${OUT_PARENT}/strong_gunc_fail_originals.sha256"

echo
echo "===== SHA256 COPIAS ====="

sha256sum "${INPUTS[@]}" \
    > "${OUT_PARENT}/strong_gunc_fail_staged.sha256"

cat "${OUT_PARENT}/strong_gunc_fail_staged.sha256"


# ============================================================
# GENERAR .fai ANTES DE BARRNAP
# ============================================================

echo
echo "===== GENERANDO ÍNDICES FASTA ====="

for F in "${INPUTS[@]}"; do

    REL="${F#${STAGE}/}"

    srun --exact \
        -n 1 \
        -c 1 \
        apptainer exec \
        --cleanenv \
        --bind "${STAGE}:/inputs" \
        "${SAMTOOLS}" \
        samtools faidx "/inputs/${REL}"

done

echo
echo "===== ÍNDICES GENERADOS ====="

find "${STAGE}" \
    -type f \
    -name '*.fai' \
    -printf '%P\t%s bytes\n' \
    | sort

N_FAI=$(find "${STAGE}" -type f -name '*.fai' | wc -l)

if [[ "${N_FAI}" -ne 8 ]]; then
    echo "ERROR: se esperaban 8 índices .fai y hay ${N_FAI}" >&2
    exit 1
fi


# ============================================================
# PRUEBA DE DEPENDENCIAS CON CLEANENV
# ============================================================

echo
echo "===== BARRNAP DEPENDENCIAS ====="

apptainer exec \
    --cleanenv \
    "${MDMCLEANER}" \
    sh -c '
        echo "which:    $(command -v which)"
        echo "nhmmer:   $(which nhmmer)"
        echo "bedtools: $(which bedtools)"
        barrnap --version
    '

echo


# ============================================================
# EJECUTAR MDMCLEANER
# ============================================================

echo "============================================================"
echo "INICIANDO MDMcleaner clean"
echo "============================================================"
echo

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    apptainer exec \
    --cleanenv \
    --bind "${DB_HOST}:/mdmdb" \
    --bind "${STAGE}:/inputs" \
    --bind "${OUT_HOST}:/out" \
    --bind "${CONFIG_HOST}:/work" \
    --bind "${PATCH}:${PATCH_TARGET}" \
    --pwd /work \
    "${MDMCLEANER}" \
    mdmcleaner clean \
    -i \
        /inputs/L1/L1_maxbin2.051.fa \
        /inputs/L1/concoct_20.fa \
        /inputs/L1/concoct_24.fa \
        /inputs/L2/L2_maxbin2.004_sub.fa \
        /inputs/L2/concoct_55.fa \
        /inputs/L2/concoct_72.fa \
        /inputs/L2/concoct_91.fa \
        /inputs/M3/concoct_37.fa \
    -o /out \
    -c /work/mdmcleaner.config \
    -t "${SLURM_CPUS_PER_TASK}" \
    --overview_files_basename strong_gunc_fail


# ============================================================
# INVENTARIO DE RESULTADOS
# ============================================================

echo
echo "============================================================"
echo "RESULTADOS"
echo "============================================================"

find "${OUT_HOST}" \
    -maxdepth 3 \
    -type f \
    -printf '%P\t%s bytes\n' \
    | sort

echo

N_FULLINFO=$(find "${OUT_HOST}" \
    -type f \
    -name 'fullcontiginfos_beforecleanup.tsv' \
    | wc -l)

echo "Informes fullcontiginfos_beforecleanup.tsv: ${N_FULLINFO}"

if [[ "${N_FULLINFO}" -ne 8 ]]; then
    echo "ERROR: se esperaban informes para 8 MAGs." >&2
    exit 1
fi

touch "${DONE}"

echo
echo "============================================================"
echo "MDMCLEANER 8 MAGS GUNC-FAIL COMPLETADO"
echo "Fin: $(date)"
echo "============================================================"
