#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00

#SBATCH --job-name=antismash_rep18
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

INPUTDIR="${ROOT}/22_final_representative_mags/fastas"
BASE="${ROOT}/26_antismash_rep18"
RESULTS="${BASE}/results"

DB="${ROOT}/databases/antismash_7.1.0"

PKGBASE="${ROOT}/databases/antismash_7.1.0_package"
PKG="${PKGBASE}/antismash"

ANTISMASH="${ANTISMASH:-/opt/ohpc/pub/containers/BIO/antismash-7.1.0.sif}"

INTERNAL_PKG="/usr/local/lib/python3.10/site-packages/antismash"

echo "============================================================"
echo "antiSMASH - 18 MAGs REPRESENTATIVOS"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "============================================================"


# ============================================================
# VALIDACIONES
# ============================================================

if [[ ! -d "${INPUTDIR}" ]]; then
    echo "ERROR: falta ${INPUTDIR}" >&2
    exit 1
fi

if [[ ! -e "${DB}/ANTISMASH_7.1.0_DB_READY" ]]; then
    echo "ERROR: base antiSMASH no preparada." >&2
    exit 1
fi

if [[ ! -e "${PKGBASE}/ANTISMASH_7.1.0_PACKAGE_PREPARED" ]]; then
    echo "ERROR: paquete antiSMASH no preparado." >&2
    exit 1
fi

BGC="${PKG}/detection/hmm_detection/data/bgc_seeds.hmm"

for EXT in "" ".h3f" ".h3i" ".h3m" ".h3p"; do
    if [[ ! -s "${BGC}${EXT}" ]]; then
        echo "ERROR: falta ${BGC}${EXT}" >&2
        exit 1
    fi
done


# ============================================================
# MAGs
# ============================================================

mapfile -t MAGS < <(
    find "${INPUTDIR}" \
    -maxdepth 1 \
    -type f \
    -name '*.fa' \
    -printf '%f\n' \
    | sed 's/\.fa$//' \
    | sort
)

NMAGS=${#MAGS[@]}

echo
echo "MAGs encontrados: ${NMAGS}"

if [[ "${NMAGS}" -ne 18 ]]; then
    echo "ERROR: se esperaban 18 MAGs y se encontraron ${NMAGS}" >&2
    exit 1
fi

printf '%s\n' "${MAGS[@]}"


# ============================================================
# PREFLIGHT
# ============================================================

echo
echo "============================================================"
echo "PREFLIGHT"
echo "============================================================"

apptainer exec \
    --cleanenv \
    --bind "${PKG}:${INTERNAL_PKG}:ro" \
    --bind "${DB}:/antismash_db:ro" \
    "${ANTISMASH}" \
    antismash \
        --databases /antismash_db \
        --check-prereqs


# ============================================================
# OUTPUT
# ============================================================

if [[ -e "${RESULTS}" ]]; then
    echo "ERROR: ya existe ${RESULTS}" >&2
    echo "No se sobrescribirá una corrida previa." >&2
    exit 1
fi

mkdir -p "${RESULTS}"


# ============================================================
# EJECUCIÓN MAG POR MAG
# ============================================================

COUNT=0

for MAG in "${MAGS[@]}"; do

    COUNT=$((COUNT + 1))

    INPUT="${INPUTDIR}/${MAG}.fa"
    OUT="${RESULTS}/${MAG}"

    echo
    echo "============================================================"
    echo "[${COUNT}/${NMAGS}] ${MAG}"
    echo "Inicio MAG: $(date)"
    echo "============================================================"

    if [[ ! -s "${INPUT}" ]]; then
        echo "ERROR: FASTA vacío o inexistente: ${INPUT}" >&2
        exit 1
    fi

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        apptainer exec \
        --cleanenv \
        --bind "${PKG}:${INTERNAL_PKG}:ro" \
        --bind "${DB}:/antismash_db:ro" \
        --bind "${INPUTDIR}:/input:ro" \
        --bind "${RESULTS}:/work" \
        "${ANTISMASH}" \
        antismash \
            --taxon bacteria \
            --genefinding-tool prodigal \
            --cpus "${SLURM_CPUS_PER_TASK}" \
            --databases /antismash_db \
            --cb-knownclusters \
            --cc-mibig \
            --rre \
            --output-dir "/work/${MAG}" \
            "/input/${MAG}.fa"

    JSON="${OUT}/${MAG}.json"

    if [[ ! -s "${JSON}" ]]; then
        echo "ERROR: ${MAG} no produjo JSON." >&2
        exit 1
    fi

    NREGIONS=$(find "${OUT}" \
        -maxdepth 1 \
        -type f \
        -name '*.region*.gbk' \
        | wc -l)

    echo "${MAG}: ${NREGIONS} regiones BGC"
    echo "Fin MAG: $(date)"

done


# ============================================================
# RESUMEN GLOBAL
# ============================================================

SUMMARY="${BASE}/antismash_rep18_run_summary.tsv"

printf "MAG\tBGC_regions\tJSON\tHTML\n" > "${SUMMARY}"

TOTAL_REGIONS=0

for MAG in "${MAGS[@]}"; do

    OUT="${RESULTS}/${MAG}"

    NREGIONS=$(find "${OUT}" \
        -maxdepth 1 \
        -type f \
        -name '*.region*.gbk' \
        | wc -l)

    TOTAL_REGIONS=$((TOTAL_REGIONS + NREGIONS))

    JSON_STATUS="NO"
    HTML_STATUS="NO"

    [[ -s "${OUT}/${MAG}.json" ]] && JSON_STATUS="YES"
    [[ -s "${OUT}/index.html" ]] && HTML_STATUS="YES"

    printf "%s\t%s\t%s\t%s\n" \
        "${MAG}" \
        "${NREGIONS}" \
        "${JSON_STATUS}" \
        "${HTML_STATUS}" \
        >> "${SUMMARY}"

done


echo
echo "===== RESUMEN ====="

cat "${SUMMARY}"

echo
echo "MAGs analizados: ${NMAGS}"
echo "Regiones BGC totales: ${TOTAL_REGIONS}"

echo
echo "============================================================"
echo "ANTISMASH REP18 COMPLETADO"
echo "Resultados: ${RESULTS}"
echo "Resumen: ${SUMMARY}"
echo "Fin: $(date)"
echo "============================================================"
