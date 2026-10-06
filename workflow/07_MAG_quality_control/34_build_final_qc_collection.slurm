#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00

#SBATCH --job-name=build_final_mags
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

PASS22="${ROOT}/14_gunc_refined/refined_candidates_90_5_GUNCpass.tsv"

DASTOOL="${ROOT}/12_dastool"

RESCUED="${ROOT}/16_mdmcleaner_kept/fastas"

FINAL="${ROOT}/19_final_qc_mags"
FASTAS="${FINAL}/fastas"

MANIFEST="${FINAL}/final_qc_manifest.tsv"
CHECKSUMS="${FINAL}/final_qc_sha256.tsv"

mkdir -p "${FINAL}"

if [[ -e "${FASTAS}" ]]; then
    echo "ERROR: ya existe ${FASTAS}" >&2
    echo "No se sobrescribirá automáticamente." >&2
    exit 1
fi

mkdir -p "${FASTAS}"

echo "============================================================"
echo "CONSTRUCCIÓN COLECCIÓN FINAL DE MAGs QC-PASS"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "============================================================"
echo


# ============================================================
# VALIDAR PASS22
# ============================================================

if [[ ! -s "${PASS22}" ]]; then
    echo "ERROR: falta ${PASS22}" >&2
    exit 1
fi

N_PASS22=$(awk 'NR>1 && NF>0 {n++} END {print n+0}' "${PASS22}")

echo "MAGs originales GUNC PASS fuertes: ${N_PASS22}"

if [[ "${N_PASS22}" -ne 22 ]]; then
    echo "ERROR: se esperaban 22 MAGs originales GUNC PASS." >&2
    exit 1
fi


# ============================================================
# MANIFIESTO
# ============================================================

printf \
"final_id\tcoassembly\tMAG\tsource_type\tsource_fasta\tfinal_fasta\n" \
> "${MANIFEST}"


# ============================================================
# COPIAR LOS 22 QC-PASS ORIGINALES
# ============================================================

echo
echo "===== COPIANDO 22 MAGs ORIGINALES QC-PASS ====="

tail -n +2 "${PASS22}" |
while IFS=$'\t' read -r \
    GROUP BIN COMPLETENESS CONTAMINATION \
    GUNC_LEVEL CSS GUNC_CONTAM RRS PASS
do

    [[ -n "${GROUP}" ]] || continue
    [[ -n "${BIN}" ]] || continue

    SRC="${DASTOOL}/${GROUP}/refined/${GROUP}_DASTool_DASTool_bins/${BIN}.fa"

    if [[ ! -s "${SRC}" ]]; then
        echo "ERROR: no existe FASTA refinado:" >&2
        echo "${SRC}" >&2
        exit 1
    fi

    FINAL_ID="${GROUP}__${BIN}"
    DEST="${FASTAS}/${FINAL_ID}.fa"

    if [[ -e "${DEST}" ]]; then
        echo "ERROR: ID final duplicado: ${FINAL_ID}" >&2
        exit 1
    fi

    cp -p "${SRC}" "${DEST}"

    printf \
    "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${FINAL_ID}" \
    "${GROUP}" \
    "${BIN}" \
    "original_refined_GUNCpass" \
    "${SRC}" \
    "${DEST}" \
    >> "${MANIFEST}"

    echo "${FINAL_ID}"

done


# ============================================================
# AÑADIR LOS 2 RESCATADOS POR MDMCLEANER
# ============================================================

echo
echo "===== AÑADIENDO 2 MAGs RESCATADOS ====="

declare -a RESCUED_IDS=(
    "L2__L2_maxbin2.004_sub"
    "L2__concoct_72"
)

declare -a RESCUED_GROUPS=(
    "L2"
    "L2"
)

declare -a RESCUED_BINS=(
    "L2_maxbin2.004_sub"
    "concoct_72"
)

for i in "${!RESCUED_IDS[@]}"; do

    FINAL_ID="${RESCUED_IDS[$i]}"
    GROUP="${RESCUED_GROUPS[$i]}"
    BIN="${RESCUED_BINS[$i]}"

    SRC="${RESCUED}/${FINAL_ID}.fa"
    DEST="${FASTAS}/${FINAL_ID}.fa"

    if [[ ! -s "${SRC}" ]]; then
        echo "ERROR: falta MAG rescatado:" >&2
        echo "${SRC}" >&2
        exit 1
    fi

    if [[ -e "${DEST}" ]]; then
        echo "ERROR: el MAG rescatado ya existe en la colección:" >&2
        echo "${FINAL_ID}" >&2
        exit 1
    fi

    cp -p "${SRC}" "${DEST}"

    printf \
    "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${FINAL_ID}" \
    "${GROUP}" \
    "${BIN}" \
    "MDMcleaner_KEEP_rescued" \
    "${SRC}" \
    "${DEST}" \
    >> "${MANIFEST}"

    echo "${FINAL_ID}"

done


# ============================================================
# VALIDACIONES
# ============================================================

echo
echo "===== VALIDACIÓN FINAL ====="

N_FASTA=$(find "${FASTAS}" \
    -maxdepth 1 \
    -type f \
    -name '*.fa' \
    | wc -l)

N_MANIFEST=$(awk 'NR>1 && NF>0 {n++} END {print n+0}' "${MANIFEST}")

N_UNIQUE=$(awk 'NR>1 {print $1}' "${MANIFEST}" \
    | sort -u \
    | wc -l)

echo "FASTA finales:       ${N_FASTA}"
echo "Filas manifiesto:    ${N_MANIFEST}"
echo "IDs únicos:          ${N_UNIQUE}"

if [[ "${N_FASTA}" -ne 24 ]]; then
    echo "ERROR: se esperaban 24 FASTA." >&2
    exit 1
fi

if [[ "${N_MANIFEST}" -ne 24 ]]; then
    echo "ERROR: se esperaban 24 filas en manifiesto." >&2
    exit 1
fi

if [[ "${N_UNIQUE}" -ne 24 ]]; then
    echo "ERROR: IDs finales no son únicos." >&2
    exit 1
fi


# ============================================================
# VERIFICAR FASTA NO VACÍOS
# ============================================================

for F in "${FASTAS}"/*.fa; do

    NSEQ=$(grep -c '^>' "${F}" || true)

    if [[ "${NSEQ}" -lt 1 ]]; then
        echo "ERROR: FASTA sin secuencias: ${F}" >&2
        exit 1
    fi

done


# ============================================================
# CHECKSUMS
# ============================================================

printf "sha256\tfinal_fasta\n" > "${CHECKSUMS}"

for F in "${FASTAS}"/*.fa; do
    SUM=$(sha256sum "${F}" | awk '{print $1}')
    printf "%s\t%s\n" "${SUM}" "${F}" >> "${CHECKSUMS}"
done


# ============================================================
# INVENTARIO
# ============================================================

echo
echo "===== COLECCIÓN FINAL ====="

find "${FASTAS}" \
    -maxdepth 1 \
    -type f \
    -name '*.fa' \
    -printf '%f\t%s bytes\n' \
    | sort

echo
echo "===== TIPO DE ORIGEN ====="

awk -F $'\t' '
NR>1 {
    n[$4]++
}
END {
    for (x in n)
        print x "\t" n[x]
}' "${MANIFEST}" | sort


echo
echo "============================================================"
echo "COLECCION FINAL QC COMPLETADA"
echo "MAGs finales: 24"
echo "Originales GUNC PASS: 22"
echo "Rescatados MDMcleaner: 2"
echo "FASTAs: ${FASTAS}"
echo "Manifiesto: ${MANIFEST}"
echo "Fin: $(date)"
echo "============================================================"
