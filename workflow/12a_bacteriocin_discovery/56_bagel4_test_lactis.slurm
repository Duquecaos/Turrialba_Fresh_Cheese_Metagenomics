#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00

#SBATCH --job-name=bagel4_test_lactis
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

ENV="${ROOT}/envs/bagel4-1.2"
BAGEL="${ROOT}/tools/BAGEL4_hpc"

INPUT="${ROOT}/34_bagel4_BAL5/input/M2__M2_maxbin2.004.fna"

BASE="${ROOT}/35_bagel4_test_lactis"
SESSION="${BASE}/session"

export PATH="${ENV}/bin:${PATH}"
export PERL5LIB="${BAGEL}/lib"

echo "============================================================"
echo "BAGEL4 TEST - Lactococcus lactis"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "============================================================"


# ============================================================
# PREFLIGHT
# ============================================================

test -s "${INPUT}" || {
    echo "ERROR: falta ${INPUT}" >&2
    exit 1
}

# ============================================================
# BASES REALMENTE NECESARIAS PARA DETECCIÓN BAGEL4
# ============================================================

for F in \
    "${BAGEL}/db_hmm/bacteriocin.hmm.descriptions" \
    "${BAGEL}/db_proteins/bagel4_class1bacteriocin_db.fasta" \
    "${BAGEL}/db_proteins/bagel4_class2bacteriocin_db.fasta" \
    "${BAGEL}/db_proteins/bagel4_class3bacteriocin_db.fasta" \
    "${BAGEL}/db_proteins/all_context.db.phr" \
    "${BAGEL}/db_proteins/all_context.db.pin" \
    "${BAGEL}/db_proteins/all_context.db.psq"
do
    test -s "${F}" || {
        echo "ERROR: falta archivo BAGEL4: ${F}" >&2
        exit 1
    }
done


# ============================================================
# HMMs USADOS POR LAS REGLAS AOI
# ============================================================

MISSING_HMM=0

while read -r HMM
do
    [[ -n "${HMM}" ]] || continue

    if [[ ! -s "${BAGEL}/db_hmm/${HMM}.hmm" ]]; then
        echo "ERROR: falta HMM requerido: ${HMM}.hmm" >&2
        MISSING_HMM=1
    fi

done < <(
    awk -F $'\t' '
    {
        for (i=2; i<=NF; i++) {
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", $i)
            if ($i != "") print $i
        }
    }' "${BAGEL}/tables/AOI_identification_rules_v2.txt" \
    | sort -u
)

if [[ "${MISSING_HMM}" -ne 0 ]]; then
    exit 1
fi

echo "Bases BAGEL4: OK"
echo "HMMs de reglas AOI: OK"

for X in \
    perl \
    blastp \
    formatdb \
    hmmsearch \
    glimmer3 \
    long-orfs \
    extract \
    build-icm
do
    command -v "${X}" >/dev/null || {
        echo "ERROR: falta ${X}" >&2
        exit 1
    }
done

echo "Perl:"
perl -v | sed -n '1,2p'

echo
echo "BLASTP:"
blastp -version | sed -n '1,2p'

echo
echo "HMMER:"
hmmsearch -h | sed -n '1,2p'

echo
echo "Glimmer:"
command -v glimmer3


# ============================================================
# SESIÓN LIMPIA
# ============================================================

if [[ -e "${BASE}" ]]; then
    echo "ERROR: ya existe ${BASE}" >&2
    echo "No se sobrescribirá la prueba." >&2
    exit 1
fi

mkdir -p \
    "${SESSION}/queryfolder"

cp -p \
    "${INPUT}" \
    "${SESSION}/queryfolder/"


# ============================================================
# BAGEL4
# ============================================================

echo
echo "============================================================"
echo "EJECUTANDO BAGEL4"
echo "============================================================"

cd "${BAGEL}"

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    perl "${BAGEL}/bagel4_wrapper.pl" \
        -s "${SESSION}" \
        -query "${SESSION}/queryfolder" \
        -r '*.fna'


# ============================================================
# VALIDACIÓN
# ============================================================

echo
echo "============================================================"
echo "VALIDANDO SALIDAS"
echo "============================================================"

N_INPUT=$(
    awk 'END{print NR-1}' \
    "${SESSION}/00.all_contigs.table"
)

N_BLAST=$(
    find "${SESSION}" \
        -maxdepth 1 \
        -type f \
        -name '*.AOI.blast.table' \
        -size +0c \
        | wc -l
)

N_HMM=$(
    find "${SESSION}" \
        -maxdepth 1 \
        -type f \
        -name '*.AOI.hmm.table' \
        -size +0c \
        | wc -l
)

N_AOI=$(
    find "${SESSION}" \
        -maxdepth 1 \
        -type f \
        -name '*.AOI.table' \
        -size +0c \
        | wc -l
)

N_GENETABLES=$(
    find "${SESSION}" \
        -maxdepth 1 \
        -type f \
        -name '*.GeneTables' \
        -size +0c \
        | wc -l
)

echo "Contigs BAGEL4 procesados: ${N_INPUT}"
echo "AOI BLAST no vacíos: ${N_BLAST}"
echo "AOI HMM no vacíos: ${N_HMM}"
echo "AOI finales no vacíos: ${N_AOI}"
echo "GeneTables no vacíos: ${N_GENETABLES}"

if [[ "${N_AOI}" -lt 1 ]]; then
    echo "ERROR: BAGEL4 no produjo ningún AOI no vacío." >&2
    exit 1
fi

if [[ "${N_GENETABLES}" -lt 1 ]]; then
    echo "ERROR: BAGEL4 produjo AOIs pero no GeneTables." >&2
    exit 1
fi


# ============================================================
# RESUMEN DE AOIs
# ============================================================

echo
echo "============================================================"
echo "AOIs DETECTADOS"
echo "============================================================"

for F in "${SESSION}"/*.AOI.table; do

    [[ -s "${F}" ]] || continue

    echo
    echo "----- $(basename "${F}") -----"
    cat "${F}"

done


echo
echo "============================================================"
echo "HITS BACTERIOCINAS EN GENETABLES"
echo "============================================================"

for F in "${SESSION}"/*.GeneTables; do

    [[ -s "${F}" ]] || continue

    echo
    echo "----- $(basename "${F}") -----"

    grep -Ei \
        'bacterioc|lactococc|lant|pediocin|enterocin|nisin|garvicin|class.?II|leader|immunity|transport' \
        "${F}" \
        || true

done


# ============================================================
# INVENTARIO DE ARCHIVOS CLAVE
# ============================================================

echo
echo "============================================================"
echo "ARCHIVOS CLAVE"
echo "============================================================"

find "${SESSION}" \
    -maxdepth 1 \
    -type f \
    \( \
        -name '*.AOI.table' \
        -o -name '*.AOI.hmm.table' \
        -o -name '*.AOI.blast.table' \
        -o -name '*.GeneTables' \
        -o -name '*.GeneTable' \
        -o -name '00.*' \
    \) \
    -printf '%f\t%s bytes\n' \
    | sort


echo
echo "============================================================"
echo "BAGEL4 TEST COMPLETADO"
echo "MAG: M2__M2_maxbin2.004"
echo "Resultados: ${SESSION}"
echo "Fin: $(date)"
echo "============================================================"
