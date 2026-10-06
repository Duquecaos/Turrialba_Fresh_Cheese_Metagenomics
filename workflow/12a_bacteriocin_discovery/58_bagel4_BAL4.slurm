#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00

#SBATCH --job-name=bagel4_BAL4
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail


# ============================================================
# CONFIGURACIÓN
# ============================================================

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

ENV="${ROOT}/envs/bagel4-1.2"
BAGEL="${ROOT}/tools/BAGEL4_hpc"

INPUT_DIR="${ROOT}/34_bagel4_BAL5/input"
BASE="${ROOT}/36_bagel4_BAL4"

PATCH="${ROOT}/tools/BAGEL4_HPC_PATCH.diff"
PROVENANCE="${ROOT}/tools/BAGEL4_HPC_provenance.txt"

export PATH="${ENV}/bin:${PATH}"
export PERL5LIB="${BAGEL}/lib"

# Evita warnings de locale del cluster
export LANG=C
export LC_ALL=C


MAG_NAMES=(
    "L1__L1_bin.86"
    "L2__L2_maxbin2.004_sub"
    "L3__L3_bin.7"
    "L3__concoct_29"
)

TAXA=(
    "Leuconostoc carnosum"
    "Lactococcus petauri"
    "Leuconostoc falkenbergense"
    "Lactococcus_A laudensis"
)


echo "============================================================"
echo "BAGEL4 - CUATRO MAGs BAL RESTANTES"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "============================================================"


# ============================================================
# PREFLIGHT GENERAL
# ============================================================

echo
echo "============================================================"
echo "PREFLIGHT"
echo "============================================================"

test -d "${BAGEL}" || {
    echo "ERROR: falta ${BAGEL}" >&2
    exit 1
}

test -d "${ENV}" || {
    echo "ERROR: falta ${ENV}" >&2
    exit 1
}

test -s "${BAGEL}/tables/codon_table.txt" || {
    echo "ERROR: falta codon_table.txt" >&2
    exit 1
}

test -s "${BAGEL}/tables/AOI_identification_rules_v2.txt" || {
    echo "ERROR: faltan reglas AOI" >&2
    exit 1
}

test -s "${PATCH}" || {
    echo "ERROR: falta BAGEL4_HPC_PATCH.diff" >&2
    exit 1
}

test -s "${PROVENANCE}" || {
    echo "ERROR: falta BAGEL4_HPC_provenance.txt" >&2
    exit 1
}


# ============================================================
# COMPROBAR QUE EL PORT HPC CORRECTO ESTÁ ACTIVO
# ============================================================

if ! grep -Eq \
    'my \$debug[[:space:]]*=[[:space:]]*0' \
    "${BAGEL}/bagel4_wrapper.pl"
then
    echo "ERROR: BAGEL4 wrapper no tiene debug=0" >&2
    exit 1
fi

if grep -q \
    '/data/molgentools/genomics/codon_table.txt' \
    "${BAGEL}/lib/anne_genomics.pm"
then
    echo "ERROR: anne_genomics conserva codon_table antiguo" >&2
    exit 1
fi

if ! grep -q \
    'HPC_SKIP_WEBBLAST' \
    "${BAGEL}/bagel4_AOI_annotation.pl"
then
    echo "ERROR: WebBlast no está desactivado" >&2
    exit 1
fi

if ! grep -q \
    'HPC_SKIP_SWISSPROT' \
    "${BAGEL}/bagel4_AOI_annotation.pl"
then
    echo "ERROR: SwissProt externo no está desactivado" >&2
    exit 1
fi


# ============================================================
# EJECUTABLES
# ============================================================

for X in \
    perl \
    blastp \
    makeblastdb \
    blastall \
    formatdb \
    hmmsearch \
    glimmer3 \
    long-orfs \
    extract \
    build-icm
do
    command -v "${X}" >/dev/null 2>&1 || {
        echo "ERROR: ejecutable no encontrado: ${X}" >&2
        exit 1
    }
done

echo "Ejecutables: OK"


# ============================================================
# BASES BAGEL4
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

echo "Bases BAGEL4: OK"


# ============================================================
# HMMs REQUERIDOS POR LAS REGLAS AOI
# ============================================================

MISSING_HMM=0

while read -r HMM
do
    [[ -n "${HMM}" ]] || continue

    if [[ ! -s "${BAGEL}/db_hmm/${HMM}.hmm" ]]; then
        echo "ERROR: falta HMM: ${HMM}.hmm" >&2
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

echo "HMMs de reglas AOI: OK"


# ============================================================
# INPUTS
# ============================================================

for MAG in "${MAG_NAMES[@]}"
do
    INPUT="${INPUT_DIR}/${MAG}.fna"

    test -s "${INPUT}" || {
        echo "ERROR: falta input ${INPUT}" >&2
        exit 1
    }
done

echo "Inputs BAL4: OK"


# ============================================================
# PREPARAR DIRECTORIO GLOBAL
# ============================================================

mkdir -p "${BASE}"

cp -p \
    "${PROVENANCE}" \
    "${BASE}/BAGEL4_HPC_provenance.txt"

sha256sum "${PATCH}" \
    > "${BASE}/BAGEL4_HPC_PATCH.sha256"


# ============================================================
# EJECUCIÓN MAG POR MAG
# ============================================================

for IDX in "${!MAG_NAMES[@]}"
do
    MAG="${MAG_NAMES[$IDX]}"
    TAXONOMY="${TAXA[$IDX]}"

    INPUT="${INPUT_DIR}/${MAG}.fna"

    FINAL="${BASE}/${MAG}"
    WORK="${BASE}/.working_${MAG}_${SLURM_JOB_ID}"
    SESSION="${WORK}/session"

    echo
    echo
    echo "############################################################"
    echo "MAG: ${MAG}"
    echo "Taxonomía: ${TAXONOMY}"
    echo "############################################################"


    # --------------------------------------------------------
    # Saltar resultados ya completados
    # --------------------------------------------------------

    if [[ -f "${FINAL}/BAGEL4_COMPLETED" ]]; then
        echo "Resultado completo previo encontrado."
        echo "Se omite: ${MAG}"
        continue
    fi

    if [[ -e "${FINAL}" ]]; then
        echo "ERROR: existe directorio final sin marcador:" >&2
        echo "${FINAL}" >&2
        exit 1
    fi

    rm -rf "${WORK}"

    mkdir -p \
        "${SESSION}/queryfolder"

    cp -p \
        "${INPUT}" \
        "${SESSION}/queryfolder/"


    # --------------------------------------------------------
    # Ejecutar BAGEL4 validado
    # --------------------------------------------------------

    echo
    echo "============================================================"
    echo "EJECUTANDO BAGEL4: ${MAG}"
    echo "============================================================"

    cd "${BAGEL}"

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        perl "${BAGEL}/bagel4_wrapper.pl" \
            -s "${SESSION}" \
            -query "${SESSION}/queryfolder" \
            -r '*.fna'


    # --------------------------------------------------------
    # Validar salida básica
    # --------------------------------------------------------

    echo
    echo "============================================================"
    echo "VALIDANDO: ${MAG}"
    echo "============================================================"

    test -s "${SESSION}/00.all_contigs.table" || {
        echo "ERROR: falta 00.all_contigs.table para ${MAG}" >&2
        exit 1
    }

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

    N_FINAL=$(
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


    echo "Contigs procesados: ${N_INPUT}"
    echo "Contigs positivos por BLAST: ${N_BLAST}"
    echo "Contigs positivos por HMM: ${N_HMM}"
    echo "Contigs con AOI final: ${N_FINAL}"
    echo "GeneTables no vacíos: ${N_GENETABLES}"


    # --------------------------------------------------------
    # Un resultado negativo es válido
    # --------------------------------------------------------

    if [[ "${N_FINAL}" -eq 0 ]]; then

        if [[ "${N_GENETABLES}" -ne 0 ]]; then
            echo "ERROR: GeneTables sin AOI final en ${MAG}" >&2
            exit 1
        fi

        echo
        echo "RESULTADO NEGATIVO VÁLIDO:"
        echo "No se detectaron AOI BAGEL4 bajo estos parámetros."

    else

        if [[ "${N_GENETABLES}" -ne "${N_FINAL}" ]]; then
            echo "ERROR: ${MAG} tiene ${N_FINAL} contigs con AOI" >&2
            echo "pero ${N_GENETABLES} GeneTables." >&2
            exit 1
        fi

        echo "Anotación de AOI: OK"
    fi


    # --------------------------------------------------------
    # Evidencia AOI
    # --------------------------------------------------------

    AOI_EVIDENCE="${WORK}/aoi_evidence.tsv"

    printf \
        'MAG\ttaxonomy\tsource\tresult\n' \
        > "${AOI_EVIDENCE}"

    for TYPE in \
        AOI.blast.table \
        AOI.hmm.table \
        AOI.table
    do
        for F in "${SESSION}"/*."${TYPE}"
        do
            [[ -s "${F}" ]] || continue

            SOURCE=$(basename "${F}")

            while IFS= read -r LINE
            do
                [[ -n "${LINE}" ]] || continue

                printf \
                    '%s\t%s\t%s\t%s\n' \
                    "${MAG}" \
                    "${TAXONOMY}" \
                    "${SOURCE}" \
                    "${LINE}"

            done < "${F}"
        done
    done


    # --------------------------------------------------------
    # Candidatos bacteriocínicos en GeneTables
    # --------------------------------------------------------

    CANDIDATES="${WORK}/bacteriocin_candidates.tsv"

    printf \
        'MAG\ttaxonomy\tsource\tGeneTable_line\n' \
        > "${CANDIDATES}"

    N_CANDIDATES=0

    for F in "${SESSION}"/*.GeneTables
    do
        [[ -s "${F}" ]] || continue

        SOURCE=$(basename "${F}")

        if grep -Eiq \
            'bacterioc|lactococc|lactocin|lacticin|leaderLanM|ComC|PF04369|PF10439' \
            "${F}"
        then
            while IFS= read -r LINE
            do
                [[ -n "${LINE}" ]] || continue

                printf \
                    '%s\t%s\t%s\t%s\n' \
                    "${MAG}" \
                    "${TAXONOMY}" \
                    "${SOURCE}" \
                    "${LINE}" \
                    >> "${CANDIDATES}"

                N_CANDIDATES=$((N_CANDIDATES + 1))

            done < <(
                grep -Ei \
                    'bacterioc|lactococc|lactocin|lacticin|leaderLanM|ComC|PF04369|PF10439' \
                    "${F}" \
                || true
            )
        fi
    done


    # --------------------------------------------------------
    # Resumen por MAG
    # --------------------------------------------------------

    SUMMARY="${WORK}/summary.tsv"

    printf \
        'MAG\ttaxonomy\tinput_contigs\tblast_positive_contigs\thmm_positive_contigs\tfinal_AOI_contigs\tGeneTables\tbacteriocin_candidate_lines\n' \
        > "${SUMMARY}"

    printf \
        '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${MAG}" \
        "${TAXONOMY}" \
        "${N_INPUT}" \
        "${N_BLAST}" \
        "${N_HMM}" \
        "${N_FINAL}" \
        "${N_GENETABLES}" \
        "${N_CANDIDATES}" \
        >> "${SUMMARY}"


    echo
    echo "Candidatos bacteriocínicos anotados: ${N_CANDIDATES}"

    if [[ "${N_CANDIDATES}" -gt 0 ]]; then

        echo
        echo "===== CANDIDATOS ${MAG} ====="

        tail -n +2 "${CANDIDATES}"
    fi


    # --------------------------------------------------------
    # Cerrar resultado de forma atómica
    # --------------------------------------------------------

    {
        echo "MAG=${MAG}"
        echo "taxonomy=${TAXONOMY}"
        echo "job_id=${SLURM_JOB_ID}"
        echo "completed=$(date)"
        echo "input=${INPUT}"
        echo "bagel4=${BAGEL}"
        echo "environment=${ENV}"
    } > "${WORK}/BAGEL4_COMPLETED"

    mv \
        "${WORK}" \
        "${FINAL}"

    echo
    echo "MAG COMPLETADO: ${MAG}"
    echo "Resultados: ${FINAL}"
done


# ============================================================
# CONSOLIDAR LOS CUATRO RESULTADOS
# ============================================================

echo
echo
echo "============================================================"
echo "CONSOLIDANDO BAL4"
echo "============================================================"

MASTER_SUMMARY="${BASE}/BAL4_summary.tsv"
MASTER_CANDIDATES="${BASE}/BAL4_bacteriocin_candidates.tsv"
MASTER_AOI="${BASE}/BAL4_AOI_evidence.tsv"


printf \
'MAG\ttaxonomy\tinput_contigs\tblast_positive_contigs\thmm_positive_contigs\tfinal_AOI_contigs\tGeneTables\tbacteriocin_candidate_lines\n' \
> "${MASTER_SUMMARY}"

printf \
'MAG\ttaxonomy\tsource\tGeneTable_line\n' \
> "${MASTER_CANDIDATES}"

printf \
'MAG\ttaxonomy\tsource\tresult\n' \
> "${MASTER_AOI}"


for MAG in "${MAG_NAMES[@]}"
do
    FINAL="${BASE}/${MAG}"

    test -f "${FINAL}/BAGEL4_COMPLETED" || {
        echo "ERROR: MAG sin completar: ${MAG}" >&2
        exit 1
    }

    tail -n +2 \
        "${FINAL}/summary.tsv" \
        >> "${MASTER_SUMMARY}"

    tail -n +2 \
        "${FINAL}/bacteriocin_candidates.tsv" \
        >> "${MASTER_CANDIDATES}"

    tail -n +2 \
        "${FINAL}/aoi_evidence.tsv" \
        >> "${MASTER_AOI}"
done


# ============================================================
# RESUMEN FINAL
# ============================================================

echo
echo "============================================================"
echo "RESUMEN FINAL BAL4"
echo "============================================================"

column -t -s $'\t' \
    "${MASTER_SUMMARY}"

echo
echo "============================================================"
echo "CANDIDATOS BACTERIOCÍNICOS BAL4"
echo "============================================================"

if [[ $(wc -l < "${MASTER_CANDIDATES}") -gt 1 ]]; then
    column -t -s $'\t' \
        "${MASTER_CANDIDATES}"
else
    echo "No se detectaron candidatos bacteriocínicos anotados."
fi


echo
echo "============================================================"
echo "BAGEL4 BAL4 COMPLETADO"
echo "Resultados: ${BASE}"
echo "Fin: $(date)"
echo "============================================================"
