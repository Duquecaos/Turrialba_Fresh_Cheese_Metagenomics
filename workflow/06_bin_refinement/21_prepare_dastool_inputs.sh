#!/bin/bash

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

COASSEMBLY_NAMES=(L1 L2 L3 M1 M2 M3)

ROOT="${SCRATCH_DIR}/12_dastool"
BINROOT="${SCRATCH_DIR}/08_binning"

mkdir -p "${ROOT}"

for G in "${COASSEMBLY_NAMES[@]}"
do

    echo "============================================================"
    echo "PREPARANDO DAS TOOL: ${G}"
    echo "============================================================"

    OUT="${ROOT}/${G}/inputs"
    mkdir -p "${OUT}"

    ASSEMBLY="${SCRATCH_DIR}/05_assembly/${G}/final.contigs.fa"

    META_DIR="${BINROOT}/${G}/metabat2"
    MAX_DIR="${BINROOT}/${G}/maxbin2"
    CON_DIR="${BINROOT}/${G}/concoct/fasta_bins"

    META_TSV="${OUT}/${G}_metabat2_contigs2bin.tsv"
    MAX_TSV="${OUT}/${G}_maxbin2_contigs2bin.tsv"
    CON_TSV="${OUT}/${G}_concoct_contigs2bin.tsv"

    # ---------------------------------------------------------
    # Verificar assembly
    # ---------------------------------------------------------

    if [[ ! -s "${ASSEMBLY}" ]]; then
        echo "ERROR: assembly no encontrado: ${ASSEMBLY}" >&2
        exit 1
    fi

    # ---------------------------------------------------------
    # MetaBAT2
    # ---------------------------------------------------------

    : > "${META_TSV}"

    shopt -s nullglob
    META_BINS=( "${META_DIR}/${G}_bin."*.fa )
    shopt -u nullglob

    if [[ "${#META_BINS[@]}" -eq 0 ]]; then
        echo "ERROR: no hay bins MetaBAT2 para ${G}" >&2
        exit 1
    fi

    for F in "${META_BINS[@]}"
    do
        BIN=$(basename "${F}" .fa)

        awk -v OFS='\t' -v bin="${BIN}" '
            /^>/ {
                id=$1
                sub(/^>/,"",id)
                print id,bin
            }
        ' "${F}" >> "${META_TSV}"
    done

    # ---------------------------------------------------------
    # MaxBin2
    # ---------------------------------------------------------

    : > "${MAX_TSV}"

    shopt -s nullglob
    MAX_BINS=( "${MAX_DIR}/${G}_maxbin2."*.fasta )
    shopt -u nullglob

    if [[ "${#MAX_BINS[@]}" -eq 0 ]]; then
        echo "ERROR: no hay bins MaxBin2 para ${G}" >&2
        exit 1
    fi

    for F in "${MAX_BINS[@]}"
    do
        BIN=$(basename "${F}" .fasta)

        awk -v OFS='\t' -v bin="${BIN}" '
            /^>/ {
                id=$1
                sub(/^>/,"",id)
                print id,bin
            }
        ' "${F}" >> "${MAX_TSV}"
    done

    # ---------------------------------------------------------
    # CONCOCT
    # ---------------------------------------------------------

    : > "${CON_TSV}"

    shopt -s nullglob
    CON_BINS=( "${CON_DIR}"/*.fa )
    shopt -u nullglob

    if [[ "${#CON_BINS[@]}" -eq 0 ]]; then
        echo "ERROR: no hay bins CONCOCT para ${G}" >&2
        exit 1
    fi

    for F in "${CON_BINS[@]}"
    do
        BIN="concoct_$(basename "${F}" .fa)"

        awk -v OFS='\t' -v bin="${BIN}" '
            /^>/ {
                id=$1
                sub(/^>/,"",id)
                print id,bin
            }
        ' "${F}" >> "${CON_TSV}"
    done

    # ---------------------------------------------------------
    # Crear lista maestra de contigs del assembly
    # ---------------------------------------------------------

    ASM_IDS="${OUT}/${G}_assembly_contigs.txt"

    awk '
        /^>/ {
            id=$1
            sub(/^>/,"",id)
            print id
        }
    ' "${ASSEMBLY}" | sort -u > "${ASM_IDS}"

    # ---------------------------------------------------------
    # Validar duplicados dentro de cada binner
    # ---------------------------------------------------------

    for TSV in "${META_TSV}" "${MAX_TSV}" "${CON_TSV}"
    do

        DUP=$(
            cut -f1 "${TSV}" |
            sort |
            uniq -d |
            wc -l
        )

        if [[ "${DUP}" -ne 0 ]]; then
            echo "ERROR: contigs asignados a más de un bin en:" >&2
            echo "${TSV}" >&2
            exit 1
        fi

    done

    # ---------------------------------------------------------
    # Validar que todos los contigs existan en el assembly
    # ---------------------------------------------------------

    for TSV in "${META_TSV}" "${MAX_TSV}" "${CON_TSV}"
    do

        TMP_IDS="${OUT}/tmp_$(basename "${TSV}").ids"

        cut -f1 "${TSV}" |
        sort -u > "${TMP_IDS}"

        MISSING=$(
            comm -23 "${TMP_IDS}" "${ASM_IDS}" |
            wc -l
        )

        if [[ "${MISSING}" -ne 0 ]]; then
            echo "ERROR: ${MISSING} contigs de $(basename "${TSV}") no existen en assembly." >&2
            exit 1
        fi

        rm -f "${TMP_IDS}"

    done

    # ---------------------------------------------------------
    # Resumen
    # ---------------------------------------------------------

    META_CONTIGS=$(wc -l < "${META_TSV}")
    MAX_CONTIGS=$(wc -l < "${MAX_TSV}")
    CON_CONTIGS=$(wc -l < "${CON_TSV}")

    META_NBINS=$(cut -f2 "${META_TSV}" | sort -u | wc -l)
    MAX_NBINS=$(cut -f2 "${MAX_TSV}" | sort -u | wc -l)
    CON_NBINS=$(cut -f2 "${CON_TSV}" | sort -u | wc -l)

    echo "MetaBAT2: ${META_NBINS} bins / ${META_CONTIGS} contigs"
    echo "MaxBin2:  ${MAX_NBINS} bins / ${MAX_CONTIGS} contigs"
    echo "CONCOCT:  ${CON_NBINS} bins / ${CON_CONTIGS} contigs"

    echo
    echo "INPUTS DAS TOOL ${G} COMPLETADOS"
    echo

done

echo "============================================================"
echo "PREPARACION DAS TOOL COMPLETADA"
echo "============================================================"
