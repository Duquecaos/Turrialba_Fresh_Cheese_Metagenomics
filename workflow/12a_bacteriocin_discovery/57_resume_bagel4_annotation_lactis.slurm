#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00

#SBATCH --job-name=bagel4_resume_lactis
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

ENV="${ROOT}/envs/bagel4-1.2"
BAGEL="${ROOT}/tools/BAGEL4_hpc"
SESSION="${ROOT}/35_bagel4_test_lactis/session"

export PATH="${ENV}/bin:${PATH}"
export PERL5LIB="${BAGEL}/lib"

export LANG=C
export LC_ALL=C

echo "============================================================"
echo "BAGEL4 - REANUDAR ANOTACION Lactococcus lactis"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "============================================================"

test -d "${SESSION}" || {
    echo "ERROR: falta ${SESSION}" >&2
    exit 1
}

test -s "${BAGEL}/tables/codon_table.txt" || {
    echo "ERROR: falta codon_table.txt" >&2
    exit 1
}

if grep -q \
'/data/molgentools/genomics/codon_table.txt' \
"${BAGEL}/lib/anne_genomics.pm"
then
    echo "ERROR: anne_genomics conserva codon_table antiguo" >&2
    exit 1
fi

mapfile -t QUERIES < <(
    find "${SESSION}" \
        -maxdepth 1 \
        -type f \
        -name '*.AOI.table' \
        -size +0c \
        -printf '%f\n' \
    | sed 's/\.AOI\.table$//' \
    | sort
)

echo "AOIs positivos: ${#QUERIES[@]}"
printf '%s\n' "${QUERIES[@]}"

if [[ "${#QUERIES[@]}" -ne 5 ]]; then
    echo "ERROR: se esperaban 5 AOI positivos." >&2
    exit 1
fi

for Q in "${QUERIES[@]}"
do
    AOI="${SESSION}/${Q}.AOI.fna"

    test -s "${AOI}" || {
        echo "ERROR: AOI vacío: ${AOI}" >&2
        exit 1
    }

    echo
    echo "============================================================"
    echo "ANOTANDO ${Q}"
    echo "============================================================"

    rm -f \
        "${SESSION}/${Q}.GeneTable" \
        "${SESSION}/${Q}.GeneTables"

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        perl "${BAGEL}/bagel4_AOI_annotation.pl" \
            -s "${SESSION}" \
            -dna "${AOI}" \
            -queryname "${Q}" \
            -cpu "${SLURM_CPUS_PER_TASK}"

    test -s "${SESSION}/${Q}.GeneTables" || {
        echo "ERROR: no se genero ${Q}.GeneTables" >&2
        exit 1
    }

    printf "GeneTables OK: "
    wc -l < "${SESSION}/${Q}.GeneTables"
done

echo
echo "============================================================"
echo "RESUMEN"
echo "============================================================"

for Q in "${QUERIES[@]}"
do
    printf "%-25s " "${Q}"
    wc -l < "${SESSION}/${Q}.GeneTables"
done

echo
echo "============================================================"
echo "CANDIDATOS PRIORITARIOS"
echo "============================================================"

for Q in \
k141_64662.8 \
k141_78158.89
do
    echo
    echo "----- ${Q} -----"

    grep -Ei \
'lactococc|bacterioc|lcn|class.?II|immunity|transport|ABC|leader' \
"${SESSION}/${Q}.GeneTables" \
|| true
done

echo
echo "============================================================"
echo "BAGEL4 ANOTACION REANUDADA COMPLETADA"
echo "Fin: $(date)"
echo "============================================================"
