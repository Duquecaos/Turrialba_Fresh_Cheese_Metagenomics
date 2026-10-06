#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=00:30:00

#SBATCH --job-name=bagel4_petauri_diag
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

ENV="${ROOT}/envs/bagel4-1.2"
BAGEL="${ROOT}/tools/BAGEL4_hpc"

SOURCE="${ROOT}/36_bagel4_BAL4/L2__L2_maxbin2.004_sub/session/k141_179388.47.fna"

TEST="${ROOT}/37_petauri_bagel4_diagnostic"
QUERY="k141_179388.47"

export PATH="${ENV}/bin:${PATH}"
export PERL5LIB="${BAGEL}/lib"
export LANG=C
export LC_ALL=C

echo "============================================================"
echo "BAGEL4 DIAGNOSTIC - L. petauri k141_179388"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "============================================================"

test -s "${SOURCE}" || {
    echo "ERROR: falta ${SOURCE}" >&2
    exit 1
}

rm -rf "${TEST}"
mkdir -p "${TEST}"

cp -p \
    "${SOURCE}" \
    "${TEST}/${QUERY}.fna"


# ============================================================
# GENERAR LAS MISMAS 6 TRADUCCIONES QUE BAGEL4
# ============================================================

perl \
    -I "${BAGEL}/lib" \
    -Mbagel4_functions \
    -Manne_genomics \
    -e '
        use strict;
        use warnings;

        my ($input, $output) = @ARGV;

        my $dna = anne_genomics::get_dna_from_fasta($input);
        my %sixframes = bagel4_functions::six_frames($dna);

        open my $fh, ">", $output
            or die "No puedo escribir $output: $!";

        for my $frame (sort {$a <=> $b} keys %sixframes) {
            print $fh ">frame_$frame\n";
            print $fh "$sixframes{$frame}\n";
        }

        close $fh;
    ' \
    "${TEST}/${QUERY}.fna" \
    "${TEST}/${QUERY}.six_frames.fna"


test -s "${TEST}/${QUERY}.six_frames.fna" || {
    echo "ERROR: no se generaron seis frames" >&2
    exit 1
}

echo
echo "Frames generados:"
grep '^>' "${TEST}/${QUERY}.six_frames.fna"


# ============================================================
# EJECUTAR EXACTAMENTE LA ETAPA BLAST AOI DE BAGEL4
# ============================================================

echo
echo "============================================================"
echo "BLAST AOI BAGEL4"
echo "============================================================"

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    perl "${BAGEL}/bagel4_AOI_identification_blast_bacteriocins.pl" \
        -s "${TEST}" \
        -dna "${TEST}/${QUERY}.fna" \
        -queryname "${QUERY}" \
    > "${TEST}/blast_AOI.stdout" \
    2> "${TEST}/blast_AOI.stderr"


# ============================================================
# MOSTRAR CLASE II ANTES/DESPUÉS DEL RECÍPROCO
# ============================================================

echo
echo "============================================================"
echo "RESULTADO BLAST CLASE II"
echo "============================================================"

CLASS2="${TEST}/${QUERY}.blast_bacteriocinII.results"

if [[ -s "${CLASS2}" ]]; then
    cat "${CLASS2}"
else
    echo "Archivo clase II vacío o inexistente"
fi


echo
echo "============================================================"
echo "HITS ELIMINADOS POR BLAST RECÍPROCO"
echo "============================================================"

if grep -F 'deleting' "${TEST}/blast_AOI.stdout"; then
    :
else
    echo "No aparece ningún mensaje 'deleting'"
fi


echo
echo "============================================================"
echo "AOI BLAST FINAL"
echo "============================================================"

AOI="${TEST}/${QUERY}.AOI.blast.table"

if [[ -s "${AOI}" ]]; then
    cat "${AOI}"
else
    echo "AOI BLAST final vacío"
fi


echo
echo "============================================================"
echo "ARCHIVOS DIAGNÓSTICOS"
echo "============================================================"

find "${TEST}" \
    -maxdepth 1 \
    -type f \
    -printf '%f\t%s bytes\n' \
| sort


echo
echo "============================================================"
echo "DIAGNÓSTICO BAGEL4 L. PETAURI COMPLETADO"
echo "Resultados: ${TEST}"
echo "Fin: $(date)"
echo "============================================================"
