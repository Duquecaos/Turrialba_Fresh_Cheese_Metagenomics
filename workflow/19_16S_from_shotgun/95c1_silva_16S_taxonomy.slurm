#!/bin/bash

#SBATCH --job-name=95c1_16S_SILVA
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=95c1_16S_SILVA_%j.out
#SBATCH --error=95c1_16S_SILVA_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"
BLAST="${BLAST:-/opt/ohpc/pub/containers/BIO/blast-2.16.0.sif}"

QUERY="${ROOT}/97_16S_shotgun/95B2_consolidated/95B2_unique_exact_representatives.fasta"

DB="${ROOT}/databases/mdmcleaner/gtdb/SILVA_138.2_SSURef_NR99_tax_silva"

OUT="${ROOT}/97_16S_shotgun/95C1_SILVA_taxonomy"

PARSER="${SCRIPT_DIR}/95c1_parse_silva_16S.py"

RAW="${OUT}/95C1_SILVA1382_blast_raw.tsv"

rm -rf "$OUT"
mkdir -p "$OUT"

export LC_ALL=C
export LANG=C
export LANGUAGE=C

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

echo "============================================================"
echo "95C1 - 16S SHOTGUN vs SILVA 138.2 SSURef NR99"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "============================================================"

[[ -s "$QUERY" ]] || {
    echo "ERROR: query faltante" >&2
    exit 2
}

[[ -s "${DB}.nsq" ]] || {
    echo "ERROR: BLAST DB no encontrada: $DB" >&2
    exit 3
}

[[ -s "$PARSER" ]] || {
    echo "ERROR: parser faltante" >&2
    exit 4
}

NQUERY=$(grep -c '^>' "$QUERY")

if [[ "$NQUERY" -ne 148 ]]
then
    echo "ERROR: query_sequences=$NQUERY expected=148" >&2
    exit 5
fi

echo
echo "Validando SILVA DB..."

"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$BLAST" \
    blastdbcmd \
        -db "$DB" \
        -info \
    > "${OUT}/95C1_blastdb_info.txt"

echo
echo "Ejecutando BLASTn..."

srun --exact \
    -n 1 \
    -c "$SLURM_CPUS_PER_TASK" \
    "$APPTAINER" exec \
        --bind "${ROOT}:${ROOT}" \
        "$BLAST" \
        blastn \
            -task blastn \
            -query "$QUERY" \
            -db "$DB" \
            -evalue 1e-10 \
            -max_target_seqs 100 \
            -max_hsps 1 \
            -num_threads "$SLURM_CPUS_PER_TASK" \
            -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs stitle' \
            -out "$RAW"

[[ -s "$RAW" ]] || {
    echo "ERROR: BLAST raw vacío" >&2
    exit 6
}

echo
echo "Parseando taxonomía..."

python3 \
    "$PARSER" \
    "$ROOT"

[[ -s "${OUT}/95C1_global_summary.tsv" ]] || {
    echo "ERROR: global summary faltante" >&2
    exit 7
}

[[ -s "${OUT}/95C1_exact_representative_taxonomy.tsv" ]] || {
    echo "ERROR: taxonomy table faltante" >&2
    exit 8
}

N=$(
    tail -n +2 \
    "${OUT}/95C1_exact_representative_taxonomy.tsv" \
    | wc -l
)

if [[ "$N" -ne 148 ]]
then
    echo "ERROR: taxonomy rows=$N expected=148" >&2
    exit 9
fi

printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"${SLURM_JOB_ID}" \
> "${OUT}/95C1_COMPLETE.ok"

echo
echo "============================================================"
echo "GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/95C1_global_summary.tsv"

echo
echo "============================================================"
echo "RESOLUTION"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/95C1_resolution_summary.tsv"

echo
echo "============================================================"
echo "95C1 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"

