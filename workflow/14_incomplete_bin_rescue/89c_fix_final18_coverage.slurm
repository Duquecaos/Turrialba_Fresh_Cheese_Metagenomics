#!/bin/bash
#SBATCH --job-name=89c_covfix
#SBATCH --partition=shared_64
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=89c_covfix_%j.out
#SBATCH --error=89c_covfix_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"
OUT="${ROOT}/75_incomplete_functional_sequence_redundancy"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"
CONTAINER="${CONTAINER:-/opt/ohpc/pub/containers/BIO/eggnog-mapper-2.1.12.sif}"

QUERY="${OUT}/selected_rescue_functional_proteins.faa"
DB="${OUT}/final18_proteins"

RAW="${OUT}/89B_vs_final18.coverage_corrected.raw.tsv"

PY="${SCRIPT_DIR}/89c_correct_final18_coverage.py"

echo "============================================================"
echo "PASO 89c - CORRECCION COBERTURA FINAL18"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

for F in \
    "$QUERY" \
    "${DB}.dmnd" \
    "$CONTAINER" \
    "$PY"
do
    [[ -s "$F" ]] || {
        echo "ERROR: falta $F" >&2
        exit 2
    }
done

rm -f "$RAW"

srun --exact -n 1 -c "${SLURM_CPUS_PER_TASK}" \
"$APPTAINER" exec \
    -B "${ROOT}:/work" \
    "$CONTAINER" \
    diamond blastp \
        --query "/work/75_incomplete_functional_sequence_redundancy/selected_rescue_functional_proteins.faa" \
        --db "/work/75_incomplete_functional_sequence_redundancy/final18_proteins" \
        --out "/work/75_incomplete_functional_sequence_redundancy/89B_vs_final18.coverage_corrected.raw.tsv" \
        --outfmt 6 \
            qseqid \
            sseqid \
            pident \
            length \
            qlen \
            slen \
            qstart \
            qend \
            sstart \
            send \
            evalue \
            bitscore \
        --very-sensitive \
        --evalue 1e-5 \
        --max-target-seqs 10 \
        --max-hsps 1 \
        --threads "${SLURM_CPUS_PER_TASK}"

[[ -s "$RAW" ]] || {
    echo "ERROR: alineamiento vacío." >&2
    exit 2
}

srun --exact -n 1 -c 1 \
    python3 "$PY"

for F in \
    89K_best_final18_protein_match_corrected.tsv \
    89L_priority_CEP_sequence_redundancy_corrected.tsv \
    89M_context_supported_BA_sequence_redundancy_corrected.tsv \
    89N_strong_BA_complete_lowcontam_corrected.tsv \
    89O_final18_match_summary_corrected.tsv
do
    [[ -s "${OUT}/${F}" ]] || {
        echo "ERROR: falta salida $F" >&2
        exit 2
    }

    echo "OK: $F"
done

echo
column -t -s $'\t' \
    "${OUT}/89O_final18_match_summary_corrected.tsv"

echo
echo "============================================================"
echo "PASO 89c FINALIZÓ CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "============================================================"
