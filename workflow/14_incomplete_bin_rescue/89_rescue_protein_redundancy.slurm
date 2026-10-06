#!/bin/bash
#SBATCH --job-name=89_rescue_redund
#SBATCH --partition=shared_64
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=89_rescue_redund_%j.out
#SBATCH --error=89_rescue_redund_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

echo "============================================================"
echo "PASO 89 - REDUNDANCIA PROTEICA DEL RESCATE FUNCIONAL"
echo "Inicio: $(date)"
echo "Host: $(hostname)"
echo "Job: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "============================================================"

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

OUT="${ROOT}/75_incomplete_functional_sequence_redundancy"

SCRIPT_A="${SCRIPT_DIR}/89a_prepare_rescue_protein_redundancy.py"
SCRIPT_B="${SCRIPT_DIR}/89b_summarize_rescue_protein_redundancy.py"

CONTAINER="${CONTAINER:-/opt/ohpc/pub/containers/BIO/eggnog-mapper-2.1.12.sif}"
APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"

FINALFAA="${ROOT}/23_eggnog_rep18/results/representative18.emapper.genepred.fasta"

QUERYFAA="${OUT}/selected_rescue_functional_proteins.faa"

DBFINAL="${OUT}/final18_proteins"
DBSELF="${OUT}/rescue_functional_proteins"

VSFINAL="${OUT}/89B_vs_final18.raw.tsv"
SELFOUT="${OUT}/89C_rescue_self_95id_90cov.raw.tsv"

mkdir -p "$OUT"

echo
echo "============================================================"
echo "VALIDANDO INPUTS"
echo "============================================================"

for F in \
    "$SCRIPT_A" \
    "$SCRIPT_B" \
    "$CONTAINER" \
    "$FINALFAA"
do
    if [[ ! -s "$F" ]]; then
        echo "ERROR: falta input requerido: $F" >&2
        exit 2
    fi
done

if [[ ! -x "$APPTAINER" ]]; then
    echo "ERROR: Apptainer no disponible: $APPTAINER" >&2
    exit 2
fi

echo "Inputs principales: OK"

echo
echo "============================================================"
echo "1. PREPARACION DEL SET FUNCIONAL"
echo "============================================================"

srun --exact -n 1 -c 1 \
    python3 "$SCRIPT_A"

if [[ ! -s "$QUERYFAA" ]]; then
    echo "ERROR: no se generó $QUERYFAA" >&2
    exit 2
fi

NQUERY=$(
    grep -c '^>' "$QUERYFAA"
)

echo "Proteínas funcionales seleccionadas: $NQUERY"

if [[ "$NQUERY" -lt 1 ]]; then
    echo "ERROR: conjunto funcional vacío." >&2
    exit 2
fi

echo
echo "============================================================"
echo "2. VALIDANDO DIAMOND"
echo "============================================================"

"$APPTAINER" exec \
    "$CONTAINER" \
    diamond version

echo
echo "============================================================"
echo "3. BASE DIAMOND FINAL18"
echo "============================================================"

rm -f "${DBFINAL}.dmnd"

srun --exact -n 1 -c "${SLURM_CPUS_PER_TASK}" \
"$APPTAINER" exec \
    -B "${ROOT}:/work" \
    "$CONTAINER" \
    diamond makedb \
        --in "/work/23_eggnog_rep18/results/representative18.emapper.genepred.fasta" \
        --db "/work/75_incomplete_functional_sequence_redundancy/final18_proteins"

test -s "${DBFINAL}.dmnd"

echo
echo "============================================================"
echo "4. RESCATE FUNCIONAL vs PROTEOMA FINAL18"
echo "============================================================"

rm -f "$VSFINAL"

srun --exact -n 1 -c "${SLURM_CPUS_PER_TASK}" \
"$APPTAINER" exec \
    -B "${ROOT}:/work" \
    "$CONTAINER" \
    diamond blastp \
        --query "/work/75_incomplete_functional_sequence_redundancy/selected_rescue_functional_proteins.faa" \
        --db "/work/75_incomplete_functional_sequence_redundancy/final18_proteins" \
        --out "/work/75_incomplete_functional_sequence_redundancy/89B_vs_final18.raw.tsv" \
        --outfmt 6 \
            qseqid \
            sseqid \
            pident \
            length \
            qlen \
            slen \
            evalue \
            bitscore \
        --very-sensitive \
        --evalue 1e-5 \
        --max-target-seqs 10 \
        --max-hsps 1 \
        --threads "${SLURM_CPUS_PER_TASK}"

test -f "$VSFINAL"

echo "Alineamientos contra final18:"
awk 'NF' "$VSFINAL" | wc -l

echo
echo "============================================================"
echo "5. BASE DIAMOND DEL PROPIO RESCATE"
echo "============================================================"

rm -f "${DBSELF}.dmnd"

srun --exact -n 1 -c "${SLURM_CPUS_PER_TASK}" \
"$APPTAINER" exec \
    -B "${ROOT}:/work" \
    "$CONTAINER" \
    diamond makedb \
        --in "/work/75_incomplete_functional_sequence_redundancy/selected_rescue_functional_proteins.faa" \
        --db "/work/75_incomplete_functional_sequence_redundancy/rescue_functional_proteins"

test -s "${DBSELF}.dmnd"

echo
echo "============================================================"
echo "6. ALL-vs-ALL PARA FAMILIAS CASI IDENTICAS"
echo "============================================================"

rm -f "$SELFOUT"

srun --exact -n 1 -c "${SLURM_CPUS_PER_TASK}" \
"$APPTAINER" exec \
    -B "${ROOT}:/work" \
    "$CONTAINER" \
    diamond blastp \
        --query "/work/75_incomplete_functional_sequence_redundancy/selected_rescue_functional_proteins.faa" \
        --db "/work/75_incomplete_functional_sequence_redundancy/rescue_functional_proteins" \
        --out "/work/75_incomplete_functional_sequence_redundancy/89C_rescue_self_95id_90cov.raw.tsv" \
        --outfmt 6 \
            qseqid \
            sseqid \
            pident \
            length \
            qlen \
            slen \
            evalue \
            bitscore \
        --very-sensitive \
        --id 95 \
        --query-cover 90 \
        --subject-cover 90 \
        --max-target-seqs 0 \
        --max-hsps 1 \
        --threads "${SLURM_CPUS_PER_TASK}"

test -f "$SELFOUT"

echo "Alineamientos all-vs-all conservados:"
awk 'NF' "$SELFOUT" | wc -l

echo
echo "============================================================"
echo "7. RESUMEN Y CLUSTERING"
echo "============================================================"

srun --exact -n 1 -c 1 \
    python3 "$SCRIPT_B"

echo
echo "============================================================"
echo "8. VALIDANDO SALIDAS"
echo "============================================================"

EXPECTED=(
    "89A_selected_rescue_functional_genes.tsv"
    "89B_vs_final18.raw.tsv"
    "89C_rescue_self_95id_90cov.raw.tsv"
    "89D_best_final18_protein_match.tsv"
    "89E_rescue_near_identity_cluster_members.tsv"
    "89F_rescue_near_identity_cluster_summary.tsv"
    "89G_match_class_by_context.tsv"
    "89H_priority_CEP_sequence_redundancy.tsv"
    "89I_context_supported_BA_sequence_redundancy.tsv"
    "89J_strong_BA_complete_lowcontam.tsv"
    "protein_redundancy_summary.tsv"
)

for F in "${EXPECTED[@]}"; do

    if [[ ! -s "${OUT}/${F}" ]]; then
        echo "ERROR: salida ausente o vacía: ${F}" >&2
        exit 2
    fi

    echo "OK: $F"
done

echo
echo "============================================================"
echo "RESUMEN"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/protein_redundancy_summary.tsv"

echo
echo "============================================================"
echo "PASO 89 FINALIZÓ CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "============================================================"
