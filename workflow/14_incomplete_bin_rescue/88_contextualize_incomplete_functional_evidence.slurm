#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00

#SBATCH --job-name=88_func_context
#SBATCH --output=88_func_context_%j.out
#SBATCH --error=88_func_context_%j.err

#SBATCH --mail-type=END,FAIL

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

echo "============================================================"
echo "PASO 88 - CONTEXTUALIZACION DE EVIDENCIA FUNCIONAL"
echo "Inicio: $(date)"
echo "Host: $(hostname)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

SCRIPT="${SCRIPT_DIR}/88_contextualize_incomplete_functional_evidence.py"

if [[ ! -s "$SCRIPT" ]]; then
    echo "ERROR: falta script Python:" >&2
    echo "$SCRIPT" >&2
    exit 2
fi

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    python3 "$SCRIPT"

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"
OUT="${ROOT}/74_incomplete_functional_evidence_context"

echo
echo "============================================================"
echo "VALIDANDO SALIDAS"
echo "============================================================"

REQUIRED=(
    "88A_functional_contig_context_best.tsv"
    "88B_functional_contig_bin_support_long.tsv"
    "88C_evidence_summary_by_table_context_integrity.tsv"
    "88D_priority_CEP_like_candidates.tsv"
    "88E_context_supported_biogenic_amine_candidates.tsv"
    "functional_evidence_context_summary.tsv"
    "README_step88.txt"
)

for F in "${REQUIRED[@]}"; do

    if [[ ! -s "${OUT}/${F}" ]]; then
        echo "ERROR: falta o está vacío ${OUT}/${F}" >&2
        exit 3
    fi

    echo "OK: ${F}"
done

NCTX=$(find "$OUT" \
    -maxdepth 1 \
    -type f \
    -name 'contextualized_87*.tsv' \
    | wc -l)

echo
echo "Tablas 87 contextualizadas: ${NCTX}/13"

if [[ "$NCTX" -ne 13 ]]; then
    echo "ERROR: se esperaban 13 tablas contextualizadas." >&2
    exit 4
fi

echo
echo "===== RESUMEN ====="

column -t -s $'\t' \
    "${OUT}/functional_evidence_context_summary.tsv"

echo
echo "============================================================"
echo "PASO 88 FINALIZÓ CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "============================================================"
