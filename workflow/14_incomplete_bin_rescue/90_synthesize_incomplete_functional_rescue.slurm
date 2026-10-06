#!/bin/bash

#SBATCH --job-name=90_func_synth
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=90_func_synth_%j.out
#SBATCH --error=90_func_synth_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

echo "============================================================"
echo "PASO 90 - SINTESIS DEL RESCATE FUNCIONAL"
echo "Inicio: $(date)"
echo "Host: $(hostname)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

SCRIPT="${SCRIPT_DIR}/90_synthesize_incomplete_functional_rescue.py"

[[ -s "$SCRIPT" ]] || {
    echo "ERROR: falta $SCRIPT" >&2
    exit 2
}

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    python3 "$SCRIPT"

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"
OUT="${ROOT}/76_incomplete_functional_synthesis"

EXPECTED=(
    "90A_global_capability_rescue_impact.tsv"
    "90B_producer_capability_rescue_impact.tsv"
    "90C_rescue_sequence_clusters_by_capability.tsv"
    "90D_high_confidence_rescue_evidence.tsv"
    "90E_priority_CEP_BA_integrated_summary.tsv"
    "functional_rescue_synthesis_summary.tsv"
    "README_step90.txt"
)

echo
echo "============================================================"
echo "VALIDANDO SALIDAS"
echo "============================================================"

for F in "${EXPECTED[@]}"; do

    [[ -s "${OUT}/${F}" ]] || {
        echo "ERROR: falta salida ${F}" >&2
        exit 3
    }

    echo "OK: ${F}"
done

echo
echo "============================================================"
echo "RESUMEN"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/functional_rescue_synthesis_summary.tsv"

echo
echo "============================================================"
echo "PASO 90 FINALIZÓ CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "============================================================"
