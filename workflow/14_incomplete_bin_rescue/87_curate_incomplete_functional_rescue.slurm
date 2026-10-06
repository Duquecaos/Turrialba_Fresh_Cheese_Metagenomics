#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=03:00:00

#SBATCH --job-name=87_func_rescue
#SBATCH --output=87_func_rescue_%j.out
#SBATCH --error=87_func_rescue_%j.err

#SBATCH --mail-type=END,FAIL

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

SCRIPT="${SCRIPT_DIR}/87_curate_incomplete_functional_rescue.py"

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"
OUT="${ROOT}/73_incomplete_functional_rescue"

echo "============================================================"
echo "PASO 87 - RESCATE FUNCIONAL CURADO"
echo "Inicio: $(date)"
echo "Host: $(hostname)"
echo "Job: ${SLURM_JOB_ID}"
echo "Usuario: ${USER}"
echo "============================================================"

test -s "${SCRIPT}" || {
    echo "ERROR: falta ${SCRIPT}" >&2
    exit 1
}

test -s \
"${ROOT}/72_incomplete_eggnog_consolidated/gene_function_inventory_all_763306.tsv" || {
    echo "ERROR: falta inventario funcional 86c." >&2
    exit 1
}

test -s \
"${ROOT}/69_incomplete_unique_contigs/unique_contig_catalog_lt90.tsv" || {
    echo "ERROR: falta catálogo de contigs del Paso 85." >&2
    exit 1
}

test -s \
"${ROOT}/70_incomplete_gene_catalog/incomplete_lt90_new.prodigal.faa" || {
    echo "ERROR: falta FASTA de proteínas del Paso 86a." >&2
    exit 1
}

mkdir -p \
"${OUT}"

echo
echo "===== EJECUTANDO ====="
echo

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    python3 \
    "${SCRIPT}"

echo
echo "============================================================"
echo "VALIDANDO SALIDAS"
echo "============================================================"

REQUIRED=(
    "87A_cheese_metabolic_marker_hits.tsv"
    "87B_cheese_metabolic_systems_on_contigs.tsv"
    "87C_biogenic_amine_candidates.tsv"
    "87I_surface_proteinase_CEP_candidates.tsv"
    "87J_curated_lipolysis_candidates.tsv"
    "87K_curated_amino_acid_aroma_markers.tsv"
    "87L_EPS_capsule_loci.tsv"
    "87M_stress_systems_on_contigs.tsv"
    "functional_rescue_summary.tsv"
    "README_step87.txt"
)

for F in "${REQUIRED[@]}"; do

    if [[ ! -s "${OUT}/${F}" ]]; then
        echo "ERROR: falta o está vacío ${OUT}/${F}" >&2
        exit 1
    fi

    echo "OK: ${F}"
done

echo
echo "===== RESUMEN ====="

column -t -s $'\t' \
"${OUT}/functional_rescue_summary.tsv"

echo
echo "============================================================"
echo "PASO 87 FINALIZÓ CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "============================================================"
