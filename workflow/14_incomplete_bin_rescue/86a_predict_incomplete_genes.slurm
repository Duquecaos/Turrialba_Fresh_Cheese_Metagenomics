#!/bin/bash
#SBATCH --job-name=86a_prodigal
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=86a_prodigal_%j.out
#SBATCH --error=86a_prodigal_%j.err
#SBATCH --mail-type=END,FAIL

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

INPUTDIR="${ROOT}/69_incomplete_unique_contigs"

INPUT="${INPUTDIR}/unique_contigs_lt90_not_in_final18.fna"

OUT="${ROOT}/70_incomplete_gene_catalog"

EGGNOG="${EGGNOG:-/opt/ohpc/pub/containers/BIO/eggnog-mapper-2.1.12.sif}"

mkdir -p "$OUT"

FAA="${OUT}/incomplete_lt90_new.prodigal.faa"
FNA="${OUT}/incomplete_lt90_new.prodigal.fna"
GFF="${OUT}/incomplete_lt90_new.prodigal.gff"

echo "============================================================"
echo "PASO 86a - PREDICCION GENICA EN CONTIGS INCOMPLETOS"
echo "Inicio: $(date)"
echo "Host: $(hostname)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"
echo

echo "Input:"
ls -lh "$INPUT"

echo
echo "Contigs de entrada:"
grep -c '^>' "$INPUT"

echo
echo "Contenedor:"
ls -lh "$EGGNOG"

echo
echo "===== PREFLIGHT PRODIGAL ====="

apptainer exec \
    --cleanenv \
    "$EGGNOG" \
    prodigal -v 2>&1 | head -n 10

echo
echo "============================================================"
echo "INICIANDO PRODIGAL"
echo "Modo: metagenómico (-p meta)"
echo "============================================================"
echo

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    apptainer exec \
        --cleanenv \
        --bind "${INPUTDIR}:/input:ro" \
        --bind "${OUT}:/work" \
        "$EGGNOG" \
        prodigal \
            -p meta \
            -i /input/unique_contigs_lt90_not_in_final18.fna \
            -a /work/incomplete_lt90_new.prodigal.faa \
            -d /work/incomplete_lt90_new.prodigal.fna \
            -f gff \
            -o /work/incomplete_lt90_new.prodigal.gff

echo
echo "============================================================"
echo "VALIDANDO OUTPUTS DE PRODIGAL"
echo "============================================================"

for F in \
    "$FAA" \
    "$FNA" \
    "$GFF"
do
    if [[ ! -s "$F" ]]; then
        echo "ERROR: output faltante/vacío: $F" >&2
        exit 1
    fi
done

NFAA=$(grep -c '^>' "$FAA")
NFNA=$(grep -c '^>' "$FNA")
NGFF=$(awk '$0 !~ /^#/ && $3=="CDS"{n++} END{print n+0}' "$GFF")

echo "Proteínas: ${NFAA}"
echo "Genes nt:  ${NFNA}"
echo "CDS GFF:   ${NGFF}"

if [[ "$NFAA" -ne "$NFNA" || "$NFAA" -ne "$NGFF" ]]; then
    echo "ERROR: conteos Prodigal inconsistentes." >&2
    exit 1
fi

echo
echo "============================================================"
echo "POSTPROCESANDO TRAZABILIDAD GEN -> CONTIG -> BIN"
echo "============================================================"

srun --exact \
    -n 1 \
    -c 1 \
    python3 \
    "${SCRIPT_DIR}/86a_build_incomplete_gene_catalog.py"

echo
echo "============================================================"
echo "RESUMEN FINAL"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/gene_catalog_summary.tsv"

echo
echo "============================================================"
echo "PASO 86a FINALIZÓ CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "============================================================"
