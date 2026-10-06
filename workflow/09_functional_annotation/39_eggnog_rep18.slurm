#!/bin/bash

#SBATCH --partition=shared_64
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=24:00:00

#SBATCH --job-name=eggnog_rep18
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

INPUTDIR="${ROOT}/22_final_representative_mags"
INPUT="${INPUTDIR}/representative18_traceable.fasta"

DB="${ROOT}/databases/eggnog_5.0.2"

BASE="${ROOT}/23_eggnog_rep18"
RESULT="${BASE}/results"
TMP="${BASE}/tmp"

EGGNOG="${EGGNOG:-/opt/ohpc/pub/containers/BIO/eggnog-mapper-2.1.12.sif}"

PREFIX="representative18"

echo "============================================================"
echo "eggNOG-mapper - 18 MAGs REPRESENTATIVOS"
echo "============================================================"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Input: ${INPUT}"
echo "eggNOG DB: 5.0.2"
echo "============================================================"
echo


# ============================================================
# VALIDACIONES
# ============================================================

if [[ ! -s "${INPUT}" ]]; then
    echo "ERROR: falta FASTA de entrada: ${INPUT}" >&2
    exit 1
fi

if [[ ! -s "${EGGNOG}" ]]; then
    echo "ERROR: falta contenedor: ${EGGNOG}" >&2
    exit 1
fi

REQUIRED_DB=(
    "${DB}/eggnog.db"
    "${DB}/eggnog_proteins.dmnd"
    "${DB}/eggnog.taxa.db"
    "${DB}/EGGNOG_DB_5.0.2_READY"
)

for F in "${REQUIRED_DB[@]}"; do
    if [[ ! -e "${F}" ]]; then
        echo "ERROR: falta componente eggNOG: ${F}" >&2
        exit 1
    fi
done

N_CONTIGS=$(grep -c '^>' "${INPUT}" || true)

echo "Contigs de entrada: ${N_CONTIGS}"

if [[ "${N_CONTIGS}" -ne 2463 ]]; then
    echo "ERROR: se esperaban 2463 contigs." >&2
    exit 1
fi


# ============================================================
# EVITAR SOBREESCRITURA
# ============================================================

if [[ -e "${RESULT}" ]]; then
    echo "ERROR: ya existe ${RESULT}" >&2
    echo "No se sobrescribirá automáticamente." >&2
    exit 1
fi

mkdir -p "${RESULT}" "${TMP}"


# ============================================================
# PREFLIGHT
# ============================================================

echo
echo "===== PREFLIGHT eggNOG ====="

apptainer exec \
    --cleanenv \
    --bind "${DB}:/eggnog_data:ro" \
    "${EGGNOG}" \
    emapper.py \
        --data_dir /eggnog_data \
        --version


# ============================================================
# ANOTACIÓN
# ============================================================

echo
echo "============================================================"
echo "INICIANDO eggNOG-mapper"
echo "Predicción génica: Prodigal"
echo "Tipo de entrada: metagenome"
echo "Búsqueda: DIAMOND"
echo "============================================================"
echo

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    apptainer exec \
    --cleanenv \
    --bind "${DB}:/eggnog_data:ro" \
    --bind "${INPUTDIR}:/input:ro" \
    --bind "${BASE}:/work" \
    "${EGGNOG}" \
    emapper.py \
        --cpu "${SLURM_CPUS_PER_TASK}" \
        -i /input/representative18_traceable.fasta \
        --itype metagenome \
        --genepred prodigal \
        -m diamond \
        --data_dir /eggnog_data \
        -o "${PREFIX}" \
        --output_dir /work/results \
        --temp_dir /work/tmp


# ============================================================
# VALIDAR OUTPUTS
# ============================================================

ANNOT="${RESULT}/${PREFIX}.emapper.annotations"
GENES="${RESULT}/${PREFIX}.emapper.genepred.fasta"
GFF="${RESULT}/${PREFIX}.emapper.genepred.gff"
HITS="${RESULT}/${PREFIX}.emapper.hits"
SEEDS="${RESULT}/${PREFIX}.emapper.seed_orthologs"

echo
echo "===== VALIDACIÓN OUTPUTS ====="

for F in \
    "${ANNOT}" \
    "${GENES}" \
    "${GFF}" \
    "${HITS}" \
    "${SEEDS}"
do
    if [[ ! -s "${F}" ]]; then
        echo "ERROR: output faltante/vacío: ${F}" >&2
        exit 1
    fi

    ls -lh "${F}"
done


# ============================================================
# RESUMEN POR MAG
# ============================================================

SUMMARY="${BASE}/eggnog_summary_by_MAG.tsv"

python3 - \
"${GENES}" \
"${ANNOT}" \
"${SUMMARY}" <<'PY'
import sys
from collections import Counter

genes_file = sys.argv[1]
annotations_file = sys.argv[2]
output_file = sys.argv[3]

predicted = Counter()
annotated = Counter()

# ------------------------------------------------------------
# Genes predichos por Prodigal
# ------------------------------------------------------------

with open(genes_file) as fh:
    for line in fh:
        if not line.startswith(">"):
            continue

        gene = line[1:].split()[0]

        if "___" not in gene:
            raise SystemExit(
                f"ERROR: gen sin MAG trazable: {gene}"
            )

        mag = gene.split("___", 1)[0]
        predicted[mag] += 1


# ------------------------------------------------------------
# Genes con registro de anotación eggNOG
# ------------------------------------------------------------

with open(annotations_file) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue

        query = line.split("\t", 1)[0]

        if "___" not in query:
            raise SystemExit(
                f"ERROR: query sin MAG trazable: {query}"
            )

        mag = query.split("___", 1)[0]
        annotated[mag] += 1


all_mags = sorted(set(predicted) | set(annotated))

if len(all_mags) != 18:
    raise SystemExit(
        f"ERROR: se detectaron {len(all_mags)} MAGs, "
        "se esperaban 18."
    )


with open(output_file, "w") as out:
    out.write(
        "MAG\tgenes_predicted\t"
        "genes_with_eggnog_record\tpercent_with_record\n"
    )

    for mag in all_mags:
        p = predicted[mag]
        a = annotated[mag]

        pct = 100 * a / p if p else 0

        out.write(
            f"{mag}\t{p}\t{a}\t{pct:.2f}\n"
        )


print(f"MAGs detectados: {len(all_mags)}")
print(f"Genes predichos: {sum(predicted.values())}")
print(
    "Genes con registro eggNOG:",
    sum(annotated.values())
)
PY


echo
echo "===== RESUMEN POR MAG ====="

column -t -s $'\t' "${SUMMARY}"


# ============================================================
# CONSERVAR TEMP SOLO SI ALGO FALLA
# ============================================================

rm -rf "${TMP}"


echo
echo "============================================================"
echo "EGGNOG REP18 COMPLETADO"
echo "MAGs: 18"
echo "Contigs: 2463"
echo "Anotaciones: ${ANNOT}"
echo "Resumen MAG: ${SUMMARY}"
echo "Fin: $(date)"
echo "============================================================"
