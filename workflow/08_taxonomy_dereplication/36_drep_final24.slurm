#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00

#SBATCH --job-name=drep_final24
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

GENOMES="${ROOT}/19_final_qc_mags/fastas"
PASS22="${ROOT}/14_gunc_refined/refined_candidates_90_5_GUNCpass.tsv"
POSTMDM="${ROOT}/17_checkm2_post_mdm/checkm2_original_vs_post_mdm.tsv"

BASE="${ROOT}/21_drep_final24"
WORK="${BASE}/drep"

DREP="${DREP:-/opt/ohpc/pub/containers/BIO/drep-3.5.0.sif}"

GENOMELIST="${BASE}/genomes.txt"
GENOMEINFO="${BASE}/genomeInfo.csv"

echo "============================================================"
echo "dRep - 24 MAGs QC-PASS"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "============================================================"

if [[ ! -s "${DREP}" ]]; then
    echo "ERROR: no existe ${DREP}" >&2
    exit 1
fi

N=$(find "${GENOMES}" -maxdepth 1 -type f -name '*.fa' | wc -l)

echo "MAGs de entrada: ${N}"

if [[ "${N}" -ne 24 ]]; then
    echo "ERROR: se esperaban 24 MAGs." >&2
    exit 1
fi

if [[ -e "${WORK}" ]]; then
    echo "ERROR: ya existe ${WORK}" >&2
    exit 1
fi

mkdir -p "${BASE}"

# ============================================================
# PREPARAR genomeInfo + listado de FASTA
# ============================================================

python3 - \
"${GENOMES}" \
"${PASS22}" \
"${POSTMDM}" \
"${GENOMELIST}" \
"${GENOMEINFO}" <<'PY'

import csv
import sys
from pathlib import Path

genomes_dir = Path(sys.argv[1])
pass22 = Path(sys.argv[2])
postmdm = Path(sys.argv[3])
genomelist = Path(sys.argv[4])
genomeinfo = Path(sys.argv[5])

fastas = sorted(genomes_dir.glob("*.fa"))

if len(fastas) != 24:
    raise SystemExit(f"ERROR: FASTA encontrados = {len(fastas)}")

quality = {}

# 22 MAGs originales
with open(pass22, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    for r in reader:
        name = f"{r['coassembly']}__{r['bin']}.fa"

        quality[name] = (
            float(r["completeness"]),
            float(r["contamination"])
        )

# 2 MAGs rescatados: usar CheckM2 POST-MDMcleaner
rescued = {
    "L2_maxbin2.004_sub": "L2__L2_maxbin2.004_sub.fa",
    "concoct_72": "L2__concoct_72.fa",
}

with open(postmdm, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    for r in reader:
        mag = r["MAG"]

        if mag in rescued:
            quality[rescued[mag]] = (
                float(r["post_completeness"]),
                float(r["post_contamination"])
            )

names = {f.name for f in fastas}

missing_quality = sorted(names - set(quality))
extra_quality = sorted(set(quality) - names)

if missing_quality:
    raise SystemExit(
        "ERROR: faltan métricas para: " +
        ", ".join(missing_quality)
    )

if extra_quality:
    raise SystemExit(
        "ERROR: métricas sin FASTA: " +
        ", ".join(extra_quality)
    )

if len(quality) != 24:
    raise SystemExit(
        f"ERROR: métricas para {len(quality)} MAGs"
    )

# dRep admite un archivo de texto con las rutas
with open(genomelist, "w") as out:
    for f in fastas:
        out.write(f"/genomes/{f.name}\n")

# dRep exige basename, completeness, contamination
with open(genomeinfo, "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow([
        "genome",
        "completeness",
        "contamination"
    ])

    for f in fastas:
        comp, cont = quality[f.name]

        writer.writerow([
            f.name,
            comp,
            cont
        ])

print(f"FASTAs: {len(fastas)}")
print(f"Métricas: {len(quality)}")
PY

echo
echo "===== genomeInfo.csv ====="
cat "${GENOMEINFO}"

echo
echo "===== DEPENDENCIAS ====="

apptainer exec \
--cleanenv \
"${DREP}" \
dRep check_dependencies

echo
echo "============================================================"
echo "INICIANDO dRep"
echo "ANI secundario: 95%"
echo "Cobertura mínima: 30%"
echo "Algoritmo: fastANI"
echo "============================================================"

srun --exact \
-n 1 \
-c "${SLURM_CPUS_PER_TASK}" \
apptainer exec \
--cleanenv \
--bind "${GENOMES}:/genomes:ro" \
--bind "${BASE}:/work" \
"${DREP}" \
dRep dereplicate \
/work/drep \
-g /work/genomes.txt \
-p "${SLURM_CPUS_PER_TASK}" \
-comp 90 \
-con 5 \
--genomeInfo /work/genomeInfo.csv \
--S_ani 0.95 \
--cov_thresh 0.30 \
--S_algorithm fastANI \
--gen_warnings

# ============================================================
# VALIDACIÓN
# ============================================================

DEREP="${WORK}/dereplicated_genomes"

if [[ ! -d "${DEREP}" ]]; then
    echo "ERROR: no se generó dereplicated_genomes." >&2
    exit 1
fi

NREP=$(find "${DEREP}" -maxdepth 1 -type f | wc -l)

echo
echo "============================================================"
echo "RESULTADO dRep"
echo "============================================================"
echo "MAGs entrada:          24"
echo "MAGs representativos:  ${NREP}"
echo

echo "===== REPRESENTANTES ====="

find "${DEREP}" \
-maxdepth 1 \
-type f \
-printf '%f\n' \
| sort

echo
echo "===== TABLAS dRep ====="

ls -lh "${WORK}/data_tables"

echo
echo "============================================================"
echo "DREP FINAL24 COMPLETADO"
echo "MAGs iniciales: 24"
echo "MAGs representativos: ${NREP}"
echo "Fin: $(date)"
echo "============================================================"
