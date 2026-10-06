#!/bin/bash

#SBATCH --partition=shared_64
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=24:00:00

#SBATCH --job-name=gtdbtk_final24
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

INPUT="${ROOT}/19_final_qc_mags/fastas"

BASE="${ROOT}/20_gtdbtk_final24"
RESULT="${BASE}/results"
SCRATCH="${BASE}/pplacer_scratch"
TMP="${BASE}/tmp"

GTDBTK="${GTDBTK:-/opt/ohpc/pub/containers/BIO/gtdbtk-2.4.0.sif}"
DB="${GTDBTK_DB:-/home/public/met-workshop/databases/release220}"

PREFIX="turrialba24"

echo "============================================================"
echo "GTDB-TK CLASIFICACIÓN FINAL 24 MAGs"
echo "============================================================"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memoria solicitada: 256G"
echo "GTDB-Tk: 2.4.0"
echo "GTDB release: r220"
echo "============================================================"
echo


# ============================================================
# VALIDACIONES
# ============================================================

if [[ ! -s "${GTDBTK}" ]]; then
    echo "ERROR: no existe contenedor GTDB-Tk:" >&2
    echo "${GTDBTK}" >&2
    exit 1
fi

if [[ ! -s "${DB}/metadata/metadata.txt" ]]; then
    echo "ERROR: base GTDB-Tk R220 no disponible:" >&2
    echo "${DB}" >&2
    exit 1
fi

N_MAG=$(find "${INPUT}" \
    -maxdepth 1 \
    -type f \
    -name '*.fa' \
    | wc -l)

echo "MAGs detectados: ${N_MAG}"

if [[ "${N_MAG}" -ne 24 ]]; then
    echo "ERROR: se esperaban exactamente 24 MAGs." >&2
    exit 1
fi

echo
echo "===== INPUT ====="

find "${INPUT}" \
    -maxdepth 1 \
    -type f \
    -name '*.fa' \
    -printf '%f\t%s bytes\n' \
    | sort


# ============================================================
# NO SOBREESCRIBIR RESULTADOS
# ============================================================

if [[ -e "${RESULT}" ]]; then
    echo "ERROR: ya existe output GTDB-Tk:" >&2
    echo "${RESULT}" >&2
    echo "No se sobrescribirá automáticamente." >&2
    exit 1
fi

mkdir -p "${BASE}" "${SCRATCH}" "${TMP}"


# ============================================================
# PREFLIGHT
# ============================================================

echo
echo "===== PREFLIGHT ====="

apptainer exec \
    --cleanenv \
    --env GTDBTK_DATA_PATH=/gtdbtk_db \
    --bind "${DB}:/gtdbtk_db:ro" \
    "${GTDBTK}" \
    gtdbtk --version

echo
echo "GTDBTK_DATA_PATH:"
echo "${DB}"

echo
echo "Metadata:"
ls -lh "${DB}/metadata/metadata.txt"


# ============================================================
# CLASSIFY_WF
# ============================================================

echo
echo "============================================================"
echo "INICIANDO GTDB-TK classify_wf"
echo "============================================================"
echo

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    apptainer exec \
    --cleanenv \
    --env GTDBTK_DATA_PATH=/gtdbtk_db \
    --bind "${DB}:/gtdbtk_db:ro" \
    --bind "${INPUT}:/genomes:ro" \
    --bind "${BASE}:/work" \
    "${GTDBTK}" \
    gtdbtk classify_wf \
        --skip_ani_screen \
        --genome_dir /genomes \
        --out_dir /work/results \
        --extension fa \
        --prefix "${PREFIX}" \
        --cpus "${SLURM_CPUS_PER_TASK}" \
        --pplacer_cpus 1 \
        --scratch_dir /work/pplacer_scratch \
        --tmpdir /work/tmp


# ============================================================
# VALIDAR RESULTADOS
# ============================================================

echo
echo "===== VALIDACIÓN RESULTADOS ====="

BAC="${RESULT}/${PREFIX}.bac120.summary.tsv"
ARC="${RESULT}/${PREFIX}.ar53.summary.tsv"

echo "Resumen bacteriano:"
if [[ -s "${BAC}" ]]; then
    ls -lh "${BAC}"
else
    echo "No generado / vacío"
fi

echo
echo "Resumen arqueano:"
if [[ -s "${ARC}" ]]; then
    ls -lh "${ARC}"
else
    echo "No generado / vacío"
fi


# ============================================================
# TABLA COMBINADA
# ============================================================

COMBINED="${BASE}/gtdbtk_final24_taxonomy.tsv"

python3 - "${BAC}" "${ARC}" "${COMBINED}" <<'PY'
import csv
import os
import sys

bac, arc, output = sys.argv[1:4]

files = []

if os.path.isfile(bac) and os.path.getsize(bac) > 0:
    files.append(("Bacteria", bac))

if os.path.isfile(arc) and os.path.getsize(arc) > 0:
    files.append(("Archaea", arc))

if not files:
    raise SystemExit(
        "ERROR: no existen summary.tsv válidos de GTDB-Tk"
    )

rows = []

for domain, path in files:

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for r in reader:

            user_genome = r.get("user_genome", "").strip()

            if not user_genome:
                continue

            rows.append({
                "user_genome": user_genome,
                "domain_marker": domain,
                "classification": r.get(
                    "classification", ""
                ),
                "closest_genome_reference": r.get(
                    "closest_genome_reference", ""
                ),
                "closest_genome_reference_radius": r.get(
                    "closest_genome_reference_radius", ""
                ),
                "closest_genome_taxonomy": r.get(
                    "closest_genome_taxonomy", ""
                ),
                "closest_genome_ani": r.get(
                    "closest_genome_ani", ""
                ),
                "closest_genome_af": r.get(
                    "closest_genome_af", ""
                ),
                "closest_placement_reference": r.get(
                    "closest_placement_reference", ""
                ),
                "classification_method": r.get(
                    "classification_method", ""
                ),
                "note": r.get("note", ""),
                "warnings": r.get("warnings", ""),
            })


unique = {r["user_genome"] for r in rows}

if len(unique) != 24:
    print(
        f"ERROR: GTDB-Tk produjo {len(unique)} "
        "genomas únicos; se esperaban 24.",
        file=sys.stderr,
    )

    print(
        "Genomas encontrados:",
        file=sys.stderr
    )

    for x in sorted(unique):
        print(x, file=sys.stderr)

    raise SystemExit(1)


fields = [
    "user_genome",
    "domain_marker",
    "classification",
    "closest_genome_reference",
    "closest_genome_reference_radius",
    "closest_genome_taxonomy",
    "closest_genome_ani",
    "closest_genome_af",
    "closest_placement_reference",
    "classification_method",
    "note",
    "warnings",
]

with open(output, "w", newline="") as out:

    writer = csv.DictWriter(
        out,
        delimiter="\t",
        fieldnames=fields
    )

    writer.writeheader()

    for r in sorted(
        rows,
        key=lambda x: x["user_genome"]
    ):
        writer.writerow(r)


print(f"Genomas únicos clasificados: {len(unique)}")
print(f"Tabla combinada: {output}")
PY


# ============================================================
# RESUMEN TAXONÓMICO
# ============================================================

echo
echo "===== TAXONOMÍA FINAL ====="

python3 - "${COMBINED}" <<'PY'
import csv
import sys

path = sys.argv[1]

with open(path, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    for r in reader:

        print(
            f"{r['user_genome']:28s}  "
            f"{r['classification']}"
        )
PY


# ============================================================
# LIMPIEZA TEMPORAL SOLO TRAS ÉXITO
# ============================================================

rm -rf "${SCRATCH}" "${TMP}"


echo
echo "============================================================"
echo "GTDB-TK FINAL24 COMPLETADO"
echo "MAGs clasificados: 24"
echo "GTDB-Tk: 2.4.0"
echo "GTDB release: R220"
echo "Tabla: ${COMBINED}"
echo "Fin: $(date)"
echo "============================================================"
