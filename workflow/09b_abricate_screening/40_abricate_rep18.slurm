#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00

#SBATCH --job-name=abricate_rep18
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

INPUTDIR="${ROOT}/22_final_representative_mags"
INPUT="${INPUTDIR}/representative18_traceable.fasta"
MASTER="${INPUTDIR}/representative_mags_master.tsv"

BASE="${ROOT}/24_abricate_rep18"
RAW="${BASE}/raw"

ABRICATE="${ABRICATE:-/opt/ohpc/pub/containers/BIO/abricate-1.0.1.sif}"

MINID=80
MINCOV=80

DATABASES=(
    resfinder
    card
    ncbi
    megares
    vfdb
    plasmidfinder
)

echo "============================================================"
echo "ABRICATE - 18 MAGs REPRESENTATIVOS"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Identidad mínima: ${MINID}%"
echo "Cobertura mínima: ${MINCOV}%"
echo "============================================================"

if [[ ! -s "${INPUT}" ]]; then
    echo "ERROR: falta ${INPUT}" >&2
    exit 1
fi

if [[ ! -s "${MASTER}" ]]; then
    echo "ERROR: falta ${MASTER}" >&2
    exit 1
fi

if [[ ! -s "${ABRICATE}" ]]; then
    echo "ERROR: falta ${ABRICATE}" >&2
    exit 1
fi

if [[ -e "${BASE}" ]]; then
    echo "ERROR: ya existe ${BASE}" >&2
    echo "No se sobrescribirá automáticamente." >&2
    exit 1
fi

mkdir -p "${RAW}"

echo
echo "===== BASES DISPONIBLES ====="

apptainer exec \
    --cleanenv \
    "${ABRICATE}" \
    abricate --list

echo
echo "===== INICIANDO SCREENING ====="

for DB in "${DATABASES[@]}"; do

    echo
    echo "----------------------------------------"
    echo "BASE: ${DB}"
    echo "----------------------------------------"

    OUT="${RAW}/${DB}.tsv"

    srun --exact \
        -n 1 \
        -c "${SLURM_CPUS_PER_TASK}" \
        apptainer exec \
        --cleanenv \
        --bind "${INPUTDIR}:/input:ro" \
        "${ABRICATE}" \
        abricate \
            --db "${DB}" \
            --minid "${MINID}" \
            --mincov "${MINCOV}" \
            --threads "${SLURM_CPUS_PER_TASK}" \
            /input/representative18_traceable.fasta \
            > "${OUT}"

    if [[ ! -s "${OUT}" ]]; then
        echo "ERROR: ${DB} no generó salida." >&2
        exit 1
    fi

    NHITS=$(awk '
        NR > 1 && NF > 0 {n++}
        END {print n+0}
    ' "${OUT}")

    echo "${DB}: ${NHITS} hits"

done


# ============================================================
# INTEGRACIÓN Y TRAZABILIDAD MAG
# ============================================================

COMBINED="${BASE}/abricate_all_hits.tsv"
HIGH="${BASE}/abricate_high_confidence_90_90.tsv"
SUMMARY="${BASE}/abricate_summary_by_MAG.tsv"

python3 - \
"${RAW}" \
"${MASTER}" \
"${COMBINED}" \
"${HIGH}" \
"${SUMMARY}" <<'PY'
import csv
import sys
from pathlib import Path
from collections import defaultdict

raw_dir = Path(sys.argv[1])
master_file = Path(sys.argv[2])
combined_file = Path(sys.argv[3])
high_file = Path(sys.argv[4])
summary_file = Path(sys.argv[5])

db_types = {
    "resfinder": "AMR",
    "card": "AMR",
    "ncbi": "AMR",
    "megares": "AMR",
    "vfdb": "virulence",
    "plasmidfinder": "plasmid_replicon",
}

mags = []

with open(master_file, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for r in reader:
        mags.append(r["representative_MAG"])

if len(mags) != 18:
    raise SystemExit(
        f"ERROR: tabla maestra contiene {len(mags)} MAGs"
    )

rows = []
all_raw_fields = []

for db, evidence_type in db_types.items():

    path = raw_dir / f"{db}.tsv"

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        if reader.fieldnames is None:
            raise SystemExit(
                f"ERROR: {path} no tiene cabecera"
            )

        for x in reader.fieldnames:
            if x not in all_raw_fields:
                all_raw_fields.append(x)

        for r in reader:

            seq = r.get("SEQUENCE", "").strip()

            if not seq:
                continue

            if "___" not in seq:
                raise SystemExit(
                    f"ERROR: contig no trazable: {seq}"
                )

            mag = seq.split("___", 1)[0]

            if mag not in mags:
                raise SystemExit(
                    f"ERROR: MAG desconocido en {db}: {mag}"
                )

            try:
                identity = float(r.get("%IDENTITY", 0) or 0)
                coverage = float(r.get("%COVERAGE", 0) or 0)
            except ValueError:
                identity = 0
                coverage = 0

            confidence = (
                "high_90_90"
                if identity >= 90 and coverage >= 90
                else "screen_80_80"
            )

            row = dict(r)
            row["MAG"] = mag
            row["source_db"] = db
            row["evidence_type"] = evidence_type
            row["screen_class"] = confidence

            rows.append(row)


leading = [
    "MAG",
    "source_db",
    "evidence_type",
    "screen_class",
]

fields = leading + [
    x for x in all_raw_fields
    if x not in leading
]

with open(combined_file, "w", newline="") as out:
    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t",
        extrasaction="ignore"
    )
    writer.writeheader()
    writer.writerows(rows)


high_rows = [
    r for r in rows
    if r["screen_class"] == "high_90_90"
]

with open(high_file, "w", newline="") as out:
    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t",
        extrasaction="ignore"
    )
    writer.writeheader()
    writer.writerows(high_rows)


# ------------------------------------------------------------
# Resumen por MAG y base
# ------------------------------------------------------------

counts = defaultdict(int)
unique_genes = defaultdict(set)

for r in rows:

    key = (r["MAG"], r["source_db"])
    counts[key] += 1

    gene = r.get("GENE", "").strip()
    if gene:
        unique_genes[key].add(gene)


with open(summary_file, "w") as out:

    out.write(
        "MAG\tdatabase\tevidence_type\t"
        "hits\tunique_gene_labels\n"
    )

    for mag in sorted(mags):
        for db, evidence_type in db_types.items():

            key = (mag, db)

            out.write(
                f"{mag}\t"
                f"{db}\t"
                f"{evidence_type}\t"
                f"{counts[key]}\t"
                f"{len(unique_genes[key])}\n"
            )


print(f"MAGs: {len(mags)}")
print(f"Hits totales 80/80: {len(rows)}")
print(f"Hits high-confidence 90/90: {len(high_rows)}")

for db in db_types:
    n = sum(
        1 for r in rows
        if r["source_db"] == db
    )
    print(f"{db}: {n}")
PY


# ============================================================
# VALIDACIONES
# ============================================================

for F in \
    "${COMBINED}" \
    "${HIGH}" \
    "${SUMMARY}"
do
    if [[ ! -s "${F}" ]]; then
        echo "ERROR: falta resultado ${F}" >&2
        exit 1
    fi
done


echo
echo "===== RESULTADOS POR BASE ====="

for DB in "${DATABASES[@]}"; do

    N=$(awk '
        NR > 1 && NF > 0 {n++}
        END {print n+0}
    ' "${RAW}/${DB}.tsv")

    printf "%-15s %s\n" "${DB}" "${N}"

done


echo
echo "===== RESUMEN NO CERO POR MAG ====="

awk -F '\t' '
NR == 1 || $4 > 0
' "${SUMMARY}" | column -t -s $'\t'


echo
echo "============================================================"
echo "ABRICATE REP18 COMPLETADO"
echo "AMR: ResFinder + CARD + NCBI + MEGARes"
echo "Virulencia: VFDB"
echo "Replicones: PlasmidFinder"
echo "Screen: >=80% identidad y >=80% cobertura"
echo "High confidence: >=90% / >=90%"
echo "Resultados: ${BASE}"
echo "Fin: $(date)"
echo "============================================================"
