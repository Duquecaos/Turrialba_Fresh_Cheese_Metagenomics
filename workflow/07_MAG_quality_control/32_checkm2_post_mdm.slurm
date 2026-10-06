#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00

#SBATCH --job-name=checkm2_post_mdm
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

INPUT="${ROOT}/16_mdmcleaner_kept/fastas"

BASE="${ROOT}/17_checkm2_post_mdm"
RESULT="${BASE}/checkm2"
TMP="${BASE}/tmp"

CHECKM2="${ROOT}/containers/checkm2-1.1.0-r1.sif"

COMPARISON="${BASE}/checkm2_original_vs_post_mdm.tsv"

mkdir -p "${BASE}" "${TMP}"

echo "============================================================"
echo "CHECKM2 POST-MDMCLEANER"
echo "============================================================"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Input: ${INPUT}"
echo "============================================================"
echo


# ============================================================
# VALIDAR CONTENEDOR E INPUT
# ============================================================

if [[ ! -s "${CHECKM2}" ]]; then
    echo "ERROR: no existe contenedor CheckM2:" >&2
    echo "${CHECKM2}" >&2
    exit 1
fi

N_MAG=$(find "${INPUT}" \
    -maxdepth 1 \
    -type f \
    -name '*.fa' \
    | wc -l)

echo "MAGs detectados: ${N_MAG}"

if [[ "${N_MAG}" -ne 8 ]]; then
    echo "ERROR: se esperaban 8 MAGs." >&2
    exit 1
fi

echo
echo "===== MAGs ====="

find "${INPUT}" \
    -maxdepth 1 \
    -type f \
    -name '*.fa' \
    -printf '%f\t%s bytes\n' \
    | sort


# ============================================================
# LOCALIZAR LA MISMA BASE CHECKM2
# ============================================================

DB=""

CANDIDATES=(
    "${ROOT}/databases/checkm2/uniref100.KO.1.dmnd"
    "${ROOT}/databases/CheckM2_database/uniref100.KO.1.dmnd"
    "$HOME/DB/CheckM2_database/uniref100.KO.1.dmnd"
    "${CHECKM2_DB:-/home/public/DB/CheckM2_database/uniref100.KO.1.dmnd}"
)

for F in "${CANDIDATES[@]}"; do
    if [[ -s "${F}" ]]; then
        DB="${F}"
        break
    fi
done

if [[ -z "${DB}" ]]; then

    SEARCH_ROOTS=()

    [[ -d "${ROOT}/databases" ]] && \
        SEARCH_ROOTS+=("${ROOT}/databases")

    [[ -d "$HOME/DB" ]] && \
        SEARCH_ROOTS+=("$HOME/DB")

    if (( ${#SEARCH_ROOTS[@]} > 0 )); then
        mapfile -t FOUND < <(
            find "${SEARCH_ROOTS[@]}" \
                -type f \
                -name 'uniref100.KO.1.dmnd' \
                2>/dev/null
        )

        if (( ${#FOUND[@]} > 0 )); then
            DB="${FOUND[0]}"
        fi
    fi
fi

if [[ -z "${DB}" || ! -s "${DB}" ]]; then
    echo "ERROR: no se localizó la base CheckM2." >&2
    exit 1
fi

echo
echo "===== BASE CHECKM2 ====="
ls -lh "${DB}"


# ============================================================
# EVITAR SOBREESCRITURA
# ============================================================

if [[ -e "${RESULT}" ]]; then
    echo "ERROR: ya existe output:" >&2
    echo "${RESULT}" >&2
    exit 1
fi


# ============================================================
# VERSIÓN
# ============================================================

echo
echo "===== CHECKM2 ====="

apptainer exec \
    --cleanenv \
    "${CHECKM2}" \
    checkm2 --version || true


# ============================================================
# EJECUTAR
# ============================================================

echo
echo "============================================================"
echo "INICIANDO CHECKM2"
echo "============================================================"

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    apptainer exec \
    --cleanenv \
    --bind "${INPUT}:/input:ro" \
    --bind "${BASE}:/results" \
    --bind "${TMP}:/cm2tmp" \
    --bind "${DB}:/checkm2db/uniref100.KO.1.dmnd:ro" \
    "${CHECKM2}" \
    checkm2 predict \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --input /input \
    --output-directory /results/checkm2 \
    --tmpdir /cm2tmp \
    --database_path /checkm2db/uniref100.KO.1.dmnd \
    -x .fa


# ============================================================
# VALIDAR QUALITY REPORT
# ============================================================

REPORT="${RESULT}/quality_report.tsv"

if [[ ! -s "${REPORT}" ]]; then
    echo "ERROR: CheckM2 no generó quality_report.tsv" >&2
    exit 1
fi

N_RESULT=$(awk 'NR>1 && NF>0 {n++} END {print n+0}' "${REPORT}")

echo
echo "Genomas evaluados por CheckM2: ${N_RESULT}"

if [[ "${N_RESULT}" -ne 8 ]]; then
    echo "ERROR: se esperaban 8 resultados CheckM2." >&2
    exit 1
fi


# ============================================================
# COMPARAR CONTRA LOS VALORES PRE-MDMCLEANER
# ============================================================

python3 - "${REPORT}" "${COMPARISON}" <<'PY'
import csv
import sys
import os

report = sys.argv[1]
outfile = sys.argv[2]

# CheckM2 inmediatamente antes de MDMcleaner
baseline = {
    "L1_maxbin2.051":       (99.80, 3.94),
    "concoct_20":           (99.99, 2.28),
    "concoct_24":           (100.00, 4.92),
    "L2_maxbin2.004_sub":   (100.00, 4.01),
    "concoct_55":           (99.97, 3.54),
    "concoct_72":           (100.00, 1.98),
    "concoct_91":           (99.99, 0.43),
    "concoct_37":           (99.16, 0.97),
}

def normalize(name):
    x = os.path.basename(name)

    for ext in (
        ".fasta.gz", ".fna.gz", ".fa.gz",
        ".fasta", ".fna", ".fa"
    ):
        if x.endswith(ext):
            x = x[:-len(ext)]
            break

    # FASTA post-MDM: GRUPO__MAG.fa
    if "__" in x:
        x = x.split("__", 1)[1]

    return x


rows = []

with open(report, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    for r in reader:

        mag = normalize(r["Name"])

        if mag not in baseline:
            raise SystemExit(
                f"ERROR: MAG CheckM2 no reconocido: {r['Name']}"
            )

        pre_comp, pre_cont = baseline[mag]

        post_comp = float(r["Completeness"])
        post_cont = float(r["Contamination"])

        rows.append({
            "MAG": mag,
            "pre_completeness": pre_comp,
            "post_completeness": post_comp,
            "delta_completeness":
                post_comp - pre_comp,
            "pre_contamination": pre_cont,
            "post_contamination": post_cont,
            "delta_contamination":
                post_cont - pre_cont,
            "post_candidate_50_10":
                post_comp >= 50 and post_cont < 10,
            "post_strong_90_5":
                post_comp >= 90 and post_cont < 5,
        })


if len(rows) != 8:
    raise SystemExit(
        f"ERROR: se esperaban 8 filas y hay {len(rows)}"
    )


fields = [
    "MAG",
    "pre_completeness",
    "post_completeness",
    "delta_completeness",
    "pre_contamination",
    "post_contamination",
    "delta_contamination",
    "post_candidate_50_10",
    "post_strong_90_5",
]

with open(outfile, "w", newline="") as out:
    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()

    for r in rows:
        writer.writerow(r)


print("===== COMPARACIÓN CHECKM2 =====")

for r in rows:
    print(
        f"{r['MAG']:22s} "
        f"Comp {r['pre_completeness']:6.2f}"
        f" -> {r['post_completeness']:6.2f} "
        f"({r['delta_completeness']:+6.2f}) | "
        f"Cont {r['pre_contamination']:5.2f}"
        f" -> {r['post_contamination']:5.2f} "
        f"({r['delta_contamination']:+5.2f})"
    )
PY


echo
echo "===== RESULTADOS CHECKM2 ====="

column -t -s $'\t' "${COMPARISON}"

echo
echo "============================================================"
echo "CHECKM2 POST-MDMCLEANER COMPLETADO"
echo "Reporte: ${REPORT}"
echo "Comparación: ${COMPARISON}"
echo "Fin: $(date)"
echo "============================================================"
