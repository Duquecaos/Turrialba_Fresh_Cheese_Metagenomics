#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

#SBATCH --job-name=consolidate_MAGs
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail


# ============================================================
# CONFIGURACIÓN
# ============================================================

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

MAPBASE="${ROOT}/39_mag_abundance"
ISBASE="${ROOT}/40_instrain_MAG_profiles"

OUT="${ROOT}/41_MAG_abundance_matrices"

MANIFEST="${MAPBASE}/samples_hostfree.tsv"
MAG_MANIFEST="${MAPBASE}/reference/MAG_manifest.tsv"

mkdir -p "${OUT}"


echo "============================================================"
echo "CONSOLIDACIÓN DE ABUNDANCIA DE 18 MAGs × 18 MUESTRAS"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "============================================================"


# ============================================================
# PREFLIGHT
# ============================================================

test -s "${MANIFEST}" || {
    echo "ERROR: falta ${MANIFEST}" >&2
    exit 1
}

test -s "${MAG_MANIFEST}" || {
    echo "ERROR: falta ${MAG_MANIFEST}" >&2
    exit 1
}

N_SAMPLES=$(
    awk 'NR>1 {n++} END {print n+0}' "${MANIFEST}"
)

N_MAGS=$(
    awk 'NR>1 {n++} END {print n+0}' "${MAG_MANIFEST}"
)

echo "Muestras en manifiesto: ${N_SAMPLES}"
echo "MAGs en manifiesto: ${N_MAGS}"

if [[ "${N_SAMPLES}" -ne 18 ]]; then
    echo "ERROR: se esperaban 18 muestras." >&2
    exit 1
fi

if [[ "${N_MAGS}" -ne 18 ]]; then
    echo "ERROR: se esperaban 18 MAGs." >&2
    exit 1
fi


N_PROFILES=$(
    find "${ISBASE}" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -name 'PROFILE_COMPLETED' \
        | wc -l
)

echo "Perfiles inStrain completados: ${N_PROFILES}"

if [[ "${N_PROFILES}" -ne 18 ]]; then
    echo "ERROR: no existen 18 perfiles completos." >&2
    exit 1
fi


# ============================================================
# CONSOLIDAR CON PYTHON ESTÁNDAR
# ============================================================

python3 - \
    "${MANIFEST}" \
    "${MAG_MANIFEST}" \
    "${ISBASE}" \
    "${MAPBASE}" \
    "${OUT}" <<'PY'

import csv
import glob
import math
import os
import statistics
import sys


manifest_file = sys.argv[1]
mag_manifest_file = sys.argv[2]
isbase = sys.argv[3]
mapbase = sys.argv[4]
outdir = sys.argv[5]


# ------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------

def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def num_or_none(x):
    if x is None:
        return None

    x = str(x).strip()

    if x == "" or x.lower() in {"nan", "na", "none"}:
        return None

    try:
        return float(x)
    except ValueError:
        return None


def fmt(x):
    if x is None:
        return "NA"
    return f"{x:.10g}"


# ------------------------------------------------------------
# MUESTRAS
# ------------------------------------------------------------

samples = read_tsv(manifest_file)

required_sample_cols = {
    "sample",
    "producer",
    "biological_unit",
    "week",
}

if not required_sample_cols.issubset(samples[0].keys()):
    raise SystemExit(
        f"ERROR: columnas faltantes en {manifest_file}"
    )

sample_order = [r["sample"] for r in samples]
sample_meta = {r["sample"]: r for r in samples}


# ------------------------------------------------------------
# MAGs
# ------------------------------------------------------------

mag_rows = read_tsv(mag_manifest_file)

if "MAG" not in mag_rows[0]:
    raise SystemExit(
        f"ERROR: no existe columna MAG en {mag_manifest_file}"
    )

mag_order = [r["MAG"] for r in mag_rows]

if len(set(mag_order)) != 18:
    raise SystemExit(
        "ERROR: los MAG IDs del manifiesto no son 18 únicos"
    )


# ------------------------------------------------------------
# LEER TODOS LOS genome_info
# ------------------------------------------------------------

profiles = {}
all_instrain_fields = None

for sample in sample_order:

    pattern = os.path.join(
        isbase,
        f"{sample}.IS",
        "**",
        "*genome_info.tsv"
    )

    matches = glob.glob(pattern, recursive=True)

    if len(matches) != 1:
        raise SystemExit(
            f"ERROR: {sample} tiene {len(matches)} archivos genome_info"
        )

    gi = matches[0]

    rows = read_tsv(gi)

    if len(rows) == 0:
        raise SystemExit(
            f"ERROR: genome_info vacío para {sample}"
        )

    if all_instrain_fields is None:
        all_instrain_fields = list(rows[0].keys())

    by_mag = {}

    for row in rows:
        mag = row["genome"]

        if mag in by_mag:
            raise SystemExit(
                f"ERROR: MAG duplicado {mag} en {sample}"
            )

        if mag not in mag_order:
            raise SystemExit(
                f"ERROR: MAG inesperado {mag} en {sample}"
            )

        by_mag[mag] = row

    profiles[sample] = {
        "path": gi,
        "rows": by_mag
    }


# ------------------------------------------------------------
# LONG TABLE 18 × 18
# ------------------------------------------------------------

long_rows = []

for sample in sample_order:

    meta = sample_meta[sample]
    observed = profiles[sample]["rows"]

    for mag in mag_order:

        reported = mag in observed

        if reported:
            src = observed[mag]
        else:
            src = {}

        coverage_raw = num_or_none(src.get("coverage"))
        breadth_raw = num_or_none(src.get("breadth"))
        breadth_min_cov_raw = num_or_none(
            src.get("breadth_minCov")
        )

        # Zero-filled exclusivamente para matrices cuantitativas.
        # reported_by_instrain permite distinguir un cero rellenado
        # de un MAG realmente listado por inStrain.
        coverage = (
            coverage_raw
            if coverage_raw is not None
            else 0.0
        )

        breadth = (
            breadth_raw
            if breadth_raw is not None
            else 0.0
        )

        breadth_min_cov = (
            breadth_min_cov_raw
            if breadth_min_cov_raw is not None
            else 0.0
        )

        row = {
            "sample": sample,
            "producer": meta["producer"],
            "biological_unit": meta["biological_unit"],
            "week": meta["week"],
            "MAG": mag,
            "reported_by_instrain": int(reported),

            "coverage": fmt(coverage),
            "breadth": fmt(breadth),
            "breadth_minCov": fmt(breadth_min_cov),

            "coverage_raw": (
                fmt(coverage_raw)
                if reported
                else "NA"
            ),
            "breadth_raw": (
                fmt(breadth_raw)
                if reported
                else "NA"
            ),

            "breadth_expected": (
                src.get("breadth_expected", "NA")
                if reported
                else "NA"
            ),

            "conANI_reference": (
                src.get("conANI_reference", "NA")
                if reported
                else "NA"
            ),

            "popANI_reference": (
                src.get("popANI_reference", "NA")
                if reported
                else "NA"
            ),

            "filtered_read_pair_count": (
                src.get("filtered_read_pair_count", "NA")
                if reported
                else "NA"
            ),

            "reads_unfiltered_pairs": (
                src.get("reads_unfiltered_pairs", "NA")
                if reported
                else "NA"
            ),

            "reads_mean_PID": (
                src.get("reads_mean_PID", "NA")
                if reported
                else "NA"
            ),

            "length": (
                src.get("length", "NA")
                if reported
                else "NA"
            ),
        }

        long_rows.append(row)


if len(long_rows) != 324:
    raise SystemExit(
        f"ERROR: tabla larga tiene {len(long_rows)} filas; esperábamos 324"
    )


long_fields = [
    "sample",
    "producer",
    "biological_unit",
    "week",
    "MAG",
    "reported_by_instrain",
    "coverage",
    "breadth",
    "breadth_minCov",
    "coverage_raw",
    "breadth_raw",
    "breadth_expected",
    "conANI_reference",
    "popANI_reference",
    "filtered_read_pair_count",
    "reads_unfiltered_pairs",
    "reads_mean_PID",
    "length",
]

write_tsv(
    os.path.join(
        outdir,
        "MAG_abundance_long_18x18.tsv"
    ),
    long_rows,
    long_fields,
)


# ------------------------------------------------------------
# MATRICES
# ------------------------------------------------------------

lookup = {
    (r["sample"], r["MAG"]): r
    for r in long_rows
}


def write_matrix(metric, filename):

    path = os.path.join(outdir, filename)

    with open(path, "w", newline="") as fh:

        writer = csv.writer(
            fh,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow(["sample"] + mag_order)

        for sample in sample_order:

            row = [sample]

            for mag in mag_order:
                row.append(
                    lookup[(sample, mag)][metric]
                )

            writer.writerow(row)


write_matrix(
    "coverage",
    "coverage_matrix_18x18.tsv"
)

write_matrix(
    "breadth",
    "breadth_matrix_18x18.tsv"
)

write_matrix(
    "breadth_minCov",
    "breadth_minCov_matrix_18x18.tsv"
)

write_matrix(
    "reported_by_instrain",
    "reported_by_instrain_matrix_18x18.tsv"
)


# ------------------------------------------------------------
# RESUMEN POR MUESTRA
# ------------------------------------------------------------

sample_summary = []

for sample in sample_order:

    vals = [
        r
        for r in long_rows
        if r["sample"] == sample
    ]

    n_reported = sum(
        int(r["reported_by_instrain"])
        for r in vals
    )

    breadth_values = [
        float(r["breadth"])
        for r in vals
    ]

    coverage_values = [
        float(r["coverage"])
        for r in vals
    ]

    # Bowtie2 overall alignment rate
    bt2log = os.path.join(
        mapbase,
        "mapping",
        sample,
        f"{sample}.bowtie2.log"
    )

    rate = "NA"

    if os.path.isfile(bt2log):
        with open(bt2log) as fh:
            for line in fh:
                if "overall alignment rate" in line:
                    rate = line.split()[0]
                    break

    sample_summary.append({
        "sample": sample,
        "producer": sample_meta[sample]["producer"],
        "biological_unit": sample_meta[sample]["biological_unit"],
        "week": sample_meta[sample]["week"],
        "bowtie2_overall_alignment_rate": rate,
        "MAGs_reported_by_instrain": n_reported,
        "max_breadth": fmt(max(breadth_values)),
        "max_coverage": fmt(max(coverage_values)),
    })


write_tsv(
    os.path.join(
        outdir,
        "sample_QC_summary.tsv"
    ),
    sample_summary,
    [
        "sample",
        "producer",
        "biological_unit",
        "week",
        "bowtie2_overall_alignment_rate",
        "MAGs_reported_by_instrain",
        "max_breadth",
        "max_coverage",
    ],
)


# ------------------------------------------------------------
# DISTRIBUCIÓN GLOBAL DE BREADTH
# ------------------------------------------------------------

breadths = sorted(
    float(r["breadth"])
    for r in long_rows
)


def percentile(values, q):

    if not values:
        return None

    if len(values) == 1:
        return values[0]

    pos = (len(values) - 1) * q
    lo = math.floor(pos)
    hi = math.ceil(pos)

    if lo == hi:
        return values[lo]

    frac = pos - lo

    return (
        values[lo] * (1 - frac)
        + values[hi] * frac
    )


quantiles = [
    ("min", 0.00),
    ("p05", 0.05),
    ("p10", 0.10),
    ("p25", 0.25),
    ("median", 0.50),
    ("p75", 0.75),
    ("p90", 0.90),
    ("p95", 0.95),
    ("p99", 0.99),
    ("max", 1.00),
]

quantile_rows = []

for label, q in quantiles:
    quantile_rows.append({
        "statistic": label,
        "breadth": fmt(
            percentile(breadths, q)
        )
    })

write_tsv(
    os.path.join(
        outdir,
        "breadth_distribution_quantiles.tsv"
    ),
    quantile_rows,
    ["statistic", "breadth"],
)


# ------------------------------------------------------------
# BINS DESCRIPTIVOS — NO SON UMBRALES DE PRESENCIA
# ------------------------------------------------------------

bins = [
    ("0", 0.0, 0.0),
    ("0-0.01", 0.0, 0.01),
    ("0.01-0.05", 0.01, 0.05),
    ("0.05-0.10", 0.05, 0.10),
    ("0.10-0.25", 0.10, 0.25),
    ("0.25-0.50", 0.25, 0.50),
    ("0.50-0.75", 0.50, 0.75),
    ("0.75-0.90", 0.75, 0.90),
    ("0.90-1.00", 0.90, 1.0000000001),
]

bin_rows = []

for label, low, high in bins:

    if label == "0":
        n = sum(v == 0 for v in breadths)
    else:
        n = sum(
            (v > low) and (v <= high)
            for v in breadths
        )

    bin_rows.append({
        "breadth_bin": label,
        "count": n,
        "percent_of_324": fmt(
            100 * n / 324
        ),
    })


write_tsv(
    os.path.join(
        outdir,
        "breadth_descriptive_bins.tsv"
    ),
    bin_rows,
    [
        "breadth_bin",
        "count",
        "percent_of_324",
    ],
)


# ------------------------------------------------------------
# TOP MAGs POR BREADTH EN CADA MUESTRA
# ------------------------------------------------------------

top_rows = []

for sample in sample_order:

    vals = [
        r
        for r in long_rows
        if r["sample"] == sample
    ]

    vals.sort(
        key=lambda r: float(r["breadth"]),
        reverse=True
    )

    for rank, r in enumerate(vals[:5], start=1):

        top_rows.append({
            "sample": sample,
            "rank": rank,
            "MAG": r["MAG"],
            "breadth": r["breadth"],
            "coverage": r["coverage"],
            "reported_by_instrain": r[
                "reported_by_instrain"
            ],
        })


write_tsv(
    os.path.join(
        outdir,
        "top5_MAGs_by_breadth_per_sample.tsv"
    ),
    top_rows,
    [
        "sample",
        "rank",
        "MAG",
        "breadth",
        "coverage",
        "reported_by_instrain",
    ],
)


# ------------------------------------------------------------
# README / INTERPRETACIÓN DE CEROS
# ------------------------------------------------------------

readme = os.path.join(
    outdir,
    "README_abundance_matrices.txt"
)

with open(readme, "w") as fh:

    fh.write(
        "MAG abundance consolidation: 18 samples x 18 dRep representative MAGs\n"
    )

    fh.write(
        "Total combinations: 324\n\n"
    )

    fh.write(
        "IMPORTANT:\n"
    )

    fh.write(
        "- reported_by_instrain=1: MAG was present in the inStrain genome_info.tsv.\n"
    )

    fh.write(
        "- reported_by_instrain=0: MAG was not reported by inStrain for that sample.\n"
    )

    fh.write(
        "- coverage and breadth are zero-filled for unreported MAGs only to create complete numeric matrices.\n"
    )

    fh.write(
        "- Zero-filled values must not be interpreted automatically as biological absence.\n"
    )

    fh.write(
        "- coverage_raw and breadth_raw are NA for unreported MAGs.\n"
    )

    fh.write(
        "- No biological presence/absence threshold has been imposed in this step.\n"
    )

    fh.write(
        "- breadth_descriptive_bins.tsv is descriptive only and does not define presence.\n"
    )


print("============================================================")
print("CONSOLIDACIÓN COMPLETADA")
print("============================================================")
print(f"Muestras: {len(sample_order)}")
print(f"MAGs: {len(mag_order)}")
print(f"Combinaciones: {len(long_rows)}")
print(f"Salida: {outdir}")

PY


# ============================================================
# VALIDACIONES
# ============================================================

LONG="${OUT}/MAG_abundance_long_18x18.tsv"

test -s "${LONG}" || {
    echo "ERROR: no se generó tabla larga." >&2
    exit 1
}

N_LONG=$(
    awk 'NR>1 {n++} END {print n+0}' \
    "${LONG}"
)

echo
echo "Filas muestra×MAG: ${N_LONG}"

if [[ "${N_LONG}" -ne 324 ]]; then
    echo "ERROR: se esperaban 324 filas." >&2
    exit 1
fi


echo
echo "============================================================"
echo "QC POR MUESTRA"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/sample_QC_summary.tsv"


echo
echo "============================================================"
echo "DISTRIBUCIÓN DE BREADTH"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/breadth_distribution_quantiles.tsv"

echo

column -t -s $'\t' \
"${OUT}/breadth_descriptive_bins.tsv"


echo
echo "============================================================"
echo "PASO 63 COMPLETADO: MATRICES DE ABUNDANCIA DE MAGs"
echo "Directorio: ${OUT}"
echo "Combinaciones muestra×MAG: ${N_LONG}"
echo "No se aplicó ningún umbral de presencia/ausencia."
echo "Fin: $(date)"
echo "============================================================"
