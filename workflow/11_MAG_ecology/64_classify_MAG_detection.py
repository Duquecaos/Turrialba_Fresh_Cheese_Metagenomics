#!/usr/bin/env python3

import csv
import os

ROOT = f"/scratch/global/{os.environ['USER']}/Shotgun_MAGs_Turrialba"

INFILE = os.path.join(
    ROOT,
    "41_MAG_abundance_matrices",
    "MAG_abundance_long_18x18.tsv"
)

OUTDIR = os.path.join(
    ROOT,
    "42_MAG_detection_classification"
)

os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# HELPERS
# ============================================================

def num(x):
    try:
        if x is None:
            return None

        x = str(x).strip()

        if x == "" or x.upper() == "NA":
            return None

        return float(x)

    except ValueError:
        return None


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(rows)


# ============================================================
# LEER
# ============================================================

with open(INFILE, newline="") as f:
    rows = list(
        csv.DictReader(f, delimiter="\t")
    )

if len(rows) != 324:
    raise SystemExit(
        f"ERROR: esperábamos 324 combinaciones; hay {len(rows)}"
    )


samples = []

mags = []

for r in rows:
    if r["sample"] not in samples:
        samples.append(r["sample"])

    if r["MAG"] not in mags:
        mags.append(r["MAG"])


if len(samples) != 18:
    raise SystemExit(
        f"ERROR: esperábamos 18 muestras; hay {len(samples)}"
    )

if len(mags) != 18:
    raise SystemExit(
        f"ERROR: esperábamos 18 MAGs; hay {len(mags)}"
    )


# ============================================================
# CLASIFICAR
# ============================================================

classified = []

for r in rows:

    breadth = num(r["breadth"])
    coverage = num(r["coverage"])
    expected = num(r["breadth_expected"])

    reported = int(
        r["reported_by_instrain"]
    )

    if (
        expected is not None
        and expected > 0
        and breadth is not None
    ):
        ratio = breadth / expected
    else:
        ratio = None


    # --------------------------------------------------------
    # CLASIFICACIÓN
    # --------------------------------------------------------

    if reported == 0 or breadth == 0:

        category = "not_reported"

    elif ratio is not None and ratio < 0.50:

        category = "partial_inconsistent"

    elif breadth < 0.10:

        category = "very_low_breadth"

    elif ratio is None:

        category = "insufficient_expected"

    elif ratio < 0.80:

        category = "intermediate"

    elif breadth < 0.50:

        category = "coherent_limited_breadth"

    else:

        category = "robust_genome_wide"


    # --------------------------------------------------------
    # DOS DEFINICIONES PARA SENSIBILIDAD
    # --------------------------------------------------------

    # Principal:
    # evidencia coherente y >=10% del genoma.
    main_detection = int(
        category in {
            "coherent_limited_breadth",
            "robust_genome_wide"
        }
    )

    # Estricta:
    # >=50% del genoma y coherencia >=80%.
    strict_detection = int(
        category == "robust_genome_wide"
    )


    new = dict(r)

    new["breadth_to_expected_ratio"] = (
        "NA"
        if ratio is None
        else f"{ratio:.10g}"
    )

    new["detection_class"] = category
    new["main_detection"] = main_detection
    new["strict_detection"] = strict_detection

    classified.append(new)


# ============================================================
# TABLA LARGA
# ============================================================

fields = list(rows[0].keys()) + [
    "breadth_to_expected_ratio",
    "detection_class",
    "main_detection",
    "strict_detection"
]

write_tsv(
    os.path.join(
        OUTDIR,
        "MAG_detection_classification_18x18.tsv"
    ),
    classified,
    fields
)


# ============================================================
# MATRICES
# ============================================================

lookup = {
    (r["sample"], r["MAG"]): r
    for r in classified
}


def matrix(field, filename):

    with open(
        os.path.join(OUTDIR, filename),
        "w",
        newline=""
    ) as f:

        w = csv.writer(
            f,
            delimiter="\t",
            lineterminator="\n"
        )

        w.writerow(["sample"] + mags)

        for sample in samples:

            w.writerow(
                [sample]
                + [
                    lookup[(sample, mag)][field]
                    for mag in mags
                ]
            )


matrix(
    "detection_class",
    "MAG_detection_class_matrix.tsv"
)

matrix(
    "main_detection",
    "MAG_detection_main_matrix.tsv"
)

matrix(
    "strict_detection",
    "MAG_detection_strict_matrix.tsv"
)


# ============================================================
# COBERTURA ENMASCARADA
# ============================================================

for r in classified:

    coverage = num(r["coverage"]) or 0.0

    r["coverage_main"] = (
        coverage
        if int(r["main_detection"]) == 1
        else 0.0
    )

    r["coverage_strict"] = (
        coverage
        if int(r["strict_detection"]) == 1
        else 0.0
    )


matrix(
    "coverage_main",
    "coverage_main_detection_matrix.tsv"
)

matrix(
    "coverage_strict",
    "coverage_strict_detection_matrix.tsv"
)


# ============================================================
# RESUMEN POR CATEGORÍA
# ============================================================

categories = [
    "robust_genome_wide",
    "coherent_limited_breadth",
    "intermediate",
    "very_low_breadth",
    "partial_inconsistent",
    "insufficient_expected",
    "not_reported"
]

summary = []

for cat in categories:

    n = sum(
        r["detection_class"] == cat
        for r in classified
    )

    summary.append({
        "detection_class": cat,
        "count": n,
        "percent_of_324": f"{100*n/324:.3f}"
    })


write_tsv(
    os.path.join(
        OUTDIR,
        "detection_class_summary.tsv"
    ),
    summary,
    [
        "detection_class",
        "count",
        "percent_of_324"
    ]
)


# ============================================================
# RESUMEN POR MUESTRA
# ============================================================

sample_summary = []

for sample in samples:

    rr = [
        r
        for r in classified
        if r["sample"] == sample
    ]

    rec = {
        "sample": sample,
        "producer": rr[0]["producer"],
        "biological_unit": rr[0]["biological_unit"],
        "week": rr[0]["week"],
    }

    for cat in categories:
        rec[cat] = sum(
            r["detection_class"] == cat
            for r in rr
        )

    rec["main_detections"] = sum(
        int(r["main_detection"])
        for r in rr
    )

    rec["strict_detections"] = sum(
        int(r["strict_detection"])
        for r in rr
    )

    sample_summary.append(rec)


write_tsv(
    os.path.join(
        OUTDIR,
        "detection_summary_by_sample.tsv"
    ),
    sample_summary,
    [
        "sample",
        "producer",
        "biological_unit",
        "week",
        *categories,
        "main_detections",
        "strict_detections"
    ]
)


# ============================================================
# CASOS DE ALTA COBERTURA PERO INCONSISTENTES
# ============================================================

discordant = []

for r in classified:

    cov = num(r["coverage"]) or 0.0
    breadth = num(r["breadth"]) or 0.0
    ratio = num(r["breadth_to_expected_ratio"])

    if (
        cov >= 1
        and r["detection_class"]
        == "partial_inconsistent"
    ):
        discordant.append({
            "sample": r["sample"],
            "MAG": r["MAG"],
            "coverage": r["coverage"],
            "breadth": r["breadth"],
            "breadth_expected":
                r["breadth_expected"],
            "breadth_to_expected_ratio":
                r["breadth_to_expected_ratio"],
            "filtered_read_pair_count":
                r["filtered_read_pair_count"],
            "reads_mean_PID":
                r["reads_mean_PID"]
        })


discordant.sort(
    key=lambda x: float(x["coverage"]),
    reverse=True
)


write_tsv(
    os.path.join(
        OUTDIR,
        "high_coverage_partial_recruitment.tsv"
    ),
    discordant,
    [
        "sample",
        "MAG",
        "coverage",
        "breadth",
        "breadth_expected",
        "breadth_to_expected_ratio",
        "filtered_read_pair_count",
        "reads_mean_PID"
    ]
)


# ============================================================
# README
# ============================================================

with open(
    os.path.join(
        OUTDIR,
        "README_detection_criteria.txt"
    ),
    "w"
) as f:

    f.write(
        "Operational MAG detection framework\n\n"
    )

    f.write(
        "breadth_consistency = observed breadth / "
        "inStrain expected breadth\n\n"
    )

    f.write(
        "robust_genome_wide:\n"
        "  breadth >= 0.50 and consistency >= 0.80\n\n"
    )

    f.write(
        "coherent_limited_breadth:\n"
        "  0.10 <= breadth < 0.50 and consistency >= 0.80\n\n"
    )

    f.write(
        "intermediate:\n"
        "  breadth >= 0.10 and "
        "0.50 <= consistency < 0.80\n\n"
    )

    f.write(
        "partial_inconsistent:\n"
        "  consistency < 0.50\n\n"
    )

    f.write(
        "very_low_breadth:\n"
        "  breadth < 0.10 after excluding inconsistent recruitment\n\n"
    )

    f.write(
        "main_detection combines robust_genome_wide + "
        "coherent_limited_breadth.\n"
    )

    f.write(
        "strict_detection includes robust_genome_wide only.\n\n"
    )

    f.write(
        "These are operational analytical criteria, not universal "
        "biological thresholds. Intermediate and excluded signals "
        "must not automatically be interpreted as true absence.\n"
    )


print("============================================================")
print("PASO 64 COMPLETADO")
print("============================================================")
print(f"Muestras: {len(samples)}")
print(f"MAGs: {len(mags)}")
print(f"Combinaciones: {len(classified)}")
print()

for r in summary:
    print(
        f'{r["detection_class"]:28s}'
        f'{r["count"]:4d} '
        f'({r["percent_of_324"]}%)'
    )

print()
print(f"Salida: {OUTDIR}")
