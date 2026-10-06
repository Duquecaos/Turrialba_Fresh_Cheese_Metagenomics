#!/usr/bin/env python3

import csv
import math
import os
import statistics

USER = os.environ["USER"]

ROOT = f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"

CLASS_FILE = os.path.join(
    ROOT,
    "42_MAG_detection_classification",
    "MAG_detection_classification_18x18.tsv"
)

MASTER_FILE = os.path.join(
    ROOT,
    "22_final_representative_mags",
    "representative_mags_master.tsv"
)

OUTDIR = os.path.join(
    ROOT,
    "43_MAG_abundance_taxonomy"
)

os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# HELPERS
# ============================================================

def read_tsv(path):
    with open(path, newline="") as f:
        return list(
            csv.DictReader(f, delimiter="\t")
        )


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def num(value, default=None):
    if value is None:
        return default

    value = str(value).strip()

    if value == "" or value.upper() in {
        "NA",
        "N/A",
        "NAN",
        "NONE"
    }:
        return default

    try:
        return float(value)
    except ValueError:
        return default


def clean_taxon(value):
    if value is None:
        return ""

    value = str(value).strip()

    if value.upper() in {
        "",
        "NA",
        "N/A",
        "NONE"
    }:
        return ""

    return value


def fmt(x):
    if x is None:
        return "NA"

    return f"{x:.10g}"


# ============================================================
# LEER TABLAS
# ============================================================

classification = read_tsv(CLASS_FILE)
master = read_tsv(MASTER_FILE)

if len(classification) != 324:
    raise SystemExit(
        f"ERROR: clasificación tiene {len(classification)} filas; "
        "se esperaban 324."
    )

if len(master) != 18:
    raise SystemExit(
        f"ERROR: master tiene {len(master)} MAGs; se esperaban 18."
    )


# ============================================================
# VALIDAR IDs DE MAGs
# ============================================================

class_mags = sorted(
    set(r["MAG"] for r in classification)
)

master_mags = sorted(
    set(r["representative_MAG"] for r in master)
)

if class_mags != master_mags:

    only_class = sorted(
        set(class_mags) - set(master_mags)
    )

    only_master = sorted(
        set(master_mags) - set(class_mags)
    )

    print("MAGs solo en clasificación:")
    for x in only_class:
        print(x)

    print()

    print("MAGs solo en master:")
    for x in only_master:
        print(x)

    raise SystemExit(
        "ERROR: los MAG IDs no coinciden entre abundancia y taxonomía."
    )

print("MAG IDs clasificación ↔ master: OK")


# ============================================================
# ÍNDICE TAXONÓMICO
# ============================================================

tax = {
    r["representative_MAG"]: r
    for r in master
}


# ============================================================
# ETIQUETA TAXONÓMICA PARA FIGURAS
# ============================================================

def taxon_label(row):

    levels = [
        ("species", "species"),
        ("genus", "genus"),
        ("family", "family"),
        ("order", "order"),
        ("class", "class"),
        ("phylum", "phylum"),
    ]

    for field, rank in levels:

        value = clean_taxon(
            row.get(field)
        )

        if value:
            return value, rank

    return row["representative_MAG"], "MAG"


# ============================================================
# UNIR CLASIFICACIÓN + TAXONOMÍA
# ============================================================

joined = []

for r in classification:

    mag = r["MAG"]
    t = tax[mag]

    label, label_rank = taxon_label(t)

    coverage = num(r["coverage"], 0.0)

    main_detection = int(
        r["main_detection"]
    )

    strict_detection = int(
        r["strict_detection"]
    )

    coverage_main = (
        coverage
        if main_detection == 1
        else 0.0
    )

    coverage_strict = (
        coverage
        if strict_detection == 1
        else 0.0
    )

    # El prefijo antes de "__" identifica el coensamblaje
    # del representante.
    origin_coassembly = mag.split("__", 1)[0]

    if origin_coassembly.startswith("L"):
        origin_producer = "L"
    elif origin_coassembly.startswith("M"):
        origin_producer = "M"
    else:
        origin_producer = "NA"

    new = dict(r)

    new.update({
        "drep_cluster":
            t["drep_cluster"],

        "n_cluster_members":
            t["n_cluster_members"],

        "cluster_members":
            t["cluster_members"],

        "coassemblies_represented":
            t["coassemblies_represented"],

        "representative_source":
            t["representative_source"],

        "representative_origin_coassembly":
            origin_coassembly,

        "representative_origin_producer":
            origin_producer,

        "completeness":
            t["completeness"],

        "contamination":
            t["contamination"],

        "domain":
            t["domain"],

        "phylum":
            t["phylum"],

        "class":
            t["class"],

        "order":
            t["order"],

        "family":
            t["family"],

        "genus":
            t["genus"],

        "species":
            t["species"],

        "classification":
            t["classification"],

        "closest_genome_reference":
            t["closest_genome_reference"],

        "closest_genome_ani":
            t["closest_genome_ani"],

        "closest_genome_af":
            t["closest_genome_af"],

        "classification_method":
            t["classification_method"],

        "warnings":
            t["warnings"],

        "taxon_label":
            label,

        "taxon_label_rank":
            label_rank,

        "coverage_main":
            coverage_main,

        "coverage_strict":
            coverage_strict,
    })

    joined.append(new)


# ============================================================
# ORDEN DE MUESTRAS Y MAGs
# ============================================================

sample_order = []

for r in joined:
    if r["sample"] not in sample_order:
        sample_order.append(r["sample"])

mag_order = master_mags


if len(sample_order) != 18:
    raise SystemExit(
        f"ERROR: se detectaron {len(sample_order)} muestras."
    )


# ============================================================
# SUMAS DE COVERAGE POR MUESTRA
# ============================================================

sample_totals = {}

for sample in sample_order:

    rr = [
        r for r in joined
        if r["sample"] == sample
    ]

    total_main = sum(
        float(r["coverage_main"])
        for r in rr
    )

    total_strict = sum(
        float(r["coverage_strict"])
        for r in rr
    )

    sample_totals[sample] = {
        "main": total_main,
        "strict": total_strict,
    }


# ============================================================
# ABUNDANCIA RELATIVA
# ============================================================

for r in joined:

    sample = r["sample"]

    total_main = sample_totals[sample]["main"]
    total_strict = sample_totals[sample]["strict"]

    cov_main = float(
        r["coverage_main"]
    )

    cov_strict = float(
        r["coverage_strict"]
    )

    if total_main > 0:
        rel_main = (
            100.0
            * cov_main
            / total_main
        )
    else:
        rel_main = 0.0

    if total_strict > 0:
        rel_strict = (
            100.0
            * cov_strict
            / total_strict
        )
    else:
        rel_strict = 0.0

    r["relative_abundance_main_pct"] = rel_main
    r["relative_abundance_strict_pct"] = rel_strict

    r["sample_total_coverage_main"] = total_main
    r["sample_total_coverage_strict"] = total_strict


# ============================================================
# VALIDAR SUMAS ~100%
# ============================================================

for sample in sample_order:

    rr = [
        r for r in joined
        if r["sample"] == sample
    ]

    s_main = sum(
        float(r["relative_abundance_main_pct"])
        for r in rr
    )

    s_strict = sum(
        float(r["relative_abundance_strict_pct"])
        for r in rr
    )

    if sample_totals[sample]["main"] > 0:
        if not math.isclose(
            s_main,
            100.0,
            rel_tol=1e-9,
            abs_tol=1e-7
        ):
            raise SystemExit(
                f"ERROR: abundancia principal de {sample} suma {s_main}"
            )

    if sample_totals[sample]["strict"] > 0:
        if not math.isclose(
            s_strict,
            100.0,
            rel_tol=1e-9,
            abs_tol=1e-7
        ):
            raise SystemExit(
                f"ERROR: abundancia estricta de {sample} suma {s_strict}"
            )


# ============================================================
# TABLA LARGA FINAL
# ============================================================

long_fields = [
    "sample",
    "producer",
    "biological_unit",
    "week",

    "MAG",

    "drep_cluster",
    "representative_origin_coassembly",
    "representative_origin_producer",
    "n_cluster_members",
    "cluster_members",
    "coassemblies_represented",
    "representative_source",

    "completeness",
    "contamination",

    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "taxon_label",
    "taxon_label_rank",
    "classification",

    "closest_genome_reference",
    "closest_genome_ani",
    "closest_genome_af",
    "classification_method",
    "warnings",

    "reported_by_instrain",
    "coverage",
    "breadth",
    "breadth_minCov",
    "breadth_expected",
    "breadth_to_expected_ratio",
    "reads_mean_PID",
    "filtered_read_pair_count",

    "detection_class",
    "main_detection",
    "strict_detection",

    "coverage_main",
    "coverage_strict",

    "sample_total_coverage_main",
    "sample_total_coverage_strict",

    "relative_abundance_main_pct",
    "relative_abundance_strict_pct",
]


# Formatear solo al escribir
long_output = []

for r in joined:

    x = dict(r)

    for field in [
        "coverage_main",
        "coverage_strict",
        "sample_total_coverage_main",
        "sample_total_coverage_strict",
        "relative_abundance_main_pct",
        "relative_abundance_strict_pct",
    ]:
        x[field] = fmt(
            float(r[field])
        )

    long_output.append(x)


LONG_OUT = os.path.join(
    OUTDIR,
    "MAG_abundance_taxonomy_long.tsv"
)

write_tsv(
    LONG_OUT,
    long_output,
    long_fields
)


# ============================================================
# MATRICES
# ============================================================

lookup = {
    (r["sample"], r["MAG"]): r
    for r in joined
}


def write_matrix(field, filename):

    path = os.path.join(
        OUTDIR,
        filename
    )

    with open(path, "w", newline="") as f:

        writer = csv.writer(
            f,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow(
            ["sample"] + mag_order
        )

        for sample in sample_order:

            values = []

            for mag in mag_order:

                value = lookup[
                    (sample, mag)
                ][field]

                values.append(
                    fmt(float(value))
                )

            writer.writerow(
                [sample] + values
            )


write_matrix(
    "relative_abundance_main_pct",
    "relative_abundance_main_18x18.tsv"
)

write_matrix(
    "relative_abundance_strict_pct",
    "relative_abundance_strict_18x18.tsv"
)

write_matrix(
    "coverage_main",
    "coverage_main_18x18.tsv"
)

write_matrix(
    "coverage_strict",
    "coverage_strict_18x18.tsv"
)


# ============================================================
# TABLA TAXONÓMICA DE LOS 18 REPRESENTANTES
# ============================================================

tax_rows = []

for mag in mag_order:

    t = tax[mag]

    label, rank = taxon_label(t)

    origin = mag.split("__", 1)[0]

    producer = (
        "L"
        if origin.startswith("L")
        else "M"
        if origin.startswith("M")
        else "NA"
    )

    tax_rows.append({
        "MAG": mag,
        "drep_cluster":
            t["drep_cluster"],
        "representative_origin_coassembly":
            origin,
        "representative_origin_producer":
            producer,
        "coassemblies_represented":
            t["coassemblies_represented"],
        "completeness":
            t["completeness"],
        "contamination":
            t["contamination"],
        "phylum":
            t["phylum"],
        "class":
            t["class"],
        "order":
            t["order"],
        "family":
            t["family"],
        "genus":
            t["genus"],
        "species":
            t["species"],
        "taxon_label":
            label,
        "taxon_label_rank":
            rank,
        "warnings":
            t["warnings"],
    })


write_tsv(
    os.path.join(
        OUTDIR,
        "representative_MAG_taxonomy.tsv"
    ),
    tax_rows,
    [
        "MAG",
        "drep_cluster",
        "representative_origin_coassembly",
        "representative_origin_producer",
        "coassemblies_represented",
        "completeness",
        "contamination",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "taxon_label",
        "taxon_label_rank",
        "warnings",
    ]
)


# ============================================================
# RESUMEN POR MUESTRA
# ============================================================

sample_summary = []

for sample in sample_order:

    rr = [
        r for r in joined
        if r["sample"] == sample
    ]

    meta = rr[0]

    main_count = sum(
        int(r["main_detection"])
        for r in rr
    )

    strict_count = sum(
        int(r["strict_detection"])
        for r in rr
    )

    main_sorted = sorted(
        rr,
        key=lambda r:
            float(
                r["relative_abundance_main_pct"]
            ),
        reverse=True
    )

    top = main_sorted[0]

    sample_summary.append({
        "sample":
            sample,

        "producer":
            meta["producer"],

        "biological_unit":
            meta["biological_unit"],

        "week":
            meta["week"],

        "main_detected_MAGs":
            main_count,

        "strict_detected_MAGs":
            strict_count,

        "total_coverage_main":
            fmt(
                sample_totals[sample]["main"]
            ),

        "total_coverage_strict":
            fmt(
                sample_totals[sample]["strict"]
            ),

        "dominant_MAG_main":
            top["MAG"],

        "dominant_taxon_main":
            top["taxon_label"],

        "dominant_relative_abundance_main_pct":
            fmt(
                float(
                    top[
                        "relative_abundance_main_pct"
                    ]
                )
            ),
    })


write_tsv(
    os.path.join(
        OUTDIR,
        "sample_MAG_abundance_summary.tsv"
    ),
    sample_summary,
    [
        "sample",
        "producer",
        "biological_unit",
        "week",
        "main_detected_MAGs",
        "strict_detected_MAGs",
        "total_coverage_main",
        "total_coverage_strict",
        "dominant_MAG_main",
        "dominant_taxon_main",
        "dominant_relative_abundance_main_pct",
    ]
)


# ============================================================
# RESUMEN DESCRIPTIVO PRODUCTOR × SEMANA
# ============================================================

groups = {}

for r in joined:

    key = (
        r["producer"],
        r["week"],
        r["MAG"]
    )

    groups.setdefault(
        key,
        []
    ).append(
        float(
            r["relative_abundance_main_pct"]
        )
    )


producer_week_rows = []

for (
    producer,
    week,
    mag
), values in sorted(groups.items()):

    t = tax[mag]

    label, rank = taxon_label(t)

    mean = statistics.mean(values)

    sd = (
        statistics.stdev(values)
        if len(values) > 1
        else 0.0
    )

    producer_week_rows.append({
        "producer":
            producer,

        "week":
            week,

        "MAG":
            mag,

        "taxon_label":
            label,

        "taxon_label_rank":
            rank,

        "n_biological_units":
            len(values),

        "mean_relative_abundance_main_pct":
            fmt(mean),

        "sd_relative_abundance_main_pct":
            fmt(sd),

        "min_relative_abundance_main_pct":
            fmt(min(values)),

        "max_relative_abundance_main_pct":
            fmt(max(values)),
    })


write_tsv(
    os.path.join(
        OUTDIR,
        "producer_week_MAG_summary_main.tsv"
    ),
    producer_week_rows,
    [
        "producer",
        "week",
        "MAG",
        "taxon_label",
        "taxon_label_rank",
        "n_biological_units",
        "mean_relative_abundance_main_pct",
        "sd_relative_abundance_main_pct",
        "min_relative_abundance_main_pct",
        "max_relative_abundance_main_pct",
    ]
)


# ============================================================
# README METODOLÓGICO
# ============================================================

with open(
    os.path.join(
        OUTDIR,
        "README_MAG_abundance_taxonomy.txt"
    ),
    "w"
) as f:

    f.write(
        "MAG abundance and GTDB taxonomy integration\n"
    )

    f.write(
        "===========================================\n\n"
    )

    f.write(
        "Reference set: 18 dRep representative MAGs.\n"
    )

    f.write(
        "Samples: 18 host-filtered shotgun metagenomes.\n\n"
    )

    f.write(
        "MAIN abundance:\n"
    )

    f.write(
        "Coverage is retained only when main_detection=1.\n"
    )

    f.write(
        "main_detection includes robust_genome_wide and "
        "coherent_limited_breadth recruitment.\n\n"
    )

    f.write(
        "STRICT abundance:\n"
    )

    f.write(
        "Coverage is retained only for robust_genome_wide "
        "recruitment.\n\n"
    )

    f.write(
        "Within-sample relative abundance is calculated as:\n"
    )

    f.write(
        "100 * retained_MAG_coverage / "
        "sum(retained_MAG_coverage in that sample).\n\n"
    )

    f.write(
        "IMPORTANT INTERPRETATION:\n"
    )

    f.write(
        "These percentages describe relative abundance within "
        "the recovered and retained MAG set. They are NOT "
        "percentages of the total microbial community.\n\n"
    )

    f.write(
        "No completeness correction was applied to mean genomic "
        "coverage. Completeness and contamination are retained "
        "as metadata.\n\n"
    )

    f.write(
        "producer_week_MAG_summary_main.tsv is descriptive only. "
        "Biological units remain the experimental units for "
        "inferential analyses.\n"
    )


# ============================================================
# FINAL QC
# ============================================================

print()
print("============================================================")
print("PASO 65 COMPLETADO")
print("============================================================")
print(f"Muestras: {len(sample_order)}")
print(f"MAGs: {len(mag_order)}")
print(f"Filas tabla larga: {len(joined)}")
print()

print("Abundancia principal por muestra:")

for row in sample_summary:

    print(
        f'{row["sample"]}\t'
        f'n={row["main_detected_MAGs"]}\t'
        f'dominante={row["dominant_taxon_main"]}\t'
        f'{row["dominant_relative_abundance_main_pct"]}%'
    )

print()
print(f"Salida: {OUTDIR}")
