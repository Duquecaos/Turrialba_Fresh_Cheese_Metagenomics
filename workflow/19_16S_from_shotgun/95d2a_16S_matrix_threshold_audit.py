#!/usr/bin/env python3

import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path(sys.argv[1])

D1 = (
    ROOT
    / "97_16S_shotgun"
    / "95D1_locus_quantification"
)

C2C = (
    ROOT
    / "97_16S_shotgun"
    / "95C2C_final_curated_bacterial_catalog"
)

INFILE = (
    D1
    / "95D1_all_samples_locus_metrics.tsv"
)

CATALOG = (
    C2C
    / "95C2C_final_bacterial_targets.tsv"
)

OUT = (
    ROOT
    / "97_16S_shotgun"
    / "95D2A_matrix_threshold_audit"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# Helpers
# ============================================================

def read_tsv(path):

    if not path.is_file():
        raise RuntimeError(
            f"Archivo faltante: {path}"
        )

    with path.open(
        newline="",
        errors="replace"
    ) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        return (
            list(reader),
            reader.fieldnames or []
        )


def write_tsv(path, rows, fields):

    with path.open(
        "w",
        newline=""
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows)


def fnum(x):

    try:
        return float(x)
    except Exception:
        return float("nan")


def quantile(values, q):

    values = sorted(
        x
        for x in values
        if not math.isnan(x)
    )

    if not values:
        return float("nan")

    if len(values) == 1:
        return values[0]

    pos = (
        (len(values) - 1)
        * q
    )

    lo = math.floor(pos)
    hi = math.ceil(pos)

    if lo == hi:
        return values[lo]

    return (
        values[lo]
        + (
            values[hi]
            - values[lo]
        )
        * (pos - lo)
    )


def fmt(x):

    if math.isnan(x):
        return "NA"

    return f"{x:.8f}"


# ============================================================
# Load
# ============================================================

rows, fields = read_tsv(
    INFILE
)

catalog, _ = read_tsv(
    CATALOG
)


if len(rows) != 363:
    raise RuntimeError(
        f"Locus-sample rows={len(rows)} expected=363"
    )

if len(catalog) != 118:
    raise RuntimeError(
        f"Final clusters={len(catalog)} expected=118"
    )


samples = sorted({
    r["sample"]
    for r in rows
})


if len(samples) != 18:
    raise RuntimeError(
        f"Samples={len(samples)} expected=18"
    )


clusters = sorted({
    r["exact_cluster_id"]
    for r in catalog
})


if len(clusters) != 118:
    raise RuntimeError(
        f"Unique final clusters={len(clusters)} expected=118"
    )


loci = sorted({
    r["locus_id"]
    for r in rows
})


if len(loci) != 121:
    raise RuntimeError(
        f"Unique loci={len(loci)} expected=121"
    )


# ============================================================
# Sample metadata
# ============================================================

sample_meta = []

for sample in samples:

    group = sample.split("_")[0]

    producer = group[0]

    biological_unit = (
        group[1:]
        if len(group) > 1
        else ""
    )

    time_code = (
        sample.rsplit(
            "_",
            1
        )[-1]
    )

    sample_meta.append({
        "sample":
            sample,

        "coassembly":
            group,

        "producer_code":
            producer,

        "biological_unit_code":
            biological_unit,

        "within_coassembly_time_code":
            time_code,

        "time_code_interpretation":
            "sample_name_code_only_not_biological_week_assignment",
    })


write_tsv(
    OUT
    / "95D2A_sample_metadata.tsv",
    sample_meta,
    list(
        sample_meta[0].keys()
    )
)


# ============================================================
# Matrix writer
# ============================================================

by_locus_sample = {
    (
        r["locus_id"],
        r["sample"]
    ): r
    for r in rows
}

by_cluster_sample = {
    (
        r["exact_cluster_id"],
        r["sample"]
    ): r
    for r in rows
}


def write_matrix(
    path,
    row_ids,
    lookup,
    metric,
    row_field
):

    matrix_rows = []

    for row_id in row_ids:

        rec = {
            row_field:
                row_id
        }

        for sample in samples:

            r = lookup.get(
                (
                    row_id,
                    sample
                )
            )

            # IMPORTANT:
            # NA = not evaluated in that native coassembly.
            # Never convert structural NA to zero.
            rec[sample] = (
                r[metric]
                if r is not None
                else "NA"
            )

        matrix_rows.append(
            rec
        )


    write_tsv(
        path,
        matrix_rows,
        [
            row_field
        ]
        + samples
    )


# ============================================================
# Locus matrices
# ============================================================

matrix_metrics = [
    (
        "all_mean_depth_per_million_mapped",
        "normalized_mean_depth_all"
    ),
    (
        "q10_mean_depth_per_million_mapped",
        "normalized_mean_depth_q10"
    ),
    (
        "all_breadth_1x_pct",
        "breadth1_all"
    ),
    (
        "q10_breadth_1x_pct",
        "breadth1_q10"
    ),
    (
        "all_breadth_5x_pct",
        "breadth5_all"
    ),
    (
        "all_breadth_10x_pct",
        "breadth10_all"
    ),
    (
        "q10_to_all_mean_depth_ratio",
        "q10_to_all_depth_ratio"
    ),
]


for metric, suffix in matrix_metrics:

    write_matrix(
        OUT
        / f"95D2A_locus_matrix_{suffix}.tsv",
        loci,
        by_locus_sample,
        metric,
        "locus_id"
    )


# ============================================================
# Exact-cluster matrices
#
# Within each coassembly every retained exact cluster corresponds
# to one primary locus. Cross-coassembly exact recurrence can
# populate additional samples. Unassessed samples remain NA.
# ============================================================

for metric, suffix in matrix_metrics:

    write_matrix(
        OUT
        / f"95D2A_cluster_matrix_{suffix}.tsv",
        clusters,
        by_cluster_sample,
        metric,
        "exact_cluster_id"
    )


# ============================================================
# Threshold sensitivity audit
#
# DESCRIPTIVE ONLY.
# This does NOT define biological presence/absence yet.
# ============================================================

all_thresholds = [
    50,
    75,
    90,
    95,
]

q10_thresholds = [
    0,
    25,
    50,
    75,
    90,
]


threshold_rows = []


for sample in samples:

    rr = [
        r
        for r in rows
        if r["sample"] == sample
    ]

    for all_t in all_thresholds:

        for q10_t in q10_thresholds:

            passed = [
                r
                for r in rr
                if (
                    fnum(
                        r[
                            "all_breadth_1x_pct"
                        ]
                    )
                    >= all_t
                    and
                    fnum(
                        r[
                            "q10_breadth_1x_pct"
                        ]
                    )
                    >= q10_t
                )
            ]


            threshold_rows.append({
                "sample":
                    sample,

                "coassembly":
                    rr[0][
                        "coassembly"
                    ],

                "total_target_loci":
                    len(rr),

                "all_breadth_threshold_pct":
                    all_t,

                "q10_breadth_threshold_pct":
                    q10_t,

                "loci_passing":
                    len(passed),

                "fraction_loci_passing":
                    f"{len(passed) / len(rr):.8f}",

                "classification_status":
                    "threshold_sensitivity_only_not_final_presence_call",
            })


write_tsv(
    OUT
    / "95D2A_threshold_sensitivity_by_sample.tsv",
    threshold_rows,
    list(
        threshold_rows[0].keys()
    )
)


# ============================================================
# Distribution audit by sample
# ============================================================

distribution_rows = []


for sample in samples:

    rr = [
        r
        for r in rows
        if r["sample"] == sample
    ]


    metrics = {
        "all_breadth_1x_pct":
            [
                fnum(
                    r[
                        "all_breadth_1x_pct"
                    ]
                )
                for r in rr
            ],

        "q10_breadth_1x_pct":
            [
                fnum(
                    r[
                        "q10_breadth_1x_pct"
                    ]
                )
                for r in rr
            ],

        "q10_to_all_mean_depth_ratio":
            [
                fnum(
                    r[
                        "q10_to_all_mean_depth_ratio"
                    ]
                )
                for r in rr
                if fnum(
                    r[
                        "all_mean_depth"
                    ]
                ) > 0
            ],

        "all_mean_depth_per_million_mapped":
            [
                fnum(
                    r[
                        "all_mean_depth_per_million_mapped"
                    ]
                )
                for r in rr
            ],
    }


    for metric, values in metrics.items():

        distribution_rows.append({
            "sample":
                sample,

            "coassembly":
                rr[0][
                    "coassembly"
                ],

            "metric":
                metric,

            "n":
                len(values),

            "min":
                fmt(
                    min(values)
                ),

            "q10":
                fmt(
                    quantile(
                        values,
                        0.10
                    )
                ),

            "q25":
                fmt(
                    quantile(
                        values,
                        0.25
                    )
                ),

            "median":
                fmt(
                    quantile(
                        values,
                        0.50
                    )
                ),

            "q75":
                fmt(
                    quantile(
                        values,
                        0.75
                    )
                ),

            "q90":
                fmt(
                    quantile(
                        values,
                        0.90
                    )
                ),

            "max":
                fmt(
                    max(values)
                ),
        })


write_tsv(
    OUT
    / "95D2A_metric_distributions_by_sample.tsv",
    distribution_rows,
    list(
        distribution_rows[0].keys()
    )
)


# ============================================================
# MAPQ ambiguity categories
# ============================================================

ambiguity_rows = []


for sample in samples:

    rr = [
        r
        for r in rows
        if r["sample"] == sample
    ]

    ratios = [
        fnum(
            r[
                "q10_to_all_mean_depth_ratio"
            ]
        )
        for r in rr
        if fnum(
            r[
                "all_mean_depth"
            ]
        ) > 0
    ]


    classes = {
        "ratio_lt_0.25":
            sum(
                x < 0.25
                for x in ratios
            ),

        "ratio_0.25_to_lt_0.50":
            sum(
                0.25 <= x < 0.50
                for x in ratios
            ),

        "ratio_0.50_to_lt_0.75":
            sum(
                0.50 <= x < 0.75
                for x in ratios
            ),

        "ratio_ge_0.75":
            sum(
                x >= 0.75
                for x in ratios
            ),
    }


    ambiguity_rows.append({
        "sample":
            sample,

        "coassembly":
            rr[0][
                "coassembly"
            ],

        "covered_loci":
            len(ratios),

        **classes,

        "interpretation":
            "low_ratio_indicates_MAPQ_sensitive_signal_not_absence",
    })


write_tsv(
    OUT
    / "95D2A_MAPQ_sensitivity_by_sample.tsv",
    ambiguity_rows,
    list(
        ambiguity_rows[0].keys()
    )
)


# ============================================================
# Taxon signal summary
#
# Sum is a cumulative 16S-locus signal.
# It is NOT organism/cell abundance because taxa may contain
# multiple rRNA operons and multiple assembled loci.
# ============================================================

taxon_groups = defaultdict(list)

for r in rows:

    key = (
        r["sample"],
        r[
            "deepest_informative_LCA_taxon"
        ]
    )

    taxon_groups[
        key
    ].append(r)


taxon_rows = []


for (
    sample,
    taxon
), rr in taxon_groups.items():

    all_values = [
        fnum(
            r[
                "all_mean_depth_per_million_mapped"
            ]
        )
        for r in rr
    ]

    q10_values = [
        fnum(
            r[
                "q10_mean_depth_per_million_mapped"
            ]
        )
        for r in rr
    ]

    taxon_rows.append({
        "sample":
            sample,

        "coassembly":
            rr[0][
                "coassembly"
            ],

        "taxon":
            taxon,

        "n_16S_loci":
            len(rr),

        "sum_normalized_locus_signal_all":
            f"{sum(all_values):.10f}",

        "max_normalized_locus_signal_all":
            f"{max(all_values):.10f}",

        "median_normalized_locus_signal_all":
            f"{statistics.median(all_values):.10f}",

        "sum_normalized_locus_signal_q10":
            f"{sum(q10_values):.10f}",

        "max_breadth_all_pct":
            f"{max(fnum(r['all_breadth_1x_pct']) for r in rr):.6f}",

        "max_breadth_q10_pct":
            f"{max(fnum(r['q10_breadth_1x_pct']) for r in rr):.6f}",

        "interpretation":
            "cumulative_16S_locus_signal_not_cell_abundance",
    })


taxon_rows.sort(
    key=lambda r: (
        r["sample"],
        -float(
            r[
                "sum_normalized_locus_signal_all"
            ]
        ),
        r["taxon"],
    )
)


write_tsv(
    OUT
    / "95D2A_taxon_locus_signal_by_sample.tsv",
    taxon_rows,
    list(
        taxon_rows[0].keys()
    )
)


# ============================================================
# Structural missingness audit for cluster matrix
# ============================================================

evaluated_cells = len(rows)

total_possible_cells = (
    len(clusters)
    * len(samples)
)

structural_na = (
    total_possible_cells
    - len(
        by_cluster_sample
    )
)


with (
    OUT
    / "95D2A_global_summary.tsv"
).open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    fh.write(
        f"samples\t{len(samples)}\n"
    )

    fh.write(
        f"primary_loci\t{len(loci)}\n"
    )

    fh.write(
        f"exact_clusters\t{len(clusters)}\n"
    )

    fh.write(
        f"locus_sample_observations\t{len(rows)}\n"
    )

    fh.write(
        f"cluster_sample_observations\t{len(by_cluster_sample)}\n"
    )

    fh.write(
        f"cluster_sample_possible_cells\t{total_possible_cells}\n"
    )

    fh.write(
        f"cluster_sample_structural_NA_cells\t{structural_na}\n"
    )

    fh.write(
        "structural_NA_interpretation\t"
        "not_evaluated_against_that_native_coassembly_not_zero\n"
    )

    fh.write(
        "final_presence_threshold_defined\tNO\n"
    )


# ============================================================
# Scope
# ============================================================

with (
    OUT
    / "95D2A_methodological_scope.tsv"
).open("w") as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        "data_source\t95D1R_existing_native_coassembly_mapping\n"
    )

    fh.write(
        "presence_absence_call\tNOT_YET_APPLIED\n"
    )

    fh.write(
        "structural_missing_values\tNA_not_zero\n"
    )

    fh.write(
        "all_mapping_metric\tBowtie2_primary_reported_alignment_depth\n"
    )

    fh.write(
        "MAPQ_sensitivity_metric\tMAPQ_ge10_depth\n"
    )

    fh.write(
        "normalized_signal\tmean_depth_per_million_mapped_alignments\n"
    )

    fh.write(
        "normalized_signal_equals_cell_abundance\tNO\n"
    )

    fh.write(
        "taxon_sum_equals_relative_abundance\tNO\n"
    )

    fh.write(
        "16S_copy_number_correction_applied\tNO\n"
    )

    fh.write(
        "cross_coassembly_unmeasured_cell_equals_absence\tNO\n"
    )


# ============================================================
# Validation
# ============================================================

required = [
    "95D2A_sample_metadata.tsv",
    "95D2A_threshold_sensitivity_by_sample.tsv",
    "95D2A_metric_distributions_by_sample.tsv",
    "95D2A_MAPQ_sensitivity_by_sample.tsv",
    "95D2A_taxon_locus_signal_by_sample.tsv",
    "95D2A_global_summary.tsv",
    "95D2A_methodological_scope.tsv",
]


for name in required:

    if not (
        OUT / name
    ).is_file():

        raise RuntimeError(
            f"Salida faltante: {name}"
        )


print(
    f"SAMPLES={len(samples)}"
)

print(
    f"LOCI={len(loci)}"
)

print(
    f"CLUSTERS={len(clusters)}"
)

print(
    f"OBSERVATIONS={len(rows)}"
)

print(
    "PRESENCE_THRESHOLD=NOT_APPLIED"
)

print(
    "95D2A=PASS"
)

