#!/usr/bin/env python3

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(sys.argv[1])

D1 = (
    ROOT
    / "97_16S_shotgun"
    / "95D1_locus_quantification"
)

D2A = (
    ROOT
    / "97_16S_shotgun"
    / "95D2A_matrix_threshold_audit"
)

INFILE = (
    D1
    / "95D1_all_samples_locus_metrics.tsv"
)

OUT = (
    ROOT
    / "97_16S_shotgun"
    / "95D2B_final_support"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# Fixed thresholds chosen after 95D2A sensitivity audit
# ============================================================

ROBUST_ALL_BREADTH = 90.0
ROBUST_Q10_BREADTH = 75.0

PARTIAL_ALL_BREADTH = 50.0


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

        return list(reader)


def write_tsv(
    path,
    rows,
    fields
):

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
        return 0.0


def support_class(
    all_breadth,
    q10_breadth
):

    if all_breadth <= 0:
        return "no_mapped_signal"

    if (
        all_breadth >= ROBUST_ALL_BREADTH
        and
        q10_breadth >= ROBUST_Q10_BREADTH
    ):
        return "robust_mapping_support"

    if (
        all_breadth >= ROBUST_ALL_BREADTH
        and
        q10_breadth < ROBUST_Q10_BREADTH
    ):
        return "high_breadth_MAPQ_sensitive"

    if (
        all_breadth
        >= PARTIAL_ALL_BREADTH
    ):
        return "partial_mapping_support"

    return "weak_fragmentary_signal"


def mapq_depth_class(
    all_depth,
    ratio
):

    if all_depth <= 0:
        return "no_depth_signal"

    if ratio >= 0.75:
        return "MAPQ_stable_ge0.75"

    if ratio >= 0.50:
        return "moderate_MAPQ_sensitivity_0.50_0.75"

    if ratio >= 0.25:
        return "high_MAPQ_sensitivity_0.25_0.50"

    return "very_high_MAPQ_sensitivity_lt0.25"


# ============================================================
# Load
# ============================================================

rows = read_tsv(
    INFILE
)

if len(rows) != 363:
    raise RuntimeError(
        f"Observations={len(rows)} expected=363"
    )


samples = sorted({
    r["sample"]
    for r in rows
})

loci = sorted({
    r["locus_id"]
    for r in rows
})

clusters = sorted({
    r["exact_cluster_id"]
    for r in rows
})


if len(samples) != 18:
    raise RuntimeError(
        f"Samples={len(samples)} expected=18"
    )

if len(loci) != 121:
    raise RuntimeError(
        f"Loci={len(loci)} expected=121"
    )


# ============================================================
# Locus × sample classification
# ============================================================

classified = []


for r in rows:

    all_breadth = fnum(
        r[
            "all_breadth_1x_pct"
        ]
    )

    q10_breadth = fnum(
        r[
            "q10_breadth_1x_pct"
        ]
    )

    all_depth = fnum(
        r[
            "all_mean_depth"
        ]
    )

    q10_depth = fnum(
        r[
            "q10_mean_depth"
        ]
    )

    depth_ratio = fnum(
        r[
            "q10_to_all_mean_depth_ratio"
        ]
    )


    sclass = support_class(
        all_breadth,
        q10_breadth
    )

    mclass = mapq_depth_class(
        all_depth,
        depth_ratio
    )


    rec = dict(r)

    rec.update({
        "final_mapping_support_class":
            sclass,

        "MAPQ_depth_sensitivity_class":
            mclass,

        "robust_all_breadth_threshold_pct":
            ROBUST_ALL_BREADTH,

        "robust_q10_breadth_threshold_pct":
            ROBUST_Q10_BREADTH,

        "partial_all_breadth_threshold_pct":
            PARTIAL_ALL_BREADTH,

        "robust_support":
            (
                "YES"
                if sclass
                == "robust_mapping_support"
                else "NO"
            ),

        "high_breadth_but_MAPQ_sensitive":
            (
                "YES"
                if sclass
                == "high_breadth_MAPQ_sensitive"
                else "NO"
            ),

        "biological_presence_call":
            "NOT_ASSIGNED",

        "interpretation_final":
            "mapping_support_for_assembled_16S_locus_not_organism_presence",
    })

    classified.append(
        rec
    )


fields = list(
    classified[0].keys()
)


write_tsv(
    OUT
    / "95D2B_locus_sample_final_support.tsv",
    classified,
    fields
)


# ============================================================
# Sample support summary
# ============================================================

sample_summary = []


for sample in samples:

    rr = [
        r
        for r in classified
        if r["sample"] == sample
    ]

    counts = Counter(
        r[
            "final_mapping_support_class"
        ]
        for r in rr
    )

    mapq_counts = Counter(
        r[
            "MAPQ_depth_sensitivity_class"
        ]
        for r in rr
    )


    sample_summary.append({
        "sample":
            sample,

        "coassembly":
            rr[0][
                "coassembly"
            ],

        "target_loci":
            len(rr),

        "robust_mapping_support":
            counts[
                "robust_mapping_support"
            ],

        "high_breadth_MAPQ_sensitive":
            counts[
                "high_breadth_MAPQ_sensitive"
            ],

        "partial_mapping_support":
            counts[
                "partial_mapping_support"
            ],

        "weak_fragmentary_signal":
            counts[
                "weak_fragmentary_signal"
            ],

        "no_mapped_signal":
            counts[
                "no_mapped_signal"
            ],

        "MAPQ_stable_ge0.75":
            mapq_counts[
                "MAPQ_stable_ge0.75"
            ],

        "very_high_MAPQ_sensitivity_lt0.25":
            mapq_counts[
                "very_high_MAPQ_sensitivity_lt0.25"
            ],

        "interpretation":
            "counts_are_mapping_support_categories_not_taxon_presence_counts",
    })


write_tsv(
    OUT
    / "95D2B_sample_support_summary.tsv",
    sample_summary,
    list(
        sample_summary[0].keys()
    )
)


# ============================================================
# Locus trajectories across its three native samples
# ============================================================

by_locus = defaultdict(list)

for r in classified:
    by_locus[
        r["locus_id"]
    ].append(r)


trajectory_rows = []


for locus in sorted(by_locus):

    rr = sorted(
        by_locus[locus],
        key=lambda x:
            x["sample"]
    )


    if len(rr) != 3:

        raise RuntimeError(
            f"{locus}: observations={len(rr)} expected=3"
        )


    robust_n = sum(
        r[
            "final_mapping_support_class"
        ]
        == "robust_mapping_support"
        for r in rr
    )

    high_sensitive_n = sum(
        r[
            "final_mapping_support_class"
        ]
        == "high_breadth_MAPQ_sensitive"
        for r in rr
    )

    any_signal_n = sum(
        r[
            "final_mapping_support_class"
        ]
        != "no_mapped_signal"
        for r in rr
    )


    if robust_n == 3:

        trajectory = (
            "robust_support_3_of_3"
        )

    elif robust_n == 2:

        trajectory = (
            "robust_support_2_of_3"
        )

    elif robust_n == 1:

        trajectory = (
            "robust_support_1_of_3"
        )

    elif any_signal_n > 0:

        trajectory = (
            "never_robust_but_signal_detected"
        )

    else:

        trajectory = (
            "no_signal_all_3"
        )


    trajectory_rows.append({
        "locus_id":
            locus,

        "coassembly":
            rr[0][
                "coassembly"
            ],

        "exact_cluster_id":
            rr[0][
                "exact_cluster_id"
            ],

        "taxon":
            rr[0][
                "deepest_informative_LCA_taxon"
            ],

        "resolution_tier":
            rr[0][
                "resolution_tier"
            ],

        "MAG_associated":
            rr[0][
                "MAG_associated"
            ],

        "final18_MAGs":
            rr[0][
                "final18_MAGs"
            ],

        "samples":
            ";".join(
                r["sample"]
                for r in rr
            ),

        "support_classes":
            ";".join(
                r[
                    "final_mapping_support_class"
                ]
                for r in rr
            ),

        "robust_samples_n":
            robust_n,

        "high_breadth_MAPQ_sensitive_samples_n":
            high_sensitive_n,

        "samples_with_any_signal_n":
            any_signal_n,

        "trajectory_class":
            trajectory,

        "interpretation":
            "within_native_coassembly_mapping_trajectory_not_prevalence",
    })


write_tsv(
    OUT
    / "95D2B_locus_trajectory_summary.tsv",
    trajectory_rows,
    list(
        trajectory_rows[0].keys()
    )
)


# ============================================================
# Category matrices
#
# Structural NA remains NA.
# ============================================================

lookup_locus = {
    (
        r["locus_id"],
        r["sample"]
    ): r[
        "final_mapping_support_class"
    ]
    for r in classified
}


# Validate unique cluster/sample observations before matrix
cluster_sample_groups = defaultdict(list)

for r in classified:

    cluster_sample_groups[
        (
            r[
                "exact_cluster_id"
            ],
            r["sample"]
        )
    ].append(r)


duplicates = {
    key: val
    for key, val
    in cluster_sample_groups.items()
    if len(val) > 1
}


if duplicates:

    raise RuntimeError(
        "Duplicate exact-cluster/sample observations "
        f"detected: {len(duplicates)}"
    )


lookup_cluster = {
    key:
        val[0][
            "final_mapping_support_class"
        ]
    for key, val
    in cluster_sample_groups.items()
}


locus_matrix = []

for locus in loci:

    rec = {
        "locus_id":
            locus
    }

    for sample in samples:

        rec[
            sample
        ] = lookup_locus.get(
            (
                locus,
                sample
            ),
            "NA"
        )

    locus_matrix.append(
        rec
    )


write_tsv(
    OUT
    / "95D2B_locus_support_class_matrix.tsv",
    locus_matrix,
    [
        "locus_id"
    ]
    + samples
)


cluster_matrix = []

for cluster in clusters:

    rec = {
        "exact_cluster_id":
            cluster
    }

    for sample in samples:

        rec[
            sample
        ] = lookup_cluster.get(
            (
                cluster,
                sample
            ),
            "NA"
        )

    cluster_matrix.append(
        rec
    )


write_tsv(
    OUT
    / "95D2B_cluster_support_class_matrix.tsv",
    cluster_matrix,
    [
        "exact_cluster_id"
    ]
    + samples
)


# ============================================================
# Taxon summary by sample
# ============================================================

taxon_groups = defaultdict(list)


for r in classified:

    taxon_groups[
        (
            r["sample"],
            r[
                "deepest_informative_LCA_taxon"
            ]
        )
    ].append(r)


taxon_rows = []


for (
    sample,
    taxon
), rr in taxon_groups.items():

    classes = Counter(
        r[
            "final_mapping_support_class"
        ]
        for r in rr
    )

    total_signal = sum(
        fnum(
            r[
                "all_mean_depth_per_million_mapped"
            ]
        )
        for r in rr
    )

    q10_signal = sum(
        fnum(
            r[
                "q10_mean_depth_per_million_mapped"
            ]
        )
        for r in rr
    )


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

        "robust_loci":
            classes[
                "robust_mapping_support"
            ],

        "high_breadth_MAPQ_sensitive_loci":
            classes[
                "high_breadth_MAPQ_sensitive"
            ],

        "partial_loci":
            classes[
                "partial_mapping_support"
            ],

        "weak_loci":
            classes[
                "weak_fragmentary_signal"
            ],

        "no_signal_loci":
            classes[
                "no_mapped_signal"
            ],

        "sum_normalized_signal_all":
            f"{total_signal:.10f}",

        "sum_normalized_signal_q10":
            f"{q10_signal:.10f}",

        "interpretation":
            "taxon_level_16S_locus_signal_not_cell_abundance",
    })


taxon_rows.sort(
    key=lambda r: (
        r["sample"],
        -int(
            r[
                "robust_loci"
            ]
        ),
        -float(
            r[
                "sum_normalized_signal_all"
            ]
        ),
        r["taxon"],
    )
)


write_tsv(
    OUT
    / "95D2B_taxon_support_by_sample.tsv",
    taxon_rows,
    list(
        taxon_rows[0].keys()
    )
)


# ============================================================
# Global summaries
# ============================================================

global_classes = Counter(
    r[
        "final_mapping_support_class"
    ]
    for r in classified
)

trajectory_counts = Counter(
    r[
        "trajectory_class"
    ]
    for r in trajectory_rows
)

mapq_classes = Counter(
    r[
        "MAPQ_depth_sensitivity_class"
    ]
    for r in classified
)


metrics = [
    (
        "samples",
        len(samples)
    ),
    (
        "primary_loci",
        len(loci)
    ),
    (
        "exact_clusters_evaluated",
        len(clusters)
    ),
    (
        "locus_sample_observations",
        len(classified)
    ),
    (
        "robust_mapping_support_observations",
        global_classes[
            "robust_mapping_support"
        ]
    ),
    (
        "high_breadth_MAPQ_sensitive_observations",
        global_classes[
            "high_breadth_MAPQ_sensitive"
        ]
    ),
    (
        "partial_mapping_support_observations",
        global_classes[
            "partial_mapping_support"
        ]
    ),
    (
        "weak_fragmentary_signal_observations",
        global_classes[
            "weak_fragmentary_signal"
        ]
    ),
    (
        "no_mapped_signal_observations",
        global_classes[
            "no_mapped_signal"
        ]
    ),
    (
        "MAPQ_stable_ge0.75_observations",
        mapq_classes[
            "MAPQ_stable_ge0.75"
        ]
    ),
    (
        "very_high_MAPQ_sensitivity_lt0.25_observations",
        mapq_classes[
            "very_high_MAPQ_sensitivity_lt0.25"
        ]
    ),
    (
        "loci_robust_3_of_3",
        trajectory_counts[
            "robust_support_3_of_3"
        ]
    ),
    (
        "loci_robust_2_of_3",
        trajectory_counts[
            "robust_support_2_of_3"
        ]
    ),
    (
        "loci_robust_1_of_3",
        trajectory_counts[
            "robust_support_1_of_3"
        ]
    ),
    (
        "loci_never_robust_but_signal",
        trajectory_counts[
            "never_robust_but_signal_detected"
        ]
    ),
    (
        "loci_no_signal_all_3",
        trajectory_counts[
            "no_signal_all_3"
        ]
    ),
]


with (
    OUT
    / "95D2B_global_summary.tsv"
).open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    for key, value in metrics:

        fh.write(
            f"{key}\t{value}\n"
        )


# ============================================================
# Methodological scope
# ============================================================

with (
    OUT
    / "95D2B_methodological_scope.tsv"
).open("w") as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        "analysis_scope\tshotgun_bacterial_16S_mapping_support\n"
    )

    fh.write(
        "robust_support_rule\t"
        "all_breadth_ge90pct_AND_MAPQ10_breadth_ge75pct\n"
    )

    fh.write(
        "high_breadth_MAPQ_sensitive_rule\t"
        "all_breadth_ge90pct_AND_MAPQ10_breadth_lt75pct\n"
    )

    fh.write(
        "partial_support_rule\t"
        "all_breadth_ge50pct_AND_lt90pct\n"
    )

    fh.write(
        "weak_signal_rule\t"
        "all_breadth_gt0pct_AND_lt50pct\n"
    )

    fh.write(
        "no_signal_rule\t"
        "all_breadth_eq0pct\n"
    )

    fh.write(
        "MAPQ_depth_sensitivity_independent_of_breadth_class\tYES\n"
    )

    fh.write(
        "robust_mapping_support_equals_organism_presence\tNO\n"
    )

    fh.write(
        "no_mapped_signal_equals_biological_absence\tNO\n"
    )

    fh.write(
        "normalized_depth_equals_cell_abundance\tNO\n"
    )

    fh.write(
        "16S_copy_number_correction_applied\tNO\n"
    )

    fh.write(
        "structural_NA_equals_zero\tNO\n"
    )

    fh.write(
        "cross_coassembly_prevalence_inference\tNO\n"
    )


# ============================================================
# Validation
# ============================================================

required = [
    "95D2B_locus_sample_final_support.tsv",
    "95D2B_sample_support_summary.tsv",
    "95D2B_locus_trajectory_summary.tsv",
    "95D2B_locus_support_class_matrix.tsv",
    "95D2B_cluster_support_class_matrix.tsv",
    "95D2B_taxon_support_by_sample.tsv",
    "95D2B_global_summary.tsv",
    "95D2B_methodological_scope.tsv",
]


for name in required:

    if not (
        OUT / name
    ).is_file():

        raise RuntimeError(
            f"Salida faltante: {name}"
        )


if len(classified) != 363:
    raise RuntimeError(
        f"Classified={len(classified)} expected=363"
    )

if len(trajectory_rows) != 121:
    raise RuntimeError(
        f"Trajectory loci={len(trajectory_rows)} expected=121"
    )


print(
    f"OBSERVATIONS={len(classified)}"
)

print(
    f"LOCI={len(trajectory_rows)}"
)

print(
    f"ROBUST={global_classes['robust_mapping_support']}"
)

print(
    f"HIGH_BREADTH_MAPQ_SENSITIVE="
    f"{global_classes['high_breadth_MAPQ_sensitive']}"
)

print(
    "BIOLOGICAL_PRESENCE_CALL=NO"
)

print(
    "95D2B=PASS"
)

