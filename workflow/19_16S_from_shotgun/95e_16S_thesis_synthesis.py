#!/usr/bin/env python3

import csv
import html
import math
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(sys.argv[1])

D2B = (
    ROOT
    / "97_16S_shotgun"
    / "95D2B_final_support"
)

D2A = (
    ROOT
    / "97_16S_shotgun"
    / "95D2A_matrix_threshold_audit"
)

C2C = (
    ROOT
    / "97_16S_shotgun"
    / "95C2C_final_curated_bacterial_catalog"
)

C2 = (
    ROOT
    / "97_16S_shotgun"
    / "95C2_bacterial_QC"
)

LOCUS_SUPPORT = (
    D2B
    / "95D2B_locus_sample_final_support.tsv"
)

SAMPLE_SUPPORT = (
    D2B
    / "95D2B_sample_support_summary.tsv"
)

LOCUS_TRAJECTORY = (
    D2B
    / "95D2B_locus_trajectory_summary.tsv"
)

TAXON_SUPPORT = (
    D2B
    / "95D2B_taxon_support_by_sample.tsv"
)

SAMPLE_META = (
    D2A
    / "95D2A_sample_metadata.tsv"
)

FINAL_CATALOG = (
    C2C
    / "95C2C_final_bacterial_targets.tsv"
)

MAG_CONSISTENCY = (
    C2
    / "95C2_MAG_16S_consistency.tsv"
)

OUT = (
    ROOT
    / "97_16S_shotgun"
    / "95E_thesis_synthesis"
)

FIG = (
    OUT
    / "figures"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)

FIG.mkdir(
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


def fnum(value):

    try:
        return float(value)

    except Exception:
        return 0.0


def inum(value):

    try:
        return int(float(value))

    except Exception:
        return 0


def mean(values):

    values = list(values)

    if not values:
        return 0.0

    return sum(values) / len(values)


def median(values):

    values = list(values)

    if not values:
        return 0.0

    return statistics.median(
        values
    )


def pct(n, d):

    if d == 0:
        return 0.0

    return 100.0 * n / d


def unique_join(values):

    vals = sorted({
        str(x).strip()
        for x in values
        if str(x).strip()
    })

    return ";".join(vals)


def split_multi(value):

    return [
        x.strip()
        for x in str(
            value or ""
        ).split(";")
        if x.strip()
    ]


def svg_text(value):

    return html.escape(
        str(value)
    )


def hex_rgb(code):

    code = code.lstrip("#")

    return tuple(
        int(
            code[i:i+2],
            16
        )
        for i in (
            0,
            2,
            4
        )
    )


def rgb_hex(rgb):

    return (
        "#"
        + "".join(
            f"{max(0, min(255, int(round(x)))):02x}"
            for x in rgb
        )
    )


def interpolate_color(
    start,
    end,
    fraction
):

    fraction = max(
        0.0,
        min(
            1.0,
            fraction
        )
    )

    a = hex_rgb(
        start
    )

    b = hex_rgb(
        end
    )

    return rgb_hex(
        tuple(
            a[i]
            + (
                b[i]
                - a[i]
            )
            * fraction
            for i in range(3)
        )
    )


# ============================================================
# Load and validate
# ============================================================

locus_rows, _ = read_tsv(
    LOCUS_SUPPORT
)

sample_rows, _ = read_tsv(
    SAMPLE_SUPPORT
)

trajectory_rows, _ = read_tsv(
    LOCUS_TRAJECTORY
)

taxon_rows, _ = read_tsv(
    TAXON_SUPPORT
)

sample_meta, _ = read_tsv(
    SAMPLE_META
)

catalog_rows, _ = read_tsv(
    FINAL_CATALOG
)

mag_consistency, _ = read_tsv(
    MAG_CONSISTENCY
)


if len(locus_rows) != 363:
    raise RuntimeError(
        f"Locus-sample rows={len(locus_rows)} expected=363"
    )

if len(sample_rows) != 18:
    raise RuntimeError(
        f"Sample rows={len(sample_rows)} expected=18"
    )

if len(trajectory_rows) != 121:
    raise RuntimeError(
        f"Trajectory rows={len(trajectory_rows)} expected=121"
    )

if len(sample_meta) != 18:
    raise RuntimeError(
        f"Sample metadata rows={len(sample_meta)} expected=18"
    )

if len(catalog_rows) != 118:
    raise RuntimeError(
        f"Final exact clusters={len(catalog_rows)} expected=118"
    )


samples = sorted({
    r["sample"]
    for r in locus_rows
})

loci = sorted({
    r["locus_id"]
    for r in locus_rows
})

clusters = sorted({
    r["exact_cluster_id"]
    for r in locus_rows
})

coassemblies = sorted({
    r["coassembly"]
    for r in locus_rows
})


if len(samples) != 18:
    raise RuntimeError(
        f"Samples={len(samples)} expected=18"
    )

if len(loci) != 121:
    raise RuntimeError(
        f"Loci={len(loci)} expected=121"
    )

if len(clusters) != 118:
    raise RuntimeError(
        f"Clusters={len(clusters)} expected=118"
    )

if len(coassemblies) != 6:
    raise RuntimeError(
        f"Coassemblies={len(coassemblies)} expected=6"
    )


meta_by_sample = {
    r["sample"]: r
    for r in sample_meta
}


# ============================================================
# 1. MASTER TABLE: locus x sample
# ============================================================

master_rows = []


for row in sorted(
    locus_rows,
    key=lambda r: (
        r["sample"],
        r["locus_id"],
    )
):

    sample = row[
        "sample"
    ]

    meta = meta_by_sample[
        sample
    ]

    rec = {
        "sample":
            sample,

        "producer_code":
            meta[
                "producer_code"
            ],

        "biological_unit_code":
            meta[
                "biological_unit_code"
            ],

        "within_coassembly_time_code":
            meta[
                "within_coassembly_time_code"
            ],

        "coassembly":
            row[
                "coassembly"
            ],

        "locus_id":
            row[
                "locus_id"
            ],

        "exact_cluster_id":
            row[
                "exact_cluster_id"
            ],

        "contig":
            row[
                "contig"
            ],

        "start":
            row[
                "start"
            ],

        "end":
            row[
                "end"
            ],

        "strand":
            row[
                "strand"
            ],

        "length_bp":
            row[
                "length_bp"
            ],

        "taxon":
            row[
                "deepest_informative_LCA_taxon"
            ],

        "SILVA_LCA_terminal":
            row[
                "SILVA_LCA_terminal"
            ],

        "resolution_tier":
            row[
                "resolution_tier"
            ],

        "MAG_associated":
            row[
                "MAG_associated"
            ],

        "final18_MAGs":
            row[
                "final18_MAGs"
            ],

        "mapping_support_class":
            row[
                "final_mapping_support_class"
            ],

        "MAPQ_depth_sensitivity_class":
            row[
                "MAPQ_depth_sensitivity_class"
            ],

        "all_mean_depth":
            row[
                "all_mean_depth"
            ],

        "q10_mean_depth":
            row[
                "q10_mean_depth"
            ],

        "all_breadth_1x_pct":
            row[
                "all_breadth_1x_pct"
            ],

        "q10_breadth_1x_pct":
            row[
                "q10_breadth_1x_pct"
            ],

        "q10_to_all_mean_depth_ratio":
            row[
                "q10_to_all_mean_depth_ratio"
            ],

        "normalized_mean_depth_all":
            row[
                "all_mean_depth_per_million_mapped"
            ],

        "normalized_mean_depth_q10":
            row[
                "q10_mean_depth_per_million_mapped"
            ],

        "mapped_alignments_library":
            row[
                "mapped_alignments_library"
            ],

        "biological_presence_call":
            "NOT_ASSIGNED",

        "interpretation":
            "16S_locus_mapping_signal_not_cell_abundance",
    }

    master_rows.append(
        rec
    )


master_fields = list(
    master_rows[0].keys()
)


write_tsv(
    OUT
    / "95E_master_16S_locus_sample_table.tsv",
    master_rows,
    master_fields
)


# ============================================================
# 2. Thesis-ready sample summary
# ============================================================

locus_by_sample = defaultdict(list)

for r in locus_rows:

    locus_by_sample[
        r["sample"]
    ].append(r)


sample_support_by_sample = {
    r["sample"]: r
    for r in sample_rows
}


thesis_sample_rows = []


for sample in samples:

    rr = locus_by_sample[
        sample
    ]

    s = sample_support_by_sample[
        sample
    ]

    meta = meta_by_sample[
        sample
    ]

    robust = inum(
        s[
            "robust_mapping_support"
        ]
    )

    targets = inum(
        s[
            "target_loci"
        ]
    )


    taxa_with_robust = {
        r[
            "deepest_informative_LCA_taxon"
        ]
        for r in rr
        if r[
            "final_mapping_support_class"
        ]
        == "robust_mapping_support"
    }


    total_all = sum(
        fnum(
            r[
                "all_mean_depth_per_million_mapped"
            ]
        )
        for r in rr
    )

    total_q10 = sum(
        fnum(
            r[
                "q10_mean_depth_per_million_mapped"
            ]
        )
        for r in rr
    )


    thesis_sample_rows.append({
        "sample":
            sample,

        "producer_code":
            meta[
                "producer_code"
            ],

        "biological_unit_code":
            meta[
                "biological_unit_code"
            ],

        "within_coassembly_time_code":
            meta[
                "within_coassembly_time_code"
            ],

        "coassembly":
            s[
                "coassembly"
            ],

        "target_16S_loci":
            targets,

        "robust_mapping_support":
            robust,

        "robust_mapping_support_pct":
            f"{pct(robust, targets):.6f}",

        "high_breadth_MAPQ_sensitive":
            s[
                "high_breadth_MAPQ_sensitive"
            ],

        "partial_mapping_support":
            s[
                "partial_mapping_support"
            ],

        "weak_fragmentary_signal":
            s[
                "weak_fragmentary_signal"
            ],

        "no_mapped_signal":
            s[
                "no_mapped_signal"
            ],

        "taxa_with_at_least_one_robust_locus":
            len(
                taxa_with_robust
            ),

        "cumulative_normalized_16S_locus_signal_all":
            f"{total_all:.10f}",

        "cumulative_normalized_16S_locus_signal_q10":
            f"{total_q10:.10f}",

        "signal_interpretation":
            "cumulative_locus_signal_not_relative_or_cell_abundance",
    })


write_tsv(
    OUT
    / "95E_table_sample_summary.tsv",
    thesis_sample_rows,
    list(
        thesis_sample_rows[0].keys()
    )
)


# ============================================================
# 3. Top taxa per sample
# ============================================================

taxon_by_sample = defaultdict(list)

for r in taxon_rows:

    taxon_by_sample[
        r["sample"]
    ].append(r)


top_taxa_rows = []


for sample in samples:

    rr = sorted(
        taxon_by_sample[
            sample
        ],
        key=lambda r:
            float(
                r[
                    "sum_normalized_signal_all"
                ]
            ),
        reverse=True
    )


    for rank, r in enumerate(
        rr[:10],
        start=1
    ):

        rec = {
            "sample":
                sample,

            "producer_code":
                meta_by_sample[
                    sample
                ][
                    "producer_code"
                ],

            "biological_unit_code":
                meta_by_sample[
                    sample
                ][
                    "biological_unit_code"
                ],

            "within_coassembly_time_code":
                meta_by_sample[
                    sample
                ][
                    "within_coassembly_time_code"
                ],

            "rank_within_sample":
                rank,

            "taxon":
                r["taxon"],

            "n_16S_loci":
                r[
                    "n_16S_loci"
                ],

            "robust_loci":
                r[
                    "robust_loci"
                ],

            "MAPQ_sensitive_high_breadth_loci":
                r[
                    "high_breadth_MAPQ_sensitive_loci"
                ],

            "partial_loci":
                r[
                    "partial_loci"
                ],

            "sum_normalized_signal_all":
                r[
                    "sum_normalized_signal_all"
                ],

            "sum_normalized_signal_q10":
                r[
                    "sum_normalized_signal_q10"
                ],

            "interpretation":
                "ranked_16S_locus_signal_not_ranked_cell_abundance",
        }

        top_taxa_rows.append(
            rec
        )


write_tsv(
    OUT
    / "95E_table_top10_taxa_per_sample.tsv",
    top_taxa_rows,
    list(
        top_taxa_rows[0].keys()
    )
)


# ============================================================
# 4. Taxon trajectory by native coassembly
# ============================================================

taxon_grouped = defaultdict(list)


for r in taxon_rows:

    taxon_grouped[
        (
            r["coassembly"],
            r["taxon"],
        )
    ].append(r)


taxon_trajectory_rows = []


for (
    group,
    taxon
), rr in sorted(
    taxon_grouped.items()
):

    rr = sorted(
        rr,
        key=lambda r:
            r["sample"]
    )


    if len(rr) != 3:

        raise RuntimeError(
            f"{group}/{taxon}: "
            f"sample rows={len(rr)} expected=3"
        )


    signals_all = [
        fnum(
            r[
                "sum_normalized_signal_all"
            ]
        )
        for r in rr
    ]

    signals_q10 = [
        fnum(
            r[
                "sum_normalized_signal_q10"
            ]
        )
        for r in rr
    ]

    robust_samples = [
        r["sample"]
        for r in rr
        if inum(
            r[
                "robust_loci"
            ]
        ) > 0
    ]


    taxon_trajectory_rows.append({
        "coassembly":
            group,

        "taxon":
            taxon,

        "n_16S_loci":
            rr[0][
                "n_16S_loci"
            ],

        "sample_code_1":
            rr[0][
                "sample"
            ],

        "signal_all_code_1":
            f"{signals_all[0]:.10f}",

        "signal_q10_code_1":
            f"{signals_q10[0]:.10f}",

        "robust_loci_code_1":
            rr[0][
                "robust_loci"
            ],

        "sample_code_2":
            rr[1][
                "sample"
            ],

        "signal_all_code_2":
            f"{signals_all[1]:.10f}",

        "signal_q10_code_2":
            f"{signals_q10[1]:.10f}",

        "robust_loci_code_2":
            rr[1][
                "robust_loci"
            ],

        "sample_code_3":
            rr[2][
                "sample"
            ],

        "signal_all_code_3":
            f"{signals_all[2]:.10f}",

        "signal_q10_code_3":
            f"{signals_q10[2]:.10f}",

        "robust_loci_code_3":
            rr[2][
                "robust_loci"
            ],

        "median_normalized_signal_all":
            f"{median(signals_all):.10f}",

        "mean_normalized_signal_all":
            f"{mean(signals_all):.10f}",

        "samples_with_robust_locus":
            len(
                robust_samples
            ),

        "robust_sample_codes":
            ";".join(
                robust_samples
            ),

        "interpretation":
            "within_native_coassembly_16S_locus_trajectory_not_prevalence",
    })


write_tsv(
    OUT
    / "95E_table_taxon_trajectory_by_coassembly.tsv",
    taxon_trajectory_rows,
    list(
        taxon_trajectory_rows[0].keys()
    )
)


# ============================================================
# 5. Top taxa by coassembly
# ============================================================

top_by_group = []


for group in coassemblies:

    rr = [
        r
        for r in taxon_trajectory_rows
        if r[
            "coassembly"
        ] == group
    ]

    rr.sort(
        key=lambda r:
            float(
                r[
                    "median_normalized_signal_all"
                ]
            ),
        reverse=True
    )


    for rank, r in enumerate(
        rr[:10],
        start=1
    ):

        rec = dict(r)

        rec[
            "rank_by_median_signal"
        ] = rank

        top_by_group.append(
            rec
        )


fields_top_group = [
    "coassembly",
    "rank_by_median_signal",
] + [
    x
    for x in taxon_trajectory_rows[0].keys()
    if x != "coassembly"
]


write_tsv(
    OUT
    / "95E_table_top10_taxa_by_coassembly.tsv",
    top_by_group,
    fields_top_group
)


# ============================================================
# 6. MAG-associated locus table and summary
# ============================================================

mag_locus_rows = [
    r
    for r in master_rows
    if r[
        "MAG_associated"
    ] == "YES"
]


write_tsv(
    OUT
    / "95E_table_MAG_associated_16S_locus_sample.tsv",
    mag_locus_rows,
    master_fields
)


mag_consistency_by_mag = {
    r["MAG"]: r
    for r in mag_consistency
}


mag_groups = defaultdict(list)


for row in mag_locus_rows:

    for mag in split_multi(
        row[
            "final18_MAGs"
        ]
    ):

        mag_groups[
            mag
        ].append(row)


mag_summary_rows = []


for mag, rr in sorted(
    mag_groups.items()
):

    classes = Counter(
        r[
            "mapping_support_class"
        ]
        for r in rr
    )

    unique_loci = sorted({
        r["locus_id"]
        for r in rr
    })

    taxa = sorted({
        r["taxon"]
        for r in rr
    })

    samples_mag = sorted({
        r["sample"]
        for r in rr
    })

    qc = mag_consistency_by_mag.get(
        mag,
        {}
    )


    mag_summary_rows.append({
        "MAG":
            mag,

        "unique_16S_loci":
            len(
                unique_loci
            ),

        "locus_ids":
            ";".join(
                unique_loci
            ),

        "taxa":
            ";".join(
                taxa
            ),

        "native_samples_evaluated":
            len(
                samples_mag
            ),

        "samples":
            ";".join(
                samples_mag
            ),

        "robust_observations":
            classes[
                "robust_mapping_support"
            ],

        "high_breadth_MAPQ_sensitive_observations":
            classes[
                "high_breadth_MAPQ_sensitive"
            ],

        "partial_observations":
            classes[
                "partial_mapping_support"
            ],

        "weak_observations":
            classes[
                "weak_fragmentary_signal"
            ],

        "no_signal_observations":
            classes[
                "no_mapped_signal"
            ],

        "cross_locus_LCA":
            qc.get(
                "cross_locus_LCA",
                ""
            ),

        "16S_consistency_flag":
            qc.get(
                "16S_consistency_flag",
                ""
            ),

        "interpretation":
            "MAG_16S_mapping_support_and_taxonomic_QC_not_contamination_proof",
    })


write_tsv(
    OUT
    / "95E_table_MAG_16S_summary.tsv",
    mag_summary_rows,
    list(
        mag_summary_rows[0].keys()
    )
)


# ============================================================
# 7. Global summary
# ============================================================

classes = Counter(
    r[
        "mapping_support_class"
    ]
    for r in master_rows
)

mapq_classes = Counter(
    r[
        "MAPQ_depth_sensitivity_class"
    ]
    for r in master_rows
)

trajectories = Counter(
    r[
        "trajectory_class"
    ]
    for r in trajectory_rows
)


global_metrics = [
    (
        "final_bacterial_exact_16S_clusters",
        len(
            catalog_rows
        )
    ),
    (
        "final_primary_bacterial_16S_loci",
        len(
            loci
        )
    ),
    (
        "shotgun_samples",
        len(
            samples
        )
    ),
    (
        "native_coassemblies",
        len(
            coassemblies
        )
    ),
    (
        "locus_sample_observations",
        len(
            master_rows
        )
    ),
    (
        "robust_mapping_support_observations",
        classes[
            "robust_mapping_support"
        ]
    ),
    (
        "robust_mapping_support_pct",
        f"{pct(classes['robust_mapping_support'], len(master_rows)):.6f}"
    ),
    (
        "high_breadth_MAPQ_sensitive_observations",
        classes[
            "high_breadth_MAPQ_sensitive"
        ]
    ),
    (
        "partial_mapping_support_observations",
        classes[
            "partial_mapping_support"
        ]
    ),
    (
        "weak_fragmentary_signal_observations",
        classes[
            "weak_fragmentary_signal"
        ]
    ),
    (
        "no_mapped_signal_observations",
        classes[
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
        trajectories[
            "robust_support_3_of_3"
        ]
    ),
    (
        "loci_robust_2_of_3",
        trajectories[
            "robust_support_2_of_3"
        ]
    ),
    (
        "loci_robust_1_of_3",
        trajectories[
            "robust_support_1_of_3"
        ]
    ),
    (
        "loci_never_robust_but_signal",
        trajectories[
            "never_robust_but_signal_detected"
        ]
    ),
    (
        "loci_no_signal_all_3",
        trajectories[
            "no_signal_all_3"
        ]
    ),
    (
        "MAG_associated_unique_16S_loci",
        len({
            r[
                "locus_id"
            ]
            for r in mag_locus_rows
        })
    ),
    (
        "MAGs_with_16S",
        len(
            mag_summary_rows
        )
    ),
    (
        "amplicon_data_used",
        "NO"
    ),
]


with (
    OUT
    / "95E_global_summary.tsv"
).open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    for key, value in global_metrics:

        fh.write(
            f"{key}\t{value}\n"
        )


# ============================================================
# 8. Thesis results snapshot
# ============================================================

with (
    OUT
    / "95E_thesis_results_snapshot.txt"
).open("w") as fh:

    fh.write(
        "95E - SYNTHESIS OF SHOTGUN 16S RESULTS\n"
    )

    fh.write(
        "======================================\n\n"
    )

    fh.write(
        f"Final curated bacterial exact 16S sequence clusters: "
        f"{len(catalog_rows)}\n"
    )

    fh.write(
        f"Final primary bacterial 16S loci: "
        f"{len(loci)}\n"
    )

    fh.write(
        f"Shotgun samples: "
        f"{len(samples)}\n"
    )

    fh.write(
        f"Locus x native-sample observations: "
        f"{len(master_rows)}\n\n"
    )

    fh.write(
        "Mapping-support observations:\n"
    )

    fh.write(
        f"  robust_mapping_support = "
        f"{classes['robust_mapping_support']} "
        f"({pct(classes['robust_mapping_support'], len(master_rows)):.2f}%)\n"
    )

    fh.write(
        f"  high_breadth_MAPQ_sensitive = "
        f"{classes['high_breadth_MAPQ_sensitive']}\n"
    )

    fh.write(
        f"  partial_mapping_support = "
        f"{classes['partial_mapping_support']}\n"
    )

    fh.write(
        f"  weak_fragmentary_signal = "
        f"{classes['weak_fragmentary_signal']}\n"
    )

    fh.write(
        f"  no_mapped_signal = "
        f"{classes['no_mapped_signal']}\n\n"
    )

    fh.write(
        "Locus trajectories across the three native sample codes:\n"
    )

    fh.write(
        f"  robust in 3/3 = "
        f"{trajectories['robust_support_3_of_3']}\n"
    )

    fh.write(
        f"  robust in 2/3 = "
        f"{trajectories['robust_support_2_of_3']}\n"
    )

    fh.write(
        f"  robust in 1/3 = "
        f"{trajectories['robust_support_1_of_3']}\n"
    )

    fh.write(
        f"  never robust but signal detected = "
        f"{trajectories['never_robust_but_signal_detected']}\n"
    )

    fh.write(
        f"  no signal in all three = "
        f"{trajectories['no_signal_all_3']}\n\n"
    )

    fh.write(
        "Interpretation guardrails:\n"
    )

    fh.write(
        "- Mapping support refers to an assembled 16S locus and "
        "does not constitute a direct organism-presence call.\n"
    )

    fh.write(
        "- Normalized mean depth is a library-normalized locus signal, "
        "not cell abundance or conventional relative abundance.\n"
    )

    fh.write(
        "- No 16S rRNA gene copy-number correction was applied.\n"
    )

    fh.write(
        "- Structural NA values across native coassemblies are not zeros.\n"
    )

    fh.write(
        "- Amplicon 16S data were not used or compared in this analysis.\n"
    )


# ============================================================
# 9. FIGURE 1:
# Mapping support categories by sample
# ============================================================

category_specs = [
    (
        "robust_mapping_support",
        "Robust",
        "#2166ac"
    ),
    (
        "high_breadth_MAPQ_sensitive",
        "High breadth / MAPQ-sensitive",
        "#fdae61"
    ),
    (
        "partial_mapping_support",
        "Partial",
        "#67a9cf"
    ),
    (
        "weak_fragmentary_signal",
        "Weak",
        "#d9d9d9"
    ),
    (
        "no_mapped_signal",
        "No mapped signal",
        "#636363"
    ),
]


sample_summary_lookup = {
    r["sample"]: r
    for r in thesis_sample_rows
}


width = 1450
height = 760

left = 90
right = 40
top = 100
bottom = 150

plot_w = (
    width
    - left
    - right
)

plot_h = (
    height
    - top
    - bottom
)

max_targets = max(
    int(
        sample_summary_lookup[s][
            "target_16S_loci"
        ]
    )
    for s in samples
)

ymax = int(
    math.ceil(
        max_targets
        / 5.0
    )
    * 5
)

bar_slot = (
    plot_w
    / len(samples)
)

bar_width = (
    bar_slot
    * 0.66
)


svg = []

svg.append(
    f'<svg xmlns="http://www.w3.org/2000/svg" '
    f'width="{width}" height="{height}" '
    f'viewBox="0 0 {width} {height}">'
)

svg.append(
    '<rect width="100%" height="100%" fill="white"/>'
)

svg.append(
    '<text x="725" y="38" text-anchor="middle" '
    'font-size="24" font-family="Arial" font-weight="bold">'
    'Shotgun 16S locus mapping-support categories by sample'
    '</text>'
)

svg.append(
    '<text x="725" y="68" text-anchor="middle" '
    'font-size="15" font-family="Arial">'
    'Counts represent assembled 16S loci, not organism abundance'
    '</text>'
)


# Y grid
for yval in range(
    0,
    ymax + 1,
    5
):

    y = (
        top
        + plot_h
        - (
            yval
            / ymax
        )
        * plot_h
    )

    svg.append(
        f'<line x1="{left}" y1="{y:.2f}" '
        f'x2="{left + plot_w}" y2="{y:.2f}" '
        'stroke="#e5e5e5" stroke-width="1"/>'
    )

    svg.append(
        f'<text x="{left - 10}" y="{y + 5:.2f}" '
        'text-anchor="end" font-size="12" '
        'font-family="Arial">'
        f'{yval}'
        '</text>'
    )


for i, sample in enumerate(
    samples
):

    x = (
        left
        + i
        * bar_slot
        + (
            bar_slot
            - bar_width
        )
        / 2
    )

    current_y = (
        top
        + plot_h
    )

    row = sample_summary_lookup[
        sample
    ]

    for field, label, color in category_specs:

        count = int(
            row[field]
        )

        h = (
            count
            / ymax
            * plot_h
        )

        current_y -= h

        if h > 0:

            svg.append(
                f'<rect x="{x:.2f}" y="{current_y:.2f}" '
                f'width="{bar_width:.2f}" height="{h:.2f}" '
                f'fill="{color}" stroke="white" stroke-width="0.5"/>'
            )

    label_x = (
        x
        + bar_width
        / 2
    )

    label_y = (
        top
        + plot_h
        + 22
    )

    svg.append(
        f'<text x="{label_x:.2f}" y="{label_y:.2f}" '
        f'transform="rotate(45 {label_x:.2f} {label_y:.2f})" '
        'font-size="12" font-family="Arial" text-anchor="start">'
        f'{svg_text(sample)}'
        '</text>'
    )


# Axis
svg.append(
    f'<line x1="{left}" y1="{top}" '
    f'x2="{left}" y2="{top + plot_h}" '
    'stroke="black" stroke-width="1.3"/>'
)

svg.append(
    f'<line x1="{left}" y1="{top + plot_h}" '
    f'x2="{left + plot_w}" y2="{top + plot_h}" '
    'stroke="black" stroke-width="1.3"/>'
)

svg.append(
    f'<text x="25" y="{top + plot_h/2}" '
    f'transform="rotate(-90 25 {top + plot_h/2})" '
    'font-size="15" font-family="Arial" text-anchor="middle">'
    'Number of 16S loci'
    '</text>'
)


# Legend
legend_x = 160
legend_y = 690

for i, (
    field,
    label,
    color
) in enumerate(
    category_specs
):

    x = (
        legend_x
        + i
        * 245
    )

    svg.append(
        f'<rect x="{x}" y="{legend_y}" '
        f'width="18" height="18" fill="{color}"/>'
    )

    svg.append(
        f'<text x="{x + 25}" y="{legend_y + 14}" '
        'font-size="12" font-family="Arial">'
        f'{svg_text(label)}'
        '</text>'
    )


svg.append(
    '</svg>'
)


with (
    FIG
    / "95E_Fig1_mapping_support_by_sample.svg"
).open("w") as fh:

    fh.write(
        "\n".join(svg)
    )


# ============================================================
# 10. FIGURE 2:
# Heatmap of top 15 taxa by cumulative normalized locus signal
#
# Missing taxon/sample combinations remain NA and are gray.
# ============================================================

taxon_total_signal = defaultdict(float)


for r in taxon_rows:

    taxon_total_signal[
        r["taxon"]
    ] += float(
        r[
            "sum_normalized_signal_all"
        ]
    )


top_taxa = [
    taxon
    for taxon, total
    in sorted(
        taxon_total_signal.items(),
        key=lambda x:
            x[1],
        reverse=True
    )[:15]
]


taxon_sample_lookup = {
    (
        r["taxon"],
        r["sample"]
    ): r
    for r in taxon_rows
}


observed_values = []


for taxon in top_taxa:

    for sample in samples:

        r = taxon_sample_lookup.get(
            (
                taxon,
                sample
            )
        )

        if r is None:
            continue

        observed_values.append(
            math.log10(
                1.0
                + float(
                    r[
                        "sum_normalized_signal_all"
                    ]
                )
            )
        )


max_heat = (
    max(
        observed_values
    )
    if observed_values
    else 1.0
)


cell_w = 55
cell_h = 32

heat_left = 285
heat_top = 110

heat_width = (
    heat_left
    + len(samples)
    * cell_w
    + 90
)

heat_height = (
    heat_top
    + len(top_taxa)
    * cell_h
    + 150
)


svg2 = []

svg2.append(
    f'<svg xmlns="http://www.w3.org/2000/svg" '
    f'width="{heat_width}" height="{heat_height}" '
    f'viewBox="0 0 {heat_width} {heat_height}">'
)

svg2.append(
    '<rect width="100%" height="100%" fill="white"/>'
)

svg2.append(
    f'<text x="{heat_width/2:.1f}" y="35" '
    'text-anchor="middle" font-size="23" '
    'font-family="Arial" font-weight="bold">'
    'Top shotgun 16S taxon-associated locus signals'
    '</text>'
)

svg2.append(
    f'<text x="{heat_width/2:.1f}" y="64" '
    'text-anchor="middle" font-size="14" font-family="Arial">'
    'Color = log10(1 + normalized cumulative 16S locus signal); '
    'gray = not evaluated in the native coassembly'
    '</text>'
)


# Sample labels
for j, sample in enumerate(
    samples
):

    x = (
        heat_left
        + j
        * cell_w
        + cell_w / 2
    )

    y = (
        heat_top
        - 12
    )

    svg2.append(
        f'<text x="{x:.2f}" y="{y:.2f}" '
        f'transform="rotate(-55 {x:.2f} {y:.2f})" '
        'font-size="11" font-family="Arial" text-anchor="start">'
        f'{svg_text(sample)}'
        '</text>'
    )


for i, taxon in enumerate(
    top_taxa
):

    y = (
        heat_top
        + i
        * cell_h
    )

    svg2.append(
        f'<text x="{heat_left - 10}" '
        f'y="{y + cell_h*0.68:.2f}" '
        'text-anchor="end" font-size="12" font-family="Arial">'
        f'{svg_text(taxon)}'
        '</text>'
    )


    for j, sample in enumerate(
        samples
    ):

        x = (
            heat_left
            + j
            * cell_w
        )

        r = taxon_sample_lookup.get(
            (
                taxon,
                sample
            )
        )


        if r is None:

            color = "#d9d9d9"

        else:

            val = math.log10(
                1.0
                + float(
                    r[
                        "sum_normalized_signal_all"
                    ]
                )
            )

            fraction = (
                val
                / max_heat
                if max_heat > 0
                else 0
            )

            color = interpolate_color(
                "#f7fbff",
                "#08306b",
                fraction
            )


        svg2.append(
            f'<rect x="{x:.2f}" y="{y:.2f}" '
            f'width="{cell_w}" height="{cell_h}" '
            f'fill="{color}" stroke="white" stroke-width="1"/>'
        )


# Scale
legend_y2 = (
    heat_top
    + len(top_taxa)
    * cell_h
    + 55
)

legend_x2 = heat_left

segments = 100

legend_width = 350


for i in range(
    segments
):

    frac = (
        i
        / (
            segments - 1
        )
    )

    color = interpolate_color(
        "#f7fbff",
        "#08306b",
        frac
    )

    x = (
        legend_x2
        + frac
        * legend_width
    )

    svg2.append(
        f'<rect x="{x:.2f}" y="{legend_y2}" '
        f'width="{legend_width/segments + 1:.2f}" '
        'height="18" '
        f'fill="{color}" stroke="none"/>'
    )


svg2.append(
    f'<text x="{legend_x2}" y="{legend_y2 + 40}" '
    'font-size="11" font-family="Arial">'
    '0'
    '</text>'
)

svg2.append(
    f'<text x="{legend_x2 + legend_width}" '
    f'y="{legend_y2 + 40}" text-anchor="end" '
    'font-size="11" font-family="Arial">'
    f'{max_heat:.2f}'
    '</text>'
)

svg2.append(
    f'<text x="{legend_x2 + legend_width/2}" '
    f'y="{legend_y2 + 40}" text-anchor="middle" '
    'font-size="11" font-family="Arial">'
    'log10(1 + normalized locus signal)'
    '</text>'
)


# NA key
na_x = (
    legend_x2
    + legend_width
    + 90
)

svg2.append(
    f'<rect x="{na_x}" y="{legend_y2}" '
    'width="20" height="18" fill="#d9d9d9"/>'
)

svg2.append(
    f'<text x="{na_x + 28}" y="{legend_y2 + 14}" '
    'font-size="11" font-family="Arial">'
    'NA = not evaluated'
    '</text>'
)

svg2.append(
    '</svg>'
)


with (
    FIG
    / "95E_Fig2_top15_taxa_signal_heatmap.svg"
).open("w") as fh:

    fh.write(
        "\n".join(
            svg2
        )
    )


# ============================================================
# 11. Figure data tables
# ============================================================

fig1_data = thesis_sample_rows

write_tsv(
    FIG
    / "95E_Fig1_data.tsv",
    fig1_data,
    list(
        fig1_data[0].keys()
    )
)


fig2_data = []


for taxon in top_taxa:

    for sample in samples:

        r = taxon_sample_lookup.get(
            (
                taxon,
                sample
            )
        )

        if r is None:

            fig2_data.append({
                "taxon":
                    taxon,

                "sample":
                    sample,

                "normalized_locus_signal_all":
                    "NA",

                "normalized_locus_signal_q10":
                    "NA",

                "robust_loci":
                    "NA",

                "structural_status":
                    "NOT_EVALUATED_NATIVE_COASSEMBLY",
            })

        else:

            fig2_data.append({
                "taxon":
                    taxon,

                "sample":
                    sample,

                "normalized_locus_signal_all":
                    r[
                        "sum_normalized_signal_all"
                    ],

                "normalized_locus_signal_q10":
                    r[
                        "sum_normalized_signal_q10"
                    ],

                "robust_loci":
                    r[
                        "robust_loci"
                    ],

                "structural_status":
                    "EVALUATED",
            })


write_tsv(
    FIG
    / "95E_Fig2_data.tsv",
    fig2_data,
    list(
        fig2_data[0].keys()
    )
)


# ============================================================
# 12. Methodological scope
# ============================================================

with (
    OUT
    / "95E_methodological_scope.tsv"
).open("w") as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        "analysis\tfinal_shotgun_16S_thesis_synthesis\n"
    )

    fh.write(
        "amplicon_16S_used\tNO\n"
    )

    fh.write(
        "Barrnap_rerun\tNO\n"
    )

    fh.write(
        "SILVA_BLAST_rerun\tNO\n"
    )

    fh.write(
        "read_mapping_rerun\tNO\n"
    )

    fh.write(
        "quantification_source\t95D1R_native_coassembly_depth\n"
    )

    fh.write(
        "final_support_source\t95D2B\n"
    )

    fh.write(
        "robust_mapping_support_rule\t"
        "all_breadth_ge90pct_AND_MAPQ10_breadth_ge75pct\n"
    )

    fh.write(
        "normalized_signal\tmean_depth_per_million_mapped_alignments\n"
    )

    fh.write(
        "normalized_signal_equals_cell_abundance\tNO\n"
    )

    fh.write(
        "normalized_signal_equals_conventional_relative_abundance\tNO\n"
    )

    fh.write(
        "16S_copy_number_correction\tNO\n"
    )

    fh.write(
        "mapping_support_equals_organism_presence\tNO\n"
    )

    fh.write(
        "no_mapped_signal_equals_biological_absence\tNO\n"
    )

    fh.write(
        "cross_coassembly_structural_NA_equals_zero\tNO\n"
    )

    fh.write(
        "cross_coassembly_structural_NA_interpretation\tNOT_EVALUATED\n"
    )

    fh.write(
        "sample_suffix_time_labels\tcodes_preserved_without_biological_week_relabeling\n"
    )

    fh.write(
        "MAG_taxonomic_discordance_equals_contamination_proof\tNO\n"
    )

    fh.write(
        "SILVA_species_label_equals_formal_species_identification\tNO\n"
    )


# ============================================================
# 13. Output inventory
# ============================================================

inventory = []


for path in sorted(
    OUT.rglob("*")
):

    if not path.is_file():
        continue

    inventory.append({
        "relative_path":
            str(
                path.relative_to(
                    OUT
                )
            ),

        "size_bytes":
            path.stat().st_size,
    })


write_tsv(
    OUT
    / "95E_output_inventory.tsv",
    inventory,
    [
        "relative_path",
        "size_bytes",
    ]
)


# ============================================================
# Final validation
# ============================================================

required = [
    "95E_master_16S_locus_sample_table.tsv",
    "95E_table_sample_summary.tsv",
    "95E_table_top10_taxa_per_sample.tsv",
    "95E_table_taxon_trajectory_by_coassembly.tsv",
    "95E_table_top10_taxa_by_coassembly.tsv",
    "95E_table_MAG_associated_16S_locus_sample.tsv",
    "95E_table_MAG_16S_summary.tsv",
    "95E_global_summary.tsv",
    "95E_thesis_results_snapshot.txt",
    "95E_methodological_scope.tsv",
    "figures/95E_Fig1_mapping_support_by_sample.svg",
    "figures/95E_Fig2_top15_taxa_signal_heatmap.svg",
    "figures/95E_Fig1_data.tsv",
    "figures/95E_Fig2_data.tsv",
]


for name in required:

    path = (
        OUT
        / name
    )

    if (
        not path.is_file()
        or path.stat().st_size == 0
    ):

        raise RuntimeError(
            f"Salida faltante/vacía: {name}"
        )


if len(master_rows) != 363:

    raise RuntimeError(
        f"Master rows={len(master_rows)} expected=363"
    )


if len(thesis_sample_rows) != 18:

    raise RuntimeError(
        f"Sample summary rows={len(thesis_sample_rows)} expected=18"
    )


if len(mag_locus_rows) != 18:

    raise RuntimeError(
        f"MAG locus-sample rows={len(mag_locus_rows)} expected=18"
    )


if len(mag_summary_rows) != 3:

    raise RuntimeError(
        f"MAG summary rows={len(mag_summary_rows)} expected=3"
    )


print(
    f"MASTER_ROWS={len(master_rows)}"
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
    f"MAG_LOCUS_SAMPLE_ROWS={len(mag_locus_rows)}"
)

print(
    f"MAGS_WITH_16S={len(mag_summary_rows)}"
)

print(
    "AMPLICON_16S_USED=NO"
)

print(
    "95E_SYNTHESIS=PASS"
)

