#!/usr/bin/env python3

import csv
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path


if len(sys.argv) != 7:
    raise SystemExit(
        "Usage: 98d1c_build_final_matrices.py "
        "<BASE> <reference_metadata.tsv> "
        "<sample_manifest.tsv> "
        "<D0_context_groups.tsv> "
        "<98C_contig_synthesis.tsv> "
        "<OUT>"
    )


BASE = Path(sys.argv[1])
REFMETA = Path(sys.argv[2])
SAMPLE_MANIFEST = Path(sys.argv[3])
D0_GROUPS = Path(sys.argv[4])
C_CONTIG = Path(sys.argv[5])
OUT = Path(sys.argv[6])

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# Helpers
# ============================================================

def read_tsv(path):

    with path.open(
        newline=""
    ) as fh:

        return list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


def write_tsv(
    path,
    rows,
    fields=None
):

    with path.open(
        "w",
        newline=""
    ) as fh:

        if fields is None:

            if rows:
                fields = list(
                    rows[0].keys()
                )

            else:
                fields = [
                    "status"
                ]


        w = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore"
        )

        w.writeheader()

        for row in rows:
            w.writerow(
                row
            )


def fnum(x):

    if x in (
        "",
        "NA",
        None,
    ):
        return None

    return float(x)


def inum(x):

    return int(
        float(x)
    )


def unique_join(values):

    return ";".join(
        dict.fromkeys(
            x
            for x in values
            if x
        )
    )


# ============================================================
# Inputs
# ============================================================

refmeta = read_tsv(
    REFMETA
)

samples = read_tsv(
    SAMPLE_MANIFEST
)

groups = read_tsv(
    D0_GROUPS
)

ccontigs = read_tsv(
    C_CONTIG
)


if len(refmeta) != 30:
    raise RuntimeError(
        f"Expected 30 references; found {len(refmeta)}"
    )


if len(samples) != 18:
    raise RuntimeError(
        f"Expected 18 samples; found {len(samples)}"
    )


ref_by_id = {
    r[
        "reference_id"
    ]:
        r
    for r in refmeta
}


sample_by_name = {
    r[
        "sample"
    ]:
        r
    for r in samples
}


group_by_id = {
    r[
        "context_group"
    ]:
        r
    for r in groups
}


ccontig_by_id = {
    r[
        "contig_id"
    ]:
        r
    for r in ccontigs
}


expected_refs = set(
    ref_by_id
)

expected_samples = set(
    sample_by_name
)


# ============================================================
# 1. Load the 18 × 30 metrics
# ============================================================

long_rows = []


for sample in [
    r[
        "sample"
    ]
    for r in samples
]:

    path = (
        BASE
        / "samples"
        / sample
        / f"{sample}.98D1B_locus_metrics.tsv"
    )


    if not path.is_file():

        raise RuntimeError(
            f"Missing sample metrics: {path}"
        )


    rows = read_tsv(
        path
    )


    if len(rows) != 30:

        raise RuntimeError(
            f"{sample}: expected 30 rows; "
            f"found {len(rows)}"
        )


    refs = {
        r[
            "reference_id"
        ]
        for r in rows
    }


    if refs != expected_refs:

        raise RuntimeError(
            f"{sample}: reference set mismatch"
        )


    smeta = sample_by_name[
        sample
    ]


    for row in rows:

        long_rows.append({
            "sample":
                sample,

            "producer_code":
                smeta[
                    "producer_code"
                ],

            "biological_unit":
                smeta[
                    "biological_unit"
                ],

            "time_code":
                smeta[
                    "time_code"
                ],

            **row,
        })


if len(long_rows) != 540:

    raise RuntimeError(
        f"Expected 540 sample-reference combinations; "
        f"found {len(long_rows)}"
    )


write_tsv(
    OUT
    / "98D1C_long_18x30.tsv",
    long_rows
)


# ============================================================
# 2. Matrix writer
# ============================================================

sample_order = [
    r[
        "sample"
    ]
    for r in samples
]


reference_order = [
    r[
        "reference_id"
    ]
    for r in refmeta
]


lookup = {
    (
        r[
            "sample"
        ],
        r[
            "reference_id"
        ]
    ):
        r
    for r in long_rows
}


def write_matrix(
    filename,
    metric,
    transform=None
):

    fields = [
        "sample",
        "producer_code",
        "biological_unit",
        "time_code",
    ] + reference_order


    rows = []


    for sample in sample_order:

        smeta = sample_by_name[
            sample
        ]

        out = {
            "sample":
                sample,

            "producer_code":
                smeta[
                    "producer_code"
                ],

            "biological_unit":
                smeta[
                    "biological_unit"
                ],

            "time_code":
                smeta[
                    "time_code"
                ],
        }


        for ref in reference_order:

            value = lookup[
                (
                    sample,
                    ref
                )
            ][
                metric
            ]


            if transform is not None:
                value = transform(
                    value
                )


            out[
                ref
            ] = value


        rows.append(
            out
        )


    write_tsv(
        OUT
        / filename,
        rows,
        fields
    )


write_matrix(
    "98D1C_matrix_detection_class.tsv",
    "operational_detection_class"
)


write_matrix(
    "98D1C_matrix_robust_binary.tsv",
    "operational_detection_class",
    lambda x:
        1
        if x
        ==
        "robust_distributed_context_signal"
        else 0
)


write_matrix(
    "98D1C_matrix_breadth_ALL.tsv",
    "breadth_all_1x"
)


write_matrix(
    "98D1C_matrix_breadth_MAPQ10.tsv",
    "breadth_MAPQ10_1x"
)


write_matrix(
    "98D1C_matrix_mean_depth_ALL.tsv",
    "mean_depth_all"
)


write_matrix(
    "98D1C_matrix_mean_depth_MAPQ10.tsv",
    "mean_depth_MAPQ10"
)


write_matrix(
    "98D1C_matrix_mean_depth_MAPQ10_per_million_input_pairs.tsv",
    "mean_depth_MAPQ10_per_million_input_pairs"
)


write_matrix(
    "98D1C_matrix_mapped_records_per_million_input_pairs.tsv",
    "mapped_records_per_million_input_pairs"
)


# ============================================================
# 3. Evidence for redundancy groups
# ============================================================

def group_member_contigs(
    ref
):

    meta = ref_by_id[
        ref
    ]


    if (
        meta[
            "reference_class"
        ]
        ==
        "existing_ATTRLOC"
    ):

        return []


    gid = meta[
        "context_group"
    ]


    grow = group_by_id.get(
        gid
    )


    if grow is None:
        return [
            meta[
                "source_id"
            ]
        ]


    members = [
        x
        for x in grow.get(
            "additional_context_members",
            ""
        ).split(";")
        if x
    ]


    if not members:

        members = [
            meta[
                "source_id"
            ]
        ]


    return members


# ============================================================
# 4. Per-reference summary
# ============================================================

reference_summary = []


for ref in reference_order:

    rows = [
        lookup[
            (
                sample,
                ref
            )
        ]
        for sample in sample_order
    ]


    classes = Counter(
        r[
            "operational_detection_class"
        ]
        for r in rows
    )


    q10_breadths = [
        fnum(
            r[
                "breadth_MAPQ10_1x"
            ]
        )
        for r in rows
    ]


    q10_depths = [
        fnum(
            r[
                "mean_depth_MAPQ10"
            ]
        )
        for r in rows
    ]


    norm_q10_depths = [
        fnum(
            r[
                "mean_depth_MAPQ10_per_million_input_pairs"
            ]
        )
        for r in rows
    ]


    mapped_q10_retention = [
        fnum(
            r[
                "MAPQ10_mean_depth_retention"
            ]
        )
        for r in rows
        if r[
            "MAPQ10_mean_depth_retention"
        ]
        not in (
            "",
            "NA",
        )
    ]


    meta = ref_by_id[
        ref
    ]


    members = group_member_contigs(
        ref
    )


    group_evidence = [
        ccontig_by_id[
            cid
        ]
        for cid in members
        if cid in ccontig_by_id
    ]


    anti_any = (
        max(
            [
                inum(
                    r[
                        "antiSMASH_RiPP_context"
                    ]
                )
                for r in group_evidence
            ],
            default=0
        )
    )


    bagel_aoi_any = (
        max(
            [
                inum(
                    r[
                        "BAGEL_AOI_support"
                    ]
                )
                for r in group_evidence
            ],
            default=0
        )
    )


    bagel_candidate_any = (
        max(
            [
                inum(
                    r[
                        "BAGEL_candidate_level_support"
                    ]
                )
                for r in group_evidence
            ],
            default=0
        )
    )


    group_patterns = unique_join(
        r[
            "final_context_pattern"
        ]
        for r in group_evidence
    )


    protein_ids = unique_join(
        r[
            "primary_protein_ids"
        ]
        for r in group_evidence
    )


    robust_samples = classes[
        "robust_distributed_context_signal"
    ]


    signal_samples = (
        18
        -
        classes[
            "no_signal"
        ]
    )


    reference_summary.append({
        "reference_id":
            ref,

        "source_id":
            meta[
                "source_id"
            ],

        "reference_class":
            meta[
                "reference_class"
            ],

        "context_group":
            meta[
                "context_group"
            ],

        "group_member_contigs":
            ";".join(
                members
            ),

        "producer_source_context":
            meta[
                "producer"
            ],

        "reference_length":
            meta[
                "reference_length"
            ],

        "candidate_protein_ids":
            (
                protein_ids
                if protein_ids
                else
                meta[
                    "candidate_protein_ids"
                ]
            ),

        "antiSMASH_RiPP_context_any_group_member":
            anti_any,

        "BAGEL_AOI_any_group_member":
            bagel_aoi_any,

        "BAGEL_candidate_level_any_group_member":
            bagel_candidate_any,

        "group_prediction_patterns":
            group_patterns,

        "samples_total":
            18,

        "samples_with_any_signal":
            signal_samples,

        "samples_robust":
            robust_samples,

        "samples_high_breadth_MAPQ_sensitive":
            classes[
                "high_breadth_MAPQ_sensitive"
            ],

        "samples_partial":
            classes[
                "partial_context_signal"
            ],

        "samples_weak_or_localized":
            classes[
                "weak_or_localized_signal"
            ],

        "samples_no_signal":
            classes[
                "no_signal"
            ],

        "robust_in_at_least_one_sample":
            int(
                robust_samples > 0
            ),

        "max_breadth_MAPQ10":
            max(
                q10_breadths
            ),

        "median_breadth_MAPQ10":
            statistics.median(
                q10_breadths
            ),

        "max_mean_depth_MAPQ10":
            max(
                q10_depths
            ),

        "median_mean_depth_MAPQ10":
            statistics.median(
                q10_depths
            ),

        "max_normalized_mean_depth_MAPQ10":
            max(
                norm_q10_depths
            ),

        "median_normalized_mean_depth_MAPQ10":
            statistics.median(
                norm_q10_depths
            ),

        "median_MAPQ10_mean_depth_retention_when_signal":
            (
                statistics.median(
                    mapped_q10_retention
                )
                if mapped_q10_retention
                else "NA"
            ),

        "DNA_recruitment_equivalent_to_expression":
            "NO",

        "functional_bacteriocin_confirmed":
            "NO",
    })


write_tsv(
    OUT
    / "98D1C_reference_summary_30.tsv",
    reference_summary
)


# ============================================================
# 5. Subsets
# ============================================================

existing_summary = [
    r
    for r in reference_summary
    if r[
        "reference_class"
    ]
    ==
    "existing_ATTRLOC"
]


specialized_summary = [
    r
    for r in reference_summary
    if r[
        "reference_class"
    ]
    ==
    "additional_specialized"
]


exploratory_summary = [
    r
    for r in reference_summary
    if r[
        "reference_class"
    ]
    ==
    "additional_exploratory"
]


write_tsv(
    OUT
    / "98D1C_existing_ATTRLOC_summary.tsv",
    existing_summary
)


write_tsv(
    OUT
    / "98D1C_additional_specialized_summary.tsv",
    specialized_summary
)


write_tsv(
    OUT
    / "98D1C_additional_exploratory_summary.tsv",
    exploratory_summary
)


# ============================================================
# 6. Per-sample summary
# ============================================================

sample_summary = []


for sample in sample_order:

    rows = [
        lookup[
            (
                sample,
                ref
            )
        ]
        for ref in reference_order
    ]


    c = Counter(
        r[
            "operational_detection_class"
        ]
        for r in rows
    )


    input_pairs = {
        r[
            "input_read_pairs"
        ]
        for r in rows
    }


    align_rates = {
        r[
            "overall_alignment_rate_pct"
        ]
        for r in rows
    }


    if len(input_pairs) != 1:
        raise RuntimeError(
            f"{sample}: inconsistent input pair counts"
        )


    if len(align_rates) != 1:
        raise RuntimeError(
            f"{sample}: inconsistent alignment rate"
        )


    smeta = sample_by_name[
        sample
    ]


    sample_summary.append({
        "sample":
            sample,

        "producer_code":
            smeta[
                "producer_code"
            ],

        "biological_unit":
            smeta[
                "biological_unit"
            ],

        "time_code":
            smeta[
                "time_code"
            ],

        "input_read_pairs":
            next(
                iter(
                    input_pairs
                )
            ),

        "overall_alignment_rate_pct":
            next(
                iter(
                    align_rates
                )
            ),

        "robust":
            c[
                "robust_distributed_context_signal"
            ],

        "high_breadth_MAPQ_sensitive":
            c[
                "high_breadth_MAPQ_sensitive"
            ],

        "partial":
            c[
                "partial_context_signal"
            ],

        "weak_or_localized":
            c[
                "weak_or_localized_signal"
            ],

        "no_signal":
            c[
                "no_signal"
            ],
    })


write_tsv(
    OUT
    / "98D1C_sample_summary_18.tsv",
    sample_summary
)


# ============================================================
# 7. ATTRLOC002 vs ATTRLOC006
# ============================================================

comparison_rows = []


for sample in sample_order:

    r2 = lookup[
        (
            sample,
            "ATTRLOC002"
        )
    ]

    r6 = lookup[
        (
            sample,
            "ATTRLOC006"
        )
    ]

    smeta = sample_by_name[
        sample
    ]


    comparison_rows.append({
        "sample":
            sample,

        "producer_code":
            smeta[
                "producer_code"
            ],

        "biological_unit":
            smeta[
                "biological_unit"
            ],

        "time_code":
            smeta[
                "time_code"
            ],

        "ATTRLOC002_breadth_ALL":
            r2[
                "breadth_all_1x"
            ],

        "ATTRLOC002_breadth_MAPQ10":
            r2[
                "breadth_MAPQ10_1x"
            ],

        "ATTRLOC002_depth_ALL":
            r2[
                "mean_depth_all"
            ],

        "ATTRLOC002_depth_MAPQ10":
            r2[
                "mean_depth_MAPQ10"
            ],

        "ATTRLOC002_class":
            r2[
                "operational_detection_class"
            ],

        "ATTRLOC006_breadth_ALL":
            r6[
                "breadth_all_1x"
            ],

        "ATTRLOC006_breadth_MAPQ10":
            r6[
                "breadth_MAPQ10_1x"
            ],

        "ATTRLOC006_depth_ALL":
            r6[
                "mean_depth_all"
            ],

        "ATTRLOC006_depth_MAPQ10":
            r6[
                "mean_depth_MAPQ10"
            ],

        "ATTRLOC006_class":
            r6[
                "operational_detection_class"
            ],
    })


write_tsv(
    OUT
    / "98D1C_ATTRLOC002_vs_ATTRLOC006.tsv",
    comparison_rows
)


# ============================================================
# 8. Global counts
# ============================================================

global_classes = Counter(
    r[
        "operational_detection_class"
    ]
    for r in long_rows
)


if sum(
    global_classes.values()
) != 540:

    raise RuntimeError(
        "Global detection-class total != 540"
    )


def count_robust_refs(rows):

    return sum(
        int(
            r[
                "samples_robust"
            ]
        ) > 0
        for r in rows
    )


def count_signal_refs(rows):

    return sum(
        int(
            r[
                "samples_with_any_signal"
            ]
        ) > 0
        for r in rows
    )


mapq_sensitive_refs = [
    r[
        "reference_id"
    ]
    for r in reference_summary
    if int(
        r[
            "samples_high_breadth_MAPQ_sensitive"
        ]
    ) > 0
]


refsum = {
    r[
        "reference_id"
    ]:
        r
    for r in reference_summary
}


with (
    OUT
    / "98D1C_global_summary.tsv"
).open(
    "w",
    newline=""
) as fh:

    w = csv.writer(
        fh,
        delimiter="\t",
        lineterminator="\n"
    )

    w.writerow([
        "metric",
        "value"
    ])


    rows = [
        (
            "samples",
            18
        ),
        (
            "competitive_contexts",
            30
        ),
        (
            "sample_context_combinations",
            540
        ),
        (
            "existing_ATTRLOC_contexts",
            7
        ),
        (
            "additional_specialized_context_groups",
            20
        ),
        (
            "additional_exploratory_context_groups",
            3
        ),
        (
            "robust_combinations",
            global_classes[
                "robust_distributed_context_signal"
            ]
        ),
        (
            "high_breadth_MAPQ_sensitive_combinations",
            global_classes[
                "high_breadth_MAPQ_sensitive"
            ]
        ),
        (
            "partial_combinations",
            global_classes[
                "partial_context_signal"
            ]
        ),
        (
            "weak_or_localized_combinations",
            global_classes[
                "weak_or_localized_signal"
            ]
        ),
        (
            "no_signal_combinations",
            global_classes[
                "no_signal"
            ]
        ),
        (
            "references_robust_in_at_least_one_sample",
            count_robust_refs(
                reference_summary
            )
        ),
        (
            "existing_ATTRLOC_robust_in_at_least_one_sample",
            count_robust_refs(
                existing_summary
            )
        ),
        (
            "additional_specialized_robust_in_at_least_one_sample",
            count_robust_refs(
                specialized_summary
            )
        ),
        (
            "additional_specialized_with_any_signal",
            count_signal_refs(
                specialized_summary
            )
        ),
        (
            "additional_exploratory_robust_in_at_least_one_sample",
            count_robust_refs(
                exploratory_summary
            )
        ),
        (
            "references_with_high_breadth_MAPQ_sensitive_calls",
            ";".join(
                mapq_sensitive_refs
            )
        ),
        (
            "ATTRLOC002_robust_samples",
            refsum[
                "ATTRLOC002"
            ][
                "samples_robust"
            ]
        ),
        (
            "ATTRLOC002_high_breadth_MAPQ_sensitive_samples",
            refsum[
                "ATTRLOC002"
            ][
                "samples_high_breadth_MAPQ_sensitive"
            ]
        ),
        (
            "ATTRLOC006_robust_samples",
            refsum[
                "ATTRLOC006"
            ][
                "samples_robust"
            ]
        ),
        (
            "read_mapping_reexecuted_in_98D1C",
            "NO"
        ),
        (
            "normalized_depth_interpretation",
            "comparative_context_recruitment_signal_not_cell_abundance"
        ),
        (
            "robust_class_equivalent_to_biological_presence",
            "NO"
        ),
        (
            "DNA_signal_equivalent_to_expression",
            "NO"
        ),
        (
            "functional_bacteriocin_confirmed",
            "NO"
        ),
        (
            "amplicon_data_used",
            "NO"
        ),
        (
            "98D1_status",
            "CLOSED"
        ),
        (
            "next_step",
            "98E_final_bacteriocin_thesis_synthesis"
        ),
    ]


    w.writerows(
        rows
    )


# ============================================================
# 9. Scope
# ============================================================

with (
    OUT
    / "98D1C_methodological_scope.tsv"
).open(
    "w",
    newline=""
) as fh:

    w = csv.writer(
        fh,
        delimiter="\t",
        lineterminator="\n"
    )

    w.writerow([
        "field",
        "value"
    ])

    rows = [
        (
            "analysis",
            "final_aggregation_of_18x30_competitive_bacteriocin_context_mapping"
        ),
        (
            "mapping_reference",
            "7_existing_ATTRLOC_plus_20_additional_specialized_plus_3_exploratory"
        ),
        (
            "mapping",
            "Bowtie2_very_sensitive_primary_alignment_only"
        ),
        (
            "conservative_signal",
            "MAPQ_ge10"
        ),
        (
            "robust_operational_rule",
            "breadth_ALL_ge0.90_AND_breadth_MAPQ10_ge0.75"
        ),
        (
            "high_breadth_MAPQ_sensitive_rule",
            "breadth_ALL_ge0.90_AND_breadth_MAPQ10_lt0.75"
        ),
        (
            "partial_rule",
            "breadth_ALL_ge0.50_AND_lt0.90"
        ),
        (
            "weak_rule",
            "breadth_ALL_gt0_AND_lt0.50"
        ),
        (
            "no_signal_rule",
            "breadth_ALL_eq0"
        ),
        (
            "normalization",
            "mean_depth_per_million_input_read_pairs"
        ),
        (
            "normalization_use",
            "comparative_signal_only"
        ),
        (
            "producer_time_inference_performed",
            "NO"
        ),
        (
            "biological_week_mapping_assumed",
            "NO"
        ),
        (
            "presence_absence_claim_from_operational_classes",
            "NO"
        ),
        (
            "expression_assessed",
            "NO"
        ),
        (
            "antimicrobial_activity_assessed",
            "NO"
        ),
        (
            "amplicon_data_used",
            "NO"
        ),
    ]

    w.writerows(
        rows
    )


# ============================================================
# 10. Validation expectations
# ============================================================

expected = {
    "robust_distributed_context_signal":
        125,

    "high_breadth_MAPQ_sensitive":
        10,

    "partial_context_signal":
        66,

    "weak_or_localized_signal":
        141,

    "no_signal":
        198,
}


for key, value in expected.items():

    observed = global_classes[
        key
    ]

    if observed != value:

        raise RuntimeError(
            f"Unexpected {key}: "
            f"{observed} != {value}"
        )


if (
    int(
        refsum[
            "ATTRLOC002"
        ][
            "samples_robust"
        ]
    )
    != 0
):

    raise RuntimeError(
        "ATTRLOC002 robust count changed unexpectedly"
    )


if (
    int(
        refsum[
            "ATTRLOC002"
        ][
            "samples_high_breadth_MAPQ_sensitive"
        ]
    )
    != 10
):

    raise RuntimeError(
        "ATTRLOC002 MAPQ-sensitive count changed unexpectedly"
    )


if (
    int(
        refsum[
            "ATTRLOC006"
        ][
            "samples_robust"
        ]
    )
    != 12
):

    raise RuntimeError(
        "ATTRLOC006 robust count changed unexpectedly"
    )


print(
    "COMBINATIONS=540"
)

print(
    f"ROBUST={global_classes['robust_distributed_context_signal']}"
)

print(
    f"MAPQ_SENSITIVE={global_classes['high_breadth_MAPQ_sensitive']}"
)

print(
    f"PARTIAL={global_classes['partial_context_signal']}"
)

print(
    f"WEAK={global_classes['weak_or_localized_signal']}"
)

print(
    f"NO_SIGNAL={global_classes['no_signal']}"
)

print(
    "98D1C=PASS"
)
