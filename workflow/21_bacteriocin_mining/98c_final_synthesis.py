#!/usr/bin/env python3

import csv
import sys
from pathlib import Path
from collections import defaultdict, Counter


if len(sys.argv) != 10:
    raise SystemExit(
        "Usage: 98c_final_synthesis.py "
        "<98B2R_global.tsv> "
        "<98B2R_loci.tsv> "
        "<98C2_positive_context.tsv> "
        "<98C3A_contig_summary.tsv> "
        "<98C3A_regions.tsv> "
        "<98C3BR2R_BAGEL_summary.tsv> "
        "<98C3C_candidate_review.tsv> "
        "<98C3C_evidence_lines.tsv> "
        "<OUT>"
    )


B2_GLOBAL = Path(sys.argv[1])
B2_LOCI = Path(sys.argv[2])
C2 = Path(sys.argv[3])
ANTI = Path(sys.argv[4])
ANTI_REGIONS = Path(sys.argv[5])
BAGEL = Path(sys.argv[6])
BAGEL_CAND = Path(sys.argv[7])
BAGEL_LINES = Path(sys.argv[8])
OUT = Path(sys.argv[9])

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
    fieldnames=None
):

    with path.open(
        "w",
        newline=""
    ) as fh:

        if fieldnames is None:

            if rows:
                fieldnames = list(
                    rows[0].keys()
                )

            else:
                fieldnames = [
                    "status"
                ]


        writer = csv.DictWriter(
            fh,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore"
        )

        writer.writeheader()

        for row in rows:
            writer.writerow(
                row
            )


def metric_table(path):

    result = {}

    with path.open(
        newline=""
    ) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        for row in reader:

            result[
                row[
                    "metric"
                ]
            ] = row[
                "value"
            ]

    return result


def uniq(values):

    return list(
        dict.fromkeys(
            x
            for x in values
            if x not in (
                "",
                None
            )
        )
    )


def as_int(value):

    try:
        return int(
            float(
                value
            )
        )
    except Exception:
        return 0


def trueish(value):

    return str(
        value
    ).strip().lower() in {
        "1",
        "true",
        "yes",
        "y",
    }


def query_root(value):

    if not value:
        return ""

    return value.split(
        "|",
        1
    )[0]


# ============================================================
# Established candidate -> ATTRLOC mapping
#
# One sequence hypothesis can be represented in >1 ATTRLOC.
# In particular, orf00015 is represented in ATTRLOC002 and
# ATTRLOC006.
# ============================================================

ATTRLOC_MAP = {
    "L2_petauri_k141_179388_lactococcin_like":
        "ATTRLOC001",

    "M2_lactis_k141_64662_orf00015":
        "ATTRLOC002;ATTRLOC006",

    "L3_laudensis_k141_46847_orf00027":
        "ATTRLOC003",

    "L3_laudensis_k141_46847_orf00030":
        "ATTRLOC003",

    "L3_laudensis_k141_84984_sORF2":
        "ATTRLOC004",

    "M2_lactis_k141_25580_orf00005":
        "ATTRLOC005",

    "M2_lactis_k141_64662_orf00009":
        "ATTRLOC006",

    "M2_lactis_k141_78158_orf00001":
        "ATTRLOC007",
}


# ============================================================
# Inputs
# ============================================================

b2_global = metric_table(
    B2_GLOBAL
)

b2_loci = read_tsv(
    B2_LOCI
)

c2 = read_tsv(
    C2
)

anti = read_tsv(
    ANTI
)

anti_regions = read_tsv(
    ANTI_REGIONS
)

bagel = read_tsv(
    BAGEL
)

bagel_cand = read_tsv(
    BAGEL_CAND
)

bagel_lines = read_tsv(
    BAGEL_LINES
)


# ============================================================
# Primary proteins only
# ============================================================

primary = [
    row
    for row in c2
    if str(
        row.get(
            "primary_for_specialized_followup",
            ""
        )
    ) == "1"
]


if len(primary) != 44:

    raise RuntimeError(
        f"Expected 44 primary proteins; observed {len(primary)}"
    )


primary_contigs = sorted({
    row[
        "contig_id"
    ]
    for row in primary
})


if len(primary_contigs) != 40:

    raise RuntimeError(
        f"Expected 40 primary contigs; observed {len(primary_contigs)}"
    )


# ============================================================
# antiSMASH lookup
# ============================================================

anti_by_contig = {
    row[
        "contig_id"
    ]:
        row
    for row in anti
}


if len(anti_by_contig) != 40:

    raise RuntimeError(
        f"Expected 40 antiSMASH contig rows; "
        f"observed {len(anti_by_contig)}"
    )


edge_by_contig = defaultdict(
    list
)


for row in anti_regions:

    cid = row.get(
        "contig_id",
        ""
    )

    if not cid:
        continue

    edge_by_contig[
        cid
    ].append(
        trueish(
            row.get(
                "contig_edge",
                ""
            )
        )
    )


# ============================================================
# BAGEL lookup
# ============================================================

bagel_by_contig = {
    row[
        "contig_id"
    ]:
        row
    for row in bagel
}


if len(bagel_by_contig) != 7:

    raise RuntimeError(
        f"Expected 7 BAGEL-screened contigs; "
        f"observed {len(bagel_by_contig)}"
    )


bagel_candidate_by_target = {
    (
        row[
            "contig_id"
        ],
        row[
            "target_protein_id"
        ]
    ):
        row
    for row in bagel_cand
}


# ============================================================
# Deduplicate strict BAGEL evidence
#
# The same ORF can occur in AOI.faa and ORF-specific .faa.
# Evidence is therefore collapsed to candidate + method type,
# not counted by repeated line occurrence.
# ============================================================

bagel_types_by_target = defaultdict(
    set
)


for row in bagel_lines:

    if row.get(
        "strict_candidate_evidence",
        ""
    ) != "1":
        continue

    key = (
        row.get(
            "contig_id",
            ""
        ),
        row.get(
            "target_protein_id",
            ""
        )
    )

    etype = row.get(
        "evidence_type",
        ""
    )

    if etype:
        bagel_types_by_target[
            key
        ].add(
            etype
        )


# ============================================================
# Protein-level synthesis
# ============================================================

protein_rows = []


for row in primary:

    pid = row[
        "protein_id"
    ]

    cid = row[
        "contig_id"
    ]


    anti_row = anti_by_contig[
        cid
    ]


    anti_support = as_int(
        anti_row.get(
            "antismash_same_contig_support",
            "0"
        )
    )


    region_edges = edge_by_contig.get(
        cid,
        []
    )


    anti_edge = int(
        anti_support == 1
        and
        bool(
            region_edges
        )
        and
        all(
            region_edges
        )
    )


    bagel_row = bagel_by_contig.get(
        cid
    )


    if bagel_row is None:

        bagel_screen_status = (
            "NOT_RUN_antiSMASH_positive"
        )

        bagel_aoi = 0

        bagel_interpretation = ""

    else:

        bagel_aoi = int(
            as_int(
                bagel_row.get(
                    "BAGEL_AOI_count",
                    "0"
                )
            ) > 0
        )

        bagel_screen_status = (
            "AOI_detected"
            if bagel_aoi
            else
            "valid_screen_no_AOI"
        )

        bagel_interpretation = (
            bagel_row.get(
                "BAGEL_interpretation",
                ""
            )
        )


    cand_row = (
        bagel_candidate_by_target
        .get(
            (
                cid,
                pid
            )
        )
    )


    if cand_row is None:

        bagel_candidate_support = 0

        bagel_candidate_interpretation = ""

    else:

        bagel_candidate_interpretation = (
            cand_row.get(
                "BAGEL_candidate_level_interpretation",
                ""
            )
        )

        bagel_candidate_support = int(
            bagel_candidate_interpretation
            ==
            "BAGEL_candidate_level_supported"
        )


    bagel_types = sorted(
        bagel_types_by_target.get(
            (
                cid,
                pid
            ),
            set()
        )
    )


    specialized_context_support = int(
        anti_support == 1
        or
        bagel_aoi == 1
    )


    specialized_context_source = (
        "antiSMASH_RiPP_like"
        if anti_support
        else
        "BAGEL_AOI"
        if bagel_aoi
        else
        "NONE"
    )


    blast_count = as_int(
        row.get(
            "BLAST_hit_count",
            "0"
        )
    )


    blast_classes = (
        row.get(
            "BLAST_classes",
            ""
        )
    )


    candidate_query = (
        row.get(
            "best_query_candidate",
            ""
        )
    )


    candidate_root = query_root(
        candidate_query
    )


    source_attrloc = (
        ATTRLOC_MAP.get(
            candidate_root,
            ""
        )
    )


    # --------------------------------------------------------
    # Evidence-pattern labels, not claims of function.
    # --------------------------------------------------------

    if blast_count > 0:

        if (
            "exact_full_length"
            in blast_classes
            or
            "exact_query_contained"
            in blast_classes
        ):

            if specialized_context_support:

                evidence_pattern = (
                    "known_sequence_exact_with_specialized_context"
                )

            else:

                evidence_pattern = (
                    "known_sequence_exact_without_specialized_context"
                )


        elif "near_exact" in blast_classes:

            if specialized_context_support:

                evidence_pattern = (
                    "known_sequence_near_exact_with_specialized_context"
                )

            else:

                evidence_pattern = (
                    "known_sequence_near_exact_without_specialized_context"
                )


        else:

            if specialized_context_support:

                evidence_pattern = (
                    "known_sequence_homology_with_specialized_context"
                )

            else:

                evidence_pattern = (
                    "known_sequence_homology_without_specialized_context"
                )


    else:

        if specialized_context_support:

            evidence_pattern = (
                "GA_profile_with_specialized_context"
            )

        else:

            evidence_pattern = (
                "GA_profile_only_after_specialized_followup"
            )


    # Explicitly preserve the known interrupted-reference caution.
    interrupted_reference_homology = int(
        "not_functional_inference"
        in candidate_query
    )


    protein_rows.append({
        "protein_id":
            pid,

        "contig_id":
            cid,

        "aa_length":
            row.get(
                "aa_length",
                ""
            ),

        "partial":
            row.get(
                "partial",
                ""
            ),

        "producer":
            row.get(
                "contig_producers",
                ""
            ),

        "incomplete_context":
            row.get(
                "strongest_incomplete_context",
                ""
            ),

        "direct_evidence_class":
            row.get(
                "direct_evidence_class",
                ""
            ),

        "GA_HMM_models":
            row.get(
                "GA_HMM_models",
                ""
            ),

        "BLAST_classes":
            blast_classes,

        "best_query_candidate":
            candidate_query,

        "candidate_sequence_root":
            candidate_root,

        "linked_existing_ATTRLOC":
            source_attrloc,

        "best_BLAST_pident":
            row.get(
                "best_BLAST_pident",
                ""
            ),

        "best_BLAST_qcov_pct":
            row.get(
                "best_BLAST_qcov_pct",
                ""
            ),

        "antiSMASH_RiPP_context":
            anti_support,

        "antiSMASH_all_detected_regions_at_contig_edge":
            anti_edge,

        "BAGEL_screen_status":
            bagel_screen_status,

        "BAGEL_AOI_support":
            bagel_aoi,

        "BAGEL_locus_interpretation":
            bagel_interpretation,

        "BAGEL_candidate_level_support":
            bagel_candidate_support,

        "BAGEL_candidate_interpretation":
            bagel_candidate_interpretation,

        "BAGEL_candidate_evidence_types":
            ";".join(
                bagel_types
            ),

        "specialized_context_support":
            specialized_context_support,

        "specialized_context_source":
            specialized_context_source,

        "reference_candidate_marked_interrupted":
            interrupted_reference_homology,

        "evidence_pattern":
            evidence_pattern,

        "functional_bacteriocin_confirmed":
            "NO",
    })


# ============================================================
# Contig-level synthesis
# ============================================================

proteins_by_contig = defaultdict(
    list
)


for row in protein_rows:

    proteins_by_contig[
        row[
            "contig_id"
        ]
    ].append(
        row
    )


contig_rows = []


for cid in sorted(
    proteins_by_contig
):

    rows = proteins_by_contig[
        cid
    ]

    anti_support = max(
        int(
            r[
                "antiSMASH_RiPP_context"
            ]
        )
        for r in rows
    )

    bagel_support = max(
        int(
            r[
                "BAGEL_AOI_support"
            ]
        )
        for r in rows
    )

    candidate_support = max(
        int(
            r[
                "BAGEL_candidate_level_support"
            ]
        )
        for r in rows
    )


    direct_homology = [
        r
        for r in rows
        if r[
            "best_query_candidate"
        ]
    ]


    if anti_support:

        final_context_pattern = (
            "antiSMASH_RiPP_like_edge_region"
        )

    elif bagel_support:

        if candidate_support:

            final_context_pattern = (
                "BAGEL_rescue_with_candidate_level_support"
            )

        else:

            final_context_pattern = (
                "BAGEL_rescue_locus_level_only"
            )

    else:

        final_context_pattern = (
            "GA_profile_only_no_specialized_region"
        )


    contig_rows.append({
        "contig_id":
            cid,

        "primary_protein_count":
            len(
                rows
            ),

        "primary_protein_ids":
            ";".join(
                r[
                    "protein_id"
                ]
                for r in rows
            ),

        "producer":
            ";".join(
                uniq([
                    r[
                        "producer"
                    ]
                    for r in rows
                ])
            ),

        "incomplete_context":
            ";".join(
                uniq([
                    r[
                        "incomplete_context"
                    ]
                    for r in rows
                ])
            ),

        "direct_homology_protein_count":
            len(
                direct_homology
            ),

        "linked_existing_ATTRLOC":
            ";".join(
                uniq([
                    r[
                        "linked_existing_ATTRLOC"
                    ]
                    for r in direct_homology
                    if r[
                        "linked_existing_ATTRLOC"
                    ]
                ])
            ),

        "antiSMASH_RiPP_context":
            anti_support,

        "antiSMASH_edge_region":
            max(
                int(
                    r[
                        "antiSMASH_all_detected_regions_at_contig_edge"
                    ]
                )
                for r in rows
            ),

        "BAGEL_AOI_support":
            bagel_support,

        "BAGEL_candidate_level_support":
            candidate_support,

        "specialized_context_support":
            int(
                anti_support
                or
                bagel_support
            ),

        "final_context_pattern":
            final_context_pattern,

        "functional_bacteriocin_confirmed":
            "NO",
    })


# ============================================================
# Direct homology recurrence summary
# ============================================================

homology_rows = [
    r
    for r in protein_rows
    if r[
        "best_query_candidate"
    ]
]


by_candidate = defaultdict(
    list
)


for row in homology_rows:

    by_candidate[
        row[
            "candidate_sequence_root"
        ]
    ].append(
        row
    )


recurrence_rows = []


for candidate in sorted(
    by_candidate
):

    rows = by_candidate[
        candidate
    ]


    classes = Counter()

    for row in rows:

        bc = row[
            "BLAST_classes"
        ]

        if (
            "exact_full_length"
            in bc
            or
            "exact_query_contained"
            in bc
        ):

            classes[
                "exact"
            ] += 1

        elif "near_exact" in bc:

            classes[
                "near_exact"
            ] += 1

        else:

            classes[
                "strong_homology"
            ] += 1


    recurrence_rows.append({
        "candidate_sequence_root":
            candidate,

        "linked_existing_ATTRLOC":
            ATTRLOC_MAP.get(
                candidate,
                ""
            ),

        "recovered_proteins":
            len(
                rows
            ),

        "recovered_contigs":
            len({
                r[
                    "contig_id"
                ]
                for r in rows
            }),

        "exact_recurrences":
            classes[
                "exact"
            ],

        "near_exact_recurrences":
            classes[
                "near_exact"
            ],

        "strong_homology_recurrences":
            classes[
                "strong_homology"
            ],

        "producers":
            ";".join(
                uniq([
                    r[
                        "producer"
                    ]
                    for r in rows
                ])
            ),

        "antiSMASH_supported_proteins":
            sum(
                int(
                    r[
                        "antiSMASH_RiPP_context"
                    ]
                )
                for r in rows
            ),

        "BAGEL_AOI_supported_proteins":
            sum(
                int(
                    r[
                        "BAGEL_AOI_support"
                    ]
                )
                for r in rows
            ),

        "BAGEL_candidate_level_supported_proteins":
            sum(
                int(
                    r[
                        "BAGEL_candidate_level_support"
                    ]
                )
                for r in rows
            ),

        "all_recurrences_have_specialized_context":
            int(
                all(
                    int(
                        r[
                            "specialized_context_support"
                        ]
                    ) == 1
                    for r in rows
                )
            ),

        "novelty_claim":
            "NO",
    })


# ============================================================
# Four specialized-followup negatives
# ============================================================

negative_contigs = [
    row
    for row in contig_rows
    if int(
        row[
            "specialized_context_support"
        ]
    ) == 0
]


# ============================================================
# BAGEL candidate-level compact table
# ============================================================

bagel_candidate_rows = []


for row in bagel_cand:

    key = (
        row[
            "contig_id"
        ],
        row[
            "target_protein_id"
        ]
    )

    method_types = sorted(
        bagel_types_by_target.get(
            key,
            set()
        )
    )


    bagel_candidate_rows.append({
        "contig_id":
            row[
                "contig_id"
            ],

        "target_protein_id":
            row[
                "target_protein_id"
            ],

        "target_aa_length":
            row[
                "target_aa_length"
            ],

        "direct_evidence_classes":
            row[
                "direct_evidence_classes"
            ],

        "BAGEL_classes":
            row[
                "BAGEL_classes"
            ],

        "BAGEL_ORF_sequence_matches":
            row[
                "BAGEL_ORF_sequence_matches"
            ],

        "exact_full_length_BAGEL_ORF_matches":
            row[
                "exact_full_length_BAGEL_ORF_matches"
            ],

        "candidate_level_supported":
            int(
                row[
                    "BAGEL_candidate_level_interpretation"
                ]
                ==
                "BAGEL_candidate_level_supported"
            ),

        "unique_BAGEL_evidence_types":
            ";".join(
                method_types
            ),

        "BAGEL_method_family_count":
            1,

        "BAGEL_candidate_level_interpretation":
            row[
                "BAGEL_candidate_level_interpretation"
            ],
    })


# ============================================================
# Evidence-pattern counts
# ============================================================

protein_pattern_counts = Counter(
    row[
        "evidence_pattern"
    ]
    for row in protein_rows
)


contig_pattern_counts = Counter(
    row[
        "final_context_pattern"
    ]
    for row in contig_rows
)


# ============================================================
# Derived totals
# ============================================================

direct_homology_contigs = {
    row[
        "contig_id"
    ]
    for row in homology_rows
}


specialized_supported_contigs = [
    row
    for row in contig_rows
    if int(
        row[
            "specialized_context_support"
        ]
    ) == 1
]


anti_positive_contigs = [
    row
    for row in contig_rows
    if int(
        row[
            "antiSMASH_RiPP_context"
        ]
    ) == 1
]


anti_edge_positive = [
    row
    for row in contig_rows
    if (
        int(
            row[
                "antiSMASH_RiPP_context"
            ]
        ) == 1
        and
        int(
            row[
                "antiSMASH_edge_region"
            ]
        ) == 1
    )
]


bagel_screened = list(
    bagel_by_contig
)


bagel_aoi_positive = [
    row
    for row in bagel
    if as_int(
        row.get(
            "BAGEL_AOI_count",
            "0"
        )
    ) > 0
]


bagel_candidate_supported = [
    row
    for row in bagel_candidate_rows
    if int(
        row[
            "candidate_level_supported"
        ]
    ) == 1
]


attrloc_recovered = set()


for row in homology_rows:

    for x in (
        row[
            "linked_existing_ATTRLOC"
        ]
        .split(";")
    ):

        if x:
            attrloc_recovered.add(
                x
            )


# ============================================================
# Outputs
# ============================================================

write_tsv(
    OUT
    / "98C_FINAL_primary_protein_synthesis.tsv",
    protein_rows
)


write_tsv(
    OUT
    / "98C_FINAL_primary_contig_synthesis.tsv",
    contig_rows
)


write_tsv(
    OUT
    / "98C_FINAL_known_candidate_recurrence.tsv",
    recurrence_rows
)


write_tsv(
    OUT
    / "98C_FINAL_BAGEL_candidate_level.tsv",
    bagel_candidate_rows
)


write_tsv(
    OUT
    / "98C_FINAL_specialized_followup_negative_contigs.tsv",
    negative_contigs
)


# ============================================================
# Global summary
# ============================================================

with (
    OUT
    / "98C_FINAL_global_summary.tsv"
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


    metrics = [
        # Existing ATTRLOC framework
        (
            "existing_ATTRLOC_loci",
            len(
                b2_loci
            )
        ),
        (
            "existing_candidate_hypotheses",
            b2_global.get(
                "candidate_hypotheses",
                "NA"
            )
        ),
        (
            "existing_loci_with_at_least_2_strict_methods",
            b2_global.get(
                "loci_with_at_least_2_strict_methods",
                "NA"
            )
        ),
        (
            "existing_strict_Comparippson_result_supported_candidates",
            b2_global.get(
                "strict_Comparippson_result_supported_candidates",
                "NA"
            )
        ),

        # Incomplete-bin positive branch
        (
            "incomplete_primary_proteins",
            len(
                protein_rows
            )
        ),
        (
            "incomplete_primary_contigs",
            len(
                contig_rows
            )
        ),
        (
            "direct_homology_proteins_to_existing_candidates",
            len(
                homology_rows
            )
        ),
        (
            "direct_homology_contigs_to_existing_candidates",
            len(
                direct_homology_contigs
            )
        ),
        (
            "existing_candidate_sequence_hypotheses_recovered",
            len(
                by_candidate
            )
        ),
        (
            "existing_ATTRLOC_labels_linked_by_recurrence",
            len(
                attrloc_recovered
            )
        ),
        (
            "antiSMASH_RiPP_positive_contigs",
            len(
                anti_positive_contigs
            )
        ),
        (
            "antiSMASH_RiPP_positive_contigs_at_edge",
            len(
                anti_edge_positive
            )
        ),
        (
            "antiSMASH_negative_contigs_screened_by_BAGEL",
            len(
                bagel_screened
            )
        ),
        (
            "BAGEL_AOI_positive_contigs",
            len(
                bagel_aoi_positive
            )
        ),
        (
            "BAGEL_candidate_level_supported_targets",
            len(
                bagel_candidate_supported
            )
        ),
        (
            "specialized_context_supported_contigs_total",
            len(
                specialized_supported_contigs
            )
        ),
        (
            "specialized_context_supported_fraction",
            (
                f"{len(specialized_supported_contigs) / len(contig_rows):.6f}"
            )
        ),
        (
            "specialized_followup_negative_contigs",
            len(
                negative_contigs
            )
        ),
        (
            "direct_homology_proteins_with_specialized_context",
            sum(
                int(
                    row[
                        "specialized_context_support"
                    ]
                )
                for row in homology_rows
            )
        ),
        (
            "all_direct_homology_proteins_have_specialized_context",
            (
                "YES"
                if all(
                    int(
                        row[
                            "specialized_context_support"
                        ]
                    ) == 1
                    for row in homology_rows
                )
                else "NO"
            )
        ),
    ]


    for key, value in sorted(
        protein_pattern_counts.items()
    ):

        metrics.append(
            (
                f"protein_evidence_pattern__{key}",
                value
            )
        )


    for key, value in sorted(
        contig_pattern_counts.items()
    ):

        metrics.append(
            (
                f"contig_pattern__{key}",
                value
            )
        )


    metrics += [
        (
            "Comparippson_counted_as_support",
            "NO"
        ),
        (
            "antiSMASH_and_98C1_HMM_counted_as_fully_independent",
            "NO"
        ),
        (
            "BAGEL_subcomponents_counted_as_multiple_method_families",
            "NO"
        ),
        (
            "additional_incomplete_context_equivalent_to_novel_bacteriocin",
            "NO"
        ),
        (
            "negative_in_incomplete_bin_equivalent_to_biological_absence",
            "NO"
        ),
        (
            "functional_bacteriocin_confirmed",
            "NO"
        ),
        (
            "expression_assessed",
            "NO"
        ),
        (
            "antimicrobial_activity_assessed_by_metagenomics",
            "NO"
        ),
        (
            "reads_remapped_in_98C",
            "NO"
        ),
        (
            "amplicon_data_used",
            "NO"
        ),
        (
            "98C_status",
            "CLOSED"
        ),
        (
            "next_step",
            "98D_deduplicate_additional_candidate_loci_before_optional_read_mapping"
        ),
    ]


    w.writerows(
        metrics
    )


# ============================================================
# Methodological scope
# ============================================================

with (
    OUT
    / "98C_FINAL_methodological_scope.tsv"
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
            "final_synthesis_of_positive_evidence_expansion_in_incomplete_bins"
        ),
        (
            "primary_screen",
            "curated_GA_HMM_or_direct_sequence_homology"
        ),
        (
            "antiSMASH_role",
            "RiPP_region_context_support"
        ),
        (
            "antiSMASH_edge_regions",
            "fragmentary_context_requires_caution"
        ),
        (
            "BAGEL_role",
            "complementary_specialized_rescue_of_antiSMASH_negative_primary_contigs"
        ),
        (
            "BAGEL_candidate_level_support",
            "requires_target_sequence_mapping_to_BAGEL_ORF_plus_specialized_BAGEL_artifact"
        ),
        (
            "BAGEL_evidence_lines",
            "deduplicated_to_method_types_for_interpretation"
        ),
        (
            "direct_homology",
            "recurrence_or_homology_to_existing_candidate_sequence_not_novelty"
        ),
        (
            "incomplete_bin_positive_evidence",
            "sequence_or_context_evidence_can_be_retained"
        ),
        (
            "incomplete_bin_negative_evidence",
            "not_biological_absence"
        ),
        (
            "bin_contamination",
            "weakens_bin_or_taxon_attribution_not_sequence_level_positive_evidence"
        ),
        (
            "DNA_evidence",
            "does_not_demonstrate_expression_production_or_antimicrobial_activity"
        ),
        (
            "amplicon_data_used",
            "NO"
        ),
    ]

    w.writerows(
        rows
    )


print(
    f"PRIMARY_PROTEINS={len(protein_rows)}"
)

print(
    f"PRIMARY_CONTIGS={len(contig_rows)}"
)

print(
    f"SPECIALIZED_SUPPORTED_CONTIGS={len(specialized_supported_contigs)}"
)

print(
    f"DIRECT_HOMOLOGY_PROTEINS={len(homology_rows)}"
)

print(
    f"SPECIALIZED_NEGATIVE_CONTIGS={len(negative_contigs)}"
)

print(
    "98C_FINAL=PASS"
)
