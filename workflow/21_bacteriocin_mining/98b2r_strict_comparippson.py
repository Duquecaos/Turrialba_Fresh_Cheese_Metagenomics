#!/usr/bin/env python3

import csv
import sys
from pathlib import Path


if len(sys.argv) != 4:
    raise SystemExit(
        "Usage: 98b2r_strict_comparippson.py "
        "<ROOT> <98B2_candidate_table.tsv> <OUT>"
    )


ROOT = Path(sys.argv[1])
INPUT = Path(sys.argv[2])
OUT = Path(sys.argv[3])

OUT.mkdir(
    parents=True,
    exist_ok=True
)


COMP = (
    ROOT
    / "29_comparippson_lactococcin"
)


RESULT_FILES = {
    "all_hits":
        COMP
        / "comparippson_lactococcin_all_hits.tsv",

    "best_hits":
        COMP
        / "comparippson_lactococcin_best_hits.tsv",

    "top10":
        COMP
        / "comparippson_lactococcin_top10.tsv",
}


# ============================================================
# Validate inputs
# ============================================================

if not INPUT.is_file():
    raise RuntimeError(
        f"Missing 98B2 input: {INPUT}"
    )


for label, path in RESULT_FILES.items():

    if not path.is_file():

        raise RuntimeError(
            f"Missing Comparippson result table "
            f"{label}: {path}"
        )


# ============================================================
# Read 98B2
# ============================================================

with INPUT.open(
    newline=""
) as fh:

    rows = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )


if not rows:
    raise RuntimeError(
        "98B2 candidate table is empty"
    )


# ============================================================
# Read genuine Comparippson OUTPUT tables as raw lines
#
# We deliberately DO NOT inspect:
#   lactococcin_full_orfs.faa
#   lactococcin_core_candidates.faa
#   lactococcin_candidate_sequences.tsv
#
# Those are query/input registries, not independent hits.
# ============================================================

result_lines = {}


for label, path in RESULT_FILES.items():

    with path.open(
        errors="replace"
    ) as fh:

        result_lines[label] = [
            (
                line_number,
                line.rstrip("\n")
            )
            for line_number, line
            in enumerate(
                fh,
                start=1
            )
            if line.strip()
        ]


# ============================================================
# Exact query-ID search in actual result tables
# ============================================================

corrected = []
detail = []


for row in rows:

    candidate = row[
        "candidate_id"
    ].strip()


    hits_by_table = {}

    total_hits = 0


    for label, lines in result_lines.items():

        hits = []

        for line_number, text in lines:

            # Exact literal candidate identifier.
            #
            # We intentionally do not use a generic token such
            # as "orf00015" because that could match unrelated
            # queries.
            if candidate in text:

                hits.append(
                    (
                        line_number,
                        text
                    )
                )


        hits_by_table[
            label
        ] = hits

        total_hits += len(
            hits
        )


        for line_number, text in hits[:20]:

            detail.append({
                "locus_id":
                    row[
                        "locus_id"
                    ],

                "candidate_id":
                    candidate,

                "comparippson_table":
                    label,

                "line_number":
                    line_number,

                "result_line":
                    text[:2000],
            })


    comp_result_support = int(
        total_hits > 0
    )


    # ========================================================
    # Strict antiSMASH candidate-level support
    # ========================================================

    antismash_candidate = int(
        row[
            "antismash_candidate_exact_translation"
        ] == "1"
        or
        row[
            "antismash_candidate_contained_translation"
        ] == "1"
    )


    # ========================================================
    # Strict BAGEL candidate-level support
    #
    # It must be from an actual same-MAG BAGEL run.
    # ========================================================

    bagel_candidate = int(
        row[
            "BAGEL_same_MAG_run_available"
        ] == "1"
        and
        (
            row[
                "BAGEL_candidate_text_hit"
            ] == "1"
            or
            row[
                "BAGEL_candidate_exact_sequence"
            ] == "1"
            or
            row[
                "BAGEL_candidate_contained_sequence"
            ] == "1"
        )
    )


    strict_families = (
        antismash_candidate
        + bagel_candidate
        + comp_result_support
    )


    if strict_families == 3:

        evidence = (
            "three_specialized_method_families"
        )

    elif strict_families == 2:

        evidence = (
            "two_specialized_method_families"
        )

    elif strict_families == 1:

        evidence = (
            "one_specialized_method_family"
        )

    else:

        evidence = (
            "no_candidate_level_specialized_support"
        )


    old_comp_sequence_signal = int(
        row[
            "Comparippson_candidate_text_hit"
        ] == "1"
        or
        row[
            "Comparippson_candidate_exact_sequence"
        ] == "1"
        or
        row[
            "Comparippson_candidate_contained_sequence"
        ] == "1"
    )


    old_comp_input_only = int(
        old_comp_sequence_signal == 1
        and
        comp_result_support == 0
    )


    new = dict(row)

    new.update({
        "strict_antismash_candidate_support":
            antismash_candidate,

        "strict_BAGEL_candidate_support":
            bagel_candidate,

        "Comparippson_all_hits_lines":
            len(
                hits_by_table[
                    "all_hits"
                ]
            ),

        "Comparippson_best_hits_lines":
            len(
                hits_by_table[
                    "best_hits"
                ]
            ),

        "Comparippson_top10_lines":
            len(
                hits_by_table[
                    "top10"
                ]
            ),

        "strict_Comparippson_result_support":
            comp_result_support,

        "previous_Comparippson_signal_input_only":
            old_comp_input_only,

        "strict_specialized_method_families":
            strict_families,

        "strict_evidence_pattern":
            evidence,
    })


    corrected.append(
        new
    )


# ============================================================
# Write candidate table
# ============================================================

candidate_out = (
    OUT
    / "98B2R_candidate_strict_evidence.tsv"
)


with candidate_out.open(
    "w",
    newline=""
) as fh:

    fields = list(
        corrected[0].keys()
    )

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(
        corrected
    )


# ============================================================
# Detail table
# ============================================================

detail_out = (
    OUT
    / "98B2R_comparippson_result_lines.tsv"
)


with detail_out.open(
    "w",
    newline=""
) as fh:

    fields = [
        "locus_id",
        "candidate_id",
        "comparippson_table",
        "line_number",
        "result_line",
    ]

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(
        detail
    )


# ============================================================
# Locus summary
# ============================================================

locus_ids = []

for row in corrected:

    if row[
        "locus_id"
    ] not in locus_ids:

        locus_ids.append(
            row[
                "locus_id"
            ]
        )


locus_rows = []


for locus_id in locus_ids:

    rr = [
        r
        for r in corrected
        if r[
            "locus_id"
        ] == locus_id
    ]


    max_methods = max(
        int(
            r[
                "strict_specialized_method_families"
            ]
        )
        for r in rr
    )


    locus_rows.append({
        "locus_id":
            locus_id,

        "MAG":
            rr[0][
                "MAG"
            ],

        "contig":
            rr[0][
                "contig"
            ],

        "candidate_count":
            len(rr),

        "candidates_with_antismash_support":
            sum(
                int(
                    r[
                        "strict_antismash_candidate_support"
                    ]
                )
                for r in rr
            ),

        "candidates_with_BAGEL_support":
            sum(
                int(
                    r[
                        "strict_BAGEL_candidate_support"
                    ]
                )
                for r in rr
            ),

        "candidates_with_Comparippson_result_support":
            sum(
                int(
                    r[
                        "strict_Comparippson_result_support"
                    ]
                )
                for r in rr
            ),

        "max_strict_method_families":
            max_methods,

        "robust_context_wide":
            rr[0].get(
                "robust_context_wide",
                ""
            ),

        "context_evidence_strict":
            rr[0].get(
                "context_evidence_strict",
                ""
            ),

        "left_edge_truncated":
            rr[0].get(
                "left_edge_truncated",
                ""
            ),

        "right_edge_truncated":
            rr[0].get(
                "right_edge_truncated",
                ""
            ),

        "contains_interrupted_candidate":
            rr[0].get(
                "contains_interrupted_candidate",
                ""
            ),

        "strict_locus_pattern":
            (
                "candidate_supported_by_3_methods"
                if max_methods == 3
                else
                "candidate_supported_by_2_methods"
                if max_methods == 2
                else
                "candidate_supported_by_1_method"
                if max_methods == 1
                else
                "no_candidate_level_specialized_support"
            ),
    })


locus_out = (
    OUT
    / "98B2R_locus_strict_summary.tsv"
)


with locus_out.open(
    "w",
    newline=""
) as fh:

    fields = list(
        locus_rows[0].keys()
    )

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(
        locus_rows
    )


# ============================================================
# Global summary
# ============================================================

global_out = (
    OUT
    / "98B2R_global_summary.tsv"
)


with global_out.open(
    "w",
    newline=""
) as fh:

    writer = csv.writer(
        fh,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow([
        "metric",
        "value"
    ])


    metrics = [
        (
            "candidate_hypotheses",
            len(
                corrected
            )
        ),
        (
            "strict_antismash_supported_candidates",
            sum(
                int(
                    r[
                        "strict_antismash_candidate_support"
                    ]
                )
                for r in corrected
            )
        ),
        (
            "strict_BAGEL_supported_candidates",
            sum(
                int(
                    r[
                        "strict_BAGEL_candidate_support"
                    ]
                )
                for r in corrected
            )
        ),
        (
            "strict_Comparippson_result_supported_candidates",
            sum(
                int(
                    r[
                        "strict_Comparippson_result_support"
                    ]
                )
                for r in corrected
            )
        ),
        (
            "previous_Comparippson_signals_that_were_input_only",
            sum(
                int(
                    r[
                        "previous_Comparippson_signal_input_only"
                    ]
                )
                for r in corrected
            )
        ),
        (
            "strict_candidates_3_methods",
            sum(
                int(
                    r[
                        "strict_specialized_method_families"
                    ]
                ) == 3
                for r in corrected
            )
        ),
        (
            "strict_candidates_2_methods",
            sum(
                int(
                    r[
                        "strict_specialized_method_families"
                    ]
                ) == 2
                for r in corrected
            )
        ),
        (
            "strict_candidates_1_method",
            sum(
                int(
                    r[
                        "strict_specialized_method_families"
                    ]
                ) == 1
                for r in corrected
            )
        ),
        (
            "strict_candidates_0_methods",
            sum(
                int(
                    r[
                        "strict_specialized_method_families"
                    ]
                ) == 0
                for r in corrected
            )
        ),
        (
            "loci_with_at_least_2_strict_methods",
            sum(
                int(
                    r[
                        "max_strict_method_families"
                    ]
                ) >= 2
                for r in locus_rows
            )
        ),
        (
            "new_prediction_performed",
            "NO"
        ),
        (
            "reads_remapped",
            "NO"
        ),
        (
            "amplicon_data_used",
            "NO"
        ),
        (
            "next_step",
            "98C_positive_evidence_expansion"
        ),
    ]


    writer.writerows(
        metrics
    )


print(
    f"CANDIDATES={len(corrected)}"
)

print(
    "STRICT_COMPARIPPSON_SUPPORTED="
    + str(
        sum(
            int(
                r[
                    "strict_Comparippson_result_support"
                ]
            )
            for r in corrected
        )
    )
)

print(
    "98B2R=PASS"
)
