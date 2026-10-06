#!/usr/bin/env python3

import csv
import sys
from pathlib import Path
from collections import defaultdict


if len(sys.argv) != 5:
    raise SystemExit(
        "Usage: 98c3c_bagel_candidate_level.py "
        "<BAGEL_summary.tsv> <positive_proteins.faa> "
        "<BAGEL_session> <OUT>"
    )


SUMMARY = Path(sys.argv[1])
PROTEINS = Path(sys.argv[2])
SESSION = Path(sys.argv[3])
OUT = Path(sys.argv[4])

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# Helpers
# ============================================================

def read_fasta(path):

    records = {}

    current = None
    seq = []

    def save():

        if current is not None:

            records[
                current
            ] = "".join(
                seq
            ).replace(
                "*",
                ""
            ).upper()

    with path.open(
        errors="replace"
    ) as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                save()

                current = (
                    line[1:]
                    .split()[0]
                )

                seq = []

            else:

                seq.append(
                    line
                )

    save()

    return records


def write_table(
    path,
    rows,
    fields=None
):

    with path.open(
        "w",
        newline=""
    ) as fh:

        if rows:

            fieldnames = (
                fields
                or list(
                    rows[0].keys()
                )
            )

        else:

            fieldnames = (
                fields
                or ["status"]
            )

        w = csv.DictWriter(
            fh,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore"
        )

        w.writeheader()

        for row in rows:
            w.writerow(
                row
            )


def classify_evidence_file(path):

    name = path.name.lower()

    if (
        "bacteriocin_hmmsearch"
        in name
        or
        (
            "bacteriocin" in name
            and "hmm" in name
        )
    ):
        return "bacteriocin_HMM"

    if (
        "bacteriocin" in name
        and "blast" in name
    ):
        return "bacteriocin_BLAST"

    if (
        "blast_bacteriocin"
        in name
    ):
        return "bacteriocin_BLAST"

    if (
        ".predict" in name
        or name.endswith(
            "predict"
        )
    ):
        return "BAGEL_predict"

    if "genetable" in name:
        return "GeneTable_context"

    if "aoi.table" in name:
        return "AOI_table_context"

    return "other_context"


STRICT_TYPES = {
    "bacteriocin_HMM",
    "bacteriocin_BLAST",
    "BAGEL_predict",
}


# ============================================================
# 1. Read validated BAGEL results
# ============================================================

with SUMMARY.open(
    newline=""
) as fh:

    rows = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )


positive = [
    r
    for r in rows
    if (
        r.get(
            "detection_pipeline_complete"
        ) == "1"
        and
        int(
            r.get(
                "BAGEL_AOI_count",
                "0"
            )
            or 0
        ) > 0
    )
]


if len(positive) != 3:

    raise RuntimeError(
        f"Expected 3 validated BAGEL AOI-positive contigs; "
        f"found {len(positive)}"
    )


# ============================================================
# 2. Target sequences
# ============================================================

protein_records = read_fasta(
    PROTEINS
)


targets = []


for row in positive:

    ids = [
        x
        for x in row[
            "primary_protein_ids"
        ].split(";")
        if x
    ]


    for pid in ids:

        if pid not in protein_records:

            raise RuntimeError(
                f"Target protein missing from FASTA: {pid}"
            )


        targets.append({
            "contig_id":
                row[
                    "contig_id"
                ],

            "target_protein_id":
                pid,

            "target_sequence":
                protein_records[
                    pid
                ],

            "target_aa_length":
                len(
                    protein_records[
                        pid
                    ]
                ),

            "direct_evidence_classes":
                row.get(
                    "direct_evidence_classes",
                    ""
                ),

            "best_query_candidates":
                row.get(
                    "best_query_candidates",
                    ""
                ),

            "GA_HMM_models":
                row.get(
                    "GA_HMM_models",
                    ""
                ),

            "BAGEL_classes":
                row.get(
                    "BAGEL_classes",
                    ""
                ),

            "BAGEL_querynames":
                row.get(
                    "BAGEL_querynames",
                    ""
                ),
        })


# ============================================================
# 3. Discover derived BAGEL protein FASTAs
#
# Exclude queryfolder: those are inputs, not BAGEL evidence.
# ============================================================

protein_fasta_files = []


for path in SESSION.rglob("*"):

    if not path.is_file():
        continue

    try:
        rel = path.relative_to(
            SESSION
        )
    except ValueError:
        continue


    if (
        rel.parts
        and rel.parts[0]
        == "queryfolder"
    ):
        continue


    if path.suffix.lower() in {
        ".faa",
        ".fa",
        ".fasta",
    }:

        protein_fasta_files.append(
            path
        )


sequence_matches = []


for path in sorted(
    protein_fasta_files
):

    try:

        records = read_fasta(
            path
        )

    except Exception:
        continue


    path_lower = path.name.lower()


    specialized_fasta = int(
        any(
            key in path_lower
            for key in [
                "bacteriocin",
                "candidate",
                "predict",
            ]
        )
    )


    for bagel_id, bagel_seq in records.items():

        if not bagel_seq:
            continue


        for target in targets:

            target_seq = target[
                "target_sequence"
            ]


            if bagel_seq == target_seq:

                relation = (
                    "exact_full_length"
                )

                target_cov = 100.0
                bagel_cov = 100.0


            elif (
                target_seq
                in bagel_seq
            ):

                relation = (
                    "target_sequence_contained_in_BAGEL_ORF"
                )

                target_cov = 100.0

                bagel_cov = (
                    100.0
                    * len(
                        target_seq
                    )
                    / len(
                        bagel_seq
                    )
                )


            elif (
                bagel_seq
                in target_seq
                and
                len(
                    bagel_seq
                )
                >= 0.80
                * len(
                    target_seq
                )
            ):

                relation = (
                    "BAGEL_ORF_contained_in_target_sequence"
                )

                target_cov = (
                    100.0
                    * len(
                        bagel_seq
                    )
                    / len(
                        target_seq
                    )
                )

                bagel_cov = 100.0


            else:

                continue


            sequence_matches.append({
                "contig_id":
                    target[
                        "contig_id"
                    ],

                "target_protein_id":
                    target[
                        "target_protein_id"
                    ],

                "target_aa_length":
                    target[
                        "target_aa_length"
                    ],

                "BAGEL_ORF_id":
                    bagel_id,

                "BAGEL_ORF_aa_length":
                    len(
                        bagel_seq
                    ),

                "sequence_relation":
                    relation,

                "target_coverage_pct":
                    round(
                        target_cov,
                        3
                    ),

                "BAGEL_ORF_coverage_pct":
                    round(
                        bagel_cov,
                        3
                    ),

                "source_fasta":
                    str(
                        path
                    ),

                "source_fasta_specialized_name":
                    specialized_fasta,
            })


# ============================================================
# 4. Text artifacts potentially containing BAGEL ORF IDs
# ============================================================

text_files = []


for path in SESSION.rglob("*"):

    if not path.is_file():
        continue

    try:
        rel = path.relative_to(
            SESSION
        )
    except ValueError:
        continue


    if (
        rel.parts
        and rel.parts[0]
        == "queryfolder"
    ):
        continue


    # Avoid binaries / huge files.
    try:
        size = path.stat().st_size
    except OSError:
        continue


    if size > 50_000_000:
        continue


    name = path.name.lower()


    if any(
        key in name
        for key in [
            "bacteriocin",
            "blast",
            "hmm",
            "predict",
            "genetable",
            "aoi.table",
        ]
    ):

        text_files.append(
            path
        )


# ============================================================
# 5. Search matched BAGEL ORF IDs in specialized evidence files
# ============================================================

evidence_lines = []


matches_by_target = defaultdict(
    list
)


for match in sequence_matches:

    matches_by_target[
        (
            match[
                "contig_id"
            ],
            match[
                "target_protein_id"
            ]
        )
    ].append(
        match
    )


for match in sequence_matches:

    bagel_id = match[
        "BAGEL_ORF_id"
    ]

    # BAGEL tables generally use first whitespace-delimited token.
    token = bagel_id.split()[0]


    for path in text_files:

        try:

            text = path.read_text(
                errors="replace"
            )

        except Exception:
            continue


        if token not in text:
            continue


        evidence_type = classify_evidence_file(
            path
        )


        for line_number, line in enumerate(
            text.splitlines(),
            start=1
        ):

            if token not in line:
                continue


            evidence_lines.append({
                "contig_id":
                    match[
                        "contig_id"
                    ],

                "target_protein_id":
                    match[
                        "target_protein_id"
                    ],

                "BAGEL_ORF_id":
                    bagel_id,

                "sequence_relation":
                    match[
                        "sequence_relation"
                    ],

                "evidence_type":
                    evidence_type,

                "strict_candidate_evidence":
                    int(
                        evidence_type
                        in STRICT_TYPES
                    ),

                "artifact":
                    str(
                        path
                    ),

                "line_number":
                    line_number,

                "evidence_line":
                    line[:4000],
            })


# ============================================================
# 6. Per-target review
# ============================================================

review_rows = []


for target in targets:

    key = (
        target[
            "contig_id"
        ],
        target[
            "target_protein_id"
        ]
    )


    seq_hits = matches_by_target.get(
        key,
        []
    )


    target_evidence = [
        e
        for e in evidence_lines
        if (
            e[
                "contig_id"
            ]
            ==
            target[
                "contig_id"
            ]
            and
            e[
                "target_protein_id"
            ]
            ==
            target[
                "target_protein_id"
            ]
        )
    ]


    exact_seq_hits = [
        x
        for x in seq_hits
        if x[
            "sequence_relation"
        ]
        == "exact_full_length"
    ]


    mapped_seq_hits = [
        x
        for x in seq_hits
        if x[
            "sequence_relation"
        ]
        in {
            "exact_full_length",
            "target_sequence_contained_in_BAGEL_ORF",
        }
    ]


    strict_evidence = [
        x
        for x in target_evidence
        if x[
            "strict_candidate_evidence"
        ] == 1
    ]


    strict_types = sorted({
        x[
            "evidence_type"
        ]
        for x in strict_evidence
    })


    specialized_fasta_hits = [
        x
        for x in seq_hits
        if x[
            "source_fasta_specialized_name"
        ] == 1
    ]


    if (
        mapped_seq_hits
        and
        (
            strict_evidence
            or specialized_fasta_hits
        )
    ):

        support = (
            "BAGEL_candidate_level_supported"
        )


    elif mapped_seq_hits:

        support = (
            "BAGEL_locus_supported_target_sequence_mapped_but_not_candidate_specific"
        )


    else:

        support = (
            "BAGEL_locus_supported_target_sequence_not_resolved_in_BAGEL_ORFs"
        )


    review_rows.append({
        "contig_id":
            target[
                "contig_id"
            ],

        "target_protein_id":
            target[
                "target_protein_id"
            ],

        "target_aa_length":
            target[
                "target_aa_length"
            ],

        "direct_evidence_classes":
            target[
                "direct_evidence_classes"
            ],

        "best_query_candidates":
            target[
                "best_query_candidates"
            ],

        "GA_HMM_models":
            target[
                "GA_HMM_models"
            ],

        "BAGEL_classes":
            target[
                "BAGEL_classes"
            ],

        "BAGEL_querynames":
            target[
                "BAGEL_querynames"
            ],

        "BAGEL_ORF_sequence_matches":
            len(
                seq_hits
            ),

        "exact_full_length_BAGEL_ORF_matches":
            len(
                exact_seq_hits
            ),

        "target_sequence_mapped_to_BAGEL_ORF":
            int(
                bool(
                    mapped_seq_hits
                )
            ),

        "specialized_FASTA_sequence_matches":
            len(
                specialized_fasta_hits
            ),

        "strict_candidate_evidence_lines":
            len(
                strict_evidence
            ),

        "strict_BAGEL_evidence_types":
            ";".join(
                strict_types
            ),

        "BAGEL_candidate_level_interpretation":
            support,
    })


# ============================================================
# 7. Artifact inventory
# ============================================================

artifact_rows = []


for path in sorted(
    text_files
):

    artifact_rows.append({
        "artifact":
            str(
                path
            ),

        "artifact_type":
            classify_evidence_file(
                path
            ),

        "bytes":
            path.stat().st_size,
    })


# ============================================================
# 8. Outputs
# ============================================================

write_table(
    OUT
    / "98C3C_target_candidate_review.tsv",
    review_rows
)


write_table(
    OUT
    / "98C3C_BAGEL_ORF_sequence_matches.tsv",
    sequence_matches,
    [
        "contig_id",
        "target_protein_id",
        "target_aa_length",
        "BAGEL_ORF_id",
        "BAGEL_ORF_aa_length",
        "sequence_relation",
        "target_coverage_pct",
        "BAGEL_ORF_coverage_pct",
        "source_fasta",
        "source_fasta_specialized_name",
    ]
)


write_table(
    OUT
    / "98C3C_candidate_evidence_lines.tsv",
    evidence_lines,
    [
        "contig_id",
        "target_protein_id",
        "BAGEL_ORF_id",
        "sequence_relation",
        "evidence_type",
        "strict_candidate_evidence",
        "artifact",
        "line_number",
        "evidence_line",
    ]
)


write_table(
    OUT
    / "98C3C_artifact_inventory.tsv",
    artifact_rows
)


# ============================================================
# 9. Global summary
# ============================================================

candidate_supported = sum(
    r[
        "BAGEL_candidate_level_interpretation"
    ]
    ==
    "BAGEL_candidate_level_supported"
    for r in review_rows
)


locus_only = sum(
    r[
        "BAGEL_candidate_level_interpretation"
    ]
    ==
    "BAGEL_locus_supported_target_sequence_mapped_but_not_candidate_specific"
    for r in review_rows
)


unresolved = sum(
    r[
        "BAGEL_candidate_level_interpretation"
    ]
    ==
    "BAGEL_locus_supported_target_sequence_not_resolved_in_BAGEL_ORFs"
    for r in review_rows
)


with (
    OUT
    / "98C3C_global_summary.tsv"
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
        (
            "BAGEL_AOI_positive_contigs_reviewed",
            3
        ),
        (
            "target_protein_hypotheses_reviewed",
            len(
                review_rows
            )
        ),
        (
            "targets_with_BAGEL_candidate_level_support",
            candidate_supported
        ),
        (
            "targets_with_locus_only_BAGEL_support",
            locus_only
        ),
        (
            "targets_with_unresolved_BAGEL_ORF_mapping",
            unresolved
        ),
        (
            "BAGEL_ORF_sequence_matches_total",
            len(
                sequence_matches
            )
        ),
        (
            "strict_candidate_evidence_lines_total",
            sum(
                int(
                    x[
                        "strict_candidate_evidence"
                    ]
                )
                for x in evidence_lines
            )
        ),
        (
            "BAGEL_counted_as_one_method_family",
            "YES"
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
            "antimicrobial_activity_assessed",
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
            "98C_final_synthesis"
        ),
    ]


    w.writerows(
        metrics
    )


print(
    f"TARGETS_REVIEWED={len(review_rows)}"
)

print(
    f"CANDIDATE_LEVEL_BAGEL_SUPPORTED={candidate_supported}"
)

print(
    f"LOCUS_ONLY_BAGEL_SUPPORTED={locus_only}"
)

print(
    "98C3C=PASS"
)
