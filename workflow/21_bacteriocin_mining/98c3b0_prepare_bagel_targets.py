#!/usr/bin/env python3

import csv
import sys
from pathlib import Path
from collections import defaultdict


if len(sys.argv) != 5:
    raise SystemExit(
        "Usage: 98c3b0_prepare_bagel_targets.py "
        "<antismash_summary.tsv> "
        "<positive_context.tsv> "
        "<source_contigs.fna> "
        "<OUT>"
    )


ANTI = Path(sys.argv[1])
CONTEXT = Path(sys.argv[2])
FASTA = Path(sys.argv[3])
OUT = Path(sys.argv[4])

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# antiSMASH contig summary
# ============================================================

with ANTI.open(
    newline=""
) as fh:

    anti_rows = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )


negative = [
    r
    for r in anti_rows
    if r[
        "antismash_same_contig_support"
    ] == "0"
]


if len(anti_rows) != 40:

    raise RuntimeError(
        f"Expected 40 antiSMASH target contigs; "
        f"observed {len(anti_rows)}"
    )


if len(negative) != 7:

    raise RuntimeError(
        f"Expected 7 antiSMASH-negative contigs; "
        f"observed {len(negative)}"
    )


negative_ids = {
    r[
        "contig_id"
    ]
    for r in negative
}


# ============================================================
# Protein/context details
# ============================================================

context_by_contig = defaultdict(
    list
)


with CONTEXT.open(
    newline=""
) as fh:

    for row in csv.DictReader(
        fh,
        delimiter="\t"
    ):

        cid = row[
            "contig_id"
        ]

        if cid in negative_ids:

            context_by_contig[
                cid
            ].append(
                row
            )


# ============================================================
# Target manifest
# ============================================================

manifest = []


for row in negative:

    cid = row[
        "contig_id"
    ]

    pp = context_by_contig.get(
        cid,
        []
    )


    evidence = sorted({
        x.get(
            "direct_evidence_class",
            ""
        )
        for x in pp
        if x.get(
            "direct_evidence_class",
            ""
        )
    })


    proteins = [
        x.get(
            "protein_id",
            ""
        )
        for x in pp
        if x.get(
            "protein_id",
            ""
        )
    ]


    candidates = sorted({
        x.get(
            "best_query_candidate",
            ""
        )
        for x in pp
        if x.get(
            "best_query_candidate",
            ""
        )
    })


    hmm = sorted({
        x.get(
            "GA_HMM_models",
            ""
        )
        for x in pp
        if x.get(
            "GA_HMM_models",
            ""
        )
    })


    context_flags = sorted({
        x.get(
            "auxiliary_context_flag",
            ""
        )
        for x in pp
        if x.get(
            "auxiliary_context_flag",
            ""
        )
    })


    if any(
        e.startswith("A_")
        or e.startswith("B_")
        for e in evidence
    ):

        priority = (
            "HIGH_direct_sequence_evidence"
        )

    else:

        priority = (
            "STANDARD_GA_HMM_screen"
        )


    manifest.append({
        "contig_id":
            cid,

        "length_bp":
            row.get(
                "length_bp",
                ""
            ),

        "primary_protein_ids":
            ";".join(
                proteins
            ),

        "direct_evidence_classes":
            ";".join(
                evidence
            ),

        "best_query_candidates":
            ";".join(
                candidates
            ),

        "GA_HMM_models":
            ";".join(
                hmm
            ),

        "producer":
            row.get(
                "producers",
                ""
            ),

        "incomplete_context":
            row.get(
                "strongest_incomplete_context",
                ""
            ),

        "auxiliary_context_flags":
            ";".join(
                context_flags
            ),

        "BAGEL_priority":
            priority,

        "antiSMASH_regions":
            row.get(
                "n_antismash_regions",
                ""
            ),

        "antiSMASH_RiPP_regions":
            row.get(
                "n_RiPP_relevant_regions",
                ""
            ),
    })


# ============================================================
# Extract FASTA
# ============================================================

found = {}

current_id = None
seq = []
keep = False


def save():

    if (
        keep
        and current_id
    ):

        found[
            current_id
        ] = "".join(
            seq
        )


with FASTA.open(
    errors="replace"
) as fh:

    for line in fh:

        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):

            save()

            current_id = (
                line[1:]
                .split()[0]
            )

            keep = (
                current_id
                in negative_ids
            )

            seq = []

        elif keep:

            seq.append(
                line
            )

save()


missing = (
    negative_ids
    -
    set(
        found
    )
)


if missing:

    raise RuntimeError(
        "Missing target contigs in FASTA: "
        + ",".join(
            sorted(
                missing
            )
        )
    )


target_fasta = (
    OUT
    / "98C3B0_antismash_negative7.fna"
)


with target_fasta.open(
    "w"
) as fh:

    for row in manifest:

        cid = row[
            "contig_id"
        ]

        sequence = found[
            cid
        ]

        fh.write(
            f">{cid}\n"
        )

        for i in range(
            0,
            len(sequence),
            80
        ):

            fh.write(
                sequence[
                    i:i+80
                ]
                + "\n"
            )


# ============================================================
# Length validation
# ============================================================

for row in manifest:

    cid = row[
        "contig_id"
    ]

    expected = int(
        row[
            "length_bp"
        ]
    )

    observed = len(
        found[
            cid
        ]
    )

    if expected != observed:

        raise RuntimeError(
            f"Length mismatch {cid}: "
            f"expected={expected}, observed={observed}"
        )


# ============================================================
# Write manifest
# ============================================================

manifest_out = (
    OUT
    / "98C3B0_negative7_manifest.tsv"
)


with manifest_out.open(
    "w",
    newline=""
) as fh:

    fields = list(
        manifest[0].keys()
    )

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(
        manifest
    )


# ============================================================
# Summary
# ============================================================

with (
    OUT
    / "98C3B0_target_summary.tsv"
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

    w.writerow([
        "primary_contigs_total",
        len(
            anti_rows
        )
    ])

    w.writerow([
        "antismash_RiPP_positive",
        len(
            anti_rows
        )
        -
        len(
            negative
        )
    ])

    w.writerow([
        "antismash_RiPP_negative",
        len(
            negative
        )
    ])

    w.writerow([
        "negative_contigs_with_direct_sequence_evidence",
        sum(
            r[
                "BAGEL_priority"
            ]
            ==
            "HIGH_direct_sequence_evidence"
            for r in manifest
        )
    ])

    w.writerow([
        "negative_contigs_GA_HMM_only",
        sum(
            r[
                "BAGEL_priority"
            ]
            ==
            "STANDARD_GA_HMM_screen"
            for r in manifest
        )
    ])

    w.writerow([
        "BAGEL_run_performed",
        "NO"
    ])

    w.writerow([
        "next_step",
        "98C3B_targeted_BAGEL_after_command_audit"
    ])


print(
    "NEGATIVE_CONTIGS=7"
)

print(
    "98C3B0_TARGET_PREPARATION=PASS"
)
