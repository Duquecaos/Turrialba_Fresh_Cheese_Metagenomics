#!/usr/bin/env python3

import sys
import statistics
from pathlib import Path
from collections import Counter


if len(sys.argv) != 5:
    raise SystemExit(
        "Usage: 98c0_audit_incomplete_catalog.py "
        "<catalog.faa> <eggnog.tsv> <master_candidates.faa> <OUT>"
    )


CATALOG = Path(sys.argv[1])
EGGNOG = Path(sys.argv[2])
MASTER = Path(sys.argv[3])
OUT = Path(sys.argv[4])

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# FASTA parser
# ============================================================

def fasta_stats(path):

    ids = []
    lengths = []

    examples = []

    current_id = None
    current_header = None
    current_len = 0

    with path.open(
        errors="replace"
    ) as fh:

        for line in fh:

            line = line.rstrip("\n")

            if line.startswith(">"):

                if current_id is not None:

                    ids.append(
                        current_id
                    )

                    lengths.append(
                        current_len
                    )

                current_header = line[1:]

                current_id = (
                    current_header
                    .split()[0]
                )

                current_len = 0

                if len(examples) < 30:
                    examples.append(
                        current_header
                    )

            else:

                current_len += len(
                    line.strip()
                    .replace("*", "")
                )


    if current_id is not None:

        ids.append(
            current_id
        )

        lengths.append(
            current_len
        )


    return (
        ids,
        lengths,
        examples
    )


catalog_ids, catalog_lengths, catalog_examples = fasta_stats(
    CATALOG
)

master_ids, master_lengths, master_examples = fasta_stats(
    MASTER
)


catalog_id_set = set(
    catalog_ids
)

duplicate_ids = (
    len(catalog_ids)
    - len(catalog_id_set)
)


# ============================================================
# Length bins
# ============================================================

bins = Counter()

for x in catalog_lengths:

    if x <= 30:
        bins["1-30"] += 1

    elif x <= 60:
        bins["31-60"] += 1

    elif x <= 100:
        bins["61-100"] += 1

    elif x <= 200:
        bins["101-200"] += 1

    elif x <= 500:
        bins["201-500"] += 1

    else:
        bins[">500"] += 1


# ============================================================
# Header-pattern audit
# ============================================================

patterns = Counter()

for x in catalog_ids[:100000]:

    if "|" in x:
        patterns["contains_pipe"] += 1

    if "__" in x:
        patterns["contains_double_underscore"] += 1

    if "_k141_" in x:
        patterns["contains__k141_"] += 1

    if "k141_" in x:
        patterns["contains_k141_"] += 1

    if "#" in x:
        patterns["contains_hash"] += 1


# ============================================================
# eggNOG schema + exact ID compatibility
# ============================================================

eggnog_noncomment = 0
eggnog_match = 0
eggnog_unique_match = set()

column_counts = Counter()

eggnog_examples = []
first_noncomment = ""


with EGGNOG.open(
    errors="replace"
) as fh:

    for line_number, line in enumerate(
        fh,
        start=1
    ):

        line = line.rstrip("\n")

        if not line:
            continue

        if line.startswith("#"):
            continue

        if not first_noncomment:
            first_noncomment = line

        fields = line.split("\t")

        column_counts[
            len(fields)
        ] += 1

        eggnog_noncomment += 1

        query_id = fields[0]

        if query_id in catalog_id_set:

            eggnog_match += 1

            eggnog_unique_match.add(
                query_id
            )


        if len(eggnog_examples) < 10:

            eggnog_examples.append(
                (
                    line_number,
                    line[:4000]
                )
            )


# ============================================================
# Write summary
# ============================================================

summary = OUT / "98C0_global_summary.tsv"


with summary.open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    values = [
        (
            "catalog_path",
            str(CATALOG)
        ),
        (
            "catalog_sequences",
            len(catalog_ids)
        ),
        (
            "catalog_unique_ids",
            len(catalog_id_set)
        ),
        (
            "catalog_duplicate_ids",
            duplicate_ids
        ),
        (
            "catalog_total_aa",
            sum(catalog_lengths)
        ),
        (
            "catalog_min_aa",
            min(catalog_lengths)
            if catalog_lengths
            else 0
        ),
        (
            "catalog_median_aa",
            round(
                statistics.median(
                    catalog_lengths
                ),
                2
            )
            if catalog_lengths
            else 0
        ),
        (
            "catalog_mean_aa",
            round(
                statistics.mean(
                    catalog_lengths
                ),
                2
            )
            if catalog_lengths
            else 0
        ),
        (
            "catalog_max_aa",
            max(catalog_lengths)
            if catalog_lengths
            else 0
        ),
        (
            "master_candidate_sequences",
            len(master_ids)
        ),
        (
            "eggnog_noncomment_rows",
            eggnog_noncomment
        ),
        (
            "eggnog_rows_exactly_matching_catalog_ID",
            eggnog_match
        ),
        (
            "eggnog_unique_catalog_ID_matches",
            len(
                eggnog_unique_match
            )
        ),
        (
            "eggnog_catalog_ID_match_fraction",
            round(
                len(
                    eggnog_unique_match
                )
                / len(
                    catalog_id_set
                ),
                6
            )
            if catalog_id_set
            else 0
        ),
        (
            "new_prediction_performed",
            "NO"
        ),
        (
            "next_step",
            "98C1_HMM_and_sequence_similarity_screen"
        ),
    ]


    for key, value in values:

        fh.write(
            f"{key}\t{value}\n"
        )


# ============================================================
# Length distribution
# ============================================================

with (
    OUT
    / "98C0_catalog_length_distribution.tsv"
).open("w") as fh:

    fh.write(
        "aa_length_class\tproteins\n"
    )

    for key in [
        "1-30",
        "31-60",
        "61-100",
        "101-200",
        "201-500",
        ">500",
    ]:

        fh.write(
            f"{key}\t{bins[key]}\n"
        )


# ============================================================
# Header examples
# ============================================================

with (
    OUT
    / "98C0_catalog_header_examples.txt"
).open("w") as fh:

    for x in catalog_examples:

        fh.write(
            x + "\n"
        )


# ============================================================
# Pattern summary
# ============================================================

with (
    OUT
    / "98C0_header_pattern_summary.tsv"
).open("w") as fh:

    fh.write(
        "pattern\tcount_in_first_100000_ids\n"
    )

    for key in sorted(patterns):

        fh.write(
            f"{key}\t{patterns[key]}\n"
        )


# ============================================================
# eggNOG schema audit
# ============================================================

with (
    OUT
    / "98C0_eggnog_schema.txt"
).open("w") as fh:

    fh.write(
        "FIRST_NONCOMMENT_LINE\n"
    )

    fh.write(
        first_noncomment[:8000]
        + "\n\n"
    )

    fh.write(
        "COLUMN_COUNTS\n"
    )

    for n, count in sorted(
        column_counts.items()
    ):

        fh.write(
            f"{n}\t{count}\n"
        )

    fh.write(
        "\nFIRST_10_NONCOMMENT_ROWS\n"
    )

    for line_number, text in eggnog_examples:

        fh.write(
            f"{line_number}\t{text}\n"
        )


print(
    f"CATALOG_SEQUENCES={len(catalog_ids)}"
)

print(
    f"UNIQUE_IDS={len(catalog_id_set)}"
)

print(
    f"EGGNOG_MATCHED_IDS={len(eggnog_unique_match)}"
)

print(
    f"MASTER_CANDIDATES={len(master_ids)}"
)

print(
    "98C0_PYTHON=PASS"
)
