#!/usr/bin/env python3

import csv
import gzip
import re
import sys
from pathlib import Path
from collections import defaultdict, Counter


if len(sys.argv) != 9:
    raise SystemExit(
        "Usage: 98d1b_summarize_mapping.py "
        "<metadata.tsv> <sample> "
        "<depth_all.tsv.gz> <depth_q10.tsv.gz> "
        "<idxstats.tsv> <bowtie2.log> "
        "<mapping_counts.tsv> <OUTDIR>"
    )


META = Path(sys.argv[1])
SAMPLE = sys.argv[2]
DEPTH_ALL = Path(sys.argv[3])
DEPTH_Q10 = Path(sys.argv[4])
IDXSTATS = Path(sys.argv[5])
BTLOG = Path(sys.argv[6])
COUNTS = Path(sys.argv[7])
OUT = Path(sys.argv[8])

OUT.mkdir(parents=True, exist_ok=True)


def read_tsv(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_metric_table(path):
    result = {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            result[r["metric"]] = r["value"]
    return result


# ============================================================
# Metadata
# ============================================================

meta = read_tsv(META)

if len(meta) != 30:
    raise RuntimeError(
        f"Expected 30 reference records; observed {len(meta)}"
    )

meta_by_ref = {
    r["reference_id"]: r
    for r in meta
}

expected_lengths = {
    r["reference_id"]:
        int(r["reference_length"])
    for r in meta
}


# ============================================================
# Bowtie2 log
# ============================================================

bt_text = BTLOG.read_text(errors="replace")

input_pairs = None
overall_rate = None


m = re.search(
    r"^\s*([\d,]+)\s+\([0-9.]+%\)\s+were paired; of these:",
    bt_text,
    flags=re.MULTILINE
)

if m:
    input_pairs = int(
        m.group(1).replace(",", "")
    )


if input_pairs is None:

    m = re.search(
        r"^\s*([\d,]+)\s+reads; of these:",
        bt_text,
        flags=re.MULTILINE
    )

    if m:
        input_pairs = int(
            m.group(1).replace(",", "")
        )


m = re.search(
    r"([0-9.]+)% overall alignment rate",
    bt_text
)

if m:
    overall_rate = float(m.group(1))


if input_pairs is None:
    raise RuntimeError(
        "Could not parse paired input count from Bowtie2 log"
    )


# ============================================================
# Mapping counts
# ============================================================

mapping_counts = read_metric_table(COUNTS)

mapped_all = int(
    mapping_counts[
        "mapped_alignment_records_all"
    ]
)

mapped_q10 = int(
    mapping_counts[
        "mapped_alignment_records_MAPQ10"
    ]
)


# ============================================================
# IDXSTATS
# ============================================================

idxstats = {}

with IDXSTATS.open() as fh:

    for line in fh:

        f = line.rstrip("\n").split("\t")

        if len(f) < 4:
            continue

        ref = f[0]

        if ref == "*":
            continue

        if ref not in expected_lengths:
            raise RuntimeError(
                f"Unexpected reference in idxstats: {ref}"
            )

        idxstats[ref] = {
            "length":
                int(f[1]),

            "mapped":
                int(f[2]),

            "unmapped":
                int(f[3]),
        }


if set(idxstats) != set(expected_lengths):
    raise RuntimeError(
        "IDXSTATS reference set differs from metadata"
    )


for ref in expected_lengths:

    if idxstats[ref]["length"] != expected_lengths[ref]:

        raise RuntimeError(
            f"Length mismatch in idxstats for {ref}: "
            f"{idxstats[ref]['length']} != "
            f"{expected_lengths[ref]}"
        )


# ============================================================
# Depth accumulator
# ============================================================

def summarize_depth(path):

    stats = {
        ref: {
            "positions": 0,
            "depth_sum": 0,
            "ge1": 0,
            "ge5": 0,
            "ge10": 0,
        }
        for ref in expected_lengths
    }


    with gzip.open(
        path,
        "rt"
    ) as fh:

        for line in fh:

            if not line.strip():
                continue

            ref, pos, depth = (
                line.rstrip("\n")
                .split("\t")[:3]
            )

            if ref not in stats:
                raise RuntimeError(
                    f"Unexpected depth reference: {ref}"
                )

            depth = int(depth)

            s = stats[ref]

            s["positions"] += 1
            s["depth_sum"] += depth

            if depth >= 1:
                s["ge1"] += 1

            if depth >= 5:
                s["ge5"] += 1

            if depth >= 10:
                s["ge10"] += 1


    for ref, length in expected_lengths.items():

        observed = stats[ref]["positions"]

        if observed != length:

            raise RuntimeError(
                f"Depth length mismatch {ref}: "
                f"{observed} positions != {length}"
            )


    return stats


all_stats = summarize_depth(
    DEPTH_ALL
)

q10_stats = summarize_depth(
    DEPTH_Q10
)


# ============================================================
# Build per-reference metrics
# ============================================================

rows = []


for ref in sorted(expected_lengths):

    length = expected_lengths[ref]

    a = all_stats[ref]
    q = q10_stats[ref]

    mean_all = (
        a["depth_sum"]
        / length
    )

    mean_q10 = (
        q["depth_sum"]
        / length
    )

    breadth_all = (
        a["ge1"]
        / length
    )

    breadth_q10 = (
        q["ge1"]
        / length
    )


    # Operational context-support classes.
    # These are descriptive screening categories only.
    if breadth_all == 0:

        detection_class = (
            "no_signal"
        )

    elif (
        breadth_all >= 0.90
        and
        breadth_q10 >= 0.75
    ):

        detection_class = (
            "robust_distributed_context_signal"
        )

    elif breadth_all >= 0.90:

        detection_class = (
            "high_breadth_MAPQ_sensitive"
        )

    elif breadth_all >= 0.50:

        detection_class = (
            "partial_context_signal"
        )

    else:

        detection_class = (
            "weak_or_localized_signal"
        )


    per_million_pairs = (
        input_pairs
        / 1_000_000.0
    )


    mapped_ref = idxstats[
        ref
    ]["mapped"]


    rows.append({
        "sample":
            SAMPLE,

        "reference_id":
            ref,

        "source_id":
            meta_by_ref[ref][
                "source_id"
            ],

        "reference_class":
            meta_by_ref[ref][
                "reference_class"
            ],

        "mapping_role":
            meta_by_ref[ref][
                "mapping_role"
            ],

        "context_group":
            meta_by_ref[ref][
                "context_group"
            ],

        "producer_source_context":
            meta_by_ref[ref][
                "producer"
            ],

        "candidate_protein_ids":
            meta_by_ref[ref][
                "candidate_protein_ids"
            ],

        "reference_length":
            length,

        "input_read_pairs":
            input_pairs,

        "overall_alignment_rate_pct":
            (
                overall_rate
                if overall_rate is not None
                else "NA"
            ),

        "mapped_alignment_records_reference":
            mapped_ref,

        "mapped_records_per_million_input_pairs":
            (
                mapped_ref
                / per_million_pairs
                if per_million_pairs > 0
                else 0
            ),

        "mean_depth_all":
            mean_all,

        "breadth_all_1x":
            breadth_all,

        "breadth_all_5x":
            a["ge5"] / length,

        "breadth_all_10x":
            a["ge10"] / length,

        "mean_depth_MAPQ10":
            mean_q10,

        "breadth_MAPQ10_1x":
            breadth_q10,

        "breadth_MAPQ10_5x":
            q["ge5"] / length,

        "breadth_MAPQ10_10x":
            q["ge10"] / length,

        "mean_depth_all_per_million_input_pairs":
            (
                mean_all
                / per_million_pairs
                if per_million_pairs > 0
                else 0
            ),

        "mean_depth_MAPQ10_per_million_input_pairs":
            (
                mean_q10
                / per_million_pairs
                if per_million_pairs > 0
                else 0
            ),

        "MAPQ10_mean_depth_retention":
            (
                mean_q10 / mean_all
                if mean_all > 0
                else "NA"
            ),

        "MAPQ10_breadth_retention":
            (
                breadth_q10 / breadth_all
                if breadth_all > 0
                else "NA"
            ),

        "operational_detection_class":
            detection_class,

        "functional_bacteriocin_confirmed":
            "NO",
    })


# ============================================================
# Write per-reference table
# ============================================================

with (
    OUT
    / f"{SAMPLE}.98D1B_locus_metrics.tsv"
).open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=list(
            rows[0].keys()
        ),
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(rows)


# ============================================================
# Sample summary
# ============================================================

class_counts = Counter(
    r[
        "operational_detection_class"
    ]
    for r in rows
)


with (
    OUT
    / f"{SAMPLE}.98D1B_sample_summary.tsv"
).open(
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


    summary = [
        (
            "sample",
            SAMPLE
        ),
        (
            "input_read_pairs",
            input_pairs
        ),
        (
            "overall_alignment_rate_pct",
            (
                overall_rate
                if overall_rate is not None
                else "NA"
            )
        ),
        (
            "mapped_alignment_records_all",
            mapped_all
        ),
        (
            "mapped_alignment_records_MAPQ10",
            mapped_q10
        ),
        (
            "MAPQ10_alignment_retention",
            (
                mapped_q10 / mapped_all
                if mapped_all > 0
                else "NA"
            )
        ),
        (
            "reference_sequences",
            len(rows)
        ),
        (
            "robust_distributed_context_signal",
            class_counts[
                "robust_distributed_context_signal"
            ]
        ),
        (
            "high_breadth_MAPQ_sensitive",
            class_counts[
                "high_breadth_MAPQ_sensitive"
            ]
        ),
        (
            "partial_context_signal",
            class_counts[
                "partial_context_signal"
            ]
        ),
        (
            "weak_or_localized_signal",
            class_counts[
                "weak_or_localized_signal"
            ]
        ),
        (
            "no_signal",
            class_counts[
                "no_signal"
            ]
        ),
        (
            "depth_all_positions",
            sum(
                r[
                    "reference_length"
                ]
                for r in rows
            )
        ),
        (
            "depth_MAPQ10_positions",
            sum(
                r[
                    "reference_length"
                ]
                for r in rows
            )
        ),
        (
            "amplicon_data_used",
            "NO"
        ),
        (
            "functional_activity_assessed",
            "NO"
        ),
    ]


    writer.writerows(summary)


print(
    f"SAMPLE={SAMPLE}"
)

print(
    f"INPUT_PAIRS={input_pairs}"
)

print(
    f"MAPPED_ALL={mapped_all}"
)

print(
    f"MAPPED_Q10={mapped_q10}"
)

print(
    "ROBUST="
    + str(
        class_counts[
            "robust_distributed_context_signal"
        ]
    )
)

print(
    "98D1B_SUMMARY=PASS"
)
