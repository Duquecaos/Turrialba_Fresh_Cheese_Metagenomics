#!/usr/bin/env python3

import csv
import math
import sys
from pathlib import Path


if len(sys.argv) != 14:
    raise SystemExit(
        "Usage: task_summary.py "
        "<task_id> <sample> <panel> <target_slug> "
        "<target_species> <priority> "
        "<classification_summary> <mapping_counts> "
        "<depth0> <depth10> <depth20> "
        "<species_q10> <output.tsv>"
    )


(
    task_id,
    sample,
    panel,
    target_slug,
    target_species,
    priority,
    classification_path,
    mapping_counts_path,
    depth0_path,
    depth10_path,
    depth20_path,
    species_q10_path,
    output_path,
) = sys.argv[1:]


def read_tsv(path):

    with Path(path).open(
        newline=""
    ) as fh:

        return list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


classification = {
    r["metric"]:
        r["value"]

    for r in read_tsv(
        classification_path
    )
}


mapping = {
    r["metric"]:
        r["value"]

    for r in read_tsv(
        mapping_counts_path
    )
}


d0 = read_tsv(
    depth0_path
)[0]

d10 = read_tsv(
    depth10_path
)[0]

d20 = read_tsv(
    depth20_path
)[0]


species_q10 = read_tsv(
    species_q10_path
)


def ivalue(
    dictionary,
    key
):

    return int(
        float(
            dictionary[
                key
            ]
        )
    )


def fvalue(
    dictionary,
    key
):

    return float(
        dictionary[
            key
        ]
    )


target_all = ivalue(
    mapping,
    "target_mapped_records_all"
)

target_q10 = ivalue(
    mapping,
    "target_mapped_records_q10"
)

target_q20 = ivalue(
    mapping,
    "target_mapped_records_q20"
)

panel_all = ivalue(
    mapping,
    "panel_mapped_records_all"
)

panel_q10 = ivalue(
    mapping,
    "panel_mapped_records_q10"
)

panel_q20 = ivalue(
    mapping,
    "panel_mapped_records_q20"
)


target_species_row = None

for r in species_q10:

    if r[
        "species_slug"
    ] == target_slug:

        target_species_row = r
        break


target_rank = None

for rank, r in enumerate(
    species_q10,
    1
):

    if r[
        "species_slug"
    ] == target_slug:

        target_rank = rank
        break


competitors = [
    r
    for r in species_q10
    if r[
        "species_slug"
    ] != target_slug
]


top_competitor = (
    competitors[0]
    if competitors
    else None
)


top_comp_count = (
    int(
        top_competitor[
            "mapped_alignment_records"
        ]
    )
    if top_competitor
    else 0
)


if top_comp_count > 0:

    target_to_comp = (
        target_q10
        / top_comp_count
    )

else:

    target_to_comp = math.inf


row = {
    "task_id":
        task_id,

    "sample":
        sample,

    "priority":
        priority,

    "panel":
        panel,

    "target_species":
        target_species,

    "target_slug":
        target_slug,

    "target_ref":
        d0[
            "target_ref"
        ],

    "target_length":
        d0[
            "target_length"
        ],

    "total_input_pairs":
        classification[
            "total_sequence_pairs"
        ],

    "panel_mapped_records_all":
        panel_all,

    "panel_mapped_records_q10":
        panel_q10,

    "panel_mapped_records_q20":
        panel_q20,

    "target_mapped_records_all":
        target_all,

    "target_mapped_records_q10":
        target_q10,

    "target_mapped_records_q20":
        target_q20,

    "target_q10_retention":
        (
            f"{target_q10/target_all:.10f}"
            if target_all
            else "0.0000000000"
        ),

    "target_q20_retention":
        (
            f"{target_q20/target_all:.10f}"
            if target_all
            else "0.0000000000"
        ),

    "target_share_of_panel_all":
        (
            f"{target_all/panel_all:.10f}"
            if panel_all
            else "0.0000000000"
        ),

    "target_share_of_panel_q10":
        (
            f"{target_q10/panel_q10:.10f}"
            if panel_q10
            else "0.0000000000"
        ),

    "target_share_of_panel_q20":
        (
            f"{target_q20/panel_q20:.10f}"
            if panel_q20
            else "0.0000000000"
        ),

    "target_proper_pair_records_all":
        mapping[
            "target_proper_pair_records_all"
        ],

    "target_proper_pair_records_q10":
        mapping[
            "target_proper_pair_records_q10"
        ],

    "target_proper_pair_records_q20":
        mapping[
            "target_proper_pair_records_q20"
        ],

    "proper_pair_fraction_target_all":
        (
            f"{ivalue(mapping, 'target_proper_pair_records_all')/target_all:.10f}"
            if target_all
            else "0.0000000000"
        ),

    "breadth_all":
        d0[
            "breadth_ge1x"
        ],

    "breadth_q10":
        d10[
            "breadth_ge1x"
        ],

    "breadth_q20":
        d20[
            "breadth_ge1x"
        ],

    "breadth5_all":
        d0[
            "breadth_ge5x"
        ],

    "breadth5_q10":
        d10[
            "breadth_ge5x"
        ],

    "breadth5_q20":
        d20[
            "breadth_ge5x"
        ],

    "mean_depth_all":
        d0[
            "mean_depth_all_positions"
        ],

    "mean_depth_q10":
        d10[
            "mean_depth_all_positions"
        ],

    "mean_depth_q20":
        d20[
            "mean_depth_all_positions"
        ],

    "window_signal_fraction_all":
        d0[
            "window_signal_fraction"
        ],

    "window_signal_fraction_q10":
        d10[
            "window_signal_fraction"
        ],

    "window_signal_fraction_q20":
        d20[
            "window_signal_fraction"
        ],

    "window_breadth10pct_fraction_q10":
        d10[
            "window_breadth_ge10pct_fraction"
        ],

    "maximum_uncovered_run_q10_bp":
        d10[
            "maximum_uncovered_run_bp"
        ],

    "target_rank_in_panel_q10":
        (
            target_rank
            if target_rank is not None
            else "NA"
        ),

    "top_competitor_q10":
        (
            top_competitor[
                "species_slug"
            ]
            if top_competitor
            else "NONE"
        ),

    "top_competitor_records_q10":
        top_comp_count,

    "target_to_top_competitor_ratio_q10":
        (
            f"{target_to_comp:.10f}"
            if math.isfinite(
                target_to_comp
            )
            else "INF"
        ),

    "interpretation":
        "competitive_chromosomal_mapping_support_not_pathogen_confirmation",
}


with Path(
    output_path
).open(
    "w",
    newline=""
) as fh:

    fields = list(
        row.keys()
    )

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerow(
        row
    )


print(
    f"TASK={task_id}"
)

print(
    f"TARGET_Q10={target_q10}"
)

print(
    f"BREADTH_Q10={d10['breadth_ge1x']}"
)

print(
    f"TARGET_RANK_Q10={target_rank}"
)

print(
    f"TOP_COMPETITOR={row['top_competitor_q10']}"
)
