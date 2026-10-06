#!/usr/bin/env python3

import csv
import math
import sys
from pathlib import Path


if len(sys.argv) != 10:
    raise SystemExit(
        "Usage: depth_metrics.py "
        "<sample> <panel> <target_species> <target_ref> "
        "<target_length> <mapq> <window_size> "
        "<summary.tsv> <windows.tsv>"
    )


sample = sys.argv[1]
panel = sys.argv[2]
target_species = sys.argv[3]
target_ref = sys.argv[4]

target_length = int(
    sys.argv[5]
)

mapq = int(
    sys.argv[6]
)

window_size = int(
    sys.argv[7]
)

summary_path = Path(
    sys.argv[8]
)

windows_path = Path(
    sys.argv[9]
)


n_windows = math.ceil(
    target_length
    / window_size
)


window_length = [
    min(
        window_size,
        target_length - i * window_size
    )
    for i in range(n_windows)
]

window_covered = [
    0
    for _ in range(n_windows)
]

window_depth_sum = [
    0
    for _ in range(n_windows)
]


positions_seen = 0
covered = 0
depth_ge5 = 0
depth_ge10 = 0
depth_sum = 0

current_zero_run = 0
max_zero_run = 0


for line in sys.stdin:

    if not line.strip():
        continue

    fields = line.rstrip(
        "\n"
    ).split(
        "\t"
    )

    if len(fields) < 3:
        continue

    ref = fields[0]
    pos = int(
        fields[1]
    )
    depth = int(
        fields[2]
    )

    if ref != target_ref:
        raise RuntimeError(
            f"Unexpected reference: {ref}"
        )

    if not (
        1 <= pos <= target_length
    ):
        raise RuntimeError(
            f"Position outside target: {pos}"
        )

    positions_seen += 1
    depth_sum += depth

    window_index = (
        pos - 1
    ) // window_size

    window_depth_sum[
        window_index
    ] += depth

    if depth > 0:

        covered += 1

        window_covered[
            window_index
        ] += 1

        current_zero_run = 0

    else:

        current_zero_run += 1

        if current_zero_run > max_zero_run:
            max_zero_run = current_zero_run

    if depth >= 5:
        depth_ge5 += 1

    if depth >= 10:
        depth_ge10 += 1


if positions_seen != target_length:

    raise RuntimeError(
        f"Expected {target_length} depth positions, "
        f"observed {positions_seen}"
    )


breadth1 = (
    covered
    / target_length
)

breadth5 = (
    depth_ge5
    / target_length
)

breadth10 = (
    depth_ge10
    / target_length
)

mean_depth_all = (
    depth_sum
    / target_length
)

mean_depth_covered = (
    depth_sum
    / covered
    if covered
    else 0.0
)


windows_with_signal = sum(
    x > 0
    for x in window_covered
)

windows_breadth_ge1pct = 0
windows_breadth_ge10pct = 0

window_rows = []


for i in range(
    n_windows
):

    start = (
        i * window_size
        + 1
    )

    end = min(
        target_length,
        (i + 1) * window_size
    )

    wlen = window_length[i]

    wbreadth = (
        window_covered[i]
        / wlen
    )

    wmean = (
        window_depth_sum[i]
        / wlen
    )

    if wbreadth >= 0.01:
        windows_breadth_ge1pct += 1

    if wbreadth >= 0.10:
        windows_breadth_ge10pct += 1

    window_rows.append({
        "sample":
            sample,

        "panel":
            panel,

        "target_species":
            target_species,

        "mapq_threshold":
            mapq,

        "window":
            i + 1,

        "start":
            start,

        "end":
            end,

        "window_length":
            wlen,

        "covered_positions":
            window_covered[i],

        "breadth":
            f"{wbreadth:.10f}",

        "mean_depth":
            f"{wmean:.10f}",
    })


summary = [{
    "sample":
        sample,

    "panel":
        panel,

    "target_species":
        target_species,

    "target_ref":
        target_ref,

    "target_length":
        target_length,

    "mapq_threshold":
        mapq,

    "positions_seen":
        positions_seen,

    "covered_positions":
        covered,

    "breadth_ge1x":
        f"{breadth1:.10f}",

    "positions_depth_ge5":
        depth_ge5,

    "breadth_ge5x":
        f"{breadth5:.10f}",

    "positions_depth_ge10":
        depth_ge10,

    "breadth_ge10x":
        f"{breadth10:.10f}",

    "mean_depth_all_positions":
        f"{mean_depth_all:.10f}",

    "mean_depth_covered_positions":
        f"{mean_depth_covered:.10f}",

    "windows_total":
        n_windows,

    "windows_with_signal":
        windows_with_signal,

    "window_signal_fraction":
        f"{windows_with_signal / n_windows:.10f}",

    "windows_breadth_ge1pct":
        windows_breadth_ge1pct,

    "window_breadth_ge1pct_fraction":
        f"{windows_breadth_ge1pct / n_windows:.10f}",

    "windows_breadth_ge10pct":
        windows_breadth_ge10pct,

    "window_breadth_ge10pct_fraction":
        f"{windows_breadth_ge10pct / n_windows:.10f}",

    "maximum_uncovered_run_bp":
        max_zero_run,
}]


with summary_path.open(
    "w",
    newline=""
) as fh:

    fields = list(
        summary[0].keys()
    )

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(
        summary
    )


with windows_path.open(
    "w",
    newline=""
) as fh:

    fields = list(
        window_rows[0].keys()
    )

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(
        window_rows
    )


print(
    f"MAPQ={mapq}"
)

print(
    f"BREADTH={breadth1:.8f}"
)

print(
    f"MEAN_DEPTH={mean_depth_all:.8f}"
)

print(
    f"WINDOW_SIGNAL_FRACTION="
    f"{windows_with_signal/n_windows:.8f}"
)
