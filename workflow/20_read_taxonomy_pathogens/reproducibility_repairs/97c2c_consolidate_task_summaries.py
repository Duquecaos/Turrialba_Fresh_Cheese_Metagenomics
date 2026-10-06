#!/usr/bin/env python3

"""
Retrospective reproducibility repair for step 97C2C.

The original retained 97C2C array/worker scripts generated one
task_summary.tsv per competitive-mapping task but did not retain the
post-array command that produced 97C2C_all_task_summary.tsv.

This script deterministically reconstructs that global table from the
completed per-task summaries.

It is a reproducibility repair, not an original pipeline script.
"""

from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path


if len(sys.argv) != 3:
    raise SystemExit(
        "Usage: 97c2c_consolidate_task_summaries.py "
        "<97C2C_competitive_mapping_dir> <output.tsv>"
    )

base = Path(sys.argv[1]).resolve()
outfile = Path(sys.argv[2]).resolve()

if not base.is_dir():
    raise RuntimeError(f"Input directory not found: {base}")

summary_files = sorted(
    base.glob("*/task_summary.tsv")
)

if len(summary_files) != 11:
    raise RuntimeError(
        f"Expected 11 task_summary.tsv files; found {len(summary_files)}"
    )

rows = []
expected_fields = None

for path in summary_files:
    complete = path.parent / "97C2C_COMPLETE.ok"

    if not complete.is_file():
        raise RuntimeError(
            f"Missing completion marker for {path}: {complete}"
        )

    with path.open(newline="") as fh:
        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        fields = reader.fieldnames

        if not fields:
            raise RuntimeError(
                f"Missing header: {path}"
            )

        if expected_fields is None:
            expected_fields = fields
        elif fields != expected_fields:
            raise RuntimeError(
                f"Header mismatch: {path}"
            )

        local_rows = list(reader)

    if len(local_rows) != 1:
        raise RuntimeError(
            f"Expected exactly one row in {path}; "
            f"found {len(local_rows)}"
        )

    row = local_rows[0]

    try:
        task_id = int(row["task_id"])
    except Exception as exc:
        raise RuntimeError(
            f"Invalid task_id in {path}"
        ) from exc

    rows.append(
        (
            task_id,
            path,
            row
        )
    )


task_ids = [x[0] for x in rows]

if len(set(task_ids)) != len(task_ids):
    raise RuntimeError(
        f"Duplicate task IDs: {task_ids}"
    )

if sorted(task_ids) != list(range(11)):
    raise RuntimeError(
        f"Expected task IDs 0-10; observed {sorted(task_ids)}"
    )


rows.sort(
    key=lambda x: x[0]
)

outfile.parent.mkdir(
    parents=True,
    exist_ok=True
)

tmp = outfile.with_suffix(
    outfile.suffix + ".tmp"
)

with tmp.open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=expected_fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="raise"
    )

    writer.writeheader()

    for _, _, row in rows:
        writer.writerow(row)

tmp.replace(outfile)


sha = hashlib.sha256(
    outfile.read_bytes()
).hexdigest()

print(f"TASK_SUMMARIES={len(rows)}")
print(f"TASK_IDS={','.join(str(x[0]) for x in rows)}")
print(f"OUTPUT_ROWS={len(rows)}")
print(f"OUTPUT={outfile}")
print(f"SHA256={sha}")
print("97C2C_CONSOLIDATION=PASS")
