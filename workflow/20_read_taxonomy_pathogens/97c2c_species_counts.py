#!/usr/bin/env python3

import csv
import sys
from collections import Counter
from pathlib import Path


if len(sys.argv) != 3:
    raise SystemExit(
        "Usage: species_counts.py "
        "<mapq_threshold> <output.tsv>"
    )


mapq = int(
    sys.argv[1]
)

outfile = Path(
    sys.argv[2]
)


counts = Counter()


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

    rname = fields[2]

    if rname == "*":
        continue

    species_slug = rname.split(
        "|",
        1
    )[0]

    counts[
        species_slug
    ] += 1


total = sum(
    counts.values()
)


rows = []


for species, count in sorted(
    counts.items(),
    key=lambda x: (
        -x[1],
        x[0]
    )
):

    rows.append({
        "mapq_threshold":
            mapq,

        "species_slug":
            species,

        "mapped_alignment_records":
            count,

        "fraction_of_panel_mapped":
            (
                f"{count/total:.10f}"
                if total
                else "0.0000000000"
            ),
    })


with outfile.open(
    "w",
    newline=""
) as fh:

    fields = [
        "mapq_threshold",
        "species_slug",
        "mapped_alignment_records",
        "fraction_of_panel_mapped",
    ]

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(
        rows
    )


print(
    f"MAPQ={mapq}"
)

print(
    f"SPECIES_WITH_ALIGNMENTS={len(rows)}"
)

print(
    f"TOTAL_ALIGNMENT_RECORDS={total}"
)
