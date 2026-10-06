#!/usr/bin/env python3

import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


if len(sys.argv) != 4:
    raise SystemExit(
        "Usage: "
        "98b0_parse_existing_antismash.py "
        "<antismash_results_dir> "
        "<attrloc_summary.tsv> "
        "<output_dir>"
    )


RESULTS = Path(sys.argv[1])
ATTRLOC = Path(sys.argv[2])
OUT = Path(sys.argv[3])

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# Helpers
# ============================================================

RIPP_RE = re.compile(
    r"bacterioc|"
    r"ripp|"
    r"lanthi|"
    r"lantipeptide|"
    r"lasso|"
    r"sacti|"
    r"ranthi|"
    r"thiopep|"
    r"cyanobactin|"
    r"linaridin|"
    r"glycocin|"
    r"microviridin|"
    r"bottromycin|"
    r"proteusin|"
    r"lipolanthine|"
    r"epipeptide|"
    r"spliceotide|"
    r"lap",
    re.I
)


def parse_quoted_qualifier(lines, start_index):

    line = lines[start_index]

    if '"' not in line:
        return "", start_index

    value = line.split(
        '"',
        1
    )[1]

    i = start_index

    while '"' not in value:

        i += 1

        if i >= len(lines):
            break

        value += lines[i].strip()

    if '"' in value:
        value = value.split(
            '"',
            1
        )[0]

    return value, i


def parse_region_feature(path):

    text = path.read_text(
        errors="replace"
    )

    lines = text.splitlines()

    region_location = ""
    products = []
    rules = []
    contig_edge = ""
    region_number = ""

    in_region = False

    i = 0

    while i < len(lines):

        line = lines[i]

        # GenBank feature key in columns 6-20.
        feature_key = ""

        if len(line) >= 21:
            feature_key = line[5:21].strip()

        if feature_key:

            if feature_key == "region":

                if in_region:
                    break

                in_region = True
                region_location = line[21:].strip()

            elif in_region:
                break

        elif in_region:

            stripped = line.strip()

            if stripped.startswith(
                "/product="
            ):

                value, i = parse_quoted_qualifier(
                    lines,
                    i
                )

                products.append(
                    value
                )

            elif stripped.startswith(
                "/rules="
            ):

                value, i = parse_quoted_qualifier(
                    lines,
                    i
                )

                rules.append(
                    value
                )

            elif stripped.startswith(
                "/contig_edge="
            ):

                contig_edge = (
                    stripped
                    .split("=", 1)[1]
                    .strip()
                    .strip('"')
                )

            elif stripped.startswith(
                "/region_number="
            ):

                region_number = (
                    stripped
                    .split("=", 1)[1]
                    .strip()
                    .strip('"')
                )

        i += 1


    # Fallback:
    # some antiSMASH region files may not expose product
    # at the region feature exactly as expected.
    if not products:

        for line in lines:

            stripped = line.strip()

            if stripped.startswith(
                "/product="
            ):

                val = (
                    stripped
                    .split("=", 1)[1]
                    .strip()
                    .strip('"')
                )

                if val:
                    products.append(
                        val
                    )


    # Deduplicate while preserving order.
    products = list(
        dict.fromkeys(products)
    )

    rules = list(
        dict.fromkeys(rules)
    )


    return {
        "region_location":
            region_location,

        "products":
            ";".join(products),

        "rules":
            ";".join(rules),

        "contig_edge":
            contig_edge,

        "region_number":
            region_number,
    }


# ============================================================
# 1. Inventory MAG directories
# ============================================================

mag_dirs = sorted(
    p
    for p in RESULTS.iterdir()
    if p.is_dir()
)


mag_inventory = []


for magdir in mag_dirs:

    mag = magdir.name

    jsons = list(
        magdir.glob("*.json")
    )

    full_gbks = [
        p
        for p in magdir.glob("*.gbk")
        if ".region" not in p.name
    ]

    regions = list(
        magdir.glob("*.region*.gbk")
    )

    html = (
        magdir
        / "index.html"
    )


    mag_inventory.append({
        "MAG":
            mag,

        "json_files":
            len(jsons),

        "main_gbk_files":
            len(full_gbks),

        "region_gbk_files":
            len(regions),

        "html_present":
            int(
                html.is_file()
            ),

        "result_complete_operational":
            int(
                len(jsons) >= 1
                and len(full_gbks) >= 1
                and html.is_file()
            ),
    })


# ============================================================
# 2. Parse regions
# ============================================================

region_rows = []


for magdir in mag_dirs:

    mag = magdir.name

    for region_file in sorted(
        magdir.glob(
            "*.region*.gbk"
        )
    ):

        basename = region_file.name

        if ".region" in basename:

            contig = basename.split(
                ".region",
                1
            )[0]

        else:

            contig = basename.rsplit(
                ".",
                1
            )[0]


        m = re.search(
            r"\.region(\d+)\.gbk$",
            basename
        )

        region_file_number = (
            m.group(1)
            if m
            else ""
        )


        parsed = parse_region_feature(
            region_file
        )


        searchable = " ".join([
            parsed[
                "products"
            ],
            parsed[
                "rules"
            ],
        ])


        relevant = bool(
            RIPP_RE.search(
                searchable
            )
        )


        region_rows.append({
            "MAG":
                mag,

            "contig":
                contig,

            "region_file_number":
                region_file_number,

            "region_number_annotation":
                parsed[
                    "region_number"
                ],

            "region_location":
                parsed[
                    "region_location"
                ],

            "products":
                parsed[
                    "products"
                ],

            "rules":
                parsed[
                    "rules"
                ],

            "contig_edge":
                parsed[
                    "contig_edge"
                ],

            "bacteriocin_RiPP_relevant":
                int(
                    relevant
                ),

            "region_file":
                str(
                    region_file
                ),
        })


# ============================================================
# 3. Product summary
# ============================================================

product_counter = Counter()


for r in region_rows:

    products = [
        x.strip()
        for x in r[
            "products"
        ].split(";")
        if x.strip()
    ]

    if not products:
        products = [
            "UNRESOLVED_PRODUCT"
        ]

    for product in products:

        product_counter[
            product
        ] += 1


# ============================================================
# 4. ATTRLOC comparison
# ============================================================

with ATTRLOC.open(
    newline=""
) as fh:

    attr_rows = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )


regions_by_mag_contig = defaultdict(
    list
)


for r in region_rows:

    regions_by_mag_contig[
        (
            r["MAG"],
            r["contig"]
        )
    ].append(
        r
    )


attr_comparison = []


for a in attr_rows:

    key = (
        a[
            "MAG"
        ],
        a[
            "contig"
        ]
    )

    exact_regions = (
        regions_by_mag_contig.get(
            key,
            []
        )
    )

    relevant_regions = [
        r
        for r in exact_regions
        if r[
            "bacteriocin_RiPP_relevant"
        ] == 1
    ]


    attr_comparison.append({
        "locus_id":
            a[
                "locus_id"
            ],

        "MAG":
            a[
                "MAG"
            ],

        "contig":
            a[
                "contig"
            ],

        "candidates":
            a[
                "candidates"
            ],

        "context_length":
            a[
                "context_length"
            ],

        "contains_interrupted_candidate":
            a[
                "contains_interrupted_candidate"
            ],

        "left_edge_truncated":
            a[
                "left_edge_truncated"
            ],

        "right_edge_truncated":
            a[
                "right_edge_truncated"
            ],

        "antismash_region_same_contig":
            int(
                len(
                    exact_regions
                ) > 0
            ),

        "antismash_bacteriocin_RiPP_region_same_contig":
            int(
                len(
                    relevant_regions
                ) > 0
            ),

        "antismash_region_products":
            " | ".join(
                sorted({
                    r[
                        "products"
                    ]
                    for r
                    in exact_regions
                    if r[
                        "products"
                    ]
                })
            ),

        "antismash_region_files":
            " | ".join(
                r[
                    "region_file"
                ]
                for r
                in exact_regions
            ),

        "interpretation":
            (
                "orthogonal_antismash_support_same_contig"
                if relevant_regions
                else
                "no_bacteriocin_RiPP_antismash_region_on_same_contig"
            ),
    })


# ============================================================
# 5. Write tables
# ============================================================

def write_rows(path, rows, fields):

    with path.open(
        "w",
        newline=""
    ) as fh:

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


write_rows(
    OUT
    / "98B0_antismash_MAG_inventory.tsv",
    mag_inventory,
    [
        "MAG",
        "json_files",
        "main_gbk_files",
        "region_gbk_files",
        "html_present",
        "result_complete_operational",
    ]
)


write_rows(
    OUT
    / "98B0_all_antismash_regions.tsv",
    region_rows,
    [
        "MAG",
        "contig",
        "region_file_number",
        "region_number_annotation",
        "region_location",
        "products",
        "rules",
        "contig_edge",
        "bacteriocin_RiPP_relevant",
        "region_file",
    ]
)


ripp_rows = [
    r
    for r in region_rows
    if r[
        "bacteriocin_RiPP_relevant"
    ] == 1
]


write_rows(
    OUT
    / "98B0_bacteriocin_RiPP_regions.tsv",
    ripp_rows,
    [
        "MAG",
        "contig",
        "region_file_number",
        "region_number_annotation",
        "region_location",
        "products",
        "rules",
        "contig_edge",
        "bacteriocin_RiPP_relevant",
        "region_file",
    ]
)


with (
    OUT
    / "98B0_product_summary.tsv"
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
        "product",
        "regions"
    ])

    for product, n in sorted(
        product_counter.items(),
        key=lambda x: (
            -x[1],
            x[0]
        )
    ):

        w.writerow([
            product,
            n
        ])


write_rows(
    OUT
    / "98B0_ATTRLOC_antismash_comparison.tsv",
    attr_comparison,
    [
        "locus_id",
        "MAG",
        "contig",
        "candidates",
        "context_length",
        "contains_interrupted_candidate",
        "left_edge_truncated",
        "right_edge_truncated",
        "antismash_region_same_contig",
        "antismash_bacteriocin_RiPP_region_same_contig",
        "antismash_region_products",
        "antismash_region_files",
        "interpretation",
    ]
)


# ============================================================
# 6. Summary
# ============================================================

complete_mags = sum(
    int(
        r[
            "result_complete_operational"
        ]
    )
    for r in mag_inventory
)

attr_same_region = sum(
    int(
        r[
            "antismash_region_same_contig"
        ]
    )
    for r in attr_comparison
)

attr_ripp_region = sum(
    int(
        r[
            "antismash_bacteriocin_RiPP_region_same_contig"
        ]
    )
    for r in attr_comparison
)


with (
    OUT
    / "98B0_global_summary.tsv"
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

    rows = [
        (
            "antismash_MAG_directories",
            len(
                mag_inventory
            )
        ),
        (
            "operationally_complete_antismash_MAG_results",
            complete_mags
        ),
        (
            "all_antismash_regions",
            len(
                region_rows
            )
        ),
        (
            "bacteriocin_RiPP_relevant_regions",
            len(
                ripp_rows
            )
        ),
        (
            "ATTRLOC_total",
            len(
                attr_comparison
            )
        ),
        (
            "ATTRLOC_on_any_antismash_region_contig",
            attr_same_region
        ),
        (
            "ATTRLOC_on_bacteriocin_RiPP_antismash_region_contig",
            attr_ripp_region
        ),
        (
            "new_antismash_run_performed",
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
            "98B1_targeted_expansion_after_existing_antismash_review"
        ),
    ]

    w.writerows(
        rows
    )


print(
    f"MAG_DIRS={len(mag_inventory)}"
)

print(
    f"COMPLETE_MAG_RESULTS={complete_mags}"
)

print(
    f"REGIONS={len(region_rows)}"
)

print(
    f"BACTERIOCIN_RIPP_REGIONS={len(ripp_rows)}"
)

print(
    f"ATTRLOC_RIPP_SAME_CONTIG={attr_ripp_region}"
)

print(
    "98B0_PARSER=PASS"
)
