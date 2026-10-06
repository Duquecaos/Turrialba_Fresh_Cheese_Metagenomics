#!/usr/bin/env python3

import csv
import sys
from pathlib import Path
from collections import defaultdict


if len(sys.argv) != 4:
    raise SystemExit(
        "Usage: 98b1_integrate_existing_bacteriocin_evidence.py "
        "<ROOT> <ATTRLOC_summary.tsv> <OUT>"
    )


ROOT = Path(sys.argv[1])
ATTRLOC_FILE = Path(sys.argv[2])
OUT = Path(sys.argv[3])

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# Existing sources
# ============================================================

SOURCES = {
    "comparippson": ROOT / "29_comparippson_lactococcin",
    "bagel35": ROOT / "35_bagel4_test_lactis",
    "bagel36": ROOT / "36_bagel4_BAL4",
    "petauri37": ROOT / "37_petauri_bagel4_diagnostic",
    "master38": ROOT / "38_bacteriocin_master",
}


ANTISMASH_TABLE = (
    ROOT
    / "98_bacteriocin_mining"
    / "98B0_existing_antismash"
    / "98B0_ATTRLOC_antismash_comparison.tsv"
)


# ============================================================
# Read ATTRLOC
# ============================================================

with ATTRLOC_FILE.open(
    newline=""
) as fh:

    loci = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )


# ============================================================
# antiSMASH evidence
# ============================================================

anti = {}

if ANTISMASH_TABLE.is_file():

    with ANTISMASH_TABLE.open(
        newline=""
    ) as fh:

        for row in csv.DictReader(
            fh,
            delimiter="\t"
        ):

            anti[
                row["locus_id"]
            ] = row


# ============================================================
# File inventory
# ============================================================

inventory = []

all_files_by_source = {}


for source, root in SOURCES.items():

    files = []

    if root.is_dir():

        for p in root.rglob("*"):

            if p.is_file():
                files.append(p)

    all_files_by_source[
        source
    ] = files


    total_bytes = sum(
        p.stat().st_size
        for p in files
        if p.exists()
    )


    priority = []

    for p in files:

        name = p.name.lower()

        if any(
            key in name
            for key in [
                "candidate",
                "evidence",
                "summary",
                "best_hit",
                "best_hits",
                "comparippson",
                "lactococcin",
                "aoi",
                "threshold",
                "diagnostic",
            ]
        ):
            priority.append(p)


    inventory.append({
        "source_family":
            source,

        "root":
            str(root),

        "directory_present":
            int(root.is_dir()),

        "files_total":
            len(files),

        "priority_files":
            len(priority),

        "total_bytes":
            total_bytes,
    })


# ============================================================
# Helpers
# ============================================================

TEXT_SUFFIXES = {
    ".tsv",
    ".txt",
    ".table",
    ".json",
    ".gbk",
    ".faa",
    ".fna",
    ".predict",
    ".log",
    ".html",
    ".csv",
}


def priority_text_file(path):

    name = path.name.lower()

    return (
        path.suffix.lower()
        in TEXT_SUFFIXES
        and
        any(
            key in name
            for key in [
                "candidate",
                "evidence",
                "summary",
                "best",
                "comparippson",
                "lactococcin",
                "aoi",
                "threshold",
                "diagnostic",
                "blast",
            ]
        )
    )


def candidate_tokens(candidate_field):

    result = []

    for candidate in candidate_field.split(";"):

        candidate = candidate.strip()

        if not candidate:
            continue

        result.append(
            candidate
        )

        # Useful exact ORF/sORF terminal token.
        pieces = candidate.split("_")

        for part in pieces:

            low = part.lower()

            if (
                low.startswith("orf")
                or low.startswith("sorf")
            ):
                result.append(
                    part
                )

    return list(
        dict.fromkeys(result)
    )


# ============================================================
# Search source families
# ============================================================

detail_rows = []
matrix_rows = []


for locus in loci:

    locus_id = locus["locus_id"]
    contig = locus["contig"]
    candidates = locus["candidates"]

    tokens = [
        contig
    ] + candidate_tokens(
        candidates
    )


    source_stats = {}


    for source, files in all_files_by_source.items():

        filename_hits = []
        content_hits = []

        total_content_matches = 0


        # --------------------------------------------
        # Filename evidence
        # --------------------------------------------

        for p in files:

            path_text = str(p)

            if contig in path_text:

                filename_hits.append(
                    path_text
                )


        # --------------------------------------------
        # Content evidence
        # only curated/summarizing files
        # --------------------------------------------

        for p in files:

            if not priority_text_file(p):
                continue

            try:

                if p.stat().st_size > 20_000_000:
                    continue

                with p.open(
                    errors="replace"
                ) as fh:

                    for line_number, line in enumerate(
                        fh,
                        start=1
                    ):

                        if any(
                            token
                            and token in line
                            for token in tokens
                        ):

                            total_content_matches += 1

                            if len(content_hits) < 30:

                                content_hits.append(
                                    (
                                        str(p),
                                        line_number,
                                        line.rstrip()[:1000]
                                    )
                                )

            except Exception:
                continue


        source_stats[
            source
        ] = {
            "filename_hits":
                len(
                    filename_hits
                ),

            "content_hits":
                total_content_matches,

            "source_hit":
                int(
                    bool(
                        filename_hits
                        or content_hits
                    )
                ),

            "aoi_file_hits":
                sum(
                    1
                    for x in filename_hits
                    if "AOI" in x
                ),
        }


        for path in filename_hits[:30]:

            detail_rows.append({
                "locus_id":
                    locus_id,

                "MAG":
                    locus["MAG"],

                "contig":
                    contig,

                "source_family":
                    source,

                "match_type":
                    "filename",

                "file":
                    path,

                "line_number":
                    "",

                "matched_text":
                    "",
            })


        for path, line_no, text in content_hits:

            detail_rows.append({
                "locus_id":
                    locus_id,

                "MAG":
                    locus["MAG"],

                "contig":
                    contig,

                "source_family":
                    source,

                "match_type":
                    "content",

                "file":
                    path,

                "line_number":
                    line_no,

                "matched_text":
                    text,
            })


    # ========================================================
    # antiSMASH
    # ========================================================

    a = anti.get(
        locus_id,
        {}
    )

    antismash_support = int(
        a.get(
            "antismash_bacteriocin_RiPP_region_same_contig",
            "0"
        )
        == "1"
    )


    # BAGEL is one methodological family,
    # regardless of how many BAGEL runs contained the locus.
    bagel_support = int(
        source_stats["bagel35"]["source_hit"]
        or
        source_stats["bagel36"]["source_hit"]
    )


    comparippson_support = int(
        source_stats[
            "comparippson"
        ][
            "source_hit"
        ]
    )


    method_family_count = (
        antismash_support
        + bagel_support
        + comparippson_support
    )


    matrix_rows.append({
        "locus_id":
            locus_id,

        "MAG":
            locus["MAG"],

        "contig":
            contig,

        "candidates":
            candidates,

        "antismash_RiPP_same_contig":
            antismash_support,

        "antismash_products":
            a.get(
                "antismash_region_products",
                ""
            ),

        "BAGEL_source_hit":
            bagel_support,

        "BAGEL35_filename_hits":
            source_stats[
                "bagel35"
            ][
                "filename_hits"
            ],

        "BAGEL35_content_hits":
            source_stats[
                "bagel35"
            ][
                "content_hits"
            ],

        "BAGEL35_AOI_file_hits":
            source_stats[
                "bagel35"
            ][
                "aoi_file_hits"
            ],

        "BAGEL36_filename_hits":
            source_stats[
                "bagel36"
            ][
                "filename_hits"
            ],

        "BAGEL36_content_hits":
            source_stats[
                "bagel36"
            ][
                "content_hits"
            ],

        "BAGEL36_AOI_file_hits":
            source_stats[
                "bagel36"
            ][
                "aoi_file_hits"
            ],

        "Comparippson_source_hit":
            comparippson_support,

        "Petauri_diagnostic_hit":
            source_stats[
                "petauri37"
            ][
                "source_hit"
            ],

        "Master38_hit":
            source_stats[
                "master38"
            ][
                "source_hit"
            ],

        "independent_prediction_method_families":
            method_family_count,

        "interpretive_scope":
            (
                "multiple_computational_method_families"
                if method_family_count >= 2
                else
                "single_computational_method_family"
                if method_family_count == 1
                else
                "no_match_in_selected_existing_prediction_sources"
            ),
    })


# ============================================================
# Write inventory
# ============================================================

def write_table(
    path,
    rows,
    fields
):

    with path.open(
        "w",
        newline=""
    ) as fh:

        w = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n"
        )

        w.writeheader()

        for row in rows:
            w.writerow(row)


write_table(
    OUT
    / "98B1_source_inventory.tsv",
    inventory,
    [
        "source_family",
        "root",
        "directory_present",
        "files_total",
        "priority_files",
        "total_bytes",
    ]
)


write_table(
    OUT
    / "98B1_ATTRLOC_evidence_matrix.tsv",
    matrix_rows,
    [
        "locus_id",
        "MAG",
        "contig",
        "candidates",
        "antismash_RiPP_same_contig",
        "antismash_products",
        "BAGEL_source_hit",
        "BAGEL35_filename_hits",
        "BAGEL35_content_hits",
        "BAGEL35_AOI_file_hits",
        "BAGEL36_filename_hits",
        "BAGEL36_content_hits",
        "BAGEL36_AOI_file_hits",
        "Comparippson_source_hit",
        "Petauri_diagnostic_hit",
        "Master38_hit",
        "independent_prediction_method_families",
        "interpretive_scope",
    ]
)


write_table(
    OUT
    / "98B1_match_details.tsv",
    detail_rows,
    [
        "locus_id",
        "MAG",
        "contig",
        "source_family",
        "match_type",
        "file",
        "line_number",
        "matched_text",
    ]
)


# ============================================================
# Summaries
# ============================================================

family_distribution = defaultdict(
    int
)

for row in matrix_rows:

    family_distribution[
        int(
            row[
                "independent_prediction_method_families"
            ]
        )
    ] += 1


with (
    OUT
    / "98B1_global_summary.tsv"
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
        "ATTRLOC_total",
        len(
            matrix_rows
        )
    ])

    w.writerow([
        "antiSMASH_supported_loci",
        sum(
            int(
                r[
                    "antismash_RiPP_same_contig"
                ]
            )
            for r in matrix_rows
        )
    ])

    w.writerow([
        "BAGEL_source_hit_loci",
        sum(
            int(
                r[
                    "BAGEL_source_hit"
                ]
            )
            for r in matrix_rows
        )
    ])

    w.writerow([
        "Comparippson_source_hit_loci",
        sum(
            int(
                r[
                    "Comparippson_source_hit"
                ]
            )
            for r in matrix_rows
        )
    ])

    w.writerow([
        "method_families_0",
        family_distribution[0]
    ])

    w.writerow([
        "method_families_1",
        family_distribution[1]
    ])

    w.writerow([
        "method_families_2",
        family_distribution[2]
    ])

    w.writerow([
        "method_families_3",
        family_distribution[3]
    ])

    w.writerow([
        "new_prediction_performed",
        "NO"
    ])

    w.writerow([
        "reads_remapped",
        "NO"
    ])

    w.writerow([
        "amplicon_data_used",
        "NO"
    ])

    w.writerow([
        "next_step",
        "98B2_exact_candidate_level_review_and_then_incomplete_bin_expansion"
    ])


print(
    f"ATTRLOC={len(matrix_rows)}"
)

print(
    "antiSMASH_supported="
    + str(
        sum(
            int(
                r[
                    "antismash_RiPP_same_contig"
                ]
            )
            for r in matrix_rows
        )
    )
)

print(
    "BAGEL_source_hits="
    + str(
        sum(
            int(
                r[
                    "BAGEL_source_hit"
                ]
            )
            for r in matrix_rows
        )
    )
)

print(
    "Comparippson_source_hits="
    + str(
        sum(
            int(
                r[
                    "Comparippson_source_hit"
                ]
            )
            for r in matrix_rows
        )
    )
)

print("98B1=PASS")
