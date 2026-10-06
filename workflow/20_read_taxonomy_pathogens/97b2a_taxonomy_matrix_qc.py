#!/usr/bin/env python3

import csv
import math
import sys
from collections import defaultdict, Counter
from pathlib import Path


ROOT = Path(sys.argv[1])

INROOT = (
    ROOT
    / "98_read_taxonomy"
    / "97B1_all_samples"
)

OUT = (
    ROOT
    / "98_read_taxonomy"
    / "97B2A_matrix_QC"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


SAMPLES = [
    "L1_1","L1_2","L1_3",
    "L2_1","L2_2","L2_3",
    "L3_1","L3_2","L3_3",
    "M1_1","M1_2","M1_3",
    "M2_1","M2_2","M2_3",
    "M3_1","M3_2","M3_3",
]


# NCBI taxonomy anchors
TAX_BACTERIA = "2"
TAX_ARCHAEA = "2157"
TAX_EUKARYOTA = "2759"
TAX_FUNGI = "4751"
TAX_VIRUSES = "10239"
TAX_HOMO = "9605"


def read_tsv(path):

    if not path.is_file():
        raise RuntimeError(
            f"Missing file: {path}"
        )

    with path.open(
        newline="",
        errors="replace"
    ) as fh:

        return list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


def write_tsv(
    path,
    rows,
    fields
):

    with path.open(
        "w",
        newline=""
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows)


def fnum(x):

    try:
        return float(x)
    except Exception:
        return 0.0


def inum(x):

    try:
        return int(float(x))
    except Exception:
        return 0


# ============================================================
# Parse lineage from Kraken report
# ============================================================

def parse_kraken_lineages(path):

    lineage_map = {}

    stack = []

    with path.open(
        errors="replace"
    ) as fh:

        for line in fh:

            parts = line.rstrip(
                "\n"
            ).split("\t")

            if len(parts) < 6:
                continue

            rank = parts[3].strip()
            taxid = parts[4].strip()

            raw_name = parts[5]

            indent = (
                len(raw_name)
                - len(
                    raw_name.lstrip(" ")
                )
            )

            name = raw_name.strip()

            while (
                stack
                and
                stack[-1][0] >= indent
            ):
                stack.pop()

            ancestors = [
                item[1]
                for item in stack
            ]

            lineage_taxids = (
                ancestors
                + [taxid]
            )

            lineage_names = [
                item[2]
                for item in stack
            ] + [name]


            # Category precedence matters.
            if TAX_HOMO in lineage_taxids:

                category = "Human"

            elif TAX_FUNGI in lineage_taxids:

                category = "Fungi"

            elif TAX_BACTERIA in lineage_taxids:

                category = "Bacteria"

            elif TAX_ARCHAEA in lineage_taxids:

                category = "Archaea"

            elif TAX_VIRUSES in lineage_taxids:

                category = "Viruses"

            elif TAX_EUKARYOTA in lineage_taxids:

                category = "Eukaryota_other"

            else:

                category = "Other_unresolved"


            lineage_map[taxid] = {
                "rank":
                    rank,

                "name":
                    name,

                "category":
                    category,

                "lineage_taxids":
                    ";".join(
                        lineage_taxids
                    ),

                "lineage_names":
                    ";".join(
                        lineage_names
                    ),
            }


            stack.append(
                (
                    indent,
                    taxid,
                    name,
                )
            )

    return lineage_map


# ============================================================
# First pass: metadata, summaries and global lineage dictionary
# ============================================================

sample_summary = {}

global_lineage = {}

lineage_conflicts = []


for sample in SAMPLES:

    sdir = (
        INROOT
        / sample
    )

    summary_file = (
        sdir
        / f"{sample}.classification_summary.tsv"
    )

    kraken_report = (
        sdir
        / f"{sample}.kraken2.report"
    )

    summary_rows = read_tsv(
        summary_file
    )

    summary = {
        r["metric"]:
            r["value"]
        for r in summary_rows
    }

    sample_summary[
        sample
    ] = summary


    local_map = parse_kraken_lineages(
        kraken_report
    )

    for taxid, info in local_map.items():

        if taxid not in global_lineage:

            global_lineage[
                taxid
            ] = info

        else:

            old = global_lineage[
                taxid
            ]

            if (
                old["category"]
                != info["category"]
            ):

                lineage_conflicts.append({
                    "taxid":
                        taxid,

                    "existing_category":
                        old[
                            "category"
                        ],

                    "new_category":
                        info[
                            "category"
                        ],

                    "sample":
                        sample,
                })


if lineage_conflicts:

    write_tsv(
        OUT
        / "97B2A_lineage_conflicts.tsv",
        lineage_conflicts,
        list(
            lineage_conflicts[0].keys()
        )
    )

    raise RuntimeError(
        f"Taxonomy lineage conflicts={len(lineage_conflicts)}"
    )


# ============================================================
# Process Bracken ranks
# ============================================================

def process_rank(rank_name):

    master = []

    unresolved = []


    for sample in SAMPLES:

        sdir = (
            INROOT
            / sample
        )

        infile = (
            sdir
            / f"{sample}.bracken.{rank_name}.augmented.tsv"
        )

        rows = read_tsv(
            infile
        )

        total_pairs = inum(
            sample_summary[
                sample
            ][
                "total_sequence_pairs"
            ]
        )


        for r in rows:

            taxid = str(
                r[
                    "taxonomy_id"
                ]
            ).strip()

            info = global_lineage.get(
                taxid
            )


            if info is None:

                category = (
                    "Other_unresolved"
                )

                lineage_names = ""

                unresolved.append({
                    "sample":
                        sample,

                    "rank":
                        rank_name,

                    "taxid":
                        taxid,

                    "name":
                        r[
                            "name"
                        ],
                })

            else:

                category = info[
                    "category"
                ]

                lineage_names = info[
                    "lineage_names"
                ]


            est = inum(
                r[
                    "new_est_reads"
                ]
            )

            frac_all = (
                est
                / total_pairs
                if total_pairs > 0
                else 0.0
            )


            rec = {
                "sample":
                    sample,

                "producer_code":
                    sample[0],

                "biological_unit_code":
                    sample.split("_")[0][1:],

                "within_unit_time_code":
                    sample.split("_")[1],

                "taxonomic_rank":
                    rank_name,

                "name":
                    r[
                        "name"
                    ].strip(),

                "taxonomy_id":
                    taxid,

                "taxonomy_category":
                    category,

                "lineage":
                    lineage_names,

                "new_est_reads":
                    est,

                "fraction_all_input_pairs":
                    f"{frac_all:.10f}",

                "bracken_fraction_total_reads":
                    r[
                        "fraction_total_reads"
                    ],

                "cellular_microbiota_included":
                    (
                        "YES"
                        if category in {
                            "Bacteria",
                            "Archaea",
                            "Fungi",
                        }
                        else "NO"
                    ),

                "bacterial_included":
                    (
                        "YES"
                        if category
                        == "Bacteria"
                        else "NO"
                    ),

                "fungal_included":
                    (
                        "YES"
                        if category
                        == "Fungi"
                        else "NO"
                    ),

                "human_excluded":
                    (
                        "YES"
                        if category
                        == "Human"
                        else "NO"
                    ),
            }

            master.append(
                rec
            )


    return master, unresolved


genus_master, genus_unresolved = process_rank(
    "genus"
)

species_master, species_unresolved = process_rank(
    "species"
)


write_tsv(
    OUT
    / "97B2A_genus_master.tsv",
    genus_master,
    list(
        genus_master[0].keys()
    )
)

write_tsv(
    OUT
    / "97B2A_species_master.tsv",
    species_master,
    list(
        species_master[0].keys()
    )
)


unresolved_all = (
    genus_unresolved
    + species_unresolved
)


if unresolved_all:

    write_tsv(
        OUT
        / "97B2A_unresolved_taxonomy.tsv",
        unresolved_all,
        list(
            unresolved_all[0].keys()
        )
    )

else:

    write_tsv(
        OUT
        / "97B2A_unresolved_taxonomy.tsv",
        [],
        [
            "sample",
            "rank",
            "taxid",
            "name",
        ]
    )


# ============================================================
# Category summaries by sample
# ============================================================

category_rows = []


for sample in SAMPLES:

    rr = [
        r
        for r in genus_master
        if r[
            "sample"
        ] == sample
    ]

    total_pairs = inum(
        sample_summary[
            sample
        ][
            "total_sequence_pairs"
        ]
    )

    classified_pairs = inum(
        sample_summary[
            sample
        ][
            "kraken_classified_pairs"
        ]
    )

    unclassified_pairs = inum(
        sample_summary[
            sample
        ][
            "kraken_unclassified_pairs"
        ]
    )


    counts = defaultdict(int)

    for r in rr:

        counts[
            r[
                "taxonomy_category"
            ]
        ] += inum(
            r[
                "new_est_reads"
            ]
        )


    genus_resolved = sum(
        counts.values()
    )


    for category in [
        "Bacteria",
        "Archaea",
        "Fungi",
        "Viruses",
        "Human",
        "Eukaryota_other",
        "Other_unresolved",
    ]:

        n = counts[
            category
        ]

        category_rows.append({
            "sample":
                sample,

            "producer_code":
                sample[0],

            "biological_unit_code":
                sample.split("_")[0][1:],

            "within_unit_time_code":
                sample.split("_")[1],

            "category":
                category,

            "estimated_reads":
                n,

            "pct_all_input_pairs":
                (
                    f"{100*n/total_pairs:.8f}"
                    if total_pairs
                    else "0"
                ),

            "pct_genus_resolved":
                (
                    f"{100*n/genus_resolved:.8f}"
                    if genus_resolved
                    else "0"
                ),
        })


    category_rows.append({
        "sample":
            sample,

        "producer_code":
            sample[0],

        "biological_unit_code":
            sample.split("_")[0][1:],

        "within_unit_time_code":
            sample.split("_")[1],

        "category":
            "Kraken_unclassified",

        "estimated_reads":
            unclassified_pairs,

        "pct_all_input_pairs":
            (
                f"{100*unclassified_pairs/total_pairs:.8f}"
            ),

        "pct_genus_resolved":
            "NA",
    })


write_tsv(
    OUT
    / "97B2A_category_summary_by_sample.tsv",
    category_rows,
    list(
        category_rows[0].keys()
    )
)


# ============================================================
# Add relative abundance within cellular microbiota
# ============================================================

cellular_totals = defaultdict(int)

bacterial_totals = defaultdict(int)

fungal_totals = defaultdict(int)


for r in genus_master:

    sample = r["sample"]
    est = inum(
        r[
            "new_est_reads"
        ]
    )

    if r[
        "cellular_microbiota_included"
    ] == "YES":

        cellular_totals[
            sample
        ] += est

    if r[
        "bacterial_included"
    ] == "YES":

        bacterial_totals[
            sample
        ] += est

    if r[
        "fungal_included"
    ] == "YES":

        fungal_totals[
            sample
        ] += est


for r in genus_master:

    sample = r[
        "sample"
    ]

    est = inum(
        r[
            "new_est_reads"
        ]
    )


    if r[
        "cellular_microbiota_included"
    ] == "YES":

        denom = cellular_totals[
            sample
        ]

        r[
            "relative_within_cellular_microbiota"
        ] = (
            f"{est/denom:.10f}"
            if denom
            else "0"
        )

    else:

        r[
            "relative_within_cellular_microbiota"
        ] = "NA"


    if r[
        "bacterial_included"
    ] == "YES":

        denom = bacterial_totals[
            sample
        ]

        r[
            "relative_within_bacteria"
        ] = (
            f"{est/denom:.10f}"
            if denom
            else "0"
        )

    else:

        r[
            "relative_within_bacteria"
        ] = "NA"


    if r[
        "fungal_included"
    ] == "YES":

        denom = fungal_totals[
            sample
        ]

        r[
            "relative_within_fungi"
        ] = (
            f"{est/denom:.10f}"
            if denom
            else "0"
        )

    else:

        r[
            "relative_within_fungi"
        ] = "NA"


# Rewrite master with relative fields
write_tsv(
    OUT
    / "97B2A_genus_master.tsv",
    genus_master,
    list(
        genus_master[0].keys()
    )
)


# ============================================================
# Matrices
#
# IMPORTANT:
# Every sample was evaluated against the same Kraken/Bracken
# database. Therefore a missing taxon in a sample can be
# represented as zero here; this is unlike the 95D native-
# coassembly structural NA situation.
# ============================================================

def make_matrix(
    rows,
    inclusion_field,
    value_field,
    outfile
):

    filtered = [
        r
        for r in rows
        if r[
            inclusion_field
        ] == "YES"
    ]

    taxa = sorted({
        (
            r[
                "taxonomy_id"
            ],
            r[
                "name"
            ]
        )
        for r in filtered
    })


    lookup = defaultdict(float)

    for r in filtered:

        key = (
            r[
                "taxonomy_id"
            ],
            r[
                "name"
            ],
            r[
                "sample"
            ],
        )

        val = r[
            value_field
        ]

        if val == "NA":
            continue

        lookup[key] += float(
            val
        )


    matrix = []

    for taxid, name in taxa:

        rec = {
            "taxonomy_id":
                taxid,

            "name":
                name,
        }

        for sample in SAMPLES:

            rec[
                sample
            ] = (
                f"{lookup[(taxid, name, sample)]:.10f}"
                if "relative" in value_field
                or "fraction" in value_field
                else str(
                    int(
                        round(
                            lookup[
                                (
                                    taxid,
                                    name,
                                    sample
                                )
                            ]
                        )
                    )
                )
            )

        matrix.append(
            rec
        )


    write_tsv(
        OUT
        / outfile,
        matrix,
        [
            "taxonomy_id",
            "name",
        ]
        + SAMPLES
    )

    return len(
        matrix
    )


cellular_taxa_n = make_matrix(
    genus_master,
    "cellular_microbiota_included",
    "new_est_reads",
    "97B2A_genus_counts_cellular_microbiota.tsv",
)

make_matrix(
    genus_master,
    "cellular_microbiota_included",
    "fraction_all_input_pairs",
    "97B2A_genus_fraction_all_pairs_cellular_microbiota.tsv",
)

make_matrix(
    genus_master,
    "cellular_microbiota_included",
    "relative_within_cellular_microbiota",
    "97B2A_genus_relative_cellular_microbiota.tsv",
)

bacterial_taxa_n = make_matrix(
    genus_master,
    "bacterial_included",
    "new_est_reads",
    "97B2A_genus_counts_bacteria.tsv",
)

make_matrix(
    genus_master,
    "bacterial_included",
    "relative_within_bacteria",
    "97B2A_genus_relative_bacteria.tsv",
)

fungal_taxa_n = make_matrix(
    genus_master,
    "fungal_included",
    "new_est_reads",
    "97B2A_genus_counts_fungi.tsv",
)

make_matrix(
    genus_master,
    "fungal_included",
    "relative_within_fungi",
    "97B2A_genus_relative_fungi.tsv",
)


# ============================================================
# Top 20 cellular microbial genera per sample
# ============================================================

top_rows = []


for sample in SAMPLES:

    rr = [
        r
        for r in genus_master
        if (
            r[
                "sample"
            ] == sample
            and
            r[
                "cellular_microbiota_included"
            ] == "YES"
        )
    ]

    rr.sort(
        key=lambda r:
            inum(
                r[
                    "new_est_reads"
                ]
            ),
        reverse=True
    )


    for rank, r in enumerate(
        rr[:20],
        1
    ):

        top_rows.append({
            "sample":
                sample,

            "rank":
                rank,

            "name":
                r[
                    "name"
                ],

            "taxonomy_id":
                r[
                    "taxonomy_id"
                ],

            "category":
                r[
                    "taxonomy_category"
                ],

            "new_est_reads":
                r[
                    "new_est_reads"
                ],

            "fraction_all_input_pairs":
                r[
                    "fraction_all_input_pairs"
                ],

            "relative_within_cellular_microbiota":
                r[
                    "relative_within_cellular_microbiota"
                ],
        })


write_tsv(
    OUT
    / "97B2A_top20_cellular_microbiota_by_sample.tsv",
    top_rows,
    list(
        top_rows[0].keys()
    )
)


# ============================================================
# Non-bacterial/high-interest taxa
# ============================================================

nonbacterial_rows = [
    r
    for r in genus_master
    if r[
        "taxonomy_category"
    ] != "Bacteria"
]

nonbacterial_rows.sort(
    key=lambda r: (
        r["sample"],
        -inum(
            r[
                "new_est_reads"
            ]
        ),
    )
)


write_tsv(
    OUT
    / "97B2A_nonbacterial_genus_signals.tsv",
    nonbacterial_rows,
    list(
        nonbacterial_rows[0].keys()
    )
)


# ============================================================
# Abundance threshold sensitivity
#
# Threshold is relative abundance WITHIN the cellular
# microbial genus composition.
# These are descriptive only; no final filter is selected yet.
# ============================================================

thresholds = [
    0.0,
    0.00001,   # 0.001%
    0.0001,    # 0.01%
    0.001,     # 0.1%
    0.01,      # 1%
]


threshold_rows = []


for sample in SAMPLES:

    rr = [
        r
        for r in genus_master
        if (
            r[
                "sample"
            ] == sample
            and
            r[
                "cellular_microbiota_included"
            ] == "YES"
        )
    ]


    for threshold in thresholds:

        kept = [
            r
            for r in rr
            if fnum(
                r[
                    "relative_within_cellular_microbiota"
                ]
            ) > threshold
        ]


        threshold_rows.append({
            "sample":
                sample,

            "threshold_fraction":
                threshold,

            "threshold_pct":
                f"{100*threshold:.6f}",

            "genera_retained":
                len(
                    kept
                ),

            "fraction_composition_retained":
                f"{sum(fnum(r['relative_within_cellular_microbiota']) for r in kept):.10f}",

            "status":
                "sensitivity_audit_not_final_filter",
        })


write_tsv(
    OUT
    / "97B2A_genus_filter_sensitivity.tsv",
    threshold_rows,
    list(
        threshold_rows[0].keys()
    )
)


# ============================================================
# Per-sample QC summary
# ============================================================

qc_rows = []


for sample in SAMPLES:

    summary = sample_summary[
        sample
    ]

    rr = [
        r
        for r in genus_master
        if r[
            "sample"
        ] == sample
    ]

    category_counts = Counter()

    for r in rr:

        category_counts[
            r[
                "taxonomy_category"
            ]
        ] += inum(
            r[
                "new_est_reads"
            ]
        )


    total_pairs = inum(
        summary[
            "total_sequence_pairs"
        ]
    )


    qc_rows.append({
        "sample":
            sample,

        "producer_code":
            sample[0],

        "biological_unit_code":
            sample.split("_")[0][1:],

        "within_unit_time_code":
            sample.split("_")[1],

        "total_pairs":
            total_pairs,

        "kraken_classified_pct":
            summary[
                "kraken_classified_pct"
            ],

        "genus_resolved_pct_all_pairs":
            summary[
                "bracken_genus_resolved_pct_all_pairs"
            ],

        "species_resolved_pct_all_pairs":
            summary[
                "bracken_species_resolved_pct_all_pairs"
            ],

        "human_pct_all_pairs":
            (
                f"{100*category_counts['Human']/total_pairs:.8f}"
            ),

        "bacteria_pct_all_pairs":
            (
                f"{100*category_counts['Bacteria']/total_pairs:.8f}"
            ),

        "archaea_pct_all_pairs":
            (
                f"{100*category_counts['Archaea']/total_pairs:.8f}"
            ),

        "fungi_pct_all_pairs":
            (
                f"{100*category_counts['Fungi']/total_pairs:.8f}"
            ),

        "viruses_pct_all_pairs":
            (
                f"{100*category_counts['Viruses']/total_pairs:.8f}"
            ),

        "cellular_microbiota_estimated_reads":
            cellular_totals[
                sample
            ],

        "bacterial_estimated_reads":
            bacterial_totals[
                sample
            ],

        "fungal_estimated_reads":
            fungal_totals[
                sample
            ],

        "interpretation":
            "Bracken_read_based_taxonomic_estimates",
    })


write_tsv(
    OUT
    / "97B2A_sample_QC_summary.tsv",
    qc_rows,
    list(
        qc_rows[0].keys()
    )
)


# ============================================================
# Sample metadata
# ============================================================

metadata_rows = []


for sample in SAMPLES:

    group, time_code = sample.split(
        "_"
    )

    metadata_rows.append({
        "sample":
            sample,

        "producer_code":
            group[0],

        "biological_unit_code":
            group[1:],

        "within_unit_time_code":
            time_code,

        "subject_id":
            group,

        "time_interpretation":
            "code_only_not_biological_week_assignment_in_97B2A",
    })


write_tsv(
    OUT
    / "97B2A_sample_metadata.tsv",
    metadata_rows,
    list(
        metadata_rows[0].keys()
    )
)


# ============================================================
# Global summary
# ============================================================

total_pairs_all = sum(
    inum(
        sample_summary[s][
            "total_sequence_pairs"
        ]
    )
    for s in SAMPLES
)


with (
    OUT
    / "97B2A_global_summary.tsv"
).open(
    "w"
) as fh:

    fh.write(
        "metric\tvalue\n"
    )

    fh.write(
        f"samples\t{len(SAMPLES)}\n"
    )

    fh.write(
        f"total_input_pairs\t{total_pairs_all}\n"
    )

    fh.write(
        f"genus_master_rows\t{len(genus_master)}\n"
    )

    fh.write(
        f"species_master_rows\t{len(species_master)}\n"
    )

    fh.write(
        f"cellular_microbiota_unique_genera\t{cellular_taxa_n}\n"
    )

    fh.write(
        f"bacterial_unique_genera\t{bacterial_taxa_n}\n"
    )

    fh.write(
        f"fungal_unique_genera\t{fungal_taxa_n}\n"
    )

    fh.write(
        f"unresolved_taxonomy_rows\t{len(unresolved_all)}\n"
    )

    fh.write(
        "human_taxa_primary_microbiota_matrix\tEXCLUDED\n"
    )

    fh.write(
        "viruses_primary_cellular_microbiota_matrix\tEXCLUDED\n"
    )

    fh.write(
        "archaea_primary_cellular_microbiota_matrix\tINCLUDED\n"
    )

    fh.write(
        "fungi_primary_cellular_microbiota_matrix\tINCLUDED\n"
    )

    fh.write(
        "zero_interpretation_in_matrices\t"
        "evaluated_same_database_no_estimated_reads_at_retained_rank\n"
    )

    fh.write(
        "final_abundance_filter_selected\tNO\n"
    )


# ============================================================
# Methodological scope
# ============================================================

with (
    OUT
    / "97B2A_methodological_scope.tsv"
).open(
    "w"
) as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        "analysis\tread_based_shotgun_taxonomy_matrix_QC\n"
    )

    fh.write(
        "classifier\tKraken2_plus_Bracken\n"
    )

    fh.write(
        "database\tPlusPF_2026-06-26\n"
    )

    fh.write(
        "primary_taxonomic_rank\tgenus\n"
    )

    fh.write(
        "species_level_role\tsecondary_support_only\n"
    )

    fh.write(
        "human_signal\tretained_in_raw_master_excluded_from_microbiota_matrix\n"
    )

    fh.write(
        "virus_signal\tretained_in_raw_master_excluded_from_cellular_microbiota_matrix\n"
    )

    fh.write(
        "cellular_microbiota_definition\tBacteria_plus_Archaea_plus_Fungi\n"
    )

    fh.write(
        "relative_abundance_denominator\tBracken_estimated_reads_within_retained_cellular_microbiota_genera\n"
    )

    fh.write(
        "fraction_all_pairs_denominator\tall_input_paired_sequences\n"
    )

    fh.write(
        "zeros_are_structural_NA\tNO\n"
    )

    fh.write(
        "all_samples_evaluated_same_reference_database\tYES\n"
    )

    fh.write(
        "final_filter_selected\tNO\n"
    )

    fh.write(
        "PERMANOVA_performed\tNO\n"
    )

    fh.write(
        "alpha_diversity_performed\tNO\n"
    )

    fh.write(
        "amplicon_data_used\tNO\n"
    )


# ============================================================
# Final validation
# ============================================================

required = [
    "97B2A_genus_master.tsv",
    "97B2A_species_master.tsv",
    "97B2A_category_summary_by_sample.tsv",
    "97B2A_genus_counts_cellular_microbiota.tsv",
    "97B2A_genus_relative_cellular_microbiota.tsv",
    "97B2A_genus_counts_bacteria.tsv",
    "97B2A_genus_relative_bacteria.tsv",
    "97B2A_genus_counts_fungi.tsv",
    "97B2A_genus_relative_fungi.tsv",
    "97B2A_top20_cellular_microbiota_by_sample.tsv",
    "97B2A_nonbacterial_genus_signals.tsv",
    "97B2A_genus_filter_sensitivity.tsv",
    "97B2A_sample_QC_summary.tsv",
    "97B2A_sample_metadata.tsv",
    "97B2A_global_summary.tsv",
    "97B2A_methodological_scope.tsv",
]


for name in required:

    path = (
        OUT
        / name
    )

    if (
        not path.is_file()
        or path.stat().st_size == 0
    ):

        raise RuntimeError(
            f"Missing/empty output: {name}"
        )


print(
    f"SAMPLES={len(SAMPLES)}"
)

print(
    f"TOTAL_PAIRS={total_pairs_all}"
)

print(
    f"CELLULAR_GENERA={cellular_taxa_n}"
)

print(
    f"BACTERIAL_GENERA={bacterial_taxa_n}"
)

print(
    f"FUNGAL_GENERA={fungal_taxa_n}"
)

print(
    f"UNRESOLVED={len(unresolved_all)}"
)

print(
    "FINAL_FILTER=NOT_SELECTED"
)

print(
    "97B2A=PASS"
)

