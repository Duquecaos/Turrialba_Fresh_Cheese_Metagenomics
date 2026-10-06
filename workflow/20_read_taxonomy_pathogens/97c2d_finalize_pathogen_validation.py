#!/usr/bin/env python3

import csv
from pathlib import Path


ROOT = Path(
    "/scratch/global"
) / Path.home().name / "Shotgun_MAGs_Turrialba"

C1 = (
    ROOT
    / "98_read_taxonomy"
    / "97C1_direct_pathogen_evidence"
    / "97C1_all_target_sample_evidence.tsv"
)

C2 = (
    ROOT
    / "98_read_taxonomy"
    / "97C2C_competitive_mapping"
    / "97C2C_all_task_summary.tsv"
)

OUT = (
    ROOT
    / "98_read_taxonomy"
    / "97C2D_final_pathogen_validation"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


def read_tsv(path):

    with path.open(
        newline=""
    ) as fh:

        return list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


def f(x):

    try:
        return float(x)
    except Exception:
        return 0.0


def i(x):

    try:
        return int(float(x))
    except Exception:
        return 0


if not C1.is_file():
    raise RuntimeError(
        f"Missing 97C1 table: {C1}"
    )

if not C2.is_file():
    raise RuntimeError(
        f"Missing 97C2C table: {C2}"
    )


c1_rows = read_tsv(C1)
c2_rows = read_tsv(C2)


c1_lookup = {
    (
        r["sample"],
        r["species"]
    ):
        r

    for r in c1_rows
}


final = []


for r in c2_rows:

    key = (
        r["sample"],
        r["target_species"]
    )

    k = c1_lookup.get(
        key
    )

    if k is None:
        raise RuntimeError(
            f"No matching 97C1 row for {key}"
        )


    q10 = i(
        r[
            "target_mapped_records_q10"
        ]
    )

    q20 = i(
        r[
            "target_mapped_records_q20"
        ]
    )

    q10_ret = f(
        r[
            "target_q10_retention"
        ]
    )

    breadth = f(
        r[
            "breadth_q10"
        ]
    )

    windows = f(
        r[
            "window_signal_fraction_q10"
        ]
    )

    q20_over_q10 = (
        q20 / q10
        if q10 > 0
        else 0.0
    )


    if (
        q10 >= 1000
        and q10_ret >= 0.50
        and breadth >= 0.50
        and windows >= 0.80
        and q20_over_q10 >= 0.50
    ):

        evidence_class = (
            "strong_distributed_species_support"
        )

        thesis_interpretation = (
            "strong_distributed_chromosomal_support_under_competitive_mapping"
        )


    elif q10 > 0:

        evidence_class = (
            "localized_species_level_support_insufficient"
        )

        thesis_interpretation = (
            "high_MAPQ_target_alignments_present_but_too_localized_for_species_level_validation"
        )


    else:

        evidence_class = (
            "no_high_MAPQ_species_specific_support"
        )

        thesis_interpretation = (
            "screening_signal_not_supported_at_species_level_by_competitive_mapping"
        )


    if (
        r["target_species"]
        == "Escherichia coli"
    ):

        caveat = (
            "E_coli_Shigella_short_read_resolution_limitation"
        )

    elif (
        r["target_species"]
        == "Bacillus cereus"
    ):

        caveat = (
            "B_cereus_group_short_read_species_resolution_limitation"
        )

    else:

        caveat = (
            "single_reference_per_species_and_short_read_limitations"
        )


    final.append({
        "task_id":
            r[
                "task_id"
            ],

        "sample":
            r[
                "sample"
            ],

        "target_species":
            r[
                "target_species"
            ],

        "priority":
            r[
                "priority"
            ],

        "kraken_species_clade_reads":
            k[
                "kraken_species_clade_reads"
            ],

        "kraken_species_direct_reads":
            k[
                "kraken_species_direct_reads"
            ],

        "bracken_species_estimated_reads":
            k[
                "bracken_species_estimated_reads"
            ],

        "bracken_fraction_all_pairs":
            k[
                "bracken_species_fraction_all_pairs"
            ],

        "competitive_target_reads_all":
            r[
                "target_mapped_records_all"
            ],

        "competitive_target_reads_q10":
            r[
                "target_mapped_records_q10"
            ],

        "competitive_target_reads_q20":
            r[
                "target_mapped_records_q20"
            ],

        "q10_retention":
            r[
                "target_q10_retention"
            ],

        "q20_over_q10":
            f"{q20_over_q10:.10f}",

        "breadth_q10":
            r[
                "breadth_q10"
            ],

        "mean_depth_q10":
            r[
                "mean_depth_q10"
            ],

        "window_signal_fraction_q10":
            r[
                "window_signal_fraction_q10"
            ],

        "proper_pair_fraction_target_all":
            r[
                "proper_pair_fraction_target_all"
            ],

        "evidence_class":
            evidence_class,

        "thesis_interpretation":
            thesis_interpretation,

        "caveat":
            caveat,
    })


final.sort(
    key=lambda r:
        int(
            r[
                "task_id"
            ]
        )
)


outfile = (
    OUT
    / "97C2D_final_validation.tsv"
)


with outfile.open(
    "w",
    newline=""
) as fh:

    fields = list(
        final[0].keys()
    )

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(
        final
    )


# ------------------------------------------------------------
# Class summary
# ------------------------------------------------------------

classes = {}

for r in final:

    classes[
        r[
            "evidence_class"
        ]
    ] = (
        classes.get(
            r[
                "evidence_class"
            ],
            0
        )
        + 1
    )


with (
    OUT
    / "97C2D_class_summary.tsv"
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
        "evidence_class",
        "n"
    ])

    for key in sorted(
        classes
    ):

        w.writerow([
            key,
            classes[key]
        ])


# ------------------------------------------------------------
# Rules
# ------------------------------------------------------------

with (
    OUT
    / "97C2D_classification_rules.tsv"
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
        "class",
        "rule"
    ])

    w.writerow([
        "strong_distributed_species_support",
        "Q10>=1000;Q10_retention>=0.50;breadth_Q10>=0.50;window_signal_Q10>=0.80;Q20_over_Q10>=0.50"
    ])

    w.writerow([
        "localized_species_level_support_insufficient",
        "Q10>0_but_strong_distributed_rule_not_met"
    ])

    w.writerow([
        "no_high_MAPQ_species_specific_support",
        "Q10=0"
    ])


# ------------------------------------------------------------
# Methodological scope
# ------------------------------------------------------------

scope = [
    (
        "workflow",
        "97C_targeted_validation_of_health_sensitive_taxonomic_screening"
    ),
    (
        "screening_source",
        "Kraken2_plus_Bracken"
    ),
    (
        "validation_source",
        "competitive_Bowtie2_mapping_against_chromosomal_RefSeq_panels"
    ),
    (
        "primary_specificity_metrics",
        "MAPQ_retention_breadth_and_spatial_distribution"
    ),
    (
        "target_rank_in_panel",
        "not_used_as_presence_or_specificity_criterion"
    ),
    (
        "target_share_of_panel",
        "not_used_as_presence_or_specificity_criterion"
    ),
    (
        "single_reference_limitation",
        "one_complete_RefSeq_genome_per_species"
    ),
    (
        "plasmids",
        "small_replicons_excluded_from_primary_validation"
    ),
    (
        "viability",
        "not_assessed"
    ),
    (
        "pathogenicity",
        "not_assessed"
    ),
    (
        "food_safety_risk",
        "not_inferred"
    ),
    (
        "negative_result",
        "lack_of_species_specific_support_not_proof_of_absence"
    ),
    (
        "amplicon_data_used",
        "NO"
    ),
]


with (
    OUT
    / "97C2D_methodological_scope.tsv"
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
        "field",
        "value"
    ])

    w.writerows(
        scope
    )


# ------------------------------------------------------------
# Human-readable snapshot
# ------------------------------------------------------------

with (
    OUT
    / "97C2D_snapshot.txt"
).open(
    "w"
) as fh:

    fh.write(
        "97C2D FINAL TARGETED PATHOGEN-TAXON VALIDATION\n"
    )

    fh.write(
        "==============================================\n\n"
    )

    for r in final:

        fh.write(
            f"{r['sample']} | "
            f"{r['target_species']} | "
            f"{r['evidence_class']} | "
            f"Q10={r['competitive_target_reads_q10']} | "
            f"breadth_Q10={100*f(r['breadth_q10']):.3f}% | "
            f"depth_Q10={f(r['mean_depth_q10']):.4f} | "
            f"windows={100*f(r['window_signal_fraction_q10']):.2f}%\n"
        )

    fh.write(
        "\nNo class implies viability, pathogenic phenotype, "
        "or food-safety risk.\n"
    )


required = [
    "97C2D_final_validation.tsv",
    "97C2D_class_summary.tsv",
    "97C2D_classification_rules.tsv",
    "97C2D_methodological_scope.tsv",
    "97C2D_snapshot.txt",
]


for name in required:

    p = OUT / name

    if (
        not p.is_file()
        or p.stat().st_size == 0
    ):

        raise RuntimeError(
            f"Missing output: {p}"
        )


print(
    f"PAIRS={len(final)}"
)

for key in sorted(
    classes
):

    print(
        f"{key}={classes[key]}"
    )

print(
    "97C2D=PASS"
)
