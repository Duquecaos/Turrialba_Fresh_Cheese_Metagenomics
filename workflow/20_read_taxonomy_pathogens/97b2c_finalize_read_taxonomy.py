#!/usr/bin/env python3

import csv
import sys
from pathlib import Path
from collections import defaultdict
from statistics import median, mean


ROOT = Path(sys.argv[1])

A = (
    ROOT
    / "98_read_taxonomy"
    / "97B2A_matrix_QC"
)

B = (
    ROOT
    / "98_read_taxonomy"
    / "97B2B_bacterial_community"
)

R = (
    ROOT
    / "98_read_taxonomy"
    / "97B2BR_statistics_correction"
)

OUT = (
    ROOT
    / "98_read_taxonomy"
    / "97B2C_final_taxonomy"
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


def read_tsv(path):

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

        w = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )

        w.writeheader()
        w.writerows(rows)


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
# 1. Load key inputs
# ============================================================

genus = read_tsv(
    A
    / "97B2A_genus_master.tsv"
)

species = read_tsv(
    A
    / "97B2A_species_master.tsv"
)

qc = read_tsv(
    A
    / "97B2A_sample_QC_summary.tsv"
)

stats_summary = {
    r["metric"]:
        r["value"]

    for r in read_tsv(
        R
        / "97B2BR_global_summary.tsv"
    )
}

alpha = read_tsv(
    B
    / "97B2B_alpha_diversity_primary.tsv"
)

axis = {
    r["axis"]:
        fnum(
            r[
                "pct_positive_eigenvalue_variation"
            ]
        )

    for r in read_tsv(
        R
        / "97B2BR_PCoA_axis_variance.tsv"
    )
}


# ============================================================
# 2. Primary genus set
# ============================================================

primary_matrix = read_tsv(
    B
    / "97B2B_bacterial_genus_matrix_primary.tsv"
)

primary_taxids = {
    r["taxonomy_id"]
    for r in primary_matrix
}


# ============================================================
# 3. Final top bacterial genera per sample
# ============================================================

top_sample = []


for sample in SAMPLES:

    rr = [
        r
        for r in genus
        if (
            r["sample"] == sample
            and
            r["taxonomy_category"] == "Bacteria"
            and
            r["taxonomy_id"] in primary_taxids
        )
    ]

    rr.sort(
        key=lambda r:
            fnum(
                r[
                    "relative_within_bacteria"
                ]
            ),
        reverse=True
    )


    for rank, r in enumerate(
        rr[:15],
        1
    ):

        top_sample.append({
            "sample":
                sample,

            "producer_code":
                r["producer_code"],

            "biological_unit_code":
                r["biological_unit_code"],

            "time_code":
                r["within_unit_time_code"],

            "rank":
                rank,

            "genus":
                r["name"],

            "taxonomy_id":
                r["taxonomy_id"],

            "relative_within_bacteria":
                r["relative_within_bacteria"],

            "fraction_all_input_pairs":
                r["fraction_all_input_pairs"],

            "interpretation":
                "Bracken_estimated_read_composition",
        })


write_tsv(
    OUT
    / "97B2C_top15_bacterial_genera_by_sample.tsv",
    top_sample,
    list(
        top_sample[0].keys()
    )
)


# ============================================================
# 4. Producer x time-code descriptive summaries
#
# Median across the 3 biological units.
# No inferential producer x time interaction implied.
# ============================================================

lookup = defaultdict(float)


for r in genus:

    if (
        r["taxonomy_category"] == "Bacteria"
        and
        r["taxonomy_id"] in primary_taxids
    ):

        key = (
            r["producer_code"],
            r["within_unit_time_code"],
            r["sample"],
            r["taxonomy_id"],
            r["name"],
        )

        lookup[key] += fnum(
            r[
                "relative_within_bacteria"
            ]
        )


group_rows = []


for producer in ["L", "M"]:

    for time_code in ["1", "2", "3"]:

        samples = [
            f"{producer}1_{time_code}",
            f"{producer}2_{time_code}",
            f"{producer}3_{time_code}",
        ]


        taxa = set()

        for key in lookup:

            p, t, sample, taxid, name = key

            if (
                p == producer
                and
                t == time_code
            ):

                taxa.add(
                    (
                        taxid,
                        name,
                    )
                )


        temp = []

        for taxid, name in taxa:

            vals = [
                lookup.get(
                    (
                        producer,
                        time_code,
                        sample,
                        taxid,
                        name,
                    ),
                    0.0
                )
                for sample in samples
            ]


            temp.append({
                "producer_code":
                    producer,

                "time_code":
                    time_code,

                "taxonomy_id":
                    taxid,

                "genus":
                    name,

                "median_relative_within_bacteria":
                    median(
                        vals
                    ),

                "mean_relative_within_bacteria":
                    mean(
                        vals
                    ),

                "minimum_relative_within_bacteria":
                    min(
                        vals
                    ),

                "maximum_relative_within_bacteria":
                    max(
                        vals
                    ),

                "samples_with_nonzero_signal":
                    sum(
                        v > 0
                        for v in vals
                    ),
            })


        temp.sort(
            key=lambda r:
                r[
                    "median_relative_within_bacteria"
                ],
            reverse=True
        )


        for rank, r in enumerate(
            temp[:15],
            1
        ):

            r["rank_by_median"] = rank

            group_rows.append(
                r
            )


write_tsv(
    OUT
    / "97B2C_top15_bacterial_genera_by_producer_timecode.tsv",
    group_rows,
    [
        "producer_code",
        "time_code",
        "rank_by_median",
        "taxonomy_id",
        "genus",
        "median_relative_within_bacteria",
        "mean_relative_within_bacteria",
        "minimum_relative_within_bacteria",
        "maximum_relative_within_bacteria",
        "samples_with_nonzero_signal",
    ]
)


# ============================================================
# 5. Fungal descriptive summary
# ============================================================

fungal = [
    r
    for r in genus
    if r[
        "taxonomy_category"
    ] == "Fungi"
]


fungal.sort(
    key=lambda r: (
        SAMPLES.index(
            r["sample"]
        ),
        -fnum(
            r[
                "fraction_all_input_pairs"
            ]
        ),
    )
)


fungal_top = []


for sample in SAMPLES:

    rr = [
        r
        for r in fungal
        if r["sample"] == sample
    ]


    for rank, r in enumerate(
        rr[:10],
        1
    ):

        fungal_top.append({
            "sample":
                sample,

            "rank":
                rank,

            "genus":
                r["name"],

            "taxonomy_id":
                r["taxonomy_id"],

            "fraction_all_input_pairs":
                r["fraction_all_input_pairs"],

            "interpretation":
                "descriptive_read_based_fungal_signal",
        })


write_tsv(
    OUT
    / "97B2C_top10_fungal_genera_by_sample.tsv",
    fungal_top,
    list(
        fungal_top[0].keys()
    )
)


# ============================================================
# 6. Food-safety-sensitive screening signals
#
# These are flags for targeted validation only.
# They are NOT confirmed pathogen detections.
# ============================================================

GENUS_FLAGS = {
    "Salmonella",
    "Listeria",
    "Staphylococcus",
    "Escherichia",
    "Cronobacter",
    "Bacillus",
    "Clostridium",
    "Campylobacter",
    "Yersinia",
}


SPECIES_FLAGS = {
    "Salmonella enterica",
    "Listeria monocytogenes",
    "Staphylococcus aureus",
    "Escherichia coli",
    "Cronobacter sakazakii",
    "Bacillus cereus",
    "Clostridium perfringens",
    "Campylobacter jejuni",
    "Yersinia enterocolitica",
}


screen = []


for r in genus:

    if (
        r["name"].strip()
        in GENUS_FLAGS
    ):

        screen.append({
            "sample":
                r["sample"],

            "rank":
                "genus",

            "taxon":
                r["name"].strip(),

            "taxonomy_id":
                r["taxonomy_id"],

            "estimated_reads":
                r["new_est_reads"],

            "fraction_all_input_pairs":
                r["fraction_all_input_pairs"],

            "status":
                "screening_signal_only",

            "required_interpretation":
                "genus_assignment_does_not_confirm_pathogenic_species_or_food_safety_risk",
        })


for r in species:

    name = r[
        "name"
    ].strip()

    if (
        name
        in SPECIES_FLAGS
    ):

        screen.append({
            "sample":
                r["sample"],

            "rank":
                "species",

            "taxon":
                name,

            "taxonomy_id":
                r["taxonomy_id"],

            "estimated_reads":
                r["new_est_reads"],

            "fraction_all_input_pairs":
                r["fraction_all_input_pairs"],

            "status":
                "screening_signal_only",

            "required_interpretation":
                "Bracken_species_assignment_requires_targeted_validation_before_pathogen_claim",
        })


screen.sort(
    key=lambda r: (
        r["taxon"],
        -fnum(
            r[
                "fraction_all_input_pairs"
            ]
        ),
        r["sample"],
    )
)


write_tsv(
    OUT
    / "97B2C_food_safety_sensitive_screening.tsv",
    screen,
    [
        "sample",
        "rank",
        "taxon",
        "taxonomy_id",
        "estimated_reads",
        "fraction_all_input_pairs",
        "status",
        "required_interpretation",
    ]
)


# ============================================================
# 7. QC summary
# ============================================================

qc_final = []


for r in qc:

    qc_final.append({
        "sample":
            r["sample"],

        "classified_pct":
            r["kraken_classified_pct"],

        "bacteria_pct_all_pairs":
            r["bacteria_pct_all_pairs"],

        "fungi_pct_all_pairs":
            r["fungi_pct_all_pairs"],

        "archaea_pct_all_pairs":
            r["archaea_pct_all_pairs"],

        "viruses_pct_all_pairs":
            r["viruses_pct_all_pairs"],

        "human_pct_all_pairs":
            r["human_pct_all_pairs"],

        "unclassified_pct":
            f"{100.0 - fnum(r['kraken_classified_pct']):.8f}",
    })


write_tsv(
    OUT
    / "97B2C_sample_classification_QC.tsv",
    qc_final,
    list(
        qc_final[0].keys()
    )
)


# ============================================================
# 8. Alpha diversity final copy
# ============================================================

write_tsv(
    OUT
    / "97B2C_alpha_diversity.tsv",
    alpha,
    list(
        alpha[0].keys()
    )
)


# ============================================================
# 9. Final frozen statistical results
# ============================================================

final_stats = [
    {
        "result":
            "time_subject_adjusted",

        "R2":
            stats_summary[
                "time_subject_adjusted_R2"
            ],

        "F":
            stats_summary[
                "time_subject_adjusted_F"
            ],

        "P":
            stats_summary[
                "time_subject_adjusted_P"
            ],

        "interpretation":
            "significant_time_associated_compositional_change_after_adjusting_for_biological_unit",
    },

    {
        "result":
            "time_within_L",

        "R2":
            stats_summary[
                "time_within_L_R2"
            ],

        "F":
            "NA",

        "P":
            stats_summary[
                "time_within_L_P"
            ],

        "interpretation":
            "no_detectable_time_effect_within_L",
    },

    {
        "result":
            "time_within_M",

        "R2":
            stats_summary[
                "time_within_M_R2"
            ],

        "F":
            "NA",

        "P":
            stats_summary[
                "time_within_M_P"
            ],

        "interpretation":
            "detectable_time_associated_compositional_change_within_M",
    },

    {
        "result":
            "producer_subject_level",

        "R2":
            stats_summary[
                "producer_subject_level_R2"
            ],

        "F":
            stats_summary[
                "producer_subject_level_F"
            ],

        "P":
            stats_summary[
                "producer_subject_level_P"
            ],

        "interpretation":
            "large_descriptive_centroid_difference_not_statistically_supported_with_n3_units_per_producer_and_unequal_dispersion",
    },
]


write_tsv(
    OUT
    / "97B2C_final_inferential_statistics.tsv",
    final_stats,
    [
        "result",
        "R2",
        "F",
        "P",
        "interpretation",
    ]
)


# ============================================================
# 10. Global summary
# ============================================================

global_rows = [
    (
        "shotgun_samples",
        "18"
    ),
    (
        "independent_biological_units",
        "6"
    ),
    (
        "units_per_producer",
        "3"
    ),
    (
        "raw_bacterial_genera",
        "1815"
    ),
    (
        "primary_bacterial_genera",
        "176"
    ),
    (
        "primary_filter",
        "genus_ge0.01pct_within_bacteria_in_ge2_samples"
    ),
    (
        "minimum_composition_preserved_pct",
        "95.895494"
    ),
    (
        "maximum_composition_preserved_pct",
        "99.951916"
    ),
    (
        "time_subject_adjusted_R2",
        stats_summary[
            "time_subject_adjusted_R2"
        ]
    ),
    (
        "time_subject_adjusted_P",
        stats_summary[
            "time_subject_adjusted_P"
        ]
    ),
    (
        "time_within_L_R2",
        stats_summary[
            "time_within_L_R2"
        ]
    ),
    (
        "time_within_L_P",
        stats_summary[
            "time_within_L_P"
        ]
    ),
    (
        "time_within_M_R2",
        stats_summary[
            "time_within_M_R2"
        ]
    ),
    (
        "time_within_M_P",
        stats_summary[
            "time_within_M_P"
        ]
    ),
    (
        "producer_R2",
        stats_summary[
            "producer_subject_level_R2"
        ]
    ),
    (
        "producer_P",
        stats_summary[
            "producer_subject_level_P"
        ]
    ),
    (
        "PERMDISP_producer_P",
        "0.001389"
    ),
    (
        "PCoA1_pct_positive_variation",
        f"{axis['PCoA1']:.8f}"
    ),
    (
        "PCoA2_pct_positive_variation",
        f"{axis['PCoA2']:.8f}"
    ),
    (
        "PCoA1_plus_PCoA2_pct",
        f"{axis['PCoA1'] + axis['PCoA2']:.8f}"
    ),
    (
        "primary_taxonomic_rank",
        "genus"
    ),
    (
        "species_level_role",
        "secondary_screening_only"
    ),
    (
        "amplicon_data_used",
        "NO"
    ),
]


with (
    OUT
    / "97B2C_global_summary.tsv"
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

    w.writerows(
        global_rows
    )


# ============================================================
# 11. Thesis-safe snapshot
# ============================================================

snapshot = f"""
97B2C - FINAL READ-BASED SHOTGUN TAXONOMY

Scope
-----
18 shotgun metagenomes.
No amplicon data used.
Primary inferential level: bacterial genus.
Classifier: Kraken2 + Bracken, PlusPF 2026-06-26.

Filtering
---------
1815 bacterial genera were detected before prevalence/abundance filtering.
The primary filter retained 176 genera:
>=0.01% of bacterial composition in >=2 samples.
This retained 95.895494-99.951916% of bacterial composition per sample.

Community structure
-------------------
Biological unit was the dominant source of compositional variation.

After adjusting for biological unit, time code explained
{stats_summary['time_subject_adjusted_R2']} of Bray-Curtis variation
(P={stats_summary['time_subject_adjusted_P']}).

Within producer L:
R2={stats_summary['time_within_L_R2']}
P={stats_summary['time_within_L_P']}.

Within producer M:
R2={stats_summary['time_within_M_R2']}
P={stats_summary['time_within_M_P']}.

Producer comparison
-------------------
At the level of the six independent biological units,
producer explained R2={stats_summary['producer_subject_level_R2']}
with P={stats_summary['producer_subject_level_P']}.

This must not be interpreted as a statistically demonstrated
producer centroid effect because only 3 independent units were
available per producer and multivariate dispersion differed strongly
between producers (PERMDISP P=0.001389).

PCoA
----
PCoA1: {axis['PCoA1']:.3f}% of positive-eigenvalue variation.
PCoA2: {axis['PCoA2']:.3f}%.
Combined: {axis['PCoA1'] + axis['PCoA2']:.3f}%.

Alpha diversity
---------------
No statistically supported global temporal change was detected
for Shannon or Simpson-family metrics using Friedman repeated-measure
tests.

Fungi
-----
Fungal read-based signals were retained as a separate descriptive
analysis and were not mixed into the primary bacterial PERMANOVA.
The strong Debaryomyces signal in L3_3 is therefore descriptive
shotgun evidence rather than part of the bacterial inference.

Food-safety-sensitive assignments
---------------------------------
Kraken2/Bracken assignments to taxa such as Salmonella, Listeria,
Staphylococcus, Escherichia or other potentially relevant groups are
screening signals only. A genus or species label is not equivalent to
confirmation of a pathogenic organism or food-safety risk.

Interpretation guardrail
------------------------
Bracken values represent estimated read composition. They are not
absolute cell abundance, viable counts, expression, or phenotype.

STATUS: READ-BASED SHOTGUN TAXONOMY CLOSED.
""".strip() + "\n"


(
    OUT
    / "97B2C_thesis_results_snapshot.txt"
).write_text(
    snapshot
)


# ============================================================
# 12. Frozen methodological scope
# ============================================================

scope = [
    (
        "analysis_status",
        "CLOSED"
    ),
    (
        "primary_data",
        "shotgun_host_removed_paired_reads"
    ),
    (
        "amplicon_data_used",
        "NO"
    ),
    (
        "classifier",
        "Kraken2_plus_Bracken"
    ),
    (
        "reference_database",
        "PlusPF_2026-06-26"
    ),
    (
        "primary_rank",
        "genus"
    ),
    (
        "primary_community",
        "bacteria"
    ),
    (
        "primary_filter",
        "ge0.01pct_within_bacteria_in_ge2_samples"
    ),
    (
        "beta_diversity",
        "Bray-Curtis"
    ),
    (
        "repeated_measure_unit",
        "L1_L2_L3_M1_M2_M3"
    ),
    (
        "time_model",
        "subject_adjusted_with_restricted_permutations"
    ),
    (
        "producer_model",
        "mean_profile_per_independent_biological_unit"
    ),
    (
        "producer_independent_n",
        "3_per_producer"
    ),
    (
        "producer_dispersion",
        "significantly_heterogeneous"
    ),
    (
        "formal_producer_time_interaction",
        "NOT_TESTED"
    ),
    (
        "fungi",
        "descriptive_separate_analysis"
    ),
    (
        "archaea",
        "descriptive_very_low_signal"
    ),
    (
        "viruses",
        "dedicated_geNomad_vOTU_workflow_remains_primary"
    ),
    (
        "human",
        "excluded_from_microbiota"
    ),
    (
        "health_sensitive_taxa",
        "screening_only_requires_targeted_validation_for_pathogen_claim"
    ),
    (
        "relative_abundance_meaning",
        "estimated_read_composition_not_cell_abundance"
    ),
]


with (
    OUT
    / "97B2C_methodological_scope.tsv"
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


# ============================================================
# 13. Validation
# ============================================================

required = [
    "97B2C_top15_bacterial_genera_by_sample.tsv",
    "97B2C_top15_bacterial_genera_by_producer_timecode.tsv",
    "97B2C_top10_fungal_genera_by_sample.tsv",
    "97B2C_food_safety_sensitive_screening.tsv",
    "97B2C_sample_classification_QC.tsv",
    "97B2C_alpha_diversity.tsv",
    "97B2C_final_inferential_statistics.tsv",
    "97B2C_global_summary.tsv",
    "97B2C_thesis_results_snapshot.txt",
    "97B2C_methodological_scope.tsv",
]


for name in required:

    p = OUT / name

    if (
        not p.exists()
        or
        p.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing/empty output: {name}"
        )


print("SAMPLES=18")
print("PRIMARY_BACTERIAL_GENERA=176")
print(
    f"HEALTH_SCREEN_ROWS={len(screen)}"
)
print(
    f"PCOA12={axis['PCoA1'] + axis['PCoA2']:.6f}"
)
print("AMPLICON_USED=NO")
print("97B2C=PASS")
