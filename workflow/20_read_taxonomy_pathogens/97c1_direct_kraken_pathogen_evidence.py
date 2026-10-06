#!/usr/bin/env python3

import csv
import sys
from pathlib import Path


ROOT = Path(sys.argv[1])

INROOT = (
    ROOT
    / "98_read_taxonomy"
    / "97B1_all_samples"
)

FINAL97 = (
    ROOT
    / "98_read_taxonomy"
    / "97B2C_final_taxonomy"
)

OUT = (
    ROOT
    / "98_read_taxonomy"
    / "97C1_direct_pathogen_evidence"
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


TARGETS = [
    {
        "priority": "tier1",
        "species": "Salmonella enterica",
        "species_taxid": "28901",
        "genus": "Salmonella",
        "genus_taxid": "590",
        "reason": "strongest_unexpected_food_safety_screening_signal",
    },
    {
        "priority": "tier1",
        "species": "Listeria monocytogenes",
        "species_taxid": "1639",
        "genus": "Listeria",
        "genus_taxid": "1637",
        "reason": "directly_related_to_thesis_antagonism_pathogen",
    },
    {
        "priority": "tier1",
        "species": "Staphylococcus aureus",
        "species_taxid": "1280",
        "genus": "Staphylococcus",
        "genus_taxid": "1279",
        "reason": "directly_related_to_thesis_antagonism_pathogen",
    },
    {
        "priority": "tier1",
        "species": "Escherichia coli",
        "species_taxid": "562",
        "genus": "Escherichia",
        "genus_taxid": "561",
        "reason": "food_safety_signal_and_enterobacteriaceae_specificity_control",
    },
    {
        "priority": "tier2",
        "species": "Bacillus cereus",
        "species_taxid": "1396",
        "genus": "Bacillus",
        "genus_taxid": "1386",
        "reason": "food_safety_screening_signal",
    },
    {
        "priority": "tier2",
        "species": "Campylobacter jejuni",
        "species_taxid": "197",
        "genus": "Campylobacter",
        "genus_taxid": "194",
        "reason": "food_safety_screening_signal",
    },
    {
        "priority": "tier2",
        "species": "Clostridium perfringens",
        "species_taxid": "1502",
        "genus": "Clostridium",
        "genus_taxid": "1485",
        "reason": "food_safety_screening_signal",
    },
    {
        "priority": "tier2",
        "species": "Cronobacter sakazakii",
        "species_taxid": "28141",
        "genus": "Cronobacter",
        "genus_taxid": "413496",
        "reason": "food_safety_screening_signal",
    },
    {
        "priority": "tier2",
        "species": "Yersinia enterocolitica",
        "species_taxid": "630",
        "genus": "Yersinia",
        "genus_taxid": "629",
        "reason": "food_safety_screening_signal",
    },
]


def read_tsv(path):

    if not path.is_file():
        raise RuntimeError(
            f"Missing input: {path}"
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


def intnum(x):

    try:
        return int(float(x))
    except Exception:
        return 0


def floatnum(x):

    try:
        return float(x)
    except Exception:
        return 0.0


def parse_kraken_report(path):

    records = {}

    with path.open(
        errors="replace"
    ) as fh:

        for line in fh:

            parts = line.rstrip(
                "\n"
            ).split("\t")

            if len(parts) < 6:
                continue

            try:
                pct = float(
                    parts[0]
                )

                clade_reads = int(
                    parts[1]
                )

                direct_reads = int(
                    parts[2]
                )

            except ValueError:
                continue


            rank_code = parts[3].strip()
            taxid = parts[4].strip()
            name = parts[5].strip()


            records[taxid] = {
                "pct_report":
                    pct,

                "clade_reads":
                    clade_reads,

                "direct_reads":
                    direct_reads,

                "rank_code":
                    rank_code,

                "name":
                    name,
            }

    return records


def bracken_lookup(
    path
):

    rows = read_tsv(
        path
    )

    return {
        str(
            r[
                "taxonomy_id"
            ]
        ).strip():
            r

        for r in rows
    }


# ============================================================
# Input validation
# ============================================================

missing = []


for sample in SAMPLES:

    sdir = INROOT / sample

    required = [
        sdir
        / f"{sample}.kraken2.report",

        sdir
        / f"{sample}.classification_summary.tsv",

        sdir
        / f"{sample}.bracken.genus.augmented.tsv",

        sdir
        / f"{sample}.bracken.species.augmented.tsv",
    ]

    for p in required:

        if (
            not p.is_file()
            or p.stat().st_size == 0
        ):

            missing.append(
                str(p)
            )


if missing:

    raise RuntimeError(
        "Missing required 97B1 files:\n"
        + "\n".join(
            missing
        )
    )


if not (
    FINAL97
    / "97B2C_food_safety_sensitive_screening.tsv"
).is_file():

    raise RuntimeError(
        "97B2C final screening table missing"
    )


# ============================================================
# Save target definition
# ============================================================

write_tsv(
    OUT
    / "97C1_target_definition.tsv",
    TARGETS,
    [
        "priority",
        "species",
        "species_taxid",
        "genus",
        "genus_taxid",
        "reason",
    ]
)


# ============================================================
# Main evidence table
# ============================================================

results = []


for sample in SAMPLES:

    sdir = INROOT / sample

    summary_rows = read_tsv(
        sdir
        / f"{sample}.classification_summary.tsv"
    )


    summary = {
        r["metric"]:
            r["value"]

        for r in summary_rows
    }


    total_pairs = intnum(
        summary[
            "total_sequence_pairs"
        ]
    )


    if total_pairs <= 0:

        raise RuntimeError(
            f"Invalid total pairs for {sample}"
        )


    report = parse_kraken_report(
        sdir
        / f"{sample}.kraken2.report"
    )


    species_bracken = bracken_lookup(
        sdir
        / f"{sample}.bracken.species.augmented.tsv"
    )


    genus_bracken = bracken_lookup(
        sdir
        / f"{sample}.bracken.genus.augmented.tsv"
    )


    for target in TARGETS:

        stax = target[
            "species_taxid"
        ]

        gtax = target[
            "genus_taxid"
        ]


        sr = report.get(
            stax,
            {
                "pct_report": 0.0,
                "clade_reads": 0,
                "direct_reads": 0,
                "rank_code": "",
                "name": "",
            }
        )


        gr = report.get(
            gtax,
            {
                "pct_report": 0.0,
                "clade_reads": 0,
                "direct_reads": 0,
                "rank_code": "",
                "name": "",
            }
        )


        sb = species_bracken.get(
            stax,
            {}
        )


        gb = genus_bracken.get(
            gtax,
            {}
        )


        species_clade = intnum(
            sr[
                "clade_reads"
            ]
        )

        species_direct = intnum(
            sr[
                "direct_reads"
            ]
        )

        genus_clade = intnum(
            gr[
                "clade_reads"
            ]
        )

        genus_direct = intnum(
            gr[
                "direct_reads"
            ]
        )


        bracken_species = intnum(
            sb.get(
                "new_est_reads",
                0
            )
        )

        bracken_genus = intnum(
            gb.get(
                "new_est_reads",
                0
            )
        )


        species_clade_fraction = (
            species_clade
            / total_pairs
        )

        species_direct_fraction = (
            species_direct
            / total_pairs
        )

        genus_clade_fraction = (
            genus_clade
            / total_pairs
        )

        bracken_species_fraction = (
            bracken_species
            / total_pairs
        )

        bracken_genus_fraction = (
            bracken_genus
            / total_pairs
        )


        if bracken_species > 0:

            kraken_clade_to_bracken_ratio = (
                species_clade
                / bracken_species
            )

        else:

            kraken_clade_to_bracken_ratio = (
                0.0
            )


        if genus_clade > 0:

            species_share_of_genus_clade = (
                species_clade
                / genus_clade
            )

        else:

            species_share_of_genus_clade = (
                0.0
            )


        results.append({
            "sample":
                sample,

            "priority":
                target[
                    "priority"
                ],

            "species":
                target[
                    "species"
                ],

            "species_taxid":
                stax,

            "genus":
                target[
                    "genus"
                ],

            "genus_taxid":
                gtax,

            "total_sequence_pairs":
                total_pairs,

            "kraken_species_row_present":
                (
                    "YES"
                    if stax in report
                    else "NO"
                ),

            "kraken_genus_row_present":
                (
                    "YES"
                    if gtax in report
                    else "NO"
                ),

            "kraken_species_clade_reads":
                species_clade,

            "kraken_species_direct_reads":
                species_direct,

            "kraken_species_clade_fraction_all_pairs":
                f"{species_clade_fraction:.10f}",

            "kraken_species_direct_fraction_all_pairs":
                f"{species_direct_fraction:.10f}",

            "kraken_genus_clade_reads":
                genus_clade,

            "kraken_genus_direct_reads":
                genus_direct,

            "kraken_genus_clade_fraction_all_pairs":
                f"{genus_clade_fraction:.10f}",

            "bracken_species_estimated_reads":
                bracken_species,

            "bracken_species_fraction_all_pairs":
                f"{bracken_species_fraction:.10f}",

            "bracken_genus_estimated_reads":
                bracken_genus,

            "bracken_genus_fraction_all_pairs":
                f"{bracken_genus_fraction:.10f}",

            "kraken_species_clade_to_bracken_species_ratio":
                f"{kraken_clade_to_bracken_ratio:.6f}",

            "kraken_species_clade_fraction_of_genus_clade":
                f"{species_share_of_genus_clade:.6f}",

            "interpretation":
                "classifier_support_only_not_pathogen_confirmation",
        })


fields = list(
    results[0].keys()
)


write_tsv(
    OUT
    / "97C1_all_target_sample_evidence.tsv",
    results,
    fields
)


# ============================================================
# Top samples per species
# ============================================================

top_rows = []


for target in TARGETS:

    species = target[
        "species"
    ]


    rr = [
        r
        for r in results
        if r[
            "species"
        ] == species
    ]


    rr.sort(
        key=lambda r:
            floatnum(
                r[
                    "kraken_species_clade_fraction_all_pairs"
                ]
            ),
        reverse=True
    )


    for rank, r in enumerate(
        rr[:6],
        1
    ):

        rec = dict(
            r
        )

        rec[
            "rank_by_kraken_species_clade_fraction"
        ] = rank

        top_rows.append(
            rec
        )


write_tsv(
    OUT
    / "97C1_top_samples_per_target.tsv",
    top_rows,
    [
        "rank_by_kraken_species_clade_fraction"
    ]
    + fields
)


# ============================================================
# Tier 1 focused table
# ============================================================

tier1 = [
    r
    for r in results
    if r[
        "priority"
    ] == "tier1"
]


tier1.sort(
    key=lambda r: (
        r["species"],
        -floatnum(
            r[
                "kraken_species_clade_fraction_all_pairs"
            ]
        ),
    )
)


write_tsv(
    OUT
    / "97C1_tier1_evidence.tsv",
    tier1,
    fields
)


# ============================================================
# Comparison: Bracken versus direct Kraken evidence
# ============================================================

comparison = []


for r in results:

    br = floatnum(
        r[
            "bracken_species_fraction_all_pairs"
        ]
    )

    kr = floatnum(
        r[
            "kraken_species_clade_fraction_all_pairs"
        ]
    )


    comparison.append({
        "sample":
            r[
                "sample"
            ],

        "priority":
            r[
                "priority"
            ],

        "species":
            r[
                "species"
            ],

        "bracken_species_fraction_all_pairs":
            f"{br:.10f}",

        "kraken_species_clade_fraction_all_pairs":
            f"{kr:.10f}",

        "absolute_difference_fraction":
            f"{(br - kr):.10f}",

        "kraken_clade_to_bracken_ratio":
            r[
                "kraken_species_clade_to_bracken_species_ratio"
            ],

        "note":
            "lower_Kraken_than_Bracken_can_reflect_Bracken_redistribution_not_false_positive_proof",
    })


comparison.sort(
    key=lambda r:
        floatnum(
            r[
                "bracken_species_fraction_all_pairs"
            ]
        ),
    reverse=True
)


write_tsv(
    OUT
    / "97C1_bracken_vs_kraken.tsv",
    comparison,
    list(
        comparison[0].keys()
    )
)


# ============================================================
# Target summary
# ============================================================

target_summary = []


for target in TARGETS:

    species = target[
        "species"
    ]


    rr = [
        r
        for r in results
        if r[
            "species"
        ] == species
    ]


    rr_sorted = sorted(
        rr,
        key=lambda r:
            floatnum(
                r[
                    "kraken_species_clade_fraction_all_pairs"
                ]
            ),
        reverse=True
    )


    top = rr_sorted[0]


    positive_rows = sum(
        intnum(
            r[
                "kraken_species_clade_reads"
            ]
        ) > 0
        for r in rr
    )


    direct_positive_rows = sum(
        intnum(
            r[
                "kraken_species_direct_reads"
            ]
        ) > 0
        for r in rr
    )


    target_summary.append({
        "priority":
            target[
                "priority"
            ],

        "species":
            species,

        "species_taxid":
            target[
                "species_taxid"
            ],

        "samples_with_kraken_species_clade_signal":
            positive_rows,

        "samples_with_kraken_species_direct_signal":
            direct_positive_rows,

        "top_sample":
            top[
                "sample"
            ],

        "top_kraken_species_clade_reads":
            top[
                "kraken_species_clade_reads"
            ],

        "top_kraken_species_clade_fraction_all_pairs":
            top[
                "kraken_species_clade_fraction_all_pairs"
            ],

        "top_bracken_species_estimated_reads":
            top[
                "bracken_species_estimated_reads"
            ],

        "top_bracken_species_fraction_all_pairs":
            top[
                "bracken_species_fraction_all_pairs"
            ],

        "top_kraken_clade_to_bracken_ratio":
            top[
                "kraken_species_clade_to_bracken_species_ratio"
            ],

        "interpretation":
            "descriptive_classifier_evidence_only",
    })


write_tsv(
    OUT
    / "97C1_target_summary.tsv",
    target_summary,
    list(
        target_summary[0].keys()
    )
)


# ============================================================
# Global summary
# ============================================================

with (
    OUT
    / "97C1_global_summary.tsv"
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
        "value",
    ])

    w.writerow([
        "samples",
        len(
            SAMPLES
        ),
    ])

    w.writerow([
        "target_species",
        len(
            TARGETS
        ),
    ])

    w.writerow([
        "target_sample_combinations",
        len(
            results
        ),
    ])

    w.writerow([
        "tier1_species",
        4,
    ])

    w.writerow([
        "tier2_species",
        5,
    ])

    w.writerow([
        "Kraken_reclassification_performed",
        "NO",
    ])

    w.writerow([
        "Bracken_reclassification_performed",
        "NO",
    ])

    w.writerow([
        "existing_97B1_kraken_reports_reused",
        "YES",
    ])

    w.writerow([
        "primary_direct_metric",
        "Kraken_species_clade_reads",
    ])

    w.writerow([
        "pathogen_confirmation",
        "NO",
    ])

    w.writerow([
        "next_step",
        "select_high_priority_sample_target_pairs_for_high_specificity_validation",
    ])


# ============================================================
# Methodological guardrails
# ============================================================

scope = [
    (
        "analysis",
        "targeted_food_safety_taxon_classifier_validation"
    ),
    (
        "input",
        "existing_Kraken2_reports_and_Bracken_outputs_from_97B1"
    ),
    (
        "Kraken_species_clade_reads",
        "reads_assigned_to_species_taxon_or_descendant_taxa_in_Kraken_report"
    ),
    (
        "Kraken_species_direct_reads",
        "reads_assigned_directly_to_species_taxon_only"
    ),
    (
        "Bracken_estimated_reads",
        "redistributed_estimate_not_direct_assignment"
    ),
    (
        "species_clade_fraction",
        "descriptive_fraction_of_all_input_sequence_pairs"
    ),
    (
        "zero_signal",
        "no_classifier_signal_not_proof_of_biological_absence"
    ),
    (
        "positive_signal",
        "classifier_support_not_confirmation_of_viable_pathogenic_organism"
    ),
    (
        "species_label",
        "requires_high_specificity_validation_before_pathogen_claim"
    ),
    (
        "amplicon_data_used",
        "NO"
    ),
]


with (
    OUT
    / "97C1_methodological_scope.tsv"
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
        "value",
    ])

    w.writerows(
        scope
    )


# ============================================================
# Validation
# ============================================================

required = [
    "97C1_target_definition.tsv",
    "97C1_all_target_sample_evidence.tsv",
    "97C1_top_samples_per_target.tsv",
    "97C1_tier1_evidence.tsv",
    "97C1_bracken_vs_kraken.tsv",
    "97C1_target_summary.tsv",
    "97C1_global_summary.tsv",
    "97C1_methodological_scope.tsv",
]


for f in required:

    p = OUT / f

    if (
        not p.is_file()
        or p.stat().st_size == 0
    ):

        raise RuntimeError(
            f"Missing/empty output: {f}"
        )


if len(results) != (
    len(SAMPLES)
    * len(TARGETS)
):

    raise RuntimeError(
        f"Expected 162 rows, got {len(results)}"
    )


print(
    f"SAMPLES={len(SAMPLES)}"
)

print(
    f"TARGETS={len(TARGETS)}"
)

print(
    f"ROWS={len(results)}"
)

print(
    "KRAKEN_RERUN=NO"
)

print(
    "BRACKEN_RERUN=NO"
)

print(
    "97C1=PASS"
)
