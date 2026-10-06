#!/usr/bin/env python3

import csv
import sys
from pathlib import Path


# script + 6 arguments
if len(sys.argv) != 7:
    raise SystemExit(
        "Usage: 98e_final_bacteriocin_synthesis.py "
        "<98D1C_reference_summary.tsv> "
        "<98D1C_global_summary.tsv> "
        "<98C_FINAL_global_summary.tsv> "
        "<98C_FINAL_recurrence.tsv> "
        "<ATTRLOC_summary.tsv> "
        "<OUT>"
    )


REFSUMMARY = Path(sys.argv[1])
D1GLOBAL = Path(sys.argv[2])
CGLOBAL = Path(sys.argv[3])
RECURRENCE = Path(sys.argv[4])
ATTRLOC = Path(sys.argv[5])
OUT = Path(sys.argv[6])

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# Helpers
# ============================================================

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


def metric_table(path):

    result = {}

    with path.open(
        newline=""
    ) as fh:

        for row in csv.DictReader(
            fh,
            delimiter="\t"
        ):

            result[
                row["metric"]
            ] = row["value"]

    return result


def write_tsv(
    path,
    rows,
    fields=None
):

    with path.open(
        "w",
        newline=""
    ) as fh:

        if fields is None:

            fields = (
                list(rows[0].keys())
                if rows
                else ["status"]
            )

        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore"
        )

        writer.writeheader()

        for row in rows:
            writer.writerow(row)


def i(value):

    try:
        return int(float(value))

    except Exception:
        return 0


# ============================================================
# Inputs
# ============================================================

refs = read_tsv(
    REFSUMMARY
)

d1 = metric_table(
    D1GLOBAL
)

cg = metric_table(
    CGLOBAL
)

recurrences = read_tsv(
    RECURRENCE
)

attrloc = read_tsv(
    ATTRLOC
)


if len(refs) != 30:

    raise RuntimeError(
        f"Expected 30 mapped contexts; found {len(refs)}"
    )


if len(attrloc) != 7:

    raise RuntimeError(
        f"Expected 7 ATTRLOC rows; found {len(attrloc)}"
    )


attr_by_id = {
    r["locus_id"]: r
    for r in attrloc
}


# ============================================================
# 1. Master catalog
# ============================================================

catalog = []


for r in refs:

    rid = r[
        "reference_id"
    ]

    rclass = r[
        "reference_class"
    ]

    robust = i(
        r[
            "samples_robust"
        ]
    )

    mapq_sensitive = i(
        r[
            "samples_high_breadth_MAPQ_sensitive"
        ]
    )

    partial = i(
        r[
            "samples_partial"
        ]
    )

    weak = i(
        r[
            "samples_weak_or_localized"
        ]
    )

    no_signal = i(
        r[
            "samples_no_signal"
        ]
    )


    # --------------------------------------------------------
    # Existing ATTRLOC
    # --------------------------------------------------------

    if rclass == "existing_ATTRLOC":

        a = attr_by_id[
            rid
        ]


        if robust > 0:

            evidence_pattern = (
                "existing_ATTRLOC_with_robust_competitive_recruitment"
            )

        elif mapq_sensitive > 0:

            evidence_pattern = (
                "existing_ATTRLOC_recruitment_limited_by_MAPQ_ambiguity"
            )

        else:

            evidence_pattern = (
                "existing_ATTRLOC_without_robust_competitive_recruitment"
            )


        prediction_support = (
            "existing_candidate_framework_98B"
        )


        context_caution_parts = []


        if i(
            a.get(
                "left_edge_truncated",
                "0"
            )
        ):

            context_caution_parts.append(
                "left_edge_truncated"
            )


        if i(
            a.get(
                "right_edge_truncated",
                "0"
            )
        ):

            context_caution_parts.append(
                "right_edge_truncated"
            )


        if i(
            a.get(
                "contains_interrupted_candidate",
                "0"
            )
        ):

            context_caution_parts.append(
                "contains_interrupted_candidate"
            )


        context_caution = (
            ";".join(
                context_caution_parts
            )
            if context_caution_parts
            else
            "none_recorded"
        )


        source_mag = a[
            "MAG"
        ]

        source_contig = a[
            "contig"
        ]

        candidate_ids = a[
            "candidates"
        ]


        antismash = (
            "see_existing_candidate_validation_98B"
        )

        bagel_aoi = (
            "see_existing_candidate_validation_98B"
        )

        bagel_candidate = (
            "see_existing_candidate_validation_98B"
        )


    # --------------------------------------------------------
    # Additional specialized
    # --------------------------------------------------------

    elif rclass == "additional_specialized":

        source_mag = ""
        source_contig = r[
            "source_id"
        ]

        candidate_ids = r[
            "candidate_protein_ids"
        ]


        anti = i(
            r[
                "antiSMASH_RiPP_context_any_group_member"
            ]
        )

        bagel = i(
            r[
                "BAGEL_AOI_any_group_member"
            ]
        )

        bagel_c = i(
            r[
                "BAGEL_candidate_level_any_group_member"
            ]
        )


        antismash = str(
            anti
        )

        bagel_aoi = str(
            bagel
        )

        bagel_candidate = str(
            bagel_c
        )


        if bagel_c:

            prediction_support = (
                "BAGEL_candidate_level_specialized_support"
            )

        elif bagel:

            prediction_support = (
                "BAGEL_locus_level_specialized_support"
            )

        elif anti:

            prediction_support = (
                "antiSMASH_RiPP_like_fragmentary_context"
            )

        else:

            prediction_support = (
                "specialized_support_unresolved"
            )


        if robust > 0:

            evidence_pattern = (
                "additional_specialized_context_with_robust_recruitment"
            )

        else:

            evidence_pattern = (
                "additional_specialized_context_without_robust_recruitment"
            )


        if anti:

            context_caution = (
                "antiSMASH_RiPP_region_fragmentary_contig_edge_context"
            )

        else:

            context_caution = (
                "incomplete_contig_context"
            )


    # --------------------------------------------------------
    # Additional exploratory
    # --------------------------------------------------------

    elif rclass == "additional_exploratory":

        source_mag = ""
        source_contig = r[
            "source_id"
        ]

        candidate_ids = r[
            "candidate_protein_ids"
        ]

        antismash = "0"
        bagel_aoi = "0"
        bagel_candidate = "0"

        prediction_support = (
            "GA_HMM_profile_only_no_specialized_region"
        )


        if robust > 0:

            evidence_pattern = (
                "exploratory_profile_only_with_robust_recruitment"
            )

        else:

            evidence_pattern = (
                "exploratory_profile_only_without_robust_recruitment"
            )


        context_caution = (
            "profile_only_incomplete_context"
        )


    else:

        raise RuntimeError(
            f"Unexpected reference class: {rclass}"
        )


    catalog.append({
        "context_id":
            rid,

        "context_category":
            rclass,

        "source_MAG":
            source_mag,

        "source_contig_or_group_representative":
            source_contig,

        "operational_group":
            r[
                "context_group"
            ],

        "group_member_contigs":
            r[
                "group_member_contigs"
            ],

        "producer_source_context":
            r[
                "producer_source_context"
            ],

        "reference_length_bp":
            r[
                "reference_length"
            ],

        "candidate_protein_ids":
            candidate_ids,

        "prediction_support":
            prediction_support,

        "antiSMASH_RiPP_context":
            antismash,

        "BAGEL_AOI_support":
            bagel_aoi,

        "BAGEL_candidate_level_support":
            bagel_candidate,

        "context_caution":
            context_caution,

        "samples_with_any_mapping_signal":
            r[
                "samples_with_any_signal"
            ],

        "samples_robust":
            robust,

        "samples_high_breadth_MAPQ_sensitive":
            mapq_sensitive,

        "samples_partial":
            partial,

        "samples_weak_or_localized":
            weak,

        "samples_no_signal":
            no_signal,

        "max_breadth_MAPQ10":
            r[
                "max_breadth_MAPQ10"
            ],

        "median_breadth_MAPQ10":
            r[
                "median_breadth_MAPQ10"
            ],

        "max_mean_depth_MAPQ10":
            r[
                "max_mean_depth_MAPQ10"
            ],

        "max_normalized_mean_depth_MAPQ10":
            r[
                "max_normalized_mean_depth_MAPQ10"
            ],

        "evidence_pattern":
            evidence_pattern,

        "novel_bacteriocin_claim":
            "NO",

        "expression_assessed":
            "NO",

        "functional_activity_assessed":
            "NO",
    })


write_tsv(
    OUT
    / "98E_FINAL_context_catalog_30.tsv",
    catalog
)


# ============================================================
# 2. Recurrence table
# ============================================================

recurrence_out = []


for r in recurrences:

    recurrence_out.append({
        **r,

        "interpretation":
            "recurrence_or_homology_to_existing_candidate_sequence",

        "count_as_independent_new_candidate_locus":
            "NO",

        "novelty_claim":
            "NO",

        "functional_confirmation":
            "NO",
    })


write_tsv(
    OUT
    / "98E_FINAL_known_candidate_recurrences.tsv",
    recurrence_out
)


# ============================================================
# 3. Useful subsets
# ============================================================

additional_specialized = [
    r
    for r in catalog
    if r[
        "context_category"
    ]
    ==
    "additional_specialized"
]


exploratory = [
    r
    for r in catalog
    if r[
        "context_category"
    ]
    ==
    "additional_exploratory"
]


existing = [
    r
    for r in catalog
    if r[
        "context_category"
    ]
    ==
    "existing_ATTRLOC"
]


additional_with_bagel_candidate = [
    r
    for r in additional_specialized
    if r[
        "BAGEL_candidate_level_support"
    ]
    ==
    "1"
]


specialized_no_robust = [
    r
    for r in additional_specialized
    if int(
        r[
            "samples_robust"
        ]
    ) == 0
]


write_tsv(
    OUT
    / "98E_FINAL_additional_specialized_20.tsv",
    additional_specialized
)


write_tsv(
    OUT
    / "98E_FINAL_exploratory_3.tsv",
    exploratory
)


write_tsv(
    OUT
    / "98E_FINAL_existing_ATTRLOC_7.tsv",
    existing
)


write_tsv(
    OUT
    / "98E_FINAL_additional_BAGEL_candidate_level.tsv",
    additional_with_bagel_candidate
)


write_tsv(
    OUT
    / "98E_FINAL_specialized_without_robust_mapping.tsv",
    specialized_no_robust
)


# ============================================================
# 4. Key findings
# ============================================================

key_findings = [
    (
        "existing_ATTRLOC_contexts",
        cg.get(
            "existing_ATTRLOC_loci",
            "7"
        )
    ),

    (
        "existing_candidate_hypotheses",
        cg.get(
            "existing_candidate_hypotheses",
            "NA"
        )
    ),

    (
        "existing_loci_with_at_least_2_strict_methods",
        cg.get(
            "existing_loci_with_at_least_2_strict_methods",
            "NA"
        )
    ),

    (
        "strict_Comparippson_supported_existing_candidates",
        cg.get(
            "existing_strict_Comparippson_result_supported_candidates",
            "NA"
        )
    ),

    (
        "direct_homology_recurrence_proteins",
        cg.get(
            "direct_homology_proteins_to_existing_candidates",
            "NA"
        )
    ),

    (
        "direct_homology_recurrence_contigs",
        cg.get(
            "direct_homology_contigs_to_existing_candidates",
            "NA"
        )
    ),

    (
        "existing_candidate_sequence_hypotheses_recovered_in_incomplete_contexts",
        cg.get(
            "existing_candidate_sequence_hypotheses_recovered",
            "NA"
        )
    ),

    (
        "mapped_operational_contexts",
        d1[
            "competitive_contexts"
        ]
    ),

    (
        "sample_context_combinations",
        d1[
            "sample_context_combinations"
        ]
    ),

    (
        "additional_specialized_operational_groups",
        d1[
            "additional_specialized_context_groups"
        ]
    ),

    (
        "additional_specialized_with_any_read_signal",
        d1[
            "additional_specialized_with_any_signal"
        ]
    ),

    (
        "additional_specialized_robust_in_at_least_one_sample",
        d1[
            "additional_specialized_robust_in_at_least_one_sample"
        ]
    ),

    (
        "additional_exploratory_operational_groups",
        d1[
            "additional_exploratory_context_groups"
        ]
    ),

    (
        "additional_exploratory_robust_in_at_least_one_sample",
        d1[
            "additional_exploratory_robust_in_at_least_one_sample"
        ]
    ),

    (
        "all_contexts_robust_in_at_least_one_sample",
        d1[
            "references_robust_in_at_least_one_sample"
        ]
    ),

    (
        "robust_sample_context_combinations",
        d1[
            "robust_combinations"
        ]
    ),

    (
        "high_breadth_MAPQ_sensitive_combinations",
        d1[
            "high_breadth_MAPQ_sensitive_combinations"
        ]
    ),

    (
        "partial_combinations",
        d1[
            "partial_combinations"
        ]
    ),

    (
        "weak_or_localized_combinations",
        d1[
            "weak_or_localized_combinations"
        ]
    ),

    (
        "no_signal_combinations",
        d1[
            "no_signal_combinations"
        ]
    ),

    (
        "reference_responsible_for_high_breadth_MAPQ_sensitive_calls",
        d1[
            "references_with_high_breadth_MAPQ_sensitive_calls"
        ]
    ),

    (
        "ATTRLOC002_robust_samples",
        d1[
            "ATTRLOC002_robust_samples"
        ]
    ),

    (
        "ATTRLOC002_MAPQ_sensitive_samples",
        d1[
            "ATTRLOC002_high_breadth_MAPQ_sensitive_samples"
        ]
    ),

    (
        "ATTRLOC006_robust_samples",
        d1[
            "ATTRLOC006_robust_samples"
        ]
    ),

    (
        "additional_BAGEL_candidate_level_operational_groups",
        len(
            additional_with_bagel_candidate
        )
    ),

    (
        "additional_specialized_without_robust_mapping",
        len(
            specialized_no_robust
        )
    ),

    (
        "functional_bacteriocin_confirmed_by_metagenomics",
        "NO"
    ),

    (
        "expression_demonstrated",
        "NO"
    ),

    (
        "antimicrobial_activity_demonstrated_by_sequence_analysis",
        "NO"
    ),

    (
        "additional_context_equivalent_to_novel_bacteriocin",
        "NO"
    ),

    (
        "98_bacteriocin_computational_branch",
        "CLOSED"
    ),
]


with (
    OUT
    / "98E_FINAL_key_findings.tsv"
).open(
    "w",
    newline=""
) as fh:

    writer = csv.writer(
        fh,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow([
        "metric",
        "value"
    ])

    writer.writerows(
        key_findings
    )


# ============================================================
# 5. Reporting guardrails
# ============================================================

guardrails = [
    (
        "ATTRLOC",
        "candidate_locus_or_context_not_functionally_confirmed_bacteriocin"
    ),

    (
        "additional_specialized",
        "specialized_prediction_context_not_equivalent_to_novel_bacteriocin"
    ),

    (
        "additional_exploratory",
        "GA_HMM_profile_signal_without_specialized_context_support"
    ),

    (
        "robust_mapping",
        "distributed_competitive_DNA_recruitment_not_absolute_presence_or_cell_abundance"
    ),

    (
        "normalized_depth",
        "comparative_context_recruitment_signal_not_relative_abundance"
    ),

    (
        "MAPQ_sensitive",
        "high_total_breadth_with_limited_discriminating_mapping_support"
    ),

    (
        "incomplete_bins",
        "negative_prediction_does_not_establish_biological_absence"
    ),

    (
        "antiSMASH_RiPP_like",
        "fragmentary_contig_edge_region_in_98C_not_complete_BGC"
    ),

    (
        "BAGEL",
        "one_specialized_method_family_even_when_HMM_BLAST_and_predict_suboutputs_agree"
    ),

    (
        "antiSMASH_plus_98C1_HMM",
        "not_fully_independent_method_families_due_shared_antiSMASH_profile_resources"
    ),

    (
        "Comparippson",
        "no_real_result_table_hits_under_search_performed"
    ),

    (
        "recurrence",
        "same_or_homologous_candidate_sequence_in_additional_context_not_novelty"
    ),

    (
        "DNA",
        "does_not_demonstrate_expression_peptide_production_or_antimicrobial_activity"
    ),

    (
        "producer_time",
        "no_inferential_producer_or_time_model_performed_in_step98"
    ),

    (
        "time_codes",
        "treated_as_codes_only_without_biological_week_mapping_inference"
    ),

    (
        "amplicon_data",
        "NOT_USED"
    ),
]


with (
    OUT
    / "98E_FINAL_reporting_guardrails.tsv"
).open(
    "w",
    newline=""
) as fh:

    writer = csv.writer(
        fh,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow([
        "topic",
        "interpretation"
    ])

    writer.writerows(
        guardrails
    )


# ============================================================
# 6. Thesis-ready Spanish summary
# ============================================================

n_robust = d1[
    "references_robust_in_at_least_one_sample"
]

n_specialized = d1[
    "additional_specialized_context_groups"
]

n_specialized_robust = d1[
    "additional_specialized_robust_in_at_least_one_sample"
]

n_exploratory = d1[
    "additional_exploratory_context_groups"
]

n_existing = d1[
    "existing_ATTRLOC_contexts"
]

n_recur_contigs = cg.get(
    "direct_homology_contigs_to_existing_candidates",
    "NA"
)

n_recur_prot = cg.get(
    "direct_homology_proteins_to_existing_candidates",
    "NA"
)


summary_text = f"""SÍNTESIS COMPUTACIONAL FINAL DE BACTERIOCINAS/RiPP

El análisis consolidó {n_existing} contextos ATTRLOC previamente establecidos.
La búsqueda positiva en bins incompletos recuperó además {n_recur_prot} proteínas
en {n_recur_contigs} contigs con homología directa a candidatos previamente
identificados. Estas señales se interpretan como recurrencias o contextos
adicionales de candidatos existentes y no como bacteriocinas nuevas.

Después de la deduplicación de los contextos adicionales se conservaron
{n_specialized} grupos operacionales con soporte especializado y
{n_exploratory} grupos exploratorios sustentados únicamente por perfiles HMM
GA. Los grupos especializados representan contextos respaldados por una región
RiPP-like de antiSMASH y/o por evidencia BAGEL. Las regiones RiPP-like
detectadas por antiSMASH en esta rama se localizaron en contigs fragmentarios
con señal de borde, por lo que no se interpretan como BGC completos.

El mapeo competitivo de las 18 bibliotecas shotgun contra 30 contextos
operacionales produjo 540 combinaciones muestra-contexto. En total,
{n_robust} de los 30 contextos alcanzaron el criterio operacional de
reclutamiento robusto en al menos una muestra. Entre los grupos adicionales
especializados, {n_specialized_robust}/{n_specialized} alcanzaron dicho
criterio y los {n_specialized} mostraron alguna señal de reclutamiento.

ATTRLOC002 constituyó el principal caso de ambigüedad de mapeo: no presentó
muestras clasificadas como robustas y mostró 10 combinaciones con amplitud
total alta pero sensible al filtrado MAPQ. En contraste, ATTRLOC006 presentó
reclutamiento robusto en 12 muestras. Este comportamiento es compatible con la
alta similitud de secuencia entre ambos contextos y limita la atribución
independiente de la señal a ATTRLOC002.

Los tres grupos exploratorios también reclutaron lecturas y alcanzaron señal
robusta en al menos una muestra. Sin embargo, el reclutamiento de ADN no
sustituye la evidencia especializada de biosíntesis y estos grupos se
mantienen separados de los candidatos con apoyo antiSMASH/BAGEL.

En conjunto, los resultados constituyen evidencia computacional de secuencias
y contextos genómicos compatibles con bacteriocinas/RiPP y de su reclutamiento
en las bibliotecas metagenómicas. No demuestran expresión génica, producción
del péptido, actividad antimicrobiana ni funcionalidad del clúster. Tampoco se
considera que un contexto adicional constituya por sí mismo una bacteriocina
nueva.
"""


(
    OUT
    / "98E_FINAL_thesis_ready_summary.txt"
).write_text(
    summary_text
)


# ============================================================
# 7. Validation
# ============================================================

if i(
    d1[
        "sample_context_combinations"
    ]
) != 540:

    raise RuntimeError(
        "Expected 540 sample-context combinations"
    )


if i(
    d1[
        "additional_specialized_context_groups"
    ]
) != 20:

    raise RuntimeError(
        "Expected 20 specialized additional groups"
    )


if i(
    d1[
        "additional_exploratory_context_groups"
    ]
) != 3:

    raise RuntimeError(
        "Expected 3 exploratory groups"
    )


if i(
    d1[
        "references_robust_in_at_least_one_sample"
    ]
) != 28:

    raise RuntimeError(
        "Expected 28 contexts robust in >=1 sample"
    )


if i(
    d1[
        "additional_specialized_robust_in_at_least_one_sample"
    ]
) != 19:

    raise RuntimeError(
        "Expected 19/20 specialized groups robust in >=1 sample"
    )


if i(
    d1[
        "ATTRLOC002_robust_samples"
    ]
) != 0:

    raise RuntimeError(
        "ATTRLOC002 robust count changed"
    )


if i(
    d1[
        "ATTRLOC006_robust_samples"
    ]
) != 12:

    raise RuntimeError(
        "ATTRLOC006 robust count changed"
    )


print(
    "FINAL_CONTEXTS=30"
)

print(
    "EXISTING_ATTRLOC=7"
)

print(
    "ADDITIONAL_SPECIALIZED=20"
)

print(
    "ADDITIONAL_EXPLORATORY=3"
)

print(
    "ROBUST_IN_AT_LEAST_ONE_SAMPLE=28"
)

print(
    "98E_FINAL=PASS"
)
