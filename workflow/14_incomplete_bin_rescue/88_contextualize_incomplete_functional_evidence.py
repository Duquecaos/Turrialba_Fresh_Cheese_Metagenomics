#!/usr/bin/env python3

from pathlib import Path
from collections import Counter, defaultdict
import csv
import os
import sys

# ============================================================
# PATHS
# ============================================================

USER = os.environ["USER"]

ROOT = Path(
    f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"
)

BIN_FILE = (
    ROOT
    / "68_incomplete_bin_inventory"
    / "raw_bin_inventory_1006.tsv"
)

CONTIG_FILE = (
    ROOT
    / "69_incomplete_unique_contigs"
    / "unique_contig_catalog_lt90.tsv"
)

IN87 = (
    ROOT
    / "73_incomplete_functional_rescue"
)

OUT = (
    ROOT
    / "74_incomplete_functional_evidence_context"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)

FILES87 = [
    "87A_cheese_metabolic_marker_hits.tsv",
    "87B_cheese_metabolic_systems_on_contigs.tsv",
    "87C_biogenic_amine_candidates.tsv",
    "87D_polyamine_arg_stress_genes.tsv",
    "87E_Opp_Dpp_components.tsv",
    "87F_Opp_Dpp_clusters.tsv",
    "87G_POT_Dtp_transporters.tsv",
    "87H_key_peptidases.tsv",
    "87I_surface_proteinase_CEP_candidates.tsv",
    "87J_curated_lipolysis_candidates.tsv",
    "87K_curated_amino_acid_aroma_markers.tsv",
    "87L_EPS_capsule_loci.tsv",
    "87M_stress_systems_on_contigs.tsv",
]


# ============================================================
# HELPERS
# ============================================================

def clean(x):
    if x is None:
        return ""

    x = str(x).strip()

    if x in {
        "",
        "-",
        "NA",
        "N/A",
        "None",
        "nan",
    }:
        return ""

    return x


def to_float(x, default=None):
    try:
        return float(clean(x))
    except Exception:
        return default


def split_semicolon(x):
    return [
        z.strip()
        for z in clean(x).split(";")
        if z.strip()
    ]


def write_tsv(path, rows, fields):
    with open(
        path,
        "w",
        newline="",
        encoding="utf-8"
    ) as fh:

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


def contamination_rank(value):
    """
    Jerarquía descriptiva ya utilizada en Paso 84.
    Menor valor = mejor contexto para atribución a bin.
    """

    x = to_float(value)

    if x is None:
        return 99

    if x <= 5:
        return 0

    if x <= 10:
        return 1

    if x <= 20:
        return 2

    return 3


def contamination_group_from_value(value):

    x = to_float(value)

    if x is None:
        return "unknown"

    if x <= 5:
        return "contam_le5"

    if x <= 10:
        return "contam_gt5_le10"

    if x <= 20:
        return "contam_gt10_le20"

    return "contam_gt20"


def prediction_integrity(row):

    if "partial" in row:

        code = clean(
            row.get("partial")
        )

        return {
            "00": "complete_predicted_CDS",
            "10": "partial_5prime_CDS",
            "01": "partial_3prime_CDS",
            "11": "partial_both_ends_CDS",
        }.get(
            code,
            "unknown_prediction_integrity"
        )

    if "contains_partial_gene" in row:

        flag = clean(
            row.get(
                "contains_partial_gene"
            )
        )

        if flag == "1":
            return "contains_partial_CDS"

        if flag == "0":
            return "all_component_CDS_complete"

        return "unknown_component_integrity"

    return "not_applicable"


# ============================================================
# VALIDAR INPUTS
# ============================================================

required_inputs = [
    BIN_FILE,
    CONTIG_FILE,
]

for fn in FILES87:
    required_inputs.append(
        IN87 / fn
    )

for path in required_inputs:

    if (
        not path.exists()
        or path.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Falta archivo requerido o está vacío: {path}"
        )


# ============================================================
# 1. INVENTARIO DE LOS 1006 BINS
# ============================================================

bins = {}

with open(
    BIN_FILE,
    newline="",
    encoding="utf-8"
) as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    required = {
        "raw_bin_uid",
        "coassembly",
        "producer",
        "binner",
        "completeness",
        "contamination",
        "quality_group",
        "contamination_group",
        "screening_mode",
    }

    missing = required - set(
        reader.fieldnames or []
    )

    if missing:
        raise RuntimeError(
            "Faltan columnas en raw_bin_inventory_1006.tsv: "
            + ", ".join(sorted(missing))
        )

    for row in reader:

        uid = clean(
            row["raw_bin_uid"]
        )

        if not uid:
            continue

        if uid in bins:
            raise RuntimeError(
                f"raw_bin_uid duplicado: {uid}"
            )

        bins[uid] = row


if len(bins) != 1006:
    raise RuntimeError(
        f"Se esperaban 1006 bins y se encontraron {len(bins)}"
    )


# ============================================================
# 2. CATÁLOGO DE CONTIGS INCOMPLETOS
# ============================================================

contigs = {}

with open(
    CONTIG_FILE,
    newline="",
    encoding="utf-8"
) as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    required = {
        "unique_contig_id",
        "member_raw_bin_uids",
        "represented_in_final18",
        "coassemblies",
        "producers",
        "strongest_incomplete_context",
    }

    missing = required - set(
        reader.fieldnames or []
    )

    if missing:
        raise RuntimeError(
            "Faltan columnas en catálogo de contigs: "
            + ", ".join(sorted(missing))
        )

    for row in reader:

        cid = clean(
            row["unique_contig_id"]
        )

        if not cid:
            continue

        if cid in contigs:
            raise RuntimeError(
                f"Contig duplicado en catálogo: {cid}"
            )

        uids = split_semicolon(
            row.get(
                "member_raw_bin_uids",
                ""
            )
        )

        if not uids:
            raise RuntimeError(
                f"Contig sin member_raw_bin_uids: {cid}"
            )

        missing_bins = [
            uid
            for uid in uids
            if uid not in bins
        ]

        if missing_bins:
            raise RuntimeError(
                f"{cid}: raw_bin_uid no encontrado: "
                + ";".join(missing_bins)
            )

        contigs[cid] = {
            "row": row,
            "uids": uids,
        }


if len(contigs) != 141173:
    raise RuntimeError(
        "Se esperaban 141173 contigs en el catálogo "
        f"y se encontraron {len(contigs)}"
    )


# ============================================================
# 3. CONTIGS QUE REALMENTE TIENEN EVIDENCIA EN PASO 87
# ============================================================

functional_contigs = set()
records_by_file = Counter()

for fn in FILES87:

    path = IN87 / fn

    with open(
        path,
        newline="",
        encoding="utf-8"
    ) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        if (
            reader.fieldnames is None
            or "contig_id"
            not in reader.fieldnames
        ):
            raise RuntimeError(
                f"{fn}: falta contig_id"
            )

        for row in reader:

            cid = clean(
                row.get(
                    "contig_id"
                )
            )

            if not cid:
                raise RuntimeError(
                    f"{fn}: fila sin contig_id"
                )

            if cid not in contigs:
                raise RuntimeError(
                    f"{fn}: contig no encontrado en catálogo: {cid}"
                )

            # Si Paso 87 ya trae los UIDs, deben concordar
            # exactamente con el catálogo del Paso 85.
            row_uids = split_semicolon(
                row.get(
                    "member_raw_bin_uids",
                    ""
                )
            )

            if row_uids:

                cat_uids = set(
                    contigs[cid]["uids"]
                )

                if set(row_uids) != cat_uids:
                    raise RuntimeError(
                        f"{fn}: member_raw_bin_uids discordantes "
                        f"para {cid}"
                    )

            functional_contigs.add(
                cid
            )

            records_by_file[
                fn
            ] += 1


# ============================================================
# 4. SELECCIÓN TRANSPARENTE DEL MEJOR CONTEXTO
# ============================================================

def choose_context(cid):

    c = contigs[cid]

    members = [
        bins[uid]
        for uid in c["uids"]
    ]

    medium = [
        x
        for x in members
        if clean(
            x.get(
                "quality_group"
            )
        ) == "B_50_89.99"
    ]

    fragmentary = [
        x
        for x in members
        if clean(
            x.get(
                "quality_group"
            )
        ) == "C_lt50"
    ]

    if medium:

        candidates = medium
        has_B = True

    else:

        candidates = fragmentary
        has_B = False

    if not candidates:
        raise RuntimeError(
            f"{cid}: no hay bins B_50_89.99 ni C_lt50"
        )

    # Primero menor categoría de contaminación.
    # Dentro de esa categoría:
    #   1. mayor completitud
    #   2. menor contaminación
    #   3. UID para desempate reproducible.
    best = sorted(
        candidates,
        key=lambda x: (
            contamination_rank(
                x.get(
                    "contamination"
                )
            ),
            -(
                to_float(
                    x.get(
                        "completeness"
                    ),
                    -1
                )
            ),
            to_float(
                x.get(
                    "contamination"
                ),
                999
            ),
            clean(
                x.get(
                    "raw_bin_uid"
                )
            ),
        )
    )[0]

    best_contam = to_float(
        best.get(
            "contamination"
        )
    )

    if has_B:

        if best_contam <= 5:

            support_class = (
                "B50_89_contam_le5"
            )

            interpretation = (
                "bin_context_positive_evidence_"
                "lower_contamination"
            )

            tax_guardrail = (
                "bin_context_permitted_but_no_species_"
                "claim_without_independent_taxonomy"
            )

        elif best_contam <= 10:

            support_class = (
                "B50_89_contam_gt5_le10"
            )

            interpretation = (
                "bin_context_positive_evidence_"
                "moderate_contamination_caution"
            )

            tax_guardrail = (
                "bin_context_cautious_no_species_claim_"
                "without_independent_taxonomy"
            )

        elif best_contam <= 20:

            support_class = (
                "B50_89_contam_gt10_le20"
            )

            interpretation = (
                "positive_gene_contig_evidence_"
                "bin_attribution_high_caution"
            )

            tax_guardrail = (
                "attribute_function_primarily_to_contig_gene_"
                "not_whole_bin_taxon"
            )

        else:

            support_class = (
                "B50_89_contam_gt20"
            )

            interpretation = (
                "positive_gene_contig_evidence_"
                "taxonomic_attribution_not_recommended"
            )

            tax_guardrail = (
                "attribute_function_to_contig_gene_"
                "not_whole_bin_taxon"
            )

    else:

        support_class = (
            "lt50_contig_gene_centric_only"
        )

        interpretation = (
            "positive_gene_contig_evidence_only_"
            "fragmentary_bins"
        )

        tax_guardrail = (
            "no_bin_or_taxon_attribution_from_"
            "fragmentary_context"
        )

    producers = sorted({
        clean(x.get("producer"))
        for x in members
        if clean(x.get("producer"))
    })

    coassemblies = sorted({
        clean(x.get("coassembly"))
        for x in members
        if clean(x.get("coassembly"))
    })

    return {
        "contig_id":
            cid,

        "context_support_class":
            support_class,

        "context_support_interpretation":
            interpretation,

        "taxonomic_attribution_guardrail":
            tax_guardrail,

        "best_support_bin_uid":
            clean(
                best.get(
                    "raw_bin_uid"
                )
            ),

        "best_support_bin_completeness":
            clean(
                best.get(
                    "completeness"
                )
            ),

        "best_support_bin_contamination":
            clean(
                best.get(
                    "contamination"
                )
            ),

        "best_support_bin_contamination_group":
            contamination_group_from_value(
                best.get(
                    "contamination"
                )
            ),

        "best_support_bin_quality_group":
            clean(
                best.get(
                    "quality_group"
                )
            ),

        "best_support_bin_coassembly":
            clean(
                best.get(
                    "coassembly"
                )
            ),

        "best_support_bin_producer":
            clean(
                best.get(
                    "producer"
                )
            ),

        "best_support_bin_binner":
            clean(
                best.get(
                    "binner"
                )
            ),

        "best_support_bin_screening_mode":
            clean(
                best.get(
                    "screening_mode"
                )
            ),

        "n_supporting_raw_bins":
            len(members),

        "n_supporting_B50_89_bins":
            len(medium),

        "n_supporting_lt50_bins":
            len(fragmentary),

        "n_supporting_coassemblies":
            len(coassemblies),

        "supporting_coassemblies":
            ";".join(coassemblies),

        "n_supporting_producers":
            len(producers),

        "supporting_producers":
            ";".join(producers),

        "cross_coassembly_exact_contig":
            int(
                len(coassemblies) > 1
            ),

        "cross_producer_exact_contig":
            int(
                len(producers) > 1
            ),

        "represented_in_final18":
            clean(
                c["row"].get(
                    "represented_in_final18"
                )
            ),
    }


context = {
    cid: choose_context(cid)
    for cid in sorted(
        functional_contigs
    )
}


# ============================================================
# 5. QC: EVIDENCIA 87 DEBE VENIR DE CONTIGS NUEVOS
# ============================================================

exact_final18_functional = [
    cid
    for cid, rec in context.items()
    if clean(
        rec[
            "represented_in_final18"
        ]
    ) == "1"
]

if exact_final18_functional:
    raise RuntimeError(
        "Hay contigs funcionales que aparecen como exactos "
        "a final18: "
        + ";".join(
            exact_final18_functional[:20]
        )
    )


# ============================================================
# 6. TABLA DE MEJOR CONTEXTO POR CONTIG FUNCIONAL
# ============================================================

context_fields = [
    "contig_id",
    "context_support_class",
    "context_support_interpretation",
    "taxonomic_attribution_guardrail",
    "best_support_bin_uid",
    "best_support_bin_completeness",
    "best_support_bin_contamination",
    "best_support_bin_contamination_group",
    "best_support_bin_quality_group",
    "best_support_bin_coassembly",
    "best_support_bin_producer",
    "best_support_bin_binner",
    "best_support_bin_screening_mode",
    "n_supporting_raw_bins",
    "n_supporting_B50_89_bins",
    "n_supporting_lt50_bins",
    "n_supporting_coassemblies",
    "supporting_coassemblies",
    "n_supporting_producers",
    "supporting_producers",
    "cross_coassembly_exact_contig",
    "cross_producer_exact_contig",
    "represented_in_final18",
]

write_tsv(
    OUT
    / "88A_functional_contig_context_best.tsv",
    [
        context[cid]
        for cid in sorted(context)
    ],
    context_fields
)


# ============================================================
# 7. TABLA LONG DE TODOS LOS BINS QUE RESPALDAN CADA CONTIG
# ============================================================

support_long = []

for cid in sorted(
    functional_contigs
):

    best_uid = (
        context[cid][
            "best_support_bin_uid"
        ]
    )

    for uid in contigs[
        cid
    ]["uids"]:

        b = bins[uid]

        support_long.append({
            "contig_id":
                cid,

            "raw_bin_uid":
                uid,

            "selected_best_support":
                int(
                    uid == best_uid
                ),

            "coassembly":
                clean(
                    b.get(
                        "coassembly"
                    )
                ),

            "producer":
                clean(
                    b.get(
                        "producer"
                    )
                ),

            "binner":
                clean(
                    b.get(
                        "binner"
                    )
                ),

            "completeness":
                clean(
                    b.get(
                        "completeness"
                    )
                ),

            "contamination":
                clean(
                    b.get(
                        "contamination"
                    )
                ),

            "quality_group":
                clean(
                    b.get(
                        "quality_group"
                    )
                ),

            "contamination_group":
                clean(
                    b.get(
                        "contamination_group"
                    )
                ),

            "screening_mode":
                clean(
                    b.get(
                        "screening_mode"
                    )
                ),

            "context_support_class_selected":
                context[
                    cid
                ][
                    "context_support_class"
                ],
        })


write_tsv(
    OUT
    / "88B_functional_contig_bin_support_long.tsv",
    support_long,
    [
        "contig_id",
        "raw_bin_uid",
        "selected_best_support",
        "coassembly",
        "producer",
        "binner",
        "completeness",
        "contamination",
        "quality_group",
        "contamination_group",
        "screening_mode",
        "context_support_class_selected",
    ]
)


# ============================================================
# 8. ENRIQUECER TODAS LAS TABLAS DEL PASO 87
# ============================================================

extra_fields = [
    "context_support_class",
    "context_support_interpretation",
    "taxonomic_attribution_guardrail",
    "best_support_bin_uid",
    "best_support_bin_completeness",
    "best_support_bin_contamination",
    "best_support_bin_contamination_group",
    "best_support_bin_quality_group",
    "best_support_bin_coassembly",
    "best_support_bin_producer",
    "best_support_bin_binner",
    "n_supporting_raw_bins",
    "n_supporting_B50_89_bins",
    "n_supporting_lt50_bins",
    "n_supporting_coassemblies",
    "supporting_coassemblies",
    "n_supporting_producers",
    "supporting_producers",
    "cross_coassembly_exact_contig",
    "cross_producer_exact_contig",
    "prediction_integrity",
]

summary_counter = Counter()
summary_contigs = defaultdict(set)

cep_priority = []
cep_priority_fields = None

ba_context_supported = []
ba_context_fields = None

total_contextualized_records = 0

for fn in FILES87:

    src = IN87 / fn

    dst = (
        OUT
        / (
            "contextualized_"
            + fn
        )
    )

    with open(
        src,
        newline="",
        encoding="utf-8"
    ) as fi:

        reader = csv.DictReader(
            fi,
            delimiter="\t"
        )

        original_fields = list(
            reader.fieldnames or []
        )

        out_fields = (
            original_fields
            + [
                x
                for x in extra_fields
                if x not in original_fields
            ]
        )

        with open(
            dst,
            "w",
            newline="",
            encoding="utf-8"
        ) as fo:

            writer = csv.DictWriter(
                fo,
                fieldnames=out_fields,
                delimiter="\t",
                lineterminator="\n",
                extrasaction="ignore"
            )

            writer.writeheader()

            for row in reader:

                cid = clean(
                    row["contig_id"]
                )

                ctx = context[cid]

                out = dict(row)

                for key in extra_fields:

                    if key == "prediction_integrity":
                        out[key] = prediction_integrity(
                            row
                        )

                    else:
                        out[key] = ctx.get(
                            key,
                            ""
                        )

                writer.writerow(
                    out
                )

                total_contextualized_records += 1

                key = (
                    fn,
                    out[
                        "context_support_class"
                    ],
                    out[
                        "prediction_integrity"
                    ],
                )

                summary_counter[
                    key
                ] += 1

                summary_contigs[
                    key
                ].add(
                    cid
                )

                # --------------------------------------------
                # CEP prioritario:
                # solamente A o T.
                # --------------------------------------------
                if fn == (
                    "87I_surface_proteinase_"
                    "CEP_candidates.tsv"
                ):

                    tier = clean(
                        out.get(
                            "rescue_tier_87"
                        )
                    )

                    if tier in {
                        "A_high_confidence_CEP_like_complete_prediction",
                        "T_truncated_CEP_like_candidate",
                    }:

                        cep_priority.append(
                            out
                        )

                        cep_priority_fields = out_fields

                # --------------------------------------------
                # Aminas con soporte de sistema/contexto
                # --------------------------------------------
                if fn == (
                    "87C_biogenic_amine_candidates.tsv"
                ):

                    level = clean(
                        out.get(
                            "evidence_level"
                        )
                    )

                    if level in {
                        "system_supported",
                        "enzyme_plus_transport_context",
                    }:

                        ba_context_supported.append(
                            out
                        )

                        ba_context_fields = out_fields


# ============================================================
# 9. RESUMEN TABLE × CONTEXT × INTEGRITY
# ============================================================

summary_rows = []

for key in sorted(
    summary_counter
):

    fn, support_class, integrity = key

    summary_rows.append({
        "input_table":
            fn,

        "context_support_class":
            support_class,

        "prediction_integrity":
            integrity,

        "n_records":
            summary_counter[key],

        "n_unique_contigs":
            len(
                summary_contigs[key]
            ),
    })


write_tsv(
    OUT
    / "88C_evidence_summary_by_table_context_integrity.tsv",
    summary_rows,
    [
        "input_table",
        "context_support_class",
        "prediction_integrity",
        "n_records",
        "n_unique_contigs",
    ]
)


# ============================================================
# 10. CEP PRIORITARIOS
# ============================================================

if cep_priority_fields is None:
    raise RuntimeError(
        "No se procesó la tabla CEP"
    )

write_tsv(
    OUT
    / "88D_priority_CEP_like_candidates.tsv",
    cep_priority,
    cep_priority_fields
)


# ============================================================
# 11. AMINAS CON SOPORTE CONTEXTUAL
# ============================================================

if ba_context_fields is None:
    raise RuntimeError(
        "No se procesó la tabla de aminas"
    )

write_tsv(
    OUT
    / "88E_context_supported_biogenic_amine_candidates.tsv",
    ba_context_supported,
    ba_context_fields
)


# ============================================================
# 12. RESUMEN GLOBAL
# ============================================================

context_counts = Counter(
    x[
        "context_support_class"
    ]
    for x in context.values()
)

cross_coassembly_n = sum(
    int(
        x[
            "cross_coassembly_exact_contig"
        ]
    )
    for x in context.values()
)

cross_producer_n = sum(
    int(
        x[
            "cross_producer_exact_contig"
        ]
    )
    for x in context.values()
)

summary_global = [
    (
        "raw_bins_inventory",
        len(bins)
    ),
    (
        "contigs_catalog_lt90",
        len(contigs)
    ),
    (
        "functional_evidence_contigs",
        len(
            functional_contigs
        )
    ),
    (
        "functional_contig_bin_support_pairs",
        len(
            support_long
        )
    ),
    (
        "contextualized_step87_records_total",
        total_contextualized_records
    ),
    (
        "functional_contigs_exact_final18",
        len(
            exact_final18_functional
        )
    ),
    (
        "cross_coassembly_exact_functional_contigs",
        cross_coassembly_n
    ),
    (
        "cross_producer_exact_functional_contigs",
        cross_producer_n
    ),
    (
        "priority_CEP_like_candidates_A_plus_T",
        len(
            cep_priority
        )
    ),
    (
        "context_supported_biogenic_amine_candidates",
        len(
            ba_context_supported
        )
    ),
]

for cls in [
    "B50_89_contam_le5",
    "B50_89_contam_gt5_le10",
    "B50_89_contam_gt10_le20",
    "B50_89_contam_gt20",
    "lt50_contig_gene_centric_only",
]:

    summary_global.append(
        (
            "functional_contigs__"
            + cls,
            context_counts.get(
                cls,
                0
            )
        )
    )


with open(
    OUT
    / "functional_evidence_context_summary.tsv",
    "w",
    newline="",
    encoding="utf-8"
) as fh:

    writer = csv.writer(
        fh,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow(
        [
            "metric",
            "value"
        ]
    )

    writer.writerows(
        summary_global
    )


# ============================================================
# 13. README
# ============================================================

readme = """\
PASO 88 - CONTEXTUALIZACION DE EVIDENCIA FUNCIONAL RESCATADA
=============================================================

OBJETIVO
--------
Agregar a la evidencia positiva del Paso 87 la pareja REAL de
completitud-contaminacion de los bins que contienen cada contig.

Este paso NO vuelve a ejecutar eggNOG, Prodigal ni las reglas funcionales.

SELECCION DEL MEJOR CONTEXTO
----------------------------
1. Si el contig aparece en uno o mas bins B_50_89.99, se utilizan estos
   para contextualizar el hallazgo.

2. Se prioriza primero la categoria de menor contaminacion:
      <=5
      >5 a <=10
      >10 a <=20
      >20

3. Dentro de la misma categoria de contaminacion se elige:
      mayor completitud,
      luego menor contaminacion,
      luego raw_bin_uid para desempate reproducible.

4. Si el contig aparece solamente en bins <50% de completitud, la
   evidencia permanece exclusivamente a nivel de gen/contig.

IMPORTANTE
----------
La contaminacion del bin modula la confiabilidad de la ATRIBUCION del
contig a un contexto genomico/taxon, no invalida automaticamente la
anotacion funcional de la secuencia.

Los bins con contaminacion >10% deben interpretarse principalmente como
soporte positivo a nivel de gen/contig.

Los bins <50% de completitud no se utilizan para inferir ausencia de
genes ni para reconstruir una capacidad negativa del organismo.

PREDICCIONES PARCIALES
----------------------
Prodigal:
    00 = CDS predicha completa
    10 = parcial 5'
    01 = parcial 3'
    11 = parcial en ambos extremos

Para sistemas o loci multigenicos se conserva contains_partial_gene.

CEP
---
88D contiene solamente:
    A_high_confidence_CEP_like_complete_prediction
    T_truncated_CEP_like_candidate

Los tiers C representan subtilasas grandes no-CEP.
Los tiers D pertenecen al screening amplio y NO deben contabilizarse
como CEP-like.

AMINAS BIOGENAS
---------------
88E contiene solamente candidatos con:
    system_supported
    enzyme_plus_transport_context

Los candidatos enzyme_only permanecen en la tabla contextualizada
completa y no se eliminan.

NO SE GENERA UN SCORE GLOBAL
----------------------------
Las dimensiones se conservan separadas:
    calidad del contexto de bin,
    contaminacion,
    completitud,
    integridad de CDS,
    contexto funcional,
    redundancia entre bins/coassemblies/productores.

Los conteos de registros de diferentes tablas NO son unidades biologicas
independientes: un mismo gen o contig puede contribuir a multiples
categorias funcionales.
"""

with open(
    OUT / "README_step88.txt",
    "w",
    encoding="utf-8"
) as fh:

    fh.write(
        readme
    )


# ============================================================
# CONSOLA
# ============================================================

print("=" * 70)
print("PASO 88 COMPLETADO")
print("=" * 70)

for metric, value in summary_global:
    print(
        f"{metric:<55} {value}"
    )

print()
print("Distribucion de contigs funcionales por contexto:")

for cls in [
    "B50_89_contam_le5",
    "B50_89_contam_gt5_le10",
    "B50_89_contam_gt10_le20",
    "B50_89_contam_gt20",
    "lt50_contig_gene_centric_only",
]:
    print(
        f"{cls:<42} "
        f"{context_counts.get(cls, 0)}"
    )

print()
print(
    f"Salida: {OUT}"
)

print()
print(
    "PASO 88 FINALIZÓ CORRECTAMENTE; "
    "ES SEGURO SALIR."
)
