#!/usr/bin/env python3

import csv
import os
from collections import defaultdict

USER = os.environ["USER"]
ROOT = f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"

GENES = os.path.join(
    ROOT,
    "53_metabolic_annotation_consolidated",
    "gene_function_inventory_all_predicted.tsv"
)

MASTER = os.path.join(
    ROOT,
    "22_final_representative_mags",
    "representative_mags_master.tsv"
)

OUTDIR = os.path.join(
    ROOT,
    "55_cheese_metabolic_capabilities"
)

os.makedirs(OUTDIR, exist_ok=True)


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(rows)


def toks(value):
    if not value:
        return set()

    return {
        x.strip().replace("ko:", "")
        for x in value.split(",")
        if x.strip() and x.strip() not in {"-", "NA", "N/A"}
    }


def ectoks(value):
    if not value:
        return set()

    return {
        x.strip()
        for x in value.split(",")
        if x.strip() and x.strip() not in {"-", "NA", "N/A"}
    }


def contig_index(gene_id):
    if "___" not in gene_id:
        return "", None

    rest = gene_id.split("___", 1)[1]

    try:
        contig, idx = rest.rsplit("_", 1)
        return contig, int(idx)
    except Exception:
        return rest, None


genes = read_tsv(GENES)
master = read_tsv(MASTER)

MAG_ORDER = [
    x["representative_MAG"]
    for x in master
]

MAG_SET = set(MAG_ORDER)

genes_by_mag = defaultdict(list)
genes_by_location = defaultdict(dict)

for row in genes:

    mag = row["MAG"]

    if mag not in MAG_SET:
        continue

    genes_by_mag[mag].append(row)

    contig, idx = contig_index(
        row["gene_id"]
    )

    if idx is not None:
        genes_by_location[
            (mag, contig)
        ][idx] = row


def neighbor_rows(mag, row, radius=5):
    contig, idx = contig_index(
        row["gene_id"]
    )

    if idx is None:
        return []

    out = []

    loc = genes_by_location[
        (mag, contig)
    ]

    for j in range(
        idx - radius,
        idx + radius + 1
    ):
        if j in loc:
            out.append(
                (
                    j - idx,
                    loc[j]
                )
            )

    return out


def name(row):
    return (
        row.get(
            "Preferred_name",
            ""
        )
        or ""
    ).strip().lower()


def desc(row):
    return (
        row.get(
            "Description",
            ""
        )
        or ""
    ).strip().lower()


# ============================================================
# PRODUCTORES POTENCIALES DE AMINAS DE INTERÉS EN ALIMENTOS
#
# Se evita:
# - contar speA como "riesgo" automáticamente
# - contar K01584 de cadA/ldcC como arginina descarboxilasa
#   cuando K01582 + contexto cadB sustentan lisina descarboxilasa
# - contar K01581 como sistema completo sin contexto
# ============================================================

amine_genes = []

for mag in MAG_ORDER:

    for row in genes_by_mag[mag]:

        k = toks(
            row.get(
                "KEGG_ko",
                ""
            )
        )

        ec = ectoks(
            row.get(
                "EC",
                ""
            )
        )

        n = name(row)
        d = desc(row)

        product = None
        marker = None

        # HISTAMINA
        if (
            "K01590" in k
            or "4.1.1.22" in ec
        ):
            product = "histamine"
            marker = "histidine_decarboxylase"

        # TIRAMINA
        elif (
            "K22330" in k
            or (
                "4.1.1.25" in ec
                and (
                    n in {
                        "tdc",
                        "tyrdc"
                    }
                    or "tyrosine decarboxylase" in d
                )
            )
        ):
            product = "tyramine"
            marker = "tyrosine_decarboxylase"

        # CADAVERINA
        # K01582 tiene prioridad frente a una asignación dual
        # K01582/K01584 cuando la anotación dice lisina/cadA/ldcC.
        elif (
            "K01582" in k
            or "4.1.1.18" in ec
        ):
            product = "cadaverine"
            marker = "lysine_decarboxylase"

        # PUTRESCINA DESDE ORNITINA
        elif (
            "K01581" in k
            or "4.1.1.17" in ec
        ):
            product = "putrescine"
            marker = "ornithine_decarboxylase"

        else:
            continue

        neighbors = neighbor_rows(
            mag,
            row,
            radius=5
        )

        neighbor_names = {
            name(x)
            for delta, x in neighbors
            if delta != 0
        }

        neighbor_kos = set()

        for delta, x in neighbors:
            if delta != 0:
                neighbor_kos.update(
                    toks(
                        x.get(
                            "KEGG_ko",
                            ""
                        )
                    )
                )

        near3 = [
            (delta, x)
            for delta, x in neighbors
            if (
                delta != 0
                and abs(delta) <= 3
            )
        ]

        level = "enzyme_only"
        context = []

        if product == "cadaverine":

            if any(
                name(x) == "cadb"
                or "K03757" in toks(
                    x.get(
                        "KEGG_ko",
                        ""
                    )
                )
                for delta, x in near3
            ):
                level = "system_supported"
                context.append(
                    "cadB_adjacent"
                )

            if any(
                name(x) == "cadc"
                for delta, x in near3
            ):
                context.append(
                    "cadC_adjacent"
                )

        elif product == "putrescine":

            if any(
                name(x) == "pote"
                or "K03756" in toks(
                    x.get(
                        "KEGG_ko",
                        ""
                    )
                )
                for delta, x in near3
            ):
                level = "system_supported"
                context.append(
                    "potE_adjacent"
                )

        elif product == "histamine":

            transporter = any(
                (
                    "transporter" in desc(x)
                    or "antiporter" in desc(x)
                    or "permease" in desc(x)
                    or "K03294" in toks(
                        x.get(
                            "KEGG_ko",
                            ""
                        )
                    )
                )
                for delta, x in near3
            )

            if transporter:
                level = "enzyme_plus_transport_context"
                context.append(
                    "adjacent_amino_acid_transport_context"
                )

        elif product == "tyramine":

            transporter = any(
                (
                    "amino acid permease" in desc(x)
                    or "amino acid transporter" in desc(x)
                    or name(x) in {
                        "tyrp",
                        "tyrt"
                    }
                )
                for delta, x in near3
            )

            if transporter:
                level = "enzyme_plus_transport_context"
                context.append(
                    "adjacent_amino_acid_transport_context"
                )

        contig, idx = contig_index(
            row["gene_id"]
        )

        amine_genes.append({
            "MAG":
                mag,

            "potential_product":
                product,

            "evidence_level":
                level,

            "marker":
                marker,

            "gene_id":
                row["gene_id"],

            "contig":
                contig,

            "gene_index":
                idx,

            "Preferred_name":
                row.get(
                    "Preferred_name",
                    ""
                ),

            "Description":
                row.get(
                    "Description",
                    ""
                ),

            "KEGG_ko":
                row.get(
                    "KEGG_ko",
                    ""
                ),

            "EC":
                row.get(
                    "EC",
                    ""
                ),

            "context_support":
                ";".join(
                    context
                )
        })


# ============================================================
# POLIAMINAS / RESISTENCIA ÁCIDA, SEPARADAS DEL RIESGO BA
# ============================================================

other_rows = []

for mag in MAG_ORDER:

    for row in genes_by_mag[mag]:

        k = toks(
            row.get(
                "KEGG_ko",
                ""
            )
        )

        n = name(row)

        category = None

        if (
            "K01585" in k
            or n.startswith("spea")
        ):
            category = (
                "biosynthetic_arginine_to_agmatine"
            )

        elif (
            "K01584" in k
            and "K01582" not in k
            and n == "adia"
        ):
            category = (
                "arginine_decarboxylase_acid_resistance"
            )

        elif "K01480" in k:
            category = (
                "agmatinase_putrescine_biosynthesis"
            )

        elif "K10536" in k:
            category = (
                "agmatine_deiminase_pathway"
            )

        else:
            continue

        contig, idx = contig_index(
            row["gene_id"]
        )

        other_rows.append({
            "MAG":
                mag,

            "category":
                category,

            "gene_id":
                row["gene_id"],

            "contig":
                contig,

            "gene_index":
                idx,

            "Preferred_name":
                row.get(
                    "Preferred_name",
                    ""
                ),

            "Description":
                row.get(
                    "Description",
                    ""
                ),

            "KEGG_ko":
                row.get(
                    "KEGG_ko",
                    ""
                ),

            "EC":
                row.get(
                    "EC",
                    ""
                )
        })


# ============================================================
# RESUMEN POR MAG
# ============================================================

by_mag = defaultdict(
    lambda: defaultdict(list)
)

for row in amine_genes:
    by_mag[
        row["MAG"]
    ][
        row["potential_product"]
    ].append(
        row
    )


summary = []

for mag in MAG_ORDER:

    rec = {
        "MAG":
            mag
    }

    n_context_supported = 0
    n_enzyme_only = 0

    for product in [
        "histamine",
        "tyramine",
        "cadaverine",
        "putrescine"
    ]:

        rows = by_mag[
            mag
        ][
            product
        ]

        if not rows:

            status = "not_detected"

        elif any(
            x["evidence_level"]
            == "system_supported"
            for x in rows
        ):

            status = "system_supported"
            n_context_supported += 1

        elif any(
            x["evidence_level"]
            == "enzyme_plus_transport_context"
            for x in rows
        ):

            status = (
                "enzyme_plus_transport_context"
            )
            n_context_supported += 1

        else:

            status = "enzyme_only"
            n_enzyme_only += 1

        rec[
            product
        ] = status

    rec[
        "n_context_or_system_supported_BA_types"
    ] = n_context_supported

    rec[
        "n_enzyme_only_BA_types"
    ] = n_enzyme_only

    summary.append(
        rec
    )


write_tsv(
    os.path.join(
        OUTDIR,
        "curated_biogenic_amine_gene_evidence.tsv"
    ),
    amine_genes,
    [
        "MAG",
        "potential_product",
        "evidence_level",
        "marker",
        "gene_id",
        "contig",
        "gene_index",
        "Preferred_name",
        "Description",
        "KEGG_ko",
        "EC",
        "context_support"
    ]
)

write_tsv(
    os.path.join(
        OUTDIR,
        "curated_biogenic_amine_matrix_18MAGs.tsv"
    ),
    summary,
    [
        "MAG",
        "histamine",
        "tyramine",
        "cadaverine",
        "putrescine",
        "n_context_or_system_supported_BA_types",
        "n_enzyme_only_BA_types"
    ]
)

write_tsv(
    os.path.join(
        OUTDIR,
        "polyamine_and_arg_acid_resistance_genes.tsv"
    ),
    other_rows,
    [
        "MAG",
        "category",
        "gene_id",
        "contig",
        "gene_index",
        "Preferred_name",
        "Description",
        "KEGG_ko",
        "EC"
    ]
)


# ============================================================
# CONSOLA
# ============================================================

print("=" * 60)
print("PASO 77b COMPLETADO")
print("=" * 60)

print(
    f"Genes específicos de BA:       {len(amine_genes)}"
)

print(
    f"Genes poliaminas/Arg-stress:   {len(other_rows)}"
)

print()

for product in [
    "histamine",
    "tyramine",
    "cadaverine",
    "putrescine"
]:

    systems = sum(
        row[product]
        in {
            "system_supported",
            "enzyme_plus_transport_context"
        }
        for row in summary
    )

    enzyme_only = sum(
        row[product]
        == "enzyme_only"
        for row in summary
    )

    print(
        f"{product:12s} "
        f"context/system={systems:2d}  "
        f"enzyme_only={enzyme_only:2d}"
    )

print()
print(f"Salida: {OUTDIR}")
print(
    "PASO 77b FINALIZÓ CORRECTAMENTE; "
    "ES SEGURO SALIR."
)
