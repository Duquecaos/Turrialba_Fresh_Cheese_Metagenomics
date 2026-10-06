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


# ============================================================
# HELPERS
# ============================================================

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


def tokens(value):
    if value is None:
        return set()

    value = value.strip()

    if not value or value in {"-", "NA", "N/A"}:
        return set()

    return {
        x.strip().replace("ko:", "")
        for x in value.split(",")
        if x.strip()
    }


def ec_tokens(value):
    if value is None:
        return set()

    value = value.strip()

    if not value or value in {"-", "NA", "N/A"}:
        return set()

    return {
        x.strip()
        for x in value.split(",")
        if x.strip()
    }


def pname(row):
    return (
        row.get("Preferred_name", "")
        or ""
    ).strip().lower()


def description(row):
    return (
        row.get("Description", "")
        or ""
    ).strip().lower()


def has_ko(row, *wanted):
    ks = tokens(row.get("KEGG_ko", ""))
    return any(k in ks for k in wanted)


def has_ec(row, *wanted):
    es = ec_tokens(row.get("EC", ""))
    return any(e in es for e in wanted)


def contig_and_index(gene_id):
    if "___" not in gene_id:
        return "", None

    rest = gene_id.split("___", 1)[1]

    try:
        contig, idx = rest.rsplit("_", 1)
        return contig, int(idx)
    except Exception:
        return rest, None


# ============================================================
# INPUT
# ============================================================

genes = read_tsv(GENES)
master = read_tsv(MASTER)

if len(master) != 18:
    raise SystemExit(
        f"ERROR: master contiene {len(master)} MAGs; esperados 18"
    )

MAG_ORDER = [
    x["representative_MAG"]
    for x in master
]

MAG_SET = set(MAG_ORDER)

genes_by_mag = defaultdict(list)

for row in genes:

    mag = row["MAG"]

    if mag not in MAG_SET:
        raise SystemExit(
            f"ERROR: MAG desconocido en inventario: {mag}"
        )

    genes_by_mag[mag].append(row)


# ============================================================
# MARCADORES ESPECÍFICOS
# ============================================================

def marker_for(row):

    n = pname(row)
    d = description(row)

    hits = set()

    # --------------------------------------------------------
    # GALACTOSA - LELOIR
    # --------------------------------------------------------

    if n == "galk" and has_ko(row, "K00849"):
        hits.add("galK")

    if n == "galt" and has_ko(row, "K00965"):
        hits.add("galT")

    if n == "gale" and has_ko(row, "K01784"):
        hits.add("galE")

    if n == "galm" and has_ko(row, "K01785"):
        hits.add("galM")


    # --------------------------------------------------------
    # LACTOSA - IMPORTACIÓN/HIDRÓLISIS
    # --------------------------------------------------------

    if n == "lacs":
        hits.add("lacS")

    if n == "lacy":
        hits.add("lacY")

    if (
        has_ec(row, "3.2.1.23")
        and "galactosidase" in d
    ):
        hits.add("beta_galactosidase")

    if n.startswith("lacz") and has_ko(row, "K01190"):
        hits.add("beta_galactosidase")

    if n in {"bgaa", "bgly"} and has_ko(row, "K12308"):
        hits.add("beta_galactosidase")


    # --------------------------------------------------------
    # LACTOSA PTS
    # --------------------------------------------------------

    if n == "lace":
        hits.add("lacE")

    if n.startswith("lacf"):
        hits.add("lacF")

    if n == "lacg":
        hits.add("lacG")


    # --------------------------------------------------------
    # LACTATO
    # --------------------------------------------------------

    if (
        n in {
            "ldh",
            "ldha",
            "ldhb",
            "ldhd",
            "ldhx"
        }
        and (
            has_ko(
                row,
                "K00016",
                "K00018",
                "K03778"
            )
            or has_ec(
                row,
                "1.1.1.27",
                "1.1.1.28",
                "1.1.1.29"
            )
        )
    ):
        hits.add("LDH")


    # --------------------------------------------------------
    # ACETATO
    # --------------------------------------------------------

    if n == "pta" and has_ko(row, "K00625"):
        hits.add("pta")

    if n == "acka" and has_ko(row, "K00925"):
        hits.add("ackA")


    # --------------------------------------------------------
    # FORMATO
    # --------------------------------------------------------

    if n == "pfla" and has_ko(row, "K04069"):
        hits.add("pflA")

    if n == "pflb" and has_ko(row, "K00656"):
        hits.add("pflB")


    # --------------------------------------------------------
    # ETANOL / ACETALDEHÍDO
    # --------------------------------------------------------

    if n == "adhe" and has_ko(row, "K04072"):
        hits.add("adhE")


    # --------------------------------------------------------
    # CITRATO: TRANSPORTE
    # --------------------------------------------------------

    if (
        n in {
            "citp",
            "citt",
            "citm",
            "citn",
            "cimh",
            "maen",
            "ysde"
        }
        and "citrate" in d
    ):
        hits.add("citrate_transport")


    # --------------------------------------------------------
    # CITRATO LIASA
    # --------------------------------------------------------

    if n == "citd" and has_ko(row, "K01646"):
        hits.add("citD")

    if n == "cite" and has_ko(row, "K01644"):
        hits.add("citE")

    if n == "citf" and has_ko(row, "K01643"):
        hits.add("citF")

    if n == "citc" and has_ko(row, "K01910"):
        hits.add("citC")

    if n == "citx" and has_ko(row, "K05964"):
        hits.add("citX")

    if n == "citg" and has_ko(
        row,
        "K05966",
        "K13927",
        "K13930"
    ):
        hits.add("citG")


    # --------------------------------------------------------
    # ACETOIN
    # --------------------------------------------------------

    if (
        n in {
            "ilvb",
            "ilvi",
            "ilvg",
            "alss"
        }
        and has_ko(row, "K01652")
    ):
        hits.add("acetolactate_synthase")

    if (
        n in {
            "aldc",
            "buda"
        }
        and has_ko(row, "K01575")
    ):
        hits.add("acetolactate_decarboxylase")


    # --------------------------------------------------------
    # 2,3-BUTANODIOL
    # --------------------------------------------------------

    if (
        n in {
            "budc",
            "bdha",
            "butb"
        }
        or "butanediol dehydrogenase" in d
    ):
        hits.add("butanediol_dehydrogenase")


    # --------------------------------------------------------
    # AMINAS BIÓGENAS ESTRICTAS
    #
    # No se acepta "cadA" por nombre solamente.
    # No se acepta tdcB/C/E/F.
    # No se acepta K06966.
    # --------------------------------------------------------

    if (
        has_ec(row, "4.1.1.22")
        or has_ko(row, "K01590")
    ):
        hits.add("histidine_decarboxylase")

    if has_ec(row, "4.1.1.25"):
        hits.add("tyrosine_decarboxylase")

    if (
        has_ec(row, "4.1.1.17")
        or has_ko(row, "K01581")
    ):
        hits.add("ornithine_decarboxylase")

    if has_ec(row, "4.1.1.18"):
        hits.add("lysine_decarboxylase")

    if has_ec(row, "4.1.1.19"):
        hits.add("arginine_decarboxylase")

    if has_ko(row, "K10536"):
        hits.add("agmatine_deiminase")


    # --------------------------------------------------------
    # TRANSPORTADORES QUE PUEDEN DAR CONTEXTO
    # A DESCARBOXILASAS
    # --------------------------------------------------------

    if n == "cadb":
        hits.add("cadB")

    if n == "pote":
        hits.add("potE")

    if n in {
        "tyrp",
        "tyrt"
    }:
        hits.add("tyramine_transporter")

    if n in {
        "hdcp",
        "hdct"
    }:
        hits.add("histamine_transporter")

    return hits


# ============================================================
# RECUPERAR HITS
# ============================================================

marker_hits = defaultdict(
    lambda: defaultdict(list)
)

hit_rows = []


for mag in MAG_ORDER:

    for row in genes_by_mag[mag]:

        markers = marker_for(row)

        if not markers:
            continue

        contig, idx = contig_and_index(
            row["gene_id"]
        )

        for marker in sorted(markers):

            marker_hits[
                mag
            ][
                marker
            ].append(
                row
            )

            hit_rows.append({
                "MAG":
                    mag,

                "marker":
                    marker,

                "gene_id":
                    row["gene_id"],

                "contig":
                    contig,

                "gene_index":
                    "" if idx is None else idx,

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


write_tsv(
    os.path.join(
        OUTDIR,
        "strict_marker_hits.tsv"
    ),
    hit_rows,
    [
        "MAG",
        "marker",
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
# EVALUAR CAPACIDADES
# ============================================================

def present(mag, marker):
    return len(
        marker_hits[
            mag
        ][
            marker
        ]
    ) > 0


def evidence(mag, markers):
    return ";".join(
        m
        for m in markers
        if present(mag, m)
    )


capability_long = []
capability_matrix = []


for mag in MAG_ORDER:

    # --------------------------------------------------------
    # LELOIR
    # --------------------------------------------------------

    leloir_required = [
        "galK",
        "galT",
        "galE"
    ]

    n_leloir = sum(
        present(mag, x)
        for x in leloir_required
    )

    if n_leloir == 3:
        leloir_status = "complete_core"
    elif n_leloir > 0:
        leloir_status = "partial"
    else:
        leloir_status = "absent"


    # --------------------------------------------------------
    # LACTOSA: TRANSPORTE + BETA-GAL
    # --------------------------------------------------------

    lactose_transport = (
        present(mag, "lacS")
        or present(mag, "lacY")
    )

    beta_gal = present(
        mag,
        "beta_galactosidase"
    )

    if lactose_transport and beta_gal:
        lactose_direct = "supported"
    elif lactose_transport or beta_gal:
        lactose_direct = "partial"
    else:
        lactose_direct = "absent"


    # --------------------------------------------------------
    # LACTOSA PTS
    # --------------------------------------------------------

    pts_markers = [
        "lacE",
        "lacF",
        "lacG"
    ]

    n_pts = sum(
        present(mag, x)
        for x in pts_markers
    )

    if n_pts == 3:
        lactose_pts = "complete_marker_set"
    elif n_pts > 0:
        lactose_pts = "partial"
    else:
        lactose_pts = "absent"


    # --------------------------------------------------------
    # LACTATO
    # --------------------------------------------------------

    lactate = (
        "supported"
        if present(mag, "LDH")
        else "not_detected"
    )


    # --------------------------------------------------------
    # ACETATO
    # --------------------------------------------------------

    acetate_markers = [
        "pta",
        "ackA"
    ]

    n_acetate = sum(
        present(mag, x)
        for x in acetate_markers
    )

    if n_acetate == 2:
        acetate = "complete_marker_set"
    elif n_acetate == 1:
        acetate = "partial"
    else:
        acetate = "absent"


    # --------------------------------------------------------
    # PFL
    # --------------------------------------------------------

    pfl_markers = [
        "pflA",
        "pflB"
    ]

    n_pfl = sum(
        present(mag, x)
        for x in pfl_markers
    )

    if n_pfl == 2:
        pfl = "complete_marker_set"
    elif n_pfl == 1:
        pfl = "partial"
    else:
        pfl = "absent"


    # --------------------------------------------------------
    # CITRATO LIASA
    # --------------------------------------------------------

    citrate_core = [
        "citD",
        "citE",
        "citF"
    ]

    n_citrate_core = sum(
        present(mag, x)
        for x in citrate_core
    )

    citrate_transport = present(
        mag,
        "citrate_transport"
    )

    citrate_activation = any(
        present(mag, x)
        for x in [
            "citC",
            "citX",
            "citG"
        ]
    )

    if (
        citrate_transport
        and n_citrate_core == 3
    ):
        citrate = "strong_system"
    elif (
        citrate_transport
        or n_citrate_core > 0
        or citrate_activation
    ):
        citrate = "partial"
    else:
        citrate = "absent"


    # --------------------------------------------------------
    # ACETOIN
    # --------------------------------------------------------

    if (
        present(
            mag,
            "acetolactate_synthase"
        )
        and present(
            mag,
            "acetolactate_decarboxylase"
        )
    ):
        acetoin = "supported_marker_pair"

    elif (
        present(
            mag,
            "acetolactate_synthase"
        )
        or present(
            mag,
            "acetolactate_decarboxylase"
        )
    ):
        acetoin = "partial"

    else:
        acetoin = "absent"


    # --------------------------------------------------------
    # 2,3-BUTANODIOL
    # --------------------------------------------------------

    butanediol = (
        "supported"
        if present(
            mag,
            "butanediol_dehydrogenase"
        )
        else "not_detected"
    )


    # --------------------------------------------------------
    # AMINAS
    # --------------------------------------------------------

    amine_markers = {
        "histamine":
            "histidine_decarboxylase",

        "tyramine":
            "tyrosine_decarboxylase",

        "putrescine_from_ornithine":
            "ornithine_decarboxylase",

        "cadaverine":
            "lysine_decarboxylase",

        "agmatine_precursor":
            "arginine_decarboxylase"
    }

    amines = [
        label
        for label, marker
        in amine_markers.items()
        if present(mag, marker)
    ]


    # --------------------------------------------------------
    # MATRIZ
    # --------------------------------------------------------

    matrix_rec = {
        "MAG":
            mag,

        "galactose_Leloir":
            leloir_status,

        "lactose_transport_betaGal":
            lactose_direct,

        "lactose_PTS_LacEFG":
            lactose_pts,

        "lactate_LDH":
            lactate,

        "acetate_PtaAckA":
            acetate,

        "formate_PFL":
            pfl,

        "citrate_fermentation":
            citrate,

        "acetoin_branch":
            acetoin,

        "butanediol_branch":
            butanediol,

        "n_strict_biogenic_amine_types":
            len(amines),

        "strict_biogenic_amine_types":
            ";".join(amines)
    }

    capability_matrix.append(
        matrix_rec
    )


    # --------------------------------------------------------
    # LONG
    # --------------------------------------------------------

    long_specs = [

        (
            "galactose_Leloir",
            leloir_status,
            [
                "galK",
                "galT",
                "galE",
                "galM"
            ]
        ),

        (
            "lactose_transport_betaGal",
            lactose_direct,
            [
                "lacS",
                "lacY",
                "beta_galactosidase"
            ]
        ),

        (
            "lactose_PTS_LacEFG",
            lactose_pts,
            pts_markers
        ),

        (
            "lactate_LDH",
            lactate,
            [
                "LDH"
            ]
        ),

        (
            "acetate_PtaAckA",
            acetate,
            acetate_markers
        ),

        (
            "formate_PFL",
            pfl,
            pfl_markers
        ),

        (
            "citrate_fermentation",
            citrate,
            [
                "citrate_transport",
                "citC",
                "citX",
                "citG",
                "citD",
                "citE",
                "citF"
            ]
        ),

        (
            "acetoin_branch",
            acetoin,
            [
                "acetolactate_synthase",
                "acetolactate_decarboxylase"
            ]
        ),

        (
            "butanediol_branch",
            butanediol,
            [
                "butanediol_dehydrogenase"
            ]
        )
    ]

    for cap, status, markers in long_specs:

        capability_long.append({
            "MAG":
                mag,

            "capability":
                cap,

            "status":
                status,

            "markers_detected":
                evidence(
                    mag,
                    markers
                )
        })


write_tsv(
    os.path.join(
        OUTDIR,
        "metabolic_capability_matrix_18MAGs.tsv"
    ),
    capability_matrix,
    list(
        capability_matrix[0].keys()
    )
)


write_tsv(
    os.path.join(
        OUTDIR,
        "metabolic_capability_evidence_long.tsv"
    ),
    capability_long,
    [
        "MAG",
        "capability",
        "status",
        "markers_detected"
    ]
)


# ============================================================
# AMINAS BIÓGENAS ESTRICTAS
# ============================================================

strict_amine_markers = {
    "histidine_decarboxylase":
        "histamine",

    "tyrosine_decarboxylase":
        "tyramine",

    "ornithine_decarboxylase":
        "putrescine_from_ornithine",

    "lysine_decarboxylase":
        "cadaverine",

    "arginine_decarboxylase":
        "agmatine_precursor"
}


amine_rows = []


for row in hit_rows:

    marker = row[
        "marker"
    ]

    if marker not in strict_amine_markers:
        continue

    out = dict(row)

    out[
        "potential_product"
    ] = strict_amine_markers[
        marker
    ]

    amine_rows.append(
        out
    )


write_tsv(
    os.path.join(
        OUTDIR,
        "strict_biogenic_amine_candidates.tsv"
    ),
    amine_rows,
    [
        "MAG",
        "potential_product",
        "marker",
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
# CONTEXTO ±5 ORFs DE DESCARBOXILASAS
# ============================================================

context_rows = []


for amine in amine_rows:

    mag = amine[
        "MAG"
    ]

    target_contig = amine[
        "contig"
    ]

    try:
        target_idx = int(
            amine["gene_index"]
        )
    except Exception:
        continue

    for row in genes_by_mag[
        mag
    ]:

        contig, idx = contig_and_index(
            row["gene_id"]
        )

        if (
            contig != target_contig
            or idx is None
        ):
            continue

        delta = (
            idx
            - target_idx
        )

        if abs(delta) > 5:
            continue

        context_rows.append({
            "MAG":
                mag,

            "target_product":
                amine[
                    "potential_product"
                ],

            "target_gene":
                amine[
                    "gene_id"
                ],

            "neighbor_gene":
                row[
                    "gene_id"
                ],

            "relative_ORF":
                delta,

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


context_rows.sort(
    key=lambda x: (
        MAG_ORDER.index(
            x["MAG"]
        ),
        x["target_gene"],
        int(
            x["relative_ORF"]
        )
    )
)


write_tsv(
    os.path.join(
        OUTDIR,
        "strict_biogenic_amine_gene_context_pm5.tsv"
    ),
    context_rows,
    [
        "MAG",
        "target_product",
        "target_gene",
        "neighbor_gene",
        "relative_ORF",
        "Preferred_name",
        "Description",
        "KEGG_ko",
        "EC"
    ]
)


# ============================================================
# INTEGRAR TAXONOMÍA + CALIDAD
# ============================================================

master_lookup = {
    x["representative_MAG"]: x
    for x in master
}


integrated_rows = []


for row in capability_matrix:

    mag = row[
        "MAG"
    ]

    m = master_lookup[
        mag
    ]

    rec = {
        "MAG":
            mag,

        "completeness":
            m.get(
                "completeness",
                ""
            ),

        "contamination":
            m.get(
                "contamination",
                ""
            ),

        "phylum":
            m.get(
                "phylum",
                ""
            ),

        "class":
            m.get(
                "class",
                ""
            ),

        "order":
            m.get(
                "order",
                ""
            ),

        "family":
            m.get(
                "family",
                ""
            ),

        "genus":
            m.get(
                "genus",
                ""
            ),

        "species":
            m.get(
                "species",
                ""
            )
    }

    rec.update(row)

    integrated_rows.append(
        rec
    )


write_tsv(
    os.path.join(
        OUTDIR,
        "MAG_taxonomy_quality_metabolism.tsv"
    ),
    integrated_rows,
    [
        "MAG",
        "completeness",
        "contamination",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        *[
            x
            for x in capability_matrix[0].keys()
            if x != "MAG"
        ]
    ]
)


# ============================================================
# README
# ============================================================

with open(
    os.path.join(
        OUTDIR,
        "README_step77.txt"
    ),
    "w"
) as f:

    f.write(
        "Paso 77 - Refinamiento de capacidades "
        "metabolicas relevantes para queso\n"
    )

    f.write(
        "====================================================\n\n"
    )

    f.write(
        "Este paso reemplaza el screening sensible "
        "del Paso 76 por marcadores especificos.\n\n"
    )

    f.write(
        "complete_core / complete_marker_set / "
        "strong_system indican que los marcadores "
        "definidos para este analisis fueron recuperados; "
        "NO prueban expresion o actividad metabolica.\n\n"
    )

    f.write(
        "Las aminas biogenas se filtran por EC/KO "
        "de descarboxilasas. No se aceptan nombres "
        "ambiguos como cadA K01534, tdcB/C/E/F "
        "ni K06966 como evidencia estricta.\n"
    )


# ============================================================
# CONSOLA
# ============================================================

print(
    "============================================================"
)

print(
    "PASO 77 COMPLETADO"
)

print(
    "============================================================"
)

print(
    f"MAGs analizados:                {len(MAG_ORDER)}"
)

print(
    f"Hits estrictos:                 {len(hit_rows)}"
)

print(
    f"Candidatos estrictos de aminas: {len(amine_rows)}"
)

print()

for cap in [
    "galactose_Leloir",
    "lactose_transport_betaGal",
    "lactose_PTS_LacEFG",
    "lactate_LDH",
    "acetate_PtaAckA",
    "formate_PFL",
    "citrate_fermentation",
    "acetoin_branch",
    "butanediol_branch"
]:

    states = defaultdict(int)

    for row in capability_long:

        if row["capability"] == cap:
            states[
                row["status"]
            ] += 1

    state_text = ", ".join(
        f"{k}={v}"
        for k, v in sorted(
            states.items()
        )
    )

    print(
        f"{cap:30s} {state_text}"
    )

print()
print(
    f"Salida: {OUTDIR}"
)

print(
    "PASO 77 FINALIZÓ CORRECTAMENTE; "
    "ES SEGURO SALIR."
)
