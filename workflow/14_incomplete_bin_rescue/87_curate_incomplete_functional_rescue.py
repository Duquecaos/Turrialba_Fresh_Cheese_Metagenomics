#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict, Counter
import csv
import os
import re
import sys

USER = os.environ["USER"]

ROOT = Path(
    f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"
)

GENES = (
    ROOT
    / "72_incomplete_eggnog_consolidated"
    / "gene_function_inventory_all_763306.tsv"
)

CONTIG_CATALOG = (
    ROOT
    / "69_incomplete_unique_contigs"
    / "unique_contig_catalog_lt90.tsv"
)

FAA = (
    ROOT
    / "70_incomplete_gene_catalog"
    / "incomplete_lt90_new.prodigal.faa"
)

OUT = (
    ROOT
    / "73_incomplete_functional_rescue"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# CONSTANTES
# ============================================================

EXPECTED_GENES = 763306
EXPECTED_CONTIGS = 141173
EXPECTED_NEW_CONTIGS = 140276


# ============================================================
# IO
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
        "nan"
    }:
        return ""

    return x


def lower(x):
    return clean(x).lower()


def read_tsv(path):
    with open(
        path,
        newline="",
        encoding="utf-8"
    ) as fh:
        return list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


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


def to_int(x, default=0):
    try:
        return int(float(clean(x)))
    except Exception:
        return default


def to_float(x, default=0.0):
    try:
        return float(clean(x))
    except Exception:
        return default


def split_values(x):
    x = clean(x)

    if not x:
        return set()

    return {
        y.strip()
        for y in re.split(
            r"[,;]",
            x
        )
        if y.strip()
    }


def ko_tokens(row):

    vals = split_values(
        row.get(
            "eggnog_KEGG_ko",
            ""
        )
    )

    return {
        x.replace(
            "ko:",
            ""
        )
        for x in vals
    }


def ec_tokens(row):

    vals = split_values(
        row.get(
            "eggnog_EC",
            ""
        )
    )

    return {
        x.replace(
            "EC:",
            ""
        )
        for x in vals
    }


def has_ko(row, *wanted):

    have = ko_tokens(row)

    return any(
        x in have
        for x in wanted
    )


def has_ec(row, *wanted):

    have = ec_tokens(row)

    return any(
        x in have
        for x in wanted
    )


def pname(row):

    return lower(
        row.get(
            "eggnog_Preferred_name",
            ""
        )
    )


def description(row):

    return lower(
        row.get(
            "eggnog_Description",
            ""
        )
    )


def pfams(row):

    return clean(
        row.get(
            "eggnog_PFAMs",
            ""
        )
    )


# ============================================================
# CONTEXTO DE PROCEDENCIA
# ============================================================

print("=" * 70)
print("PASO 87 - RESCATE FUNCIONAL CURADO")
print("=" * 70)
print()
print("Cargando catálogo de contigs...")


contig_catalog_rows = read_tsv(
    CONTIG_CATALOG
)


if len(contig_catalog_rows) != EXPECTED_CONTIGS:
    raise RuntimeError(
        f"Se esperaban {EXPECTED_CONTIGS} contigs "
        f"en el catálogo y hay "
        f"{len(contig_catalog_rows)}"
    )


contig_meta = {
    row["unique_contig_id"]: row
    for row in contig_catalog_rows
}


new_contigs = {
    row["unique_contig_id"]
    for row in contig_catalog_rows
    if clean(
        row.get(
            "represented_in_final18",
            ""
        )
    ) in {
        "0",
        "False",
        "false"
    }
}


if len(new_contigs) != EXPECTED_NEW_CONTIGS:
    raise RuntimeError(
        f"Se esperaban {EXPECTED_NEW_CONTIGS} "
        f"contigs nuevos respecto a final18 y "
        f"hay {len(new_contigs)}"
    )


# ============================================================
# CARGA DEL INVENTARIO FUNCIONAL
# ============================================================

print("Cargando inventario funcional de genes...")


genes = read_tsv(
    GENES
)


if len(genes) != EXPECTED_GENES:
    raise RuntimeError(
        f"Se esperaban {EXPECTED_GENES} genes y "
        f"se encontraron {len(genes)}"
    )


gene_ids = [
    row["gene_id"]
    for row in genes
]


if len(set(gene_ids)) != EXPECTED_GENES:
    raise RuntimeError(
        "Los gene_id no son únicos."
    )


for row in genes:

    cid = row["contig_id"]

    if cid not in contig_meta:
        raise RuntimeError(
            f"Contig sin metadata: {cid}"
        )

    meta = contig_meta[cid]

    row["member_raw_bin_uids"] = meta.get(
        "member_raw_bin_uids",
        ""
    )

    row["annotation_priority"] = meta.get(
        "annotation_priority",
        ""
    )

    row["_name"] = pname(row)
    row["_desc"] = description(row)
    row["_pfam"] = pfams(row)
    row["_orf"] = to_int(
        row.get(
            "gene_number_on_contig",
            0
        )
    )


genes_by_contig = defaultdict(list)

for row in genes:
    genes_by_contig[
        row["contig_id"]
    ].append(row)


for cid in genes_by_contig:

    genes_by_contig[cid].sort(
        key=lambda x: x["_orf"]
    )


# ============================================================
# CAMPOS COMUNES DE SALIDA
# ============================================================

COMMON_FIELDS = [
    "block",
    "marker",
    "gene_id",
    "contig_id",
    "gene_number_on_contig",
    "start",
    "end",
    "strand",
    "aa_length",
    "partial",
    "gene_evidence_scope",
    "strongest_incomplete_context",
    "n_bin_memberships",
    "coassemblies",
    "producers",
    "binners",
    "quality_groups",
    "min_member_bin_completeness",
    "max_member_bin_completeness",
    "min_member_bin_contamination",
    "max_member_bin_contamination",
    "member_raw_bin_uids",
    "eggnog_Preferred_name",
    "eggnog_Description",
    "eggnog_EC",
    "eggnog_KEGG_ko",
    "eggnog_PFAMs",
]


def gene_output(
    row,
    block,
    marker
):

    return {
        "block":
            block,

        "marker":
            marker,

        "gene_id":
            row["gene_id"],

        "contig_id":
            row["contig_id"],

        "gene_number_on_contig":
            row.get(
                "gene_number_on_contig",
                ""
            ),

        "start":
            row.get(
                "start",
                ""
            ),

        "end":
            row.get(
                "end",
                ""
            ),

        "strand":
            row.get(
                "strand",
                ""
            ),

        "aa_length":
            row.get(
                "aa_length",
                ""
            ),

        "partial":
            row.get(
                "partial",
                ""
            ),

        "gene_evidence_scope":
            row.get(
                "gene_evidence_scope",
                ""
            ),

        "strongest_incomplete_context":
            row.get(
                "strongest_incomplete_context",
                ""
            ),

        "n_bin_memberships":
            row.get(
                "n_bin_memberships",
                ""
            ),

        "coassemblies":
            row.get(
                "coassemblies",
                ""
            ),

        "producers":
            row.get(
                "producers",
                ""
            ),

        "binners":
            row.get(
                "binners",
                ""
            ),

        "quality_groups":
            row.get(
                "quality_groups",
                ""
            ),

        "min_member_bin_completeness":
            row.get(
                "min_member_bin_completeness",
                ""
            ),

        "max_member_bin_completeness":
            row.get(
                "max_member_bin_completeness",
                ""
            ),

        "min_member_bin_contamination":
            row.get(
                "min_member_bin_contamination",
                ""
            ),

        "max_member_bin_contamination":
            row.get(
                "max_member_bin_contamination",
                ""
            ),

        "member_raw_bin_uids":
            row.get(
                "member_raw_bin_uids",
                ""
            ),

        "eggnog_Preferred_name":
            row.get(
                "eggnog_Preferred_name",
                ""
            ),

        "eggnog_Description":
            row.get(
                "eggnog_Description",
                ""
            ),

        "eggnog_EC":
            row.get(
                "eggnog_EC",
                ""
            ),

        "eggnog_KEGG_ko":
            row.get(
                "eggnog_KEGG_ko",
                ""
            ),

        "eggnog_PFAMs":
            row.get(
                "eggnog_PFAMs",
                ""
            ),
    }


# ============================================================
# ÍNDICE DE MARCADORES
# ============================================================

marker_hits = []
marker_by_contig = defaultdict(
    lambda: defaultdict(list)
)


def add_marker(
    row,
    block,
    marker
):

    rec = gene_output(
        row,
        block,
        marker
    )

    marker_hits.append(rec)

    marker_by_contig[
        row["contig_id"]
    ][
        marker
    ].append(row)


# ============================================================
# 1. METABOLISMO RELEVANTE PARA QUESO
# ============================================================

print(
    "Detectando marcadores metabólicos..."
)


for row in genes:

    if row.get(
        "has_eggnog_annotation"
    ) != "1":
        continue

    n = row["_name"]
    d = row["_desc"]


    # --------------------------------------------------------
    # GALACTOSA - LELOIR
    # --------------------------------------------------------

    if (
        n == "galk"
        and has_ko(
            row,
            "K00849"
        )
    ):
        add_marker(
            row,
            "cheese_metabolism",
            "galK"
        )

    if (
        n == "galt"
        and has_ko(
            row,
            "K00965"
        )
    ):
        add_marker(
            row,
            "cheese_metabolism",
            "galT"
        )

    if (
        n == "gale"
        and has_ko(
            row,
            "K01784"
        )
    ):
        add_marker(
            row,
            "cheese_metabolism",
            "galE"
        )

    if (
        n == "galm"
        and has_ko(
            row,
            "K01785"
        )
    ):
        add_marker(
            row,
            "cheese_metabolism",
            "galM"
        )


    # --------------------------------------------------------
    # LACTOSA
    # --------------------------------------------------------

    if n == "lacs":
        add_marker(
            row,
            "cheese_metabolism",
            "lacS"
        )

    if n == "lacy":
        add_marker(
            row,
            "cheese_metabolism",
            "lacY"
        )

    beta_gal = False

    if (
        has_ec(
            row,
            "3.2.1.23"
        )
        and "galactosidase" in d
    ):
        beta_gal = True

    if (
        n.startswith("lacz")
        and has_ko(
            row,
            "K01190"
        )
    ):
        beta_gal = True

    if (
        n in {
            "bgaa",
            "bgly"
        }
        and has_ko(
            row,
            "K12308"
        )
    ):
        beta_gal = True

    if beta_gal:
        add_marker(
            row,
            "cheese_metabolism",
            "beta_galactosidase"
        )


    # --------------------------------------------------------
    # LACTOSA PTS
    # --------------------------------------------------------

    if n == "lace":
        add_marker(
            row,
            "cheese_metabolism",
            "lacE"
        )

    if n.startswith("lacf"):
        add_marker(
            row,
            "cheese_metabolism",
            "lacF"
        )

    if n == "lacg":
        add_marker(
            row,
            "cheese_metabolism",
            "lacG"
        )


    # --------------------------------------------------------
    # LDH
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
        add_marker(
            row,
            "cheese_metabolism",
            "LDH"
        )


    # --------------------------------------------------------
    # PTA / ACKA
    # --------------------------------------------------------

    if (
        n == "pta"
        and has_ko(
            row,
            "K00625"
        )
    ):
        add_marker(
            row,
            "cheese_metabolism",
            "pta"
        )

    if (
        n == "acka"
        and has_ko(
            row,
            "K00925"
        )
    ):
        add_marker(
            row,
            "cheese_metabolism",
            "ackA"
        )


    # --------------------------------------------------------
    # PFL
    # --------------------------------------------------------

    if (
        n == "pfla"
        and has_ko(
            row,
            "K04069"
        )
    ):
        add_marker(
            row,
            "cheese_metabolism",
            "pflA"
        )

    if (
        n == "pflb"
        and has_ko(
            row,
            "K00656"
        )
    ):
        add_marker(
            row,
            "cheese_metabolism",
            "pflB"
        )


    # --------------------------------------------------------
    # adhE
    # --------------------------------------------------------

    if (
        n == "adhe"
        and has_ko(
            row,
            "K04072"
        )
    ):
        add_marker(
            row,
            "cheese_metabolism",
            "adhE"
        )


    # --------------------------------------------------------
    # CITRATO
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
        add_marker(
            row,
            "cheese_metabolism",
            "citrate_transport"
        )

    if (
        n == "citd"
        and has_ko(
            row,
            "K01646"
        )
    ):
        add_marker(
            row,
            "cheese_metabolism",
            "citD"
        )

    if (
        n == "cite"
        and has_ko(
            row,
            "K01644"
        )
    ):
        add_marker(
            row,
            "cheese_metabolism",
            "citE"
        )

    if (
        n == "citf"
        and has_ko(
            row,
            "K01643"
        )
    ):
        add_marker(
            row,
            "cheese_metabolism",
            "citF"
        )

    if (
        n == "citc"
        and has_ko(
            row,
            "K01910"
        )
    ):
        add_marker(
            row,
            "cheese_metabolism",
            "citC"
        )

    if (
        n == "citx"
        and has_ko(
            row,
            "K05964"
        )
    ):
        add_marker(
            row,
            "cheese_metabolism",
            "citX"
        )

    if (
        n == "citg"
        and has_ko(
            row,
            "K05966",
            "K13927",
            "K13930"
        )
    ):
        add_marker(
            row,
            "cheese_metabolism",
            "citG"
        )


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
        and has_ko(
            row,
            "K01652"
        )
    ):
        add_marker(
            row,
            "cheese_metabolism",
            "acetolactate_synthase"
        )

    if (
        n in {
            "aldc",
            "buda"
        }
        and has_ko(
            row,
            "K01575"
        )
    ):
        add_marker(
            row,
            "cheese_metabolism",
            "acetolactate_decarboxylase"
        )


    # --------------------------------------------------------
    # BUTANODIOL
    # --------------------------------------------------------

    if (
        n in {
            "budc",
            "bdha",
            "butb"
        }
        or "butanediol dehydrogenase" in d
    ):
        add_marker(
            row,
            "cheese_metabolism",
            "butanediol_dehydrogenase"
        )


# ============================================================
# RECONSTRUCCIÓN CONTIG-CÉNTRICA DE SISTEMAS
# ============================================================

def nearest_span(
    marker_dict,
    required
):
    """
    Devuelve el menor span de ORFs que contiene al menos
    un gen de cada marcador requerido.
    """

    pools = []

    for marker in required:

        rows = marker_dict.get(
            marker,
            []
        )

        if not rows:
            return None, []

        pools.append(rows)

    flattened = []

    for marker, rows in zip(
        required,
        pools
    ):
        for row in rows:
            flattened.append(
                (
                    row["_orf"],
                    marker,
                    row
                )
            )

    flattened.sort()

    best_span = None
    best_rows = []

    for i in range(len(flattened)):

        seen = {}
        start = flattened[i][0]

        for j in range(
            i,
            len(flattened)
        ):

            pos, marker, row = (
                flattened[j]
            )

            if marker not in seen:
                seen[marker] = row

            if len(seen) == len(required):

                span = pos - start

                if (
                    best_span is None
                    or span < best_span
                ):
                    best_span = span
                    best_rows = list(
                        seen.values()
                    )

                break

    return best_span, best_rows


metabolic_systems = []


def add_system(
    cid,
    system,
    markers_required,
    markers_detected,
    status,
    rows,
    orf_span=""
):

    meta = contig_meta[cid]

    metabolic_systems.append({
        "contig_id":
            cid,

        "system":
            system,

        "status":
            status,

        "markers_required":
            ";".join(
                markers_required
            ),

        "markers_detected":
            ";".join(
                sorted(
                    markers_detected
                )
            ),

        "n_required_markers":
            len(
                markers_required
            ),

        "n_detected_markers":
            len(
                set(
                    markers_detected
                )
            ),

        "ORF_span":
            orf_span,

        "gene_ids":
            ";".join(
                sorted({
                    x["gene_id"]
                    for x in rows
                })
            ),

        "contains_partial_gene":
            int(
                any(
                    clean(
                        x.get(
                            "partial",
                            ""
                        )
                    ) != "00"
                    for x in rows
                )
            ),

        "strongest_incomplete_context":
            meta.get(
                "strongest_incomplete_context",
                ""
            ),

        "n_bin_memberships":
            meta.get(
                "n_bin_memberships",
                ""
            ),

        "coassemblies":
            meta.get(
                "coassemblies",
                ""
            ),

        "producers":
            meta.get(
                "producers",
                ""
            ),

        "binners":
            meta.get(
                "binners",
                ""
            ),

        "quality_groups":
            meta.get(
                "quality_groups",
                ""
            ),

        "min_member_bin_completeness":
            meta.get(
                "min_member_bin_completeness",
                ""
            ),

        "max_member_bin_completeness":
            meta.get(
                "max_member_bin_completeness",
                ""
            ),

        "min_member_bin_contamination":
            meta.get(
                "min_member_bin_contamination",
                ""
            ),

        "max_member_bin_contamination":
            meta.get(
                "max_member_bin_contamination",
                ""
            ),

        "member_raw_bin_uids":
            meta.get(
                "member_raw_bin_uids",
                ""
            ),
    })


for cid, mdict in marker_by_contig.items():

    # Leloir
    req = [
        "galK",
        "galT",
        "galE"
    ]

    detected = [
        x
        for x in req
        if mdict.get(x)
    ]

    if detected:

        span, rr = nearest_span(
            mdict,
            req
        )

        if (
            span is not None
            and span <= 30
        ):
            status = (
                "complete_marker_set_"
                "same_contig"
            )
        else:
            status = (
                "partial_marker_evidence"
            )

            rr = [
                r
                for x in detected
                for r in mdict[x]
            ]

        add_system(
            cid,
            "galactose_Leloir",
            req,
            detected,
            status,
            rr,
            "" if span is None else span
        )


    # Lactosa transporte + betaGal
    transport = []

    for x in [
        "lacS",
        "lacY"
    ]:
        if mdict.get(x):
            transport.append(x)

    bg = bool(
        mdict.get(
            "beta_galactosidase"
        )
    )

    if transport or bg:

        rows = []

        for x in transport:
            rows.extend(
                mdict[x]
            )

        rows.extend(
            mdict.get(
                "beta_galactosidase",
                []
            )
        )

        if transport and bg:
            status = (
                "transport_plus_betaGal_"
                "same_contig"
            )
        else:
            status = (
                "partial_marker_evidence"
            )

        add_system(
            cid,
            "lactose_transport_betaGal",
            [
                "lacS_or_lacY",
                "beta_galactosidase"
            ],
            (
                transport
                + (
                    ["beta_galactosidase"]
                    if bg else []
                )
            ),
            status,
            rows
        )


    # Lactose PTS
    req = [
        "lacE",
        "lacF",
        "lacG"
    ]

    detected = [
        x
        for x in req
        if mdict.get(x)
    ]

    if detected:

        span, rr = nearest_span(
            mdict,
            req
        )

        if (
            span is not None
            and span <= 20
        ):
            status = (
                "complete_marker_set_"
                "same_contig"
            )
        else:
            status = (
                "partial_marker_evidence"
            )

            rr = [
                r
                for x in detected
                for r in mdict[x]
            ]

        add_system(
            cid,
            "lactose_PTS_LacEFG",
            req,
            detected,
            status,
            rr,
            "" if span is None else span
        )


    # LDH
    if mdict.get("LDH"):

        rr = mdict["LDH"]

        add_system(
            cid,
            "lactate_LDH",
            ["LDH"],
            ["LDH"],
            "single_marker_supported",
            rr,
            0
        )


    # Pta AckA
    req = [
        "pta",
        "ackA"
    ]

    detected = [
        x
        for x in req
        if mdict.get(x)
    ]

    if detected:

        span, rr = nearest_span(
            mdict,
            req
        )

        if (
            span is not None
            and span <= 20
        ):
            status = (
                "complete_marker_pair_"
                "same_contig"
            )
        else:
            status = (
                "partial_marker_evidence"
            )

            rr = [
                r
                for x in detected
                for r in mdict[x]
            ]

        add_system(
            cid,
            "acetate_PtaAckA",
            req,
            detected,
            status,
            rr,
            "" if span is None else span
        )


    # PFL
    req = [
        "pflA",
        "pflB"
    ]

    detected = [
        x
        for x in req
        if mdict.get(x)
    ]

    if detected:

        span, rr = nearest_span(
            mdict,
            req
        )

        if (
            span is not None
            and span <= 20
        ):
            status = (
                "complete_marker_pair_"
                "same_contig"
            )
        else:
            status = (
                "partial_marker_evidence"
            )

            rr = [
                r
                for x in detected
                for r in mdict[x]
            ]

        add_system(
            cid,
            "formate_PFL",
            req,
            detected,
            status,
            rr,
            "" if span is None else span
        )


    # Citrato
    core = [
        "citD",
        "citE",
        "citF"
    ]

    core_detected = [
        x
        for x in core
        if mdict.get(x)
    ]

    trans = bool(
        mdict.get(
            "citrate_transport"
        )
    )

    activation = [
        x
        for x in [
            "citC",
            "citX",
            "citG"
        ]
        if mdict.get(x)
    ]

    citrate_markers = (
        core_detected
        + activation
        + (
            ["citrate_transport"]
            if trans
            else []
        )
    )

    if citrate_markers:

        req_strong = [
            "citrate_transport",
            "citD",
            "citE",
            "citF"
        ]

        span, rr = nearest_span(
            mdict,
            req_strong
        )

        if (
            span is not None
            and span <= 30
        ):
            status = (
                "strong_same_contig_"
                "transport_plus_core"
            )

        elif len(core_detected) == 3:

            span_core, rr = nearest_span(
                mdict,
                core
            )

            status = (
                "citrate_lyase_core_"
                "same_contig"
                if (
                    span_core is not None
                    and span_core <= 20
                )
                else
                "partial_marker_evidence"
            )

            span = span_core

        else:

            status = (
                "partial_marker_evidence"
            )

            rr = [
                r
                for x in citrate_markers
                for r in mdict[x]
            ]

        add_system(
            cid,
            "citrate_fermentation",
            req_strong,
            citrate_markers,
            status,
            rr,
            "" if span is None else span
        )


    # Acetoin
    req = [
        "acetolactate_synthase",
        "acetolactate_decarboxylase"
    ]

    detected = [
        x
        for x in req
        if mdict.get(x)
    ]

    if detected:

        span, rr = nearest_span(
            mdict,
            req
        )

        if (
            span is not None
            and span <= 20
        ):
            status = (
                "supported_marker_pair_"
                "same_contig"
            )
        else:
            status = (
                "partial_marker_evidence"
            )

            rr = [
                r
                for x in detected
                for r in mdict[x]
            ]

        add_system(
            cid,
            "acetoin_branch",
            req,
            detected,
            status,
            rr,
            "" if span is None else span
        )


    # Butanediol
    if mdict.get(
        "butanediol_dehydrogenase"
    ):

        rr = mdict[
            "butanediol_dehydrogenase"
        ]

        add_system(
            cid,
            "butanediol_branch",
            [
                "butanediol_dehydrogenase"
            ],
            [
                "butanediol_dehydrogenase"
            ],
            "single_marker_supported",
            rr,
            0
        )


# ============================================================
# 2. AMINAS BIÓGENAS
# ============================================================

print(
    "Curando aminas biógenas..."
)


ba_candidates = []


def nearby_rows(
    row,
    radius=3
):

    cid = row["contig_id"]
    idx = row["_orf"]

    return [
        x
        for x in genes_by_contig[cid]
        if abs(
            x["_orf"] - idx
        ) <= radius
    ]


for row in genes:

    if row.get(
        "has_eggnog_annotation"
    ) != "1":
        continue

    k = ko_tokens(row)
    ec = ec_tokens(row)
    n = row["_name"]
    d = row["_desc"]

    product = ""
    marker = ""

    # Histamina
    if (
        "K01590" in k
        or "4.1.1.22" in ec
    ):
        product = "histamine"
        marker = (
            "histidine_decarboxylase"
        )

    # Tiramina
    elif (
        "K22330" in k
        or (
            "4.1.1.25" in ec
            and (
                n in {
                    "tdc",
                    "tyrdc"
                }
                or (
                    "tyrosine "
                    "decarboxylase"
                ) in d
            )
        )
    ):
        product = "tyramine"
        marker = (
            "tyrosine_decarboxylase"
        )

    # Cadaverina
    elif (
        "K01582" in k
        or "4.1.1.18" in ec
    ):
        product = "cadaverine"
        marker = (
            "lysine_decarboxylase"
        )

    # Putrescina
    elif (
        "K01581" in k
        or "4.1.1.17" in ec
    ):
        product = "putrescine"
        marker = (
            "ornithine_decarboxylase"
        )

    else:
        continue


    near = nearby_rows(
        row,
        radius=3
    )

    level = "enzyme_only"
    context_support = []


    if product == "cadaverine":

        if any(
            (
                x["_name"] == "cadb"
                or has_ko(
                    x,
                    "K03757"
                )
            )
            for x in near
            if x["gene_id"]
            != row["gene_id"]
        ):
            level = "system_supported"
            context_support.append(
                "cadB_adjacent"
            )


    elif product == "putrescine":

        if any(
            (
                x["_name"] == "pote"
                or has_ko(
                    x,
                    "K03756"
                )
            )
            for x in near
            if x["gene_id"]
            != row["gene_id"]
        ):
            level = "system_supported"
            context_support.append(
                "potE_adjacent"
            )


    elif product == "histamine":

        transporter = any(
            (
                "transporter" in x["_desc"]
                or "antiporter" in x["_desc"]
                or "permease" in x["_desc"]
                or has_ko(
                    x,
                    "K03294"
                )
            )
            for x in near
            if x["gene_id"]
            != row["gene_id"]
        )

        if transporter:

            level = (
                "enzyme_plus_transport_"
                "context"
            )

            context_support.append(
                "adjacent_transport_context"
            )


    elif product == "tyramine":

        transporter = any(
            (
                (
                    "amino acid permease"
                    in x["_desc"]
                )
                or (
                    "amino acid transporter"
                    in x["_desc"]
                )
                or x["_name"] in {
                    "tyrp",
                    "tyrt"
                }
            )
            for x in near
            if x["gene_id"]
            != row["gene_id"]
        )

        if transporter:

            level = (
                "enzyme_plus_transport_"
                "context"
            )

            context_support.append(
                "adjacent_amino_acid_"
                "transport_context"
            )


    rec = gene_output(
        row,
        "biogenic_amines",
        marker
    )

    rec.update({
        "potential_product":
            product,

        "evidence_level":
            level,

        "context_support":
            ";".join(
                context_support
            ),

        "neighbor_gene_ids_pm3":
            ";".join(
                x["gene_id"]
                for x in near
                if (
                    x["gene_id"]
                    != row["gene_id"]
                )
            )
    })

    ba_candidates.append(rec)


# Polyamine / acid resistance contextual genes

polyamine_rows = []

for row in genes:

    if row.get(
        "has_eggnog_annotation"
    ) != "1":
        continue

    k = ko_tokens(row)
    n = row["_name"]

    category = ""

    if (
        "K01585" in k
        or n.startswith("spea")
    ):
        category = (
            "biosynthetic_arginine_"
            "to_agmatine"
        )

    elif (
        "K01584" in k
        and "K01582" not in k
        and n == "adia"
    ):
        category = (
            "arginine_decarboxylase_"
            "acid_resistance"
        )

    elif "K01480" in k:
        category = (
            "agmatinase_putrescine_"
            "biosynthesis"
        )

    elif "K10536" in k:
        category = (
            "agmatine_deiminase_pathway"
        )

    else:
        continue

    rec = gene_output(
        row,
        "polyamine_arg_stress",
        category
    )

    polyamine_rows.append(
        rec
    )


# ============================================================
# 3. PROTEÓLISIS: OPP / DPP / POT / PEPTIDASAS
# ============================================================

print(
    "Reconstruyendo utilización de péptidos..."
)


def normalize_named_component(
    name,
    prefix
):

    n = lower(name)

    m = re.fullmatch(
        rf"{prefix}([abcdf])"
        rf"(?:[_-]?[0-9]+)?",
        n
    )

    if not m:
        return None

    return m.group(1).upper()


peptide_components = []
peptide_by_contig_system = defaultdict(
    lambda: defaultdict(list)
)


for row in genes:

    if row.get(
        "has_eggnog_annotation"
    ) != "1":
        continue

    for system in [
        "opp",
        "dpp"
    ]:

        comp = normalize_named_component(
            row.get(
                "eggnog_Preferred_name",
                ""
            ),
            system
        )

        if comp is None:
            continue

        rec = gene_output(
            row,
            "peptide_transport",
            f"{system.upper()}_{comp}"
        )

        rec["system"] = (
            system.upper()
        )

        rec["component"] = comp

        peptide_components.append(rec)

        peptide_by_contig_system[
            row["contig_id"]
        ][
            system.upper()
        ].append(
            (
                comp,
                row
            )
        )


peptide_clusters = []

required_components = {
    "A",
    "B",
    "C",
    "D",
    "F"
}


for cid, systems in (
    peptide_by_contig_system.items()
):

    for system, pairs in (
        systems.items()
    ):

        comps = {
            comp
            for comp, row
            in pairs
        }

        rows = [
            row
            for comp, row
            in pairs
        ]

        indices = [
            row["_orf"]
            for row in rows
        ]

        span = (
            max(indices)
            - min(indices)
            if indices
            else None
        )

        n_req = len(
            required_components
            .intersection(comps)
        )

        if (
            required_components
            .issubset(comps)
            and span is not None
            and span <= 20
        ):
            status = (
                "complete_named_cluster"
            )

        elif (
            n_req >= 3
            and span is not None
            and span <= 20
        ):
            status = (
                "partial_named_cluster"
            )

        else:
            status = (
                "dispersed_or_insufficient"
            )

        meta = contig_meta[cid]

        peptide_clusters.append({
            "contig_id":
                cid,

            "system":
                system,

            "components":
                ",".join(
                    sorted(comps)
                ),

            "n_required_components":
                n_req,

            "ORF_span":
                (
                    ""
                    if span is None
                    else span
                ),

            "cluster_status":
                status,

            "genes":
                ";".join(
                    sorted(
                        x["gene_id"]
                        for x in rows
                    )
                ),

            "contains_partial_gene":
                int(
                    any(
                        clean(
                            x.get(
                                "partial",
                                ""
                            )
                        ) != "00"
                        for x in rows
                    )
                ),

            "strongest_incomplete_context":
                meta.get(
                    "strongest_incomplete_context",
                    ""
                ),

            "producers":
                meta.get(
                    "producers",
                    ""
                ),

            "coassemblies":
                meta.get(
                    "coassemblies",
                    ""
                ),

            "member_raw_bin_uids":
                meta.get(
                    "member_raw_bin_uids",
                    ""
                ),
        })


# POT / DTP

pot_rows = []

for row in genes:

    if row.get(
        "has_eggnog_annotation"
    ) != "1":
        continue

    n = row["_name"]

    name_match = bool(
        re.fullmatch(
            r"dtp[abdt]"
            r"(?:[_-]?[0-9]+)?",
            n
        )
    )

    yjdl_match = (
        n == "yjdl"
        and has_ko(
            row,
            "K03305"
        )
    )

    if (
        name_match
        or yjdl_match
    ):

        pot_rows.append(
            gene_output(
                row,
                "peptide_transport",
                "POT_Dtp"
            )
        )


# Peptidasas

PEPTIDASE_TYPES = {
    "pepn": "PepN",
    "pepc": "PepC",
    "pepx": "PepX",
    "pepo": "PepO",
    "pepf": "PepF",
    "pepq": "PepQ",
    "pepv": "PepV",
    "pept": "PepT",
    "pepd": "PepD",
    "pepp": "PepP",
    "pepa": "PepA",
    "pepda": "PepDA",
    "pepb": "PepB",
    "pepe": "PepE",
    "peps": "PepS",
}


peptidase_rows = []

for row in genes:

    if row.get(
        "has_eggnog_annotation"
    ) != "1":
        continue

    n = row["_name"]

    base = re.sub(
        r"[_-]?[0-9]+$",
        "",
        n
    )

    if base not in PEPTIDASE_TYPES:
        continue

    rec = gene_output(
        row,
        "proteolysis",
        PEPTIDASE_TYPES[base]
    )

    rec["peptidase_type"] = (
        PEPTIDASE_TYPES[base]
    )

    peptidase_rows.append(rec)


# ============================================================
# 4. CEP-LIKE / SUBTILASAS
# ============================================================

print(
    "Curando proteinasa superficial CEP-like..."
)


explicit_prt_names = {
    "prtp",
    "prtb",
    "prth",
    "prts",
    "prtr",
}


accessory_patterns = {
    "PA":
        r"(^|,)PA($|,)",

    "SLAP":
        r"(^|,)SLAP($|,)",

    "fn3":
        r"fn3",

    "FIVAR":
        r"FIVAR",

    "Gram_pos_anchor":
        r"Gram_pos_anchor",

    "PKD":
        r"(^|,)PKD($|,)",

    "MAM":
        r"(^|,)MAM($|,)",

    "Big":
        r"(^|,)Big[_0-9]",

    "CARDB":
        r"CARDB",

    "Laminin_G":
        r"Laminin_G",

    "NHL":
        r"(^|,)NHL($|,)",
}


discordant_description_patterns = [
    r"gluconolactonase",
    r"xylan catabolic",
]


cep_pre = []
cep_gene_ids = set()


for row in genes:

    if row.get(
        "has_eggnog_annotation"
    ) != "1":
        continue

    n = row["_name"]
    d = row["_desc"]
    p = row["_pfam"]

    explicit_prt = (
        n in explicit_prt_names
        or bool(
            re.fullmatch(
                r"prt[a-z0-9]+",
                n
            )
        )
        or (
            "cell envelope proteinase"
            in d
        )
        or "surface proteinase" in d
        or "lactocepin" in d
        or "caseinase" in d
        or (
            "caseinolytic proteinase"
            in d
        )
        or (
            "extracellular proteinase"
            in d
        )
        or (
            "extracellular protease"
            in d
        )
    )

    has_s8_s53 = (
        "peptidase_s8"
        in p.lower()
        or "peptidase_s53"
        in p.lower()
        or "subtilase" in d
        or "subtilisin" in d
        or "peptidase s8" in d
        or "peptidase s53" in d
    )

    if not (
        explicit_prt
        or has_s8_s53
    ):
        continue

    cep_gene_ids.add(
        row["gene_id"]
    )

    cep_pre.append({
        "row":
            row,

        "explicit_prt":
            explicit_prt,

        "has_s8_s53":
            has_s8_s53
    })


# ------------------------------------------------------------
# Cargar sólo las secuencias candidatas CEP
# ------------------------------------------------------------

cep_sequences = {}

current = None
capture = False
seq = []


with open(
    FAA,
    encoding="utf-8"
) as fh:

    for line in fh:

        line = line.rstrip()

        if line.startswith(">"):

            if (
                current is not None
                and capture
            ):
                cep_sequences[
                    current
                ] = "".join(seq)

            current = (
                line[1:]
                .split()[0]
            )

            capture = (
                current
                in cep_gene_ids
            )

            seq = []

        elif capture:
            seq.append(
                line.strip()
            )

    if (
        current is not None
        and capture
    ):
        cep_sequences[
            current
        ] = "".join(seq)


cep_rows = []


for item in cep_pre:

    row = item["row"]

    gid = row["gene_id"]

    seq = cep_sequences.get(
        gid,
        ""
    )

    aa = to_int(
        row.get(
            "aa_length",
            0
        )
    )

    n = row["_name"]
    d = row["_desc"]
    p = row["_pfam"]

    has_s8_s53 = (
        item["has_s8_s53"]
    )

    explicit_prt = (
        item["explicit_prt"]
    )

    tail = (
        seq[-120:]
        if seq
        else ""
    )

    lpxtg = bool(
        re.search(
            r"LP.TG",
            tail,
            flags=re.IGNORECASE
        )
    )

    gram_anchor = (
        "gram_pos_anchor"
        in p.lower()
    )

    anchor_support = (
        lpxtg
        or gram_anchor
    )

    accessory_hits = []

    for label, pattern in (
        accessory_patterns.items()
    ):

        if re.search(
            pattern,
            p,
            flags=re.IGNORECASE
        ):
            accessory_hits.append(
                label
            )

    discordant = any(
        re.search(
            pattern,
            d
        )
        for pattern
        in discordant_description_patterns
    )


    # --------------------------------------------------------
    # Tier equivalente a 78c antes de considerar truncamiento
    # --------------------------------------------------------

    if (
        has_s8_s53
        and explicit_prt
        and aa >= 1200
        and anchor_support
        and not discordant
    ):
        base_tier = (
            "A_high_confidence_CEP_like"
        )

    elif (
        has_s8_s53
        and explicit_prt
        and aa >= 1200
        and len(
            accessory_hits
        ) >= 2
        and not discordant
    ):
        base_tier = (
            "B_probable_CEP_like"
        )

    elif (
        has_s8_s53
        and aa >= 1200
        and len(
            accessory_hits
        ) >= 1
    ):
        base_tier = (
            "C_large_multidomain_"
            "subtilase_non_CEP"
        )

    else:
        base_tier = (
            "D_other_subtilase"
        )


    # --------------------------------------------------------
    # Curación específica para rescate incompleto
    # --------------------------------------------------------

    partial = clean(
        row.get(
            "partial",
            ""
        )
    )

    if (
        partial != "00"
        and base_tier in {
            "A_high_confidence_CEP_like",
            "B_probable_CEP_like",
        }
    ):

        rescue_tier = (
            "T_truncated_CEP_like_"
            "candidate"
        )

    elif (
        partial == "00"
        and base_tier
        == "A_high_confidence_CEP_like"
    ):

        rescue_tier = (
            "A_high_confidence_CEP_like_"
            "complete_prediction"
        )

    elif (
        partial == "00"
        and base_tier
        == "B_probable_CEP_like"
    ):

        rescue_tier = (
            "B_probable_CEP_like_"
            "complete_prediction"
        )

    else:
        rescue_tier = base_tier


    rec = gene_output(
        row,
        "surface_proteinase",
        "CEP_subtilase_candidate"
    )

    rec.update({
        "has_S8_S53":
            int(
                has_s8_s53
            ),

        "explicit_Prt_annotation":
            int(
                explicit_prt
            ),

        "large_ge_1200aa":
            int(
                aa >= 1200
            ),

        "LPXTG_like_Cterminal":
            int(
                lpxtg
            ),

        "Gram_pos_anchor_support":
            int(
                gram_anchor
            ),

        "anchor_support_any":
            int(
                anchor_support
            ),

        "n_accessory_architecture_domains":
            len(
                accessory_hits
            ),

        "accessory_architecture_domains":
            ";".join(
                accessory_hits
            ),

        "discordant_primary_annotation":
            int(
                discordant
            ),

        "base_tier_78c_equivalent":
            base_tier,

        "rescue_tier_87":
            rescue_tier,
    })

    cep_rows.append(rec)


# ============================================================
# 5. LIPÓLISIS
# ============================================================

print(
    "Curando candidatos de lipólisis..."
)


lipolysis_rows = []


for row in genes:

    if row.get(
        "has_eggnog_annotation"
    ) != "1":
        continue

    pref = row["_name"]
    desc = row["_desc"]
    pfam_l = row["_pfam"].lower()
    ko_l = lower(
        row.get(
            "eggnog_KEGG_ko",
            ""
        )
    )
    ec_l = lower(
        row.get(
            "eggnog_EC",
            ""
        )
    )

    lipoyl_related = any(
        x in desc
        for x in (
            "lipoyl",
            "lipoate",
            "lipoic acid",
            "lipoylation",
        )
    )

    phospholipase = (
        "phospholipase" in desc
        or (
            re.match(
                r"^pl[acd]"
                r"[a-z0-9_-]*$",
                pref or ""
            )
            is not None
        )
    )

    explicit_lipase = (
        (
            re.search(
                r"\blipase\b",
                desc
            )
            is not None
            or (
                "triacylglycerol lipase"
                in desc
            )
            or (
                "lipase family"
                in desc
            )
            or (
                re.search(
                    r"\blipase\b",
                    pfam_l
                )
                is not None
            )
        )
        and not lipoyl_related
        and not phospholipase
    )

    if not explicit_lipase:
        continue


    if (
        pref == "lifo"
        or "lipase_chap" in pfam_l
        or (
            "folding of the extracellular "
            "lipase"
        ) in desc
    ):

        tier = (
            "D_lipase_accessory_not_enzyme"
        )

        reason = (
            "lipase_chaperone_or_"
            "folding_accessory"
        )


    elif (
        "lipase (class 3)" in desc
        or "k01046" in ko_l
        or "3.1.1.3" in ec_l
        or "secretory lipase" in desc
    ):

        tier = (
            "A_specific_lipase_candidate"
        )

        reasons = []

        if (
            "lipase (class 3)"
            in desc
        ):
            reasons.append(
                "class_3_lipase_annotation"
            )

        if "k01046" in ko_l:
            reasons.append("K01046")

        if "3.1.1.3" in ec_l:
            reasons.append("EC_3.1.1.3")

        if "secretory lipase" in desc:
            reasons.append(
                "explicit_secretory_lipase"
            )

        reason = ";".join(
            reasons
        )


    elif (
        "k12686" in ko_l
        or (
            "autotransporter"
            in pfam_l
            and "lipase_gdsl"
            in pfam_l
        )
        or (
            "cog0657 esterase lipase"
            in desc
        )
        or (
            "esterase lipase"
            in desc
        )
        or "k03928" in ko_l
        or "k03929" in ko_l
        or "k14731" in ko_l
    ):

        tier = (
            "B_esterase_lipase_candidate"
        )

        reason = (
            "specific_esterase_lipase_"
            "annotation_or_architecture"
        )


    elif (
        "gdsl" in desc
        or "lipase_gdsl" in pfam_l
    ):

        tier = (
            "C_broad_GDSL_acylhydrolase"
        )

        reason = (
            "broad_GDSL_lipase_"
            "acylhydrolase_family"
        )


    else:

        tier = (
            "E_unresolved_lipolysis_"
            "candidate"
        )

        reason = (
            "insufficient_specificity"
        )


    rec = gene_output(
        row,
        "lipolysis",
        "explicit_lipase_candidate"
    )

    rec.update({
        "curated_lipolysis_tier":
            tier,

        "curation_reason":
            reason,
    })

    lipolysis_rows.append(rec)


# ============================================================
# 6. AMINOÁCIDOS / AROMA
# ============================================================

print(
    "Curando marcadores de aminoácidos/aroma..."
)


aroma_rows = []


for row in genes:

    if row.get(
        "has_eggnog_annotation"
    ) != "1":
        continue

    n = row["_name"]
    d = row["_desc"]
    ko = lower(
        row.get(
            "eggnog_KEGG_ko",
            ""
        )
    )
    ec = lower(
        row.get(
            "eggnog_EC",
            ""
        )
    )

    marker = ""

    if (
        n in {
            "ilve",
            "bcat",
            "bca_t",
            "bcaat"
        }
        or (
            "branched-chain amino acid "
            "aminotransferase"
        ) in d
        or (
            "branched chain amino acid "
            "aminotransferase"
        ) in d
    ):
        marker = (
            "branched_chain_"
            "aminotransferase"
        )

    elif (
        n in {
            "arat",
            "tyrb"
        }
        or (
            "aromatic amino acid "
            "aminotransferase"
        ) in d
        or (
            "aromatic-amino-acid "
            "aminotransferase"
        ) in d
        or (
            "aromatic-amino-acid "
            "transaminase"
        ) in d
    ):
        marker = (
            "aromatic_amino_acid_"
            "aminotransferase"
        )

    elif (
        n in {
            "mdea",
            "mgl"
        }
        or (
            "methionine gamma-lyase"
            in d
        )
        or (
            "methionine gamma lyase"
            in d
        )
    ):
        marker = (
            "methionine_gamma_lyase"
        )

    elif (
        "cystathionine beta-lyase"
        in d
        or (
            "cystathionine gamma-lyase"
            in d
        )
        or (
            "cystathionine beta lyase"
            in d
        )
        or (
            "cystathionine gamma lyase"
            in d
        )
    ):
        marker = (
            "sulfur_amino_acid_"
            "lyase_candidate"
        )

    elif (
        n == "kivd"
        or (
            "alpha-keto acid "
            "decarboxylase"
        ) in d
        or (
            "alpha keto acid "
            "decarboxylase"
        ) in d
        or (
            "2-ketoacid decarboxylase"
        ) in d
    ):
        marker = (
            "ketoacid_decarboxylase"
        )

    else:
        continue


    tier = ""
    interpretation = ""
    discordance = 0


    if (
        marker
        == "branched_chain_"
           "aminotransferase"
    ):

        if "k00826" in ko:

            tier = (
                "A_BCAT_IlvE_supported"
            )

            interpretation = (
                "branched_chain_amino_"
                "acid_transamination"
            )

        elif (
            "k02619" in ko
            or n == "pabc"
        ):

            tier = (
                "X_excluded_PabC_like"
            )

            interpretation = (
                "not_counted_as_BCAT_"
                "aroma_marker"
            )

        else:

            tier = (
                "C_BCAT_annotation_only"
            )

            interpretation = (
                "branched_chain_"
                "transaminase_candidate"
            )


    elif (
        marker
        == "aromatic_amino_acid_"
           "aminotransferase"
    ):

        if (
            n.startswith("arat")
            or "k00841" in ko
        ):

            tier = (
                "A_AraT_supported"
            )

            interpretation = (
                "aromatic_amino_acid_"
                "transamination"
            )

        elif (
            n.startswith("tyrb")
            or "k00832" in ko
        ):

            tier = (
                "B_TyrB_broad_aromatic_"
                "aminotransferase"
            )

            interpretation = (
                "broad_aromatic_amino_"
                "acid_transamination"
            )

        else:

            tier = (
                "C_aromatic_"
                "aminotransferase_candidate"
            )

            interpretation = (
                "aromatic_transamination_"
                "candidate"
            )


    elif (
        marker
        == "methionine_gamma_lyase"
    ):

        tier = (
            "A_methionine_gamma_"
            "lyase_supported"
        )

        interpretation = (
            "direct_sulfur_volatile_"
            "precursor_candidate"
        )


    elif (
        marker
        == "sulfur_amino_acid_"
           "lyase_candidate"
    ):

        if (
            "k01761" in ko
            or "4.4.1.11" in ec
        ):

            tier = (
                "A_K01761_methionine_"
                "gamma_lyase_like"
            )

            interpretation = (
                "direct_sulfur_volatile_"
                "precursor_candidate"
            )

            if (
                "methionine gamma-lyase"
                not in d
                and not n.startswith("mgl")
            ):
                discordance = 1

        elif (
            n.startswith("mccb")
            or "k17217" in ko
            or (
                "volatile sulfur compounds"
                in d
            )
        ):

            tier = (
                "B_cheese_relevant_"
                "CBL_like"
            )

            interpretation = (
                "sulfur_aroma_candidate"
            )

        elif (
            n.startswith("metc")
            or "k01760" in ko
        ):

            tier = (
                "C_MetC_cystathionine_"
                "beta_lyase"
            )

            interpretation = (
                "sulfur_amino_acid_"
                "metabolism_possible_"
                "aroma_role"
            )

        elif (
            n.startswith("metb")
            or "k01739" in ko
        ):

            tier = (
                "D_MetB_sulfur_amino_"
                "acid_metabolism"
            )

            interpretation = (
                "sulfur_amino_acid_"
                "metabolism_not_direct_"
                "aroma_proof"
            )

        else:

            tier = (
                "E_other_sulfur_"
                "lyase_candidate"
            )

            interpretation = (
                "uncertain_sulfur_"
                "aroma_relevance"
            )


    elif (
        marker
        == "ketoacid_decarboxylase"
    ):

        tier = (
            "A_ketoacid_decarboxylase_"
            "candidate"
        )

        interpretation = (
            "Ehrlich_route_marker_"
            "candidate"
        )


    rec = gene_output(
        row,
        "amino_acid_aroma",
        marker
    )

    rec.update({
        "curated_aroma_tier":
            tier,

        "curated_interpretation":
            interpretation,

        "annotation_discordance":
            discordance,
    })

    aroma_rows.append(rec)


# ============================================================
# 7. EPS / CÁPSULA
# ============================================================

print(
    "Reconstruyendo loci EPS/cápsula..."
)


def eps_anchor_classes(row):

    n = row["_name"]
    d = row["_desc"]

    classes = set()

    if re.match(
        r"^eps[a-z0-9_-]+$",
        n or ""
    ):
        classes.add(
            "eps_named"
        )

    if re.match(
        r"^cps[a-z0-9_-]+$",
        n or ""
    ):
        classes.add(
            "cps_named"
        )

    if n in {
        "wzx",
        "wzy",
        "wzz"
    }:
        classes.add(
            "flippase_or_polymerase"
        )

    if n in {
        "wza",
        "wzb",
        "wzc"
    }:
        classes.add(
            "capsule_export_regulation"
        )

    if (
        "exopolysaccharide"
        in d
    ):
        classes.add(
            "explicit_exopolysaccharide"
        )

    if (
        "capsular polysaccharide"
        in d
    ):
        classes.add(
            "capsular_polysaccharide"
        )

    if (
        "polysaccharide export"
        in d
        or (
            "polysaccharide polymerase"
            in d
        )
        or (
            "polysaccharide flippase"
            in d
        )
    ):
        classes.add(
            "polysaccharide_export_"
            "polymerization"
        )

    return classes


def is_glycosyltransferase(row):

    d = row["_desc"]
    p = row["_pfam"].lower()

    return (
        "glycosyltransferase"
        in d
        or "glycosyl transferase"
        in d
        or "glycosyltransf"
        in p
    )


eps_rows = []
eps_counter = 0


for cid, contig_genes in (
    genes_by_contig.items()
):

    anchors = []

    for row in contig_genes:

        if row.get(
            "has_eggnog_annotation"
        ) != "1":
            continue

        classes = (
            eps_anchor_classes(row)
        )

        if classes:

            anchors.append(
                (
                    row["_orf"],
                    row,
                    classes
                )
            )

    if not anchors:
        continue

    anchors.sort(
        key=lambda x: x[0]
    )

    groups = []
    current = [
        anchors[0]
    ]

    for anchor in anchors[1:]:

        if (
            anchor[0]
            - current[-1][0]
            <= 20
        ):
            current.append(anchor)

        else:

            groups.append(current)

            current = [anchor]

    groups.append(current)


    for group in groups:

        min_anchor = min(
            x[0]
            for x in group
        )

        max_anchor = max(
            x[0]
            for x in group
        )

        start = max(
            1,
            min_anchor - 10
        )

        end = (
            max_anchor + 10
        )

        window = [
            row
            for row in contig_genes
            if (
                start
                <= row["_orf"]
                <= end
            )
        ]

        classes = set()
        marker_gene_ids = []
        gt_gene_ids = []

        for row in window:

            cs = (
                eps_anchor_classes(row)
            )

            if cs:

                classes.update(cs)

                marker_gene_ids.append(
                    row["gene_id"]
                )

            if is_glycosyltransferase(
                row
            ):

                classes.add(
                    "glycosyltransferase"
                )

                gt_gene_ids.append(
                    row["gene_id"]
                )


        has_named = bool(
            classes.intersection({
                "eps_named",
                "cps_named",
                "explicit_exopolysaccharide",
                "capsular_polysaccharide",
            })
        )

        has_export = bool(
            classes.intersection({
                "flippase_or_polymerase",
                "capsule_export_regulation",
                "polysaccharide_export_"
                "polymerization",
            })
        )

        has_gt = (
            "glycosyltransferase"
            in classes
        )

        n_anchor_genes = len(
            set(
                marker_gene_ids
            )
        )

        if (
            has_named
            and has_export
            and has_gt
        ):

            tier = (
                "A_strong_EPS_capsule_"
                "like_locus"
            )

        elif (
            (
                has_named
                and has_gt
            )
            or (
                has_export
                and has_gt
            )
            or n_anchor_genes >= 2
        ):

            tier = (
                "B_probable_EPS_capsule_"
                "like_locus"
            )

        else:

            tier = (
                "C_single_marker_context"
            )


        has_eps = bool(
            {
                "explicit_exopolysaccharide",
                "eps_named",
            }
            & classes
        )

        has_capsule = bool(
            {
                "capsular_polysaccharide",
                "cps_named",
                "capsule_export_regulation",
            }
            & classes
        )


        if (
            has_eps
            and has_capsule
        ):

            interpretation = (
                "mixed_EPS_capsule_"
                "like_locus"
            )

        elif has_eps:

            interpretation = (
                "EPS_biosynthesis_"
                "like_locus"
            )

        elif has_capsule:

            interpretation = (
                "capsule_biosynthesis_"
                "like_locus"
            )

        elif (
            has_export
            and has_gt
        ):

            interpretation = (
                "generic_surface_"
                "polysaccharide_locus"
            )

        elif len(
            set(gt_gene_ids)
        ) >= 2:

            interpretation = (
                "glycosyltransferase_"
                "rich_context"
            )

        else:

            interpretation = (
                "single_or_weak_"
                "polysaccharide_context"
            )


        eps_counter += 1

        meta = contig_meta[cid]

        eps_rows.append({
            "locus_id":
                f"IEPSLOC{eps_counter:05d}",

            "contig_id":
                cid,

            "window_start_ORF":
                start,

            "window_end_ORF":
                end,

            "anchor_min_ORF":
                min_anchor,

            "anchor_max_ORF":
                max_anchor,

            "n_anchor_genes":
                n_anchor_genes,

            "n_glycosyltransferase_genes":
                len(
                    set(gt_gene_ids)
                ),

            "marker_classes":
                ";".join(
                    sorted(classes)
                ),

            "marker_gene_ids":
                ";".join(
                    sorted(
                        set(
                            marker_gene_ids
                        )
                    )
                ),

            "glycosyltransferase_gene_ids":
                ";".join(
                    sorted(
                        set(
                            gt_gene_ids
                        )
                    )
                ),

            "curated_tier":
                tier,

            "curated_interpretation":
                interpretation,

            "contains_partial_gene":
                int(
                    any(
                        clean(
                            x.get(
                                "partial",
                                ""
                            )
                        ) != "00"
                        for x in window
                        if (
                            x["gene_id"]
                            in set(
                                marker_gene_ids
                                + gt_gene_ids
                            )
                        )
                    )
                ),

            "strongest_incomplete_context":
                meta.get(
                    "strongest_incomplete_context",
                    ""
                ),

            "n_bin_memberships":
                meta.get(
                    "n_bin_memberships",
                    ""
                ),

            "coassemblies":
                meta.get(
                    "coassemblies",
                    ""
                ),

            "producers":
                meta.get(
                    "producers",
                    ""
                ),

            "binners":
                meta.get(
                    "binners",
                    ""
                ),

            "quality_groups":
                meta.get(
                    "quality_groups",
                    ""
                ),

            "min_member_bin_completeness":
                meta.get(
                    "min_member_bin_completeness",
                    ""
                ),

            "max_member_bin_completeness":
                meta.get(
                    "max_member_bin_completeness",
                    ""
                ),

            "min_member_bin_contamination":
                meta.get(
                    "min_member_bin_contamination",
                    ""
                ),

            "max_member_bin_contamination":
                meta.get(
                    "max_member_bin_contamination",
                    ""
                ),

            "member_raw_bin_uids":
                meta.get(
                    "member_raw_bin_uids",
                    ""
                ),
        })


# ============================================================
# 8. ESTRÉS ÁCIDO / OSMÓTICO / OXIDATIVO
# ============================================================

print(
    "Reconstruyendo sistemas de estrés..."
)


stress_markers = defaultdict(
    lambda: defaultdict(list)
)


def add_stress(
    row,
    marker
):

    stress_markers[
        row["contig_id"]
    ][
        marker
    ].append(row)


for row in genes:

    if row.get(
        "has_eggnog_annotation"
    ) != "1":
        continue

    n = row["_name"]
    d = row["_desc"]


    # Ácido
    if (
        n in {
            "gada",
            "gadb"
        }
        or (
            "glutamate decarboxylase"
            in d
        )
    ):
        add_stress(
            row,
            "glutamate_decarboxylase"
        )

    elif n == "gadc":
        add_stress(
            row,
            "glutamate_GABA_antiporter"
        )

    elif n == "adia":
        add_stress(
            row,
            "arginine_decarboxylase_AR"
        )

    elif n == "adic":
        add_stress(
            row,
            "arginine_agmatine_antiporter"
        )

    elif (
        n == "arca"
        and "arginine deiminase"
        in d
    ):
        add_stress(
            row,
            "ADI_arcA"
        )

    elif (
        n == "arcb"
        and (
            "ornithine carbamoyltransferase"
            in d
            or (
                "ornithine transcarbamylase"
                in d
            )
        )
    ):
        add_stress(
            row,
            "ADI_arcB"
        )

    elif (
        n == "arcc"
        and "carbamate kinase"
        in d
    ):
        add_stress(
            row,
            "ADI_arcC"
        )

    elif (
        n == "arcd"
        and "arginine" in d
        and "antiporter" in d
    ):
        add_stress(
            row,
            "ADI_arcD"
        )

    elif (
        n == "urea"
        and "urease" in d
    ):
        add_stress(
            row,
            "urease_A"
        )

    elif (
        n == "ureb"
        and "urease" in d
    ):
        add_stress(
            row,
            "urease_B"
        )

    elif (
        n == "urec"
        and "urease" in d
    ):
        add_stress(
            row,
            "urease_C"
        )


    # Osmótico
    if re.match(
        r"^opu[a-z0-9_-]+$",
        n or ""
    ):
        add_stress(
            row,
            "Opu_compatible_solute_transport"
        )

    if n == "bett":
        add_stress(
            row,
            "BetT_choline_transport"
        )

    if n == "prop":
        add_stress(
            row,
            "ProP_osmolyte_transport"
        )

    if n in {
        "prov",
        "prow",
        "prox"
    }:
        add_stress(
            row,
            f"ProU_{n[-1].upper()}"
        )

    if n == "beta":
        add_stress(
            row,
            "BetA"
        )

    if n == "betb":
        add_stress(
            row,
            "BetB"
        )

    if n in {
        "ecta",
        "ectb",
        "ectc",
        "ectd"
    }:
        add_stress(
            row,
            f"Ect_{n[-1].upper()}"
        )

    if n == "otsa":
        add_stress(
            row,
            "OtsA"
        )

    if n == "otsb":
        add_stress(
            row,
            "OtsB"
        )

    if (
        "glycine betaine transporter"
        in d
        or "betaine transporter"
        in d
    ):
        add_stress(
            row,
            "explicit_betaine_transport"
        )

    if "ectoine synthase" in d:
        add_stress(
            row,
            "Ect_C"
        )


    # Oxidativo
    if (
        n in {
            "kata",
            "kate",
            "katg"
        }
        or (
            re.search(
                r"\bcatalase\b",
                d
            )
            is not None
        )
    ):
        add_stress(
            row,
            "catalase"
        )

    if (
        n in {
            "soda",
            "sodb",
            "sodc"
        }
        or (
            "superoxide dismutase"
            in d
        )
    ):
        add_stress(
            row,
            "superoxide_dismutase"
        )

    if (
        n in {
            "ahpc",
            "ahpf"
        }
        or (
            "alkyl hydroperoxide reductase"
            in d
        )
    ):
        add_stress(
            row,
            "alkyl_hydroperoxide_reductase"
        )

    if (
        n == "tpx"
        or "thiol peroxidase"
        in d
    ):
        add_stress(
            row,
            "thiol_peroxidase"
        )

    if (
        n == "bcp"
        or (
            "bacterioferritin "
            "comigratory protein"
        ) in d
    ):
        add_stress(
            row,
            "peroxiredoxin_Bcp"
        )

    if n in {
        "msra",
        "msrb"
    }:
        add_stress(
            row,
            "methionine_sulfoxide_repair"
        )


stress_system_rows = []


def stress_system(
    cid,
    block,
    system,
    required,
    min_required=None,
    max_span=None
):

    mdict = stress_markers[cid]

    detected = [
        x
        for x in required
        if mdict.get(x)
    ]

    if not detected:
        return

    if min_required is None:
        min_required = len(
            required
        )

    rows = [
        row
        for marker in detected
        for row in mdict[marker]
    ]

    indices = [
        row["_orf"]
        for row in rows
    ]

    span = (
        max(indices)
        - min(indices)
        if indices
        else None
    )

    supported = (
        len(
            set(detected)
        ) >= min_required
        and (
            max_span is None
            or (
                span is not None
                and span <= max_span
            )
        )
    )

    status = (
        "system_supported_same_contig"
        if supported
        else
        "marker_evidence_only"
    )

    meta = contig_meta[cid]

    stress_system_rows.append({
        "contig_id":
            cid,

        "block":
            block,

        "system":
            system,

        "status":
            status,

        "markers_required":
            ";".join(
                required
            ),

        "markers_detected":
            ";".join(
                sorted(
                    set(detected)
                )
            ),

        "n_marker_types_detected":
            len(
                set(detected)
            ),

        "ORF_span":
            (
                ""
                if span is None
                else span
            ),

        "gene_ids":
            ";".join(
                sorted({
                    x["gene_id"]
                    for x in rows
                })
            ),

        "contains_partial_gene":
            int(
                any(
                    clean(
                        x.get(
                            "partial",
                            ""
                        )
                    ) != "00"
                    for x in rows
                )
            ),

        "strongest_incomplete_context":
            meta.get(
                "strongest_incomplete_context",
                ""
            ),

        "producers":
            meta.get(
                "producers",
                ""
            ),

        "coassemblies":
            meta.get(
                "coassemblies",
                ""
            ),

        "member_raw_bin_uids":
            meta.get(
                "member_raw_bin_uids",
                ""
            ),
    })


for cid in stress_markers:

    # Ácido
    stress_system(
        cid,
        "acid_stress",
        "glutamate_GAD_system",
        [
            "glutamate_decarboxylase",
            "glutamate_GABA_antiporter",
        ],
        min_required=2,
        max_span=10
    )

    stress_system(
        cid,
        "acid_stress",
        "arginine_decarboxylase_AR_system",
        [
            "arginine_decarboxylase_AR",
            "arginine_agmatine_antiporter",
        ],
        min_required=2,
        max_span=10
    )

    stress_system(
        cid,
        "acid_stress",
        "arginine_deiminase_ADI_system",
        [
            "ADI_arcA",
            "ADI_arcB",
            "ADI_arcC",
            "ADI_arcD",
        ],
        min_required=3,
        max_span=15
    )

    stress_system(
        cid,
        "acid_stress",
        "urease_core_system",
        [
            "urease_A",
            "urease_B",
            "urease_C",
        ],
        min_required=3,
        max_span=10
    )


    # Osmótico
    stress_system(
        cid,
        "osmotic_stress",
        "ProU_transport_system",
        [
            "ProU_V",
            "ProU_W",
            "ProU_X",
        ],
        min_required=2,
        max_span=12
    )

    stress_system(
        cid,
        "osmotic_stress",
        "glycine_betaine_BetAB",
        [
            "BetA",
            "BetB",
        ],
        min_required=2,
        max_span=10
    )

    stress_system(
        cid,
        "osmotic_stress",
        "ectoine_EctABC",
        [
            "Ect_A",
            "Ect_B",
            "Ect_C",
        ],
        min_required=3,
        max_span=12
    )

    stress_system(
        cid,
        "osmotic_stress",
        "trehalose_OtsAB",
        [
            "OtsA",
            "OtsB",
        ],
        min_required=2,
        max_span=10
    )


    # Opu cluster
    opu = stress_markers[cid].get(
        "Opu_compatible_solute_transport",
        []
    )

    if opu:

        idx = sorted(
            x["_orf"]
            for x in opu
        )

        span = (
            max(idx) - min(idx)
            if idx
            else None
        )

        meta = contig_meta[cid]

        status = (
            "Opu_transport_cluster"
            if (
                len(opu) >= 2
                and span is not None
                and span <= 20
            )
            else
            "Opu_marker"
        )

        stress_system_rows.append({
            "contig_id":
                cid,

            "block":
                "osmotic_stress",

            "system":
                "Opu_compatible_solute_transport",

            "status":
                status,

            "markers_required":
                "multiple_Opu_components",

            "markers_detected":
                "Opu_compatible_solute_transport",

            "n_marker_types_detected":
                1,

            "ORF_span":
                (
                    ""
                    if span is None
                    else span
                ),

            "gene_ids":
                ";".join(
                    sorted(
                        x["gene_id"]
                        for x in opu
                    )
                ),

            "contains_partial_gene":
                int(
                    any(
                        clean(
                            x.get(
                                "partial",
                                ""
                            )
                        ) != "00"
                        for x in opu
                    )
                ),

            "strongest_incomplete_context":
                meta.get(
                    "strongest_incomplete_context",
                    ""
                ),

            "producers":
                meta.get(
                    "producers",
                    ""
                ),

            "coassemblies":
                meta.get(
                    "coassemblies",
                    ""
                ),

            "member_raw_bin_uids":
                meta.get(
                    "member_raw_bin_uids",
                    ""
                ),
        })


    # Transportadores osmóticos individuales
    for marker, system in [
        (
            "ProP_osmolyte_transport",
            "ProP_transport"
        ),
        (
            "BetT_choline_transport",
            "BetT_precursor_transport"
        ),
        (
            "explicit_betaine_transport",
            "explicit_betaine_transport"
        ),
    ]:

        rr = stress_markers[cid].get(
            marker,
            []
        )

        if not rr:
            continue

        meta = contig_meta[cid]

        stress_system_rows.append({
            "contig_id":
                cid,

            "block":
                "osmotic_stress",

            "system":
                system,

            "status":
                "single_marker_supported",

            "markers_required":
                marker,

            "markers_detected":
                marker,

            "n_marker_types_detected":
                1,

            "ORF_span":
                0,

            "gene_ids":
                ";".join(
                    sorted(
                        x["gene_id"]
                        for x in rr
                    )
                ),

            "contains_partial_gene":
                int(
                    any(
                        clean(
                            x.get(
                                "partial",
                                ""
                            )
                        ) != "00"
                        for x in rr
                    )
                ),

            "strongest_incomplete_context":
                meta.get(
                    "strongest_incomplete_context",
                    ""
                ),

            "producers":
                meta.get(
                    "producers",
                    ""
                ),

            "coassemblies":
                meta.get(
                    "coassemblies",
                    ""
                ),

            "member_raw_bin_uids":
                meta.get(
                    "member_raw_bin_uids",
                    ""
                ),
        })


    # Oxidativo: repertorio en el mismo contig
    oxidative = [
        "catalase",
        "superoxide_dismutase",
        "alkyl_hydroperoxide_reductase",
        "thiol_peroxidase",
        "peroxiredoxin_Bcp",
        "methionine_sulfoxide_repair",
    ]

    detected = [
        x
        for x in oxidative
        if stress_markers[cid].get(x)
    ]

    if detected:

        rows = [
            r
            for x in detected
            for r in stress_markers[cid][x]
        ]

        meta = contig_meta[cid]

        status = (
            "broad_antioxidant_repertoire_"
            "same_contig"
            if len(
                set(detected)
            ) >= 3
            else
            "oxidative_marker_evidence"
        )

        stress_system_rows.append({
            "contig_id":
                cid,

            "block":
                "oxidative_stress",

            "system":
                "antioxidant_repertoire",

            "status":
                status,

            "markers_required":
                ";".join(
                    oxidative
                ),

            "markers_detected":
                ";".join(
                    sorted(
                        set(detected)
                    )
                ),

            "n_marker_types_detected":
                len(
                    set(detected)
                ),

            "ORF_span":
                (
                    max(
                        x["_orf"]
                        for x in rows
                    )
                    - min(
                        x["_orf"]
                        for x in rows
                    )
                    if rows
                    else ""
                ),

            "gene_ids":
                ";".join(
                    sorted({
                        x["gene_id"]
                        for x in rows
                    })
                ),

            "contains_partial_gene":
                int(
                    any(
                        clean(
                            x.get(
                                "partial",
                                ""
                            )
                        ) != "00"
                        for x in rows
                    )
                ),

            "strongest_incomplete_context":
                meta.get(
                    "strongest_incomplete_context",
                    ""
                ),

            "producers":
                meta.get(
                    "producers",
                    ""
                ),

            "coassemblies":
                meta.get(
                    "coassemblies",
                    ""
                ),

            "member_raw_bin_uids":
                meta.get(
                    "member_raw_bin_uids",
                    ""
                ),
        })


# ============================================================
# SALIDAS
# ============================================================

print(
    "Escribiendo resultados..."
)


write_tsv(
    OUT
    / "87A_cheese_metabolic_marker_hits.tsv",
    [
        x
        for x in marker_hits
        if x["block"]
        == "cheese_metabolism"
    ],
    COMMON_FIELDS
)


SYSTEM_FIELDS = [
    "contig_id",
    "system",
    "status",
    "markers_required",
    "markers_detected",
    "n_required_markers",
    "n_detected_markers",
    "ORF_span",
    "gene_ids",
    "contains_partial_gene",
    "strongest_incomplete_context",
    "n_bin_memberships",
    "coassemblies",
    "producers",
    "binners",
    "quality_groups",
    "min_member_bin_completeness",
    "max_member_bin_completeness",
    "min_member_bin_contamination",
    "max_member_bin_contamination",
    "member_raw_bin_uids",
]


write_tsv(
    OUT
    / "87B_cheese_metabolic_systems_on_contigs.tsv",
    metabolic_systems,
    SYSTEM_FIELDS
)


BA_FIELDS = (
    COMMON_FIELDS
    + [
        "potential_product",
        "evidence_level",
        "context_support",
        "neighbor_gene_ids_pm3",
    ]
)


write_tsv(
    OUT
    / "87C_biogenic_amine_candidates.tsv",
    ba_candidates,
    BA_FIELDS
)


write_tsv(
    OUT
    / "87D_polyamine_arg_stress_genes.tsv",
    polyamine_rows,
    COMMON_FIELDS
)


PEPTIDE_COMPONENT_FIELDS = (
    COMMON_FIELDS
    + [
        "system",
        "component",
    ]
)


write_tsv(
    OUT
    / "87E_Opp_Dpp_components.tsv",
    peptide_components,
    PEPTIDE_COMPONENT_FIELDS
)


write_tsv(
    OUT
    / "87F_Opp_Dpp_clusters.tsv",
    peptide_clusters,
    [
        "contig_id",
        "system",
        "components",
        "n_required_components",
        "ORF_span",
        "cluster_status",
        "genes",
        "contains_partial_gene",
        "strongest_incomplete_context",
        "producers",
        "coassemblies",
        "member_raw_bin_uids",
    ]
)


write_tsv(
    OUT
    / "87G_POT_Dtp_transporters.tsv",
    pot_rows,
    COMMON_FIELDS
)


PEPTIDASE_FIELDS = (
    COMMON_FIELDS
    + [
        "peptidase_type",
    ]
)


write_tsv(
    OUT
    / "87H_key_peptidases.tsv",
    peptidase_rows,
    PEPTIDASE_FIELDS
)


CEP_FIELDS = (
    COMMON_FIELDS
    + [
        "has_S8_S53",
        "explicit_Prt_annotation",
        "large_ge_1200aa",
        "LPXTG_like_Cterminal",
        "Gram_pos_anchor_support",
        "anchor_support_any",
        "n_accessory_architecture_domains",
        "accessory_architecture_domains",
        "discordant_primary_annotation",
        "base_tier_78c_equivalent",
        "rescue_tier_87",
    ]
)


write_tsv(
    OUT
    / "87I_surface_proteinase_CEP_candidates.tsv",
    cep_rows,
    CEP_FIELDS
)


LIP_FIELDS = (
    COMMON_FIELDS
    + [
        "curated_lipolysis_tier",
        "curation_reason",
    ]
)


write_tsv(
    OUT
    / "87J_curated_lipolysis_candidates.tsv",
    lipolysis_rows,
    LIP_FIELDS
)


AROMA_FIELDS = (
    COMMON_FIELDS
    + [
        "curated_aroma_tier",
        "curated_interpretation",
        "annotation_discordance",
    ]
)


write_tsv(
    OUT
    / "87K_curated_amino_acid_aroma_markers.tsv",
    aroma_rows,
    AROMA_FIELDS
)


write_tsv(
    OUT
    / "87L_EPS_capsule_loci.tsv",
    eps_rows,
    [
        "locus_id",
        "contig_id",
        "window_start_ORF",
        "window_end_ORF",
        "anchor_min_ORF",
        "anchor_max_ORF",
        "n_anchor_genes",
        "n_glycosyltransferase_genes",
        "marker_classes",
        "marker_gene_ids",
        "glycosyltransferase_gene_ids",
        "curated_tier",
        "curated_interpretation",
        "contains_partial_gene",
        "strongest_incomplete_context",
        "n_bin_memberships",
        "coassemblies",
        "producers",
        "binners",
        "quality_groups",
        "min_member_bin_completeness",
        "max_member_bin_completeness",
        "min_member_bin_contamination",
        "max_member_bin_contamination",
        "member_raw_bin_uids",
    ]
)


write_tsv(
    OUT
    / "87M_stress_systems_on_contigs.tsv",
    stress_system_rows,
    [
        "contig_id",
        "block",
        "system",
        "status",
        "markers_required",
        "markers_detected",
        "n_marker_types_detected",
        "ORF_span",
        "gene_ids",
        "contains_partial_gene",
        "strongest_incomplete_context",
        "producers",
        "coassemblies",
        "member_raw_bin_uids",
    ]
)


# ============================================================
# RESUMEN GLOBAL
# ============================================================

summary = []


def summary_add(
    metric,
    value
):
    summary.append({
        "metric": metric,
        "value": value
    })


summary_add(
    "input_genes",
    len(genes)
)

summary_add(
    "input_contigs_catalog",
    len(contig_catalog_rows)
)

summary_add(
    "new_contigs_not_exact_final18",
    len(new_contigs)
)

summary_add(
    "cheese_metabolic_marker_hits",
    sum(
        x["block"]
        == "cheese_metabolism"
        for x in marker_hits
    )
)

summary_add(
    "cheese_metabolic_system_records",
    len(metabolic_systems)
)

summary_add(
    "biogenic_amine_candidate_genes",
    len(ba_candidates)
)

summary_add(
    "polyamine_arg_stress_genes",
    len(polyamine_rows)
)

summary_add(
    "Opp_Dpp_component_genes",
    len(peptide_components)
)

summary_add(
    "Opp_Dpp_cluster_records",
    len(peptide_clusters)
)

summary_add(
    "POT_Dtp_transporter_genes",
    len(pot_rows)
)

summary_add(
    "key_peptidase_genes",
    len(peptidase_rows)
)

summary_add(
    "surface_proteinase_candidates",
    len(cep_rows)
)

summary_add(
    "specific_or_broad_lipolysis_candidates",
    len(lipolysis_rows)
)

summary_add(
    "amino_acid_aroma_marker_genes",
    len(aroma_rows)
)

summary_add(
    "EPS_capsule_loci",
    len(eps_rows)
)

summary_add(
    "stress_system_records",
    len(stress_system_rows)
)


# BA por producto
for product, n in sorted(
    Counter(
        x["potential_product"]
        for x in ba_candidates
    ).items()
):
    summary_add(
        f"BA_{product}_genes",
        n
    )


# CEP tiers
for tier, n in sorted(
    Counter(
        x["rescue_tier_87"]
        for x in cep_rows
    ).items()
):
    summary_add(
        f"CEP_{tier}",
        n
    )


# EPS tiers
for tier, n in sorted(
    Counter(
        x["curated_tier"]
        for x in eps_rows
    ).items()
):
    summary_add(
        f"EPS_{tier}",
        n
    )


write_tsv(
    OUT
    / "functional_rescue_summary.tsv",
    summary,
    [
        "metric",
        "value"
    ]
)


# ============================================================
# QC POR NIVEL DE EVIDENCIA
# ============================================================

qc_scope = Counter()


def register_scope(
    rows
):

    for row in rows:

        scope = clean(
            row.get(
                "gene_evidence_scope",
                ""
            )
        )

        if scope:
            qc_scope[scope] += 1


register_scope(
    [
        x
        for x in marker_hits
        if x["block"]
        == "cheese_metabolism"
    ]
)

register_scope(
    ba_candidates
)

register_scope(
    polyamine_rows
)

register_scope(
    peptide_components
)

register_scope(
    pot_rows
)

register_scope(
    peptidase_rows
)

register_scope(
    cep_rows
)

register_scope(
    lipolysis_rows
)

register_scope(
    aroma_rows
)


write_tsv(
    OUT
    / "QC_positive_gene_hits_by_evidence_scope.tsv",
    [
        {
            "gene_evidence_scope":
                k,
            "n_hit_records":
                v
        }
        for k, v
        in sorted(
            qc_scope.items()
        )
    ],
    [
        "gene_evidence_scope",
        "n_hit_records",
    ]
)


# ============================================================
# README
# ============================================================

with open(
    OUT / "README_step87.txt",
    "w",
    encoding="utf-8"
) as fh:

    fh.write(
"""PASO 87 - RESCATE FUNCIONAL CURADO DE CONTIGS/BINS INCOMPLETOS
================================================================

UNIDAD DE ANALISIS
------------------
La unidad primaria es el contig deduplicado ICONTIG y sus genes
predichos. Un ICONTIG puede pertenecer a múltiples bins originales.

Los bins incompletos NO se incorporan como MAGs ecológicos
independientes ni se suman como organismos diferentes.

NIVELES DE EVIDENCIA
--------------------
bin_context_50_89_positive_evidence:
    Evidencia positiva recuperada en un contig presente al menos en
    un bin con completitud de 50-89.99%.

lt50_contig_gene_centric_only:
    Evidencia positiva gene/contig-céntrica. No permite inferir
    ausencia ni reconstrucción genómica completa.

AUSENCIA
--------
Este paso NO genera conclusiones de ausencia funcional.
No recuperar un marcador en un bin incompleto no implica que
el organismo carezca de ese gen o sistema.

SISTEMAS MULTIGENICOS
---------------------
Las reconstrucciones fuertes requieren colocalización de los
componentes en el mismo ICONTIG y dentro de la ventana de ORFs
especificada.

Esto es deliberadamente más conservador que reconstruir sistemas
mediante la simple coocurrencia de genes en bins incompletos.

CEP-LIKE
--------
Se reutiliza la lógica de curación del Paso 78c.

Una predicción parcial de Prodigal que cumpliría criterios A/B no
se promueve automáticamente a la misma interpretación que una
proteína completa y se etiqueta:
T_truncated_CEP_like_candidate.

LIPOLISIS
---------
Las familias GDSL y esterasas no equivalen automáticamente a
hidrólisis de grasa láctea.

AROMA
-----
Los marcadores de transaminación y metabolismo de aminoácidos
representan potencial genómico. No demuestran producción de
compuestos volátiles.

EPS/CAPSULA
-----------
Los loci se clasifican como EPS/capsule-like y no demuestran
fenotipo de producción de EPS ni cambios de textura.

AMINAS BIOGENAS
----------------
La presencia de una descarboxilasa representa potencial genético.
Los niveles system_supported o enzyme_plus_transport_context
requieren evidencia contextual en el mismo contig.

CONTAMINACION
-------------
La procedencia de cada ICONTIG conserva los bins originales y sus
rangos de contaminación. En contextos contaminados, la atribución
debe mantenerse a nivel de contig y no extrapolarse automáticamente
al taxón/bin completo.

BACTERIOCINAS / RiPP
--------------------
NO se consideran cerradas en este paso.

Los precursores pequeños y ORFs no llamados por Prodigal requieren
un análisis nucleotídico/sORF separado.

NOVEDAD
-------
"Nuevo" significa solamente que la secuencia del ICONTIG no era
exactamente idéntica a un contig de los 18 MAGs finales.

No significa especie nueva, población nueva ni función novedosa.

INTERPRETACION GENERAL
----------------------
Todos los resultados representan potencial genómico recuperado.
No demuestran expresión, actividad enzimática ni fenotipo.
"""
    )


# ============================================================
# CONSOLA
# ============================================================

print()
print("=" * 70)
print("PASO 87 COMPLETADO")
print("=" * 70)

for rec in summary:
    print(
        f"{rec['metric']:<55} "
        f"{rec['value']}"
    )

print()
print(
    f"Salida: {OUT}"
)
print()
print(
    "PASO 87 FINALIZÓ CORRECTAMENTE; "
    "ES SEGURO SALIR."
)
