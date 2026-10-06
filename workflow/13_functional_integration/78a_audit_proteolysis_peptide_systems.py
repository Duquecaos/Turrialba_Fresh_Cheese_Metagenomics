#!/usr/bin/env python3

import csv
import os
import re
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
    "56_proteolysis_peptide_audit"
)

os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# IO
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


genes = read_tsv(GENES)
master = read_tsv(MASTER)

MAG_ORDER = [
    x["representative_MAG"]
    for x in master
]

MAG_SET = set(MAG_ORDER)

taxonomy = {
    x["representative_MAG"]: x
    for x in master
}


# ============================================================
# HELPERS
# ============================================================

def norm(x):
    return (
        x or ""
    ).strip().lower()


def pname(row):
    return norm(
        row.get(
            "Preferred_name",
            ""
        )
    )


def desc(row):
    return norm(
        row.get(
            "Description",
            ""
        )
    )


def split_ko(x):
    if not x:
        return set()

    return {
        y.strip()
        for y in x.split(",")
        if y.strip()
    }


# ============================================================
# CLASIFICACIÓN EXPLORATORIA
#
# IMPORTANTE:
# esto identifica candidatos.
# NO declara todavía capacidad proteolítica completa.
# ============================================================

def classify(row):

    n = pname(row)
    d = desc(row)

    categories = set()

    # --------------------------------------------------------
    # 1. PROTEASAS EXTRACELULARES / ENVOLTURA
    # --------------------------------------------------------

    strong_surface_names = {
        "prtp",
        "prth",
        "prt",
        "prtbh"
    }

    if n in strong_surface_names:
        categories.add(
            "surface_or_extracellular_proteinase"
        )

    if (
        "cell envelope proteinase" in d
        or "cell-envelope proteinase" in d
        or "caseinolytic extracellular" in d
        or "caseinase" in d
        or "caseinolytic proteinase" in d
        or "extracellular proteinase" in d
        or "extracellular protease" in d
    ):
        categories.add(
            "surface_or_extracellular_proteinase"
        )

    # Subtilisinas secretadas pueden ser relevantes,
    # pero se mantienen separadas porque no equivalen a PrtP.
    if (
        "subtilisin" in d
        and (
            "secret" in d
            or "extracellular" in d
            or "cell wall" in d
        )
    ):
        categories.add(
            "other_extracellular_protease"
        )


    # --------------------------------------------------------
    # 2. TRANSPORTE DE OLIGOPÉPTIDOS
    # --------------------------------------------------------

    if re.fullmatch(
        r"opp[a-f][0-9_-]*",
        n
    ):
        categories.add(
            "oligopeptide_transport_Opp"
        )

    if (
        "oligopeptide" in d
        and (
            "transport" in d
            or "permease" in d
            or "abc" in d
            or "binding protein" in d
        )
    ):
        categories.add(
            "oligopeptide_transport_Opp"
        )


    # --------------------------------------------------------
    # 3. DI/TRIPÉPTIDOS
    # --------------------------------------------------------

    if re.fullmatch(
        r"dpp[a-f][0-9_-]*",
        n
    ):
        categories.add(
            "dipeptide_transport_Dpp"
        )

    if n in {
        "dtpt",
        "dtpa",
        "dtpb",
        "pept1"
    }:
        categories.add(
            "di_tripeptide_transport"
        )

    if (
        (
            "dipeptide" in d
            or "tripeptide" in d
            or "di- and tripeptide" in d
            or "di and tripeptide" in d
        )
        and (
            "transport" in d
            or "permease" in d
            or "symporter" in d
            or "binding protein" in d
        )
    ):
        categories.add(
            "di_tripeptide_transport"
        )


    # --------------------------------------------------------
    # 4. PEPTIDASAS INTRACELULARES CLÁSICAS
    # --------------------------------------------------------

    pep_names = {
        "pepn",
        "pepc",
        "pepx",
        "pepo",
        "pepf",
        "pepq",
        "pepv",
        "pept",
        "pepd",
        "pepa",
        "pepm",
        "pepb",
        "pepe",
        "pepi",
        "pepr",
        "pepz",
        "pepda"
    }

    if n in pep_names:
        categories.add(
            "intracellular_peptidase_named"
        )

    if (
        "aminopeptidase" in d
        or "dipeptidase" in d
        or "tripeptidase" in d
        or "endopeptidase pep" in d
        or "proline iminopeptidase" in d
        or "prolidase" in d
        or "proline dipeptidase" in d
        or "x-prolyl" in d
        or "x-prolyl dipeptidyl" in d
    ):
        categories.add(
            "intracellular_peptidase_description"
        )


    # --------------------------------------------------------
    # 5. EXCLUSIONES HOUSEKEEPING
    #
    # Se marcan para auditoría, pero NO se incluyen como
    # evidencia tecnológica de proteólisis láctea.
    # --------------------------------------------------------

    housekeeping = False

    if n in {
        "clpp",
        "clpx",
        "clpa",
        "clpb",
        "clpc",
        "clpe",
        "lon",
        "ftsH".lower()
    }:
        housekeeping = True

    if (
        "atp-dependent protease" in d
        or "protein quality control" in d
        or "proteasome" in d
    ):
        housekeeping = True

    if housekeeping:
        categories.add(
            "housekeeping_protease_excluded"
        )


    return categories


# ============================================================
# RECUPERAR CANDIDATOS
# ============================================================

hits = []

for row in genes:

    mag = row["MAG"]

    if mag not in MAG_SET:
        continue

    cats = classify(row)

    if not cats:
        continue

    m = taxonomy[mag]

    for cat in sorted(cats):

        hits.append({
            "MAG":
                mag,

            "genus":
                m.get(
                    "genus",
                    ""
                ),

            "species":
                m.get(
                    "species",
                    ""
                ),

            "category":
                cat,

            "gene_id":
                row["gene_id"],

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

            "PFAMs":
                row.get(
                    "PFAMs",
                    ""
                )
        })


write_tsv(
    os.path.join(
        OUTDIR,
        "proteolysis_candidate_gene_audit.tsv"
    ),
    hits,
    [
        "MAG",
        "genus",
        "species",
        "category",
        "gene_id",
        "Preferred_name",
        "Description",
        "KEGG_ko",
        "EC",
        "PFAMs"
    ]
)


# ============================================================
# MATRIZ DE CONTEOS
# ============================================================

categories = [
    "surface_or_extracellular_proteinase",
    "other_extracellular_protease",
    "oligopeptide_transport_Opp",
    "dipeptide_transport_Dpp",
    "di_tripeptide_transport",
    "intracellular_peptidase_named",
    "intracellular_peptidase_description",
    "housekeeping_protease_excluded"
]

counts = defaultdict(
    lambda: defaultdict(set)
)

for row in hits:

    counts[
        row["MAG"]
    ][
        row["category"]
    ].add(
        row["gene_id"]
    )


matrix = []

for mag in MAG_ORDER:

    m = taxonomy[mag]

    rec = {
        "MAG":
            mag,

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

    for cat in categories:

        rec[cat] = len(
            counts[
                mag
            ][
                cat
            ]
        )

    matrix.append(
        rec
    )


write_tsv(
    os.path.join(
        OUTDIR,
        "proteolysis_candidate_counts_18MAGs.tsv"
    ),
    matrix,
    [
        "MAG",
        "genus",
        "species",
        *categories
    ]
)


# ============================================================
# PREVALENCIA POR CATEGORÍA
# ============================================================

prev = []

for cat in categories:

    mags = [
        mag
        for mag in MAG_ORDER
        if len(
            counts[
                mag
            ][
                cat
            ]
        ) > 0
    ]

    unique_genes = {
        row["gene_id"]
        for row in hits
        if row["category"] == cat
    }

    prev.append({
        "category":
            cat,

        "n_MAGs":
            len(mags),

        "n_genes":
            len(unique_genes),

        "MAGs":
            ";".join(
                mags
            )
    })


write_tsv(
    os.path.join(
        OUTDIR,
        "proteolysis_candidate_prevalence.tsv"
    ),
    prev,
    [
        "category",
        "n_MAGs",
        "n_genes",
        "MAGs"
    ]
)


# ============================================================
# NOMBRES / ANOTACIONES ÚNICAS PARA AUDITAR
# ============================================================

unique_rows = []
seen = set()

for row in hits:

    if row[
        "category"
    ] == "housekeeping_protease_excluded":
        continue

    key = (
        row["category"],
        row["Preferred_name"],
        row["Description"],
        row["KEGG_ko"],
        row["EC"]
    )

    if key in seen:
        continue

    seen.add(key)

    unique_rows.append({
        "category":
            row["category"],

        "Preferred_name":
            row["Preferred_name"],

        "Description":
            row["Description"],

        "KEGG_ko":
            row["KEGG_ko"],

        "EC":
            row["EC"]
    })


unique_rows.sort(
    key=lambda x: (
        x["category"],
        x["Preferred_name"],
        x["Description"]
    )
)


write_tsv(
    os.path.join(
        OUTDIR,
        "proteolysis_unique_annotations.tsv"
    ),
    unique_rows,
    [
        "category",
        "Preferred_name",
        "Description",
        "KEGG_ko",
        "EC"
    ]
)


# ============================================================
# MAGs VINCULADOS A LOCI BACTERIOCÍNICOS
# ============================================================

focus_mags = {
    "L2__L2_maxbin2.004_sub",
    "L3__concoct_29",
    "M2__M2_maxbin2.004",
    "L2__L2_maxbin2.011_sub"
}

focus_rows = [
    x
    for x in hits
    if (
        x["MAG"] in focus_mags
        and x["category"]
        != "housekeeping_protease_excluded"
    )
]

write_tsv(
    os.path.join(
        OUTDIR,
        "proteolysis_focus_bacteriocin_MAGs.tsv"
    ),
    focus_rows,
    [
        "MAG",
        "genus",
        "species",
        "category",
        "gene_id",
        "Preferred_name",
        "Description",
        "KEGG_ko",
        "EC",
        "PFAMs"
    ]
)


# ============================================================
# CONSOLA
# ============================================================

print("=" * 60)
print("PASO 78a COMPLETADO")
print("=" * 60)

print(
    f"MAGs analizados: {len(MAG_ORDER)}"
)

print(
    f"Hits totales de auditoría: {len(hits)}"
)

print()

for row in prev:

    print(
        f"{row['category']:40s} "
        f"{int(row['n_MAGs']):2d}/18 MAGs  "
        f"{int(row['n_genes']):4d} genes"
    )

print()
print(
    f"Salida: {OUTDIR}"
)

print(
    "PASO 78a FINALIZÓ CORRECTAMENTE; "
    "ES SEGURO SALIR."
)
