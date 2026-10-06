#!/usr/bin/env python3

import csv
import os
from collections import Counter, defaultdict

USER = os.environ["USER"]

ROOT = f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"

EGGDIR = os.path.join(
    ROOT,
    "23_eggnog_rep18"
)

ANN = os.path.join(
    EGGDIR,
    "results",
    "representative18.emapper.annotations"
)

FAA = os.path.join(
    EGGDIR,
    "results",
    "representative18.emapper.genepred.fasta"
)

SUMMARY = os.path.join(
    EGGDIR,
    "eggnog_summary_by_MAG.tsv"
)

OUTDIR = os.path.join(
    ROOT,
    "53_metabolic_annotation_consolidated"
)

os.makedirs(
    OUTDIR,
    exist_ok=True
)


# ============================================================
# HELPERS
# ============================================================

def read_tsv(path):
    with open(path, newline="") as f:
        return list(
            csv.DictReader(
                f,
                delimiter="\t"
            )
        )


def write_tsv(path, rows, fields):
    with open(
        path,
        "w",
        newline=""
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(rows)


def clean(value):
    if value is None:
        return ""

    value = str(value).strip()

    if value in {
        "",
        "-",
        "NA",
        "N/A"
    }:
        return ""

    return value


def split_tokens(value):
    value = clean(value)

    if not value:
        return []

    tokens = []

    for token in value.split(","):

        token = token.strip()

        if (
            token
            and token != "-"
        ):
            tokens.append(token)

    return sorted(
        set(tokens)
    )


def normalize_ko(token):
    token = token.strip()

    if token.startswith("ko:"):
        token = token[3:]

    return token


def normalize_reaction(token):
    token = token.strip()

    if token.startswith("rn:"):
        token = token[3:]

    return token


def normalize_pathway(token):
    """
    eggNOG puede reportar, por ejemplo:
    ko00010, map00010.

    Para evitar duplicar la misma vía se conserva
    solo el código numérico canónico cuando aplica.
    """

    token = token.strip()

    for prefix in (
        "ko",
        "map"
    ):
        if (
            token.startswith(prefix)
            and len(token) == 7
            and token[2:].isdigit()
        ):
            return token[2:]

    return token


def mag_from_gene_id(gene_id):
    if "___" not in gene_id:
        raise ValueError(
            "ID sin separador ___: "
            f"{gene_id}"
        )

    return gene_id.split(
        "___",
        1
    )[0]


# ============================================================
# VALIDAR ARCHIVOS
# ============================================================

for path in [
    ANN,
    FAA,
    SUMMARY
]:
    if (
        not os.path.isfile(path)
        or os.path.getsize(path) == 0
    ):
        raise SystemExit(
            f"ERROR: archivo ausente/vacío: {path}"
        )


# ============================================================
# 1. RESUMEN eggNOG POR MAG
# ============================================================

summary_rows = read_tsv(
    SUMMARY
)

if len(summary_rows) != 18:
    raise SystemExit(
        "ERROR: eggnog_summary_by_MAG.tsv "
        f"tiene {len(summary_rows)} MAGs; "
        "se esperaban 18."
    )


summary_lookup = {
    row["MAG"]: row
    for row in summary_rows
}

MAG_ORDER = [
    row["MAG"]
    for row in summary_rows
]

MAG_SET = set(
    MAG_ORDER
)

if len(MAG_SET) != 18:
    raise SystemExit(
        "ERROR: existen MAGs duplicados "
        "en eggnog_summary_by_MAG.tsv"
    )


# ============================================================
# 2. PROTEÍNAS PREDICHAS
# ============================================================

predicted_genes = {}
predicted_count = Counter()


with open(FAA) as f:

    for line in f:

        if not line.startswith(">"):
            continue

        gene_id = (
            line[1:]
            .split()[0]
            .strip()
        )

        mag = mag_from_gene_id(
            gene_id
        )

        if mag not in MAG_SET:
            raise SystemExit(
                "ERROR: proteína asignada a "
                f"MAG desconocido: {gene_id}"
            )

        if gene_id in predicted_genes:
            raise SystemExit(
                "ERROR: proteína duplicada: "
                f"{gene_id}"
            )

        predicted_genes[
            gene_id
        ] = mag

        predicted_count[
            mag
        ] += 1


# ============================================================
# 3. LEER .emapper.annotations
# ============================================================

header = None
annotation_by_gene = {}
annotated_count = Counter()


with open(ANN, newline="") as f:

    for raw in f:

        raw = raw.rstrip("\n")

        if not raw:
            continue

        if raw.startswith("##"):
            continue

        if raw.startswith("#query"):

            header = raw.split("\t")

            header[0] = (
                header[0]
                .lstrip("#")
            )

            continue

        if raw.startswith("#"):
            continue

        if header is None:
            raise SystemExit(
                "ERROR: no se encontró cabecera "
                "#query en annotations."
            )

        fields = raw.split("\t")

        if len(fields) != len(header):
            raise SystemExit(
                "ERROR: número inesperado de "
                "columnas en annotations."
            )

        row = dict(
            zip(
                header,
                fields
            )
        )

        gene_id = row["query"]

        if gene_id not in predicted_genes:
            raise SystemExit(
                "ERROR: gen anotado no existe "
                f"en FAA: {gene_id}"
            )

        if gene_id in annotation_by_gene:
            raise SystemExit(
                "ERROR: anotación duplicada: "
                f"{gene_id}"
            )

        mag = predicted_genes[
            gene_id
        ]

        row["MAG"] = mag

        annotation_by_gene[
            gene_id
        ] = row

        annotated_count[
            mag
        ] += 1


# ============================================================
# 4. QC CONTRA eggnog_summary_by_MAG.tsv
# ============================================================

qc_rows = []


for mag in MAG_ORDER:

    expected_pred = int(
        summary_lookup[
            mag
        ][
            "genes_predicted"
        ]
    )

    expected_ann = int(
        summary_lookup[
            mag
        ][
            "genes_with_eggnog_record"
        ]
    )

    observed_pred = (
        predicted_count[
            mag
        ]
    )

    observed_ann = (
        annotated_count[
            mag
        ]
    )

    qc_rows.append({
        "MAG":
            mag,

        "summary_genes_predicted":
            expected_pred,

        "observed_genes_predicted":
            observed_pred,

        "predicted_match":
            int(
                expected_pred
                == observed_pred
            ),

        "summary_eggnog_records":
            expected_ann,

        "observed_eggnog_records":
            observed_ann,

        "annotations_match":
            int(
                expected_ann
                == observed_ann
            ),

        "percent_with_record":
            summary_lookup[
                mag
            ][
                "percent_with_record"
            ],
    })


if not all(
    row["predicted_match"] == 1
    and row["annotations_match"] == 1
    for row in qc_rows
):
    raise SystemExit(
        "ERROR: conteos por MAG no coinciden "
        "con eggnog_summary_by_MAG.tsv"
    )


write_tsv(
    os.path.join(
        OUTDIR,
        "QC_eggnog_assignment_by_MAG.tsv"
    ),
    qc_rows,
    [
        "MAG",
        "summary_genes_predicted",
        "observed_genes_predicted",
        "predicted_match",
        "summary_eggnog_records",
        "observed_eggnog_records",
        "annotations_match",
        "percent_with_record",
    ]
)


# ============================================================
# 5. INVENTARIO DE TODOS LOS GENES
#
# Incluye también proteínas SIN registro eggNOG.
# ============================================================

SELECTED_ANNOTATION_FIELDS = [
    "seed_ortholog",
    "evalue",
    "score",
    "eggNOG_OGs",
    "max_annot_lvl",
    "COG_category",
    "Description",
    "Preferred_name",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "PFAMs",
]


gene_inventory = []


for gene_id, mag in predicted_genes.items():

    ann = annotation_by_gene.get(
        gene_id
    )

    rec = {
        "MAG":
            mag,

        "gene_id":
            gene_id,

        "has_eggnog_record":
            int(
                ann is not None
            ),
    }

    for field in SELECTED_ANNOTATION_FIELDS:

        if ann is None:
            rec[field] = ""
        else:
            rec[field] = ann.get(
                field,
                ""
            )

    gene_inventory.append(
        rec
    )


gene_inventory.sort(
    key=lambda x: (
        MAG_ORDER.index(
            x["MAG"]
        ),
        x["gene_id"]
    )
)


write_tsv(
    os.path.join(
        OUTDIR,
        "gene_function_inventory_all_predicted.tsv"
    ),
    gene_inventory,
    [
        "MAG",
        "gene_id",
        "has_eggnog_record",
        *SELECTED_ANNOTATION_FIELDS,
    ]
)


# ============================================================
# 6. EXTRAER ANOTACIONES FUNCIONALES NORMALIZADAS
# ============================================================

FEATURE_FIELDS = {
    "KEGG_KO":
        "KEGG_ko",

    "KEGG_Module":
        "KEGG_Module",

    "KEGG_Pathway":
        "KEGG_Pathway",

    "KEGG_Reaction":
        "KEGG_Reaction",

    "EC":
        "EC",

    "CAZy":
        "CAZy",
}


feature_links = {
    feature_type: set()
    for feature_type in FEATURE_FIELDS
}


for gene_id, ann in annotation_by_gene.items():

    mag = ann["MAG"]

    # ----------------------------
    # KO
    # ----------------------------

    for token in split_tokens(
        ann.get(
            "KEGG_ko"
        )
    ):

        token = normalize_ko(
            token
        )

        if token:
            feature_links[
                "KEGG_KO"
            ].add(
                (
                    mag,
                    token,
                    gene_id
                )
            )

    # ----------------------------
    # MODULE
    # ----------------------------

    for token in split_tokens(
        ann.get(
            "KEGG_Module"
        )
    ):

        feature_links[
            "KEGG_Module"
        ].add(
            (
                mag,
                token,
                gene_id
            )
        )

    # ----------------------------
    # PATHWAY
    # ----------------------------

    for token in split_tokens(
        ann.get(
            "KEGG_Pathway"
        )
    ):

        token = normalize_pathway(
            token
        )

        if token:
            feature_links[
                "KEGG_Pathway"
            ].add(
                (
                    mag,
                    token,
                    gene_id
                )
            )

    # ----------------------------
    # REACTION
    # ----------------------------

    for token in split_tokens(
        ann.get(
            "KEGG_Reaction"
        )
    ):

        token = normalize_reaction(
            token
        )

        if token:
            feature_links[
                "KEGG_Reaction"
            ].add(
                (
                    mag,
                    token,
                    gene_id
                )
            )

    # ----------------------------
    # EC
    # ----------------------------

    for token in split_tokens(
        ann.get(
            "EC"
        )
    ):

        feature_links[
            "EC"
        ].add(
            (
                mag,
                token,
                gene_id
            )
        )

    # ----------------------------
    # CAZy
    # ----------------------------

    for token in split_tokens(
        ann.get(
            "CAZy"
        )
    ):

        feature_links[
            "CAZy"
        ].add(
            (
                mag,
                token,
                gene_id
            )
        )


# ============================================================
# 7. TABLAS LONG POR TIPO DE FUNCIÓN
# ============================================================

feature_summary_global = []
feature_counts_by_mag = []


for feature_type, links in feature_links.items():

    links = sorted(
        links,
        key=lambda x: (
            MAG_ORDER.index(
                x[0]
            ),
            x[1],
            x[2]
        )
    )

    outfile = os.path.join(
        OUTDIR,
        f"{feature_type}_gene_links.tsv"
    )

    rows = [
        {
            "MAG":
                mag,
            "feature":
                feature,
            "gene_id":
                gene_id,
        }
        for mag, feature, gene_id in links
    ]

    write_tsv(
        outfile,
        rows,
        [
            "MAG",
            "feature",
            "gene_id",
        ]
    )


    # ----------------------------
    # Conteos MAG × feature
    # ----------------------------

    counts = Counter(
        (
            mag,
            feature
        )
        for mag, feature, gene_id
        in links
    )

    count_rows = []

    for (
        mag,
        feature
    ), n_genes in sorted(
        counts.items(),
        key=lambda x: (
            MAG_ORDER.index(
                x[0][0]
            ),
            x[0][1]
        )
    ):

        count_rows.append({
            "MAG":
                mag,
            "feature":
                feature,
            "n_genes":
                n_genes,
        })

    write_tsv(
        os.path.join(
            OUTDIR,
            f"{feature_type}_counts_by_MAG.tsv"
        ),
        count_rows,
        [
            "MAG",
            "feature",
            "n_genes",
        ]
    )


    unique_features = set(
        feature
        for mag, feature, gene_id
        in links
    )

    genes_with_feature = set(
        gene_id
        for mag, feature, gene_id
        in links
    )

    feature_summary_global.append({
        "feature_type":
            feature_type,

        "unique_features":
            len(
                unique_features
            ),

        "genes_with_feature":
            len(
                genes_with_feature
            ),

        "gene_feature_links":
            len(
                links
            ),
    })


    # ----------------------------
    # Número de features por MAG
    # ----------------------------

    by_mag = defaultdict(
        lambda: {
            "features": set(),
            "genes": set(),
        }
    )

    for (
        mag,
        feature,
        gene_id
    ) in links:

        by_mag[
            mag
        ][
            "features"
        ].add(
            feature
        )

        by_mag[
            mag
        ][
            "genes"
        ].add(
            gene_id
        )

    for mag in MAG_ORDER:

        feature_counts_by_mag.append({
            "MAG":
                mag,

            "feature_type":
                feature_type,

            "unique_features":
                len(
                    by_mag[
                        mag
                    ][
                        "features"
                    ]
                ),

            "genes_with_feature":
                len(
                    by_mag[
                        mag
                    ][
                        "genes"
                    ]
                ),
        })


# ============================================================
# 8. RESUMEN GLOBAL
# ============================================================

write_tsv(
    os.path.join(
        OUTDIR,
        "functional_annotation_global_summary.tsv"
    ),
    feature_summary_global,
    [
        "feature_type",
        "unique_features",
        "genes_with_feature",
        "gene_feature_links",
    ]
)


write_tsv(
    os.path.join(
        OUTDIR,
        "functional_feature_counts_by_MAG.tsv"
    ),
    feature_counts_by_mag,
    [
        "MAG",
        "feature_type",
        "unique_features",
        "genes_with_feature",
    ]
)


# ============================================================
# 9. PREVALENCIA DE CADA FEATURE ENTRE LOS 18 MAGs
# ============================================================

for feature_type, links in feature_links.items():

    feature_to_mags = defaultdict(
        set
    )

    feature_to_genes = defaultdict(
        set
    )

    for (
        mag,
        feature,
        gene_id
    ) in links:

        feature_to_mags[
            feature
        ].add(
            mag
        )

        feature_to_genes[
            feature
        ].add(
            gene_id
        )

    prevalence_rows = []

    for feature in feature_to_mags:

        prevalence_rows.append({
            "feature":
                feature,

            "n_MAGs":
                len(
                    feature_to_mags[
                        feature
                    ]
                ),

            "n_genes":
                len(
                    feature_to_genes[
                        feature
                    ]
                ),

            "MAGs":
                ";".join(
                    sorted(
                        feature_to_mags[
                            feature
                        ]
                    )
                ),
        })

    prevalence_rows.sort(
        key=lambda x: (
            -x["n_MAGs"],
            x["feature"]
        )
    )

    write_tsv(
        os.path.join(
            OUTDIR,
            f"{feature_type}_prevalence_across_MAGs.tsv"
        ),
        prevalence_rows,
        [
            "feature",
            "n_MAGs",
            "n_genes",
            "MAGs",
        ]
    )


# ============================================================
# 10. COG
#
# COG_category puede contener varias letras juntas.
# ============================================================

cog_links = set()


for gene_id, ann in annotation_by_gene.items():

    mag = ann["MAG"]

    value = clean(
        ann.get(
            "COG_category"
        )
    )

    if not value:
        continue

    for char in value:

        if char.isalpha():

            cog_links.add(
                (
                    mag,
                    char,
                    gene_id
                )
            )


cog_rows = [
    {
        "MAG":
            mag,

        "COG_category":
            cog,

        "gene_id":
            gene_id,
    }
    for mag, cog, gene_id
    in sorted(
        cog_links,
        key=lambda x: (
            MAG_ORDER.index(
                x[0]
            ),
            x[1],
            x[2]
        )
    )
]


write_tsv(
    os.path.join(
        OUTDIR,
        "COG_gene_links.tsv"
    ),
    cog_rows,
    [
        "MAG",
        "COG_category",
        "gene_id",
    ]
)


# ============================================================
# 11. README
# ============================================================

with open(
    os.path.join(
        OUTDIR,
        "README_step75.txt"
    ),
    "w"
) as f:

    f.write(
        "Paso 75 - Consolidacion de la anotacion "
        "funcional eggNOG existente\n"
    )

    f.write(
        "====================================================\n\n"
    )

    f.write(
        "No se reejecuto Prodigal ni eggNOG-mapper.\n"
    )

    f.write(
        "Los 18 MAGs se recuperan desde el prefijo "
        "del gene_id antes de '___'.\n\n"
    )

    f.write(
        "gene_function_inventory_all_predicted.tsv "
        "incluye todas las proteinas predichas, "
        "incluyendo aquellas sin registro eggNOG.\n\n"
    )

    f.write(
        "Las tablas KEGG_Module/Pathway/KO indican "
        "anotaciones asociadas a genes. "
        "NO representan por si solas rutas completas.\n\n"
    )

    f.write(
        "La completitud de rutas se evaluara "
        "en el paso siguiente utilizando los "
        "componentes requeridos de cada ruta/modulo "
        "y considerando la completitud de los MAGs.\n"
    )


# ============================================================
# 12. CONSOLA
# ============================================================

print(
    "============================================================"
)

print(
    "PASO 75 COMPLETADO"
)

print(
    "============================================================"
)

print(
    f"MAGs:                        {len(MAG_ORDER)}"
)

print(
    f"Proteínas predichas:          {len(predicted_genes)}"
)

print(
    f"Genes con registro eggNOG:    {len(annotation_by_gene)}"
)

print()

print(
    "QC conteos por MAG:           OK"
)

print()

for row in feature_summary_global:

    print(
        f"{row['feature_type']:18s}"
        f"features={row['unique_features']:6d}  "
        f"genes={row['genes_with_feature']:6d}"
    )

print()

print(
    f"Salida: {OUTDIR}"
)
