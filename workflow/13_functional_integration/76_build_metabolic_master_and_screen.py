#!/usr/bin/env python3

import csv
import os
import re
from collections import defaultdict

USER = os.environ["USER"]

ROOT = f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"

MASTER = os.path.join(
    ROOT,
    "22_final_representative_mags",
    "representative_mags_master.tsv"
)

GENES = os.path.join(
    ROOT,
    "53_metabolic_annotation_consolidated",
    "gene_function_inventory_all_predicted.tsv"
)

FEATURE_COUNTS = os.path.join(
    ROOT,
    "53_metabolic_annotation_consolidated",
    "functional_feature_counts_by_MAG.tsv"
)

BACT = os.path.join(
    ROOT,
    "51_bacteriocin_final_evidence",
    "final_summary_by_locus.tsv"
)

OUTDIR = os.path.join(
    ROOT,
    "54_metabolic_functional_integration"
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


def fnum(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


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


# ============================================================
# INPUT
# ============================================================

master = read_tsv(
    MASTER
)

genes = read_tsv(
    GENES
)

feature_counts = read_tsv(
    FEATURE_COUNTS
)

bact = read_tsv(
    BACT
)


if len(master) != 18:
    raise SystemExit(
        f"ERROR: master tiene {len(master)} MAGs; esperados 18."
    )


MAG_ORDER = [
    row["representative_MAG"]
    for row in master
]

MAG_SET = set(
    MAG_ORDER
)

if len(MAG_SET) != 18:
    raise SystemExit(
        "ERROR: MAGs duplicados en master."
    )


# ============================================================
# QC CALIDAD DE MAGs
# ============================================================

quality_rows = []

for row in master:

    mag = row[
        "representative_MAG"
    ]

    comp = fnum(
        row["completeness"]
    )

    contam = fnum(
        row["contamination"]
    )

    if comp is None or contam is None:
        raise SystemExit(
            f"ERROR: calidad ausente para {mag}"
        )

    passes = int(
        comp >= 90
        and contam <= 5
    )

    quality_rows.append({
        "MAG":
            mag,

        "completeness":
            comp,

        "contamination":
            contam,

        "passes_90_5":
            passes
    })


if not all(
    x["passes_90_5"] == 1
    for x in quality_rows
):
    raise SystemExit(
        "ERROR: algún representante no pasa 90/5."
    )


write_tsv(
    os.path.join(
        OUTDIR,
        "QC_MAG_quality.tsv"
    ),
    quality_rows,
    [
        "MAG",
        "completeness",
        "contamination",
        "passes_90_5"
    ]
)


# ============================================================
# COUNTS FUNCIONALES POR MAG
# ============================================================

functional = defaultdict(dict)

for row in feature_counts:

    mag = row["MAG"]

    if mag not in MAG_SET:
        raise SystemExit(
            f"ERROR: MAG desconocido en funciones: {mag}"
        )

    ftype = row[
        "feature_type"
    ]

    functional[
        mag
    ][
        f"{ftype}_unique"
    ] = row[
        "unique_features"
    ]

    functional[
        mag
    ][
        f"{ftype}_genes"
    ] = row[
        "genes_with_feature"
    ]


# ============================================================
# LOCI BACTERIOCÍNICOS INDEPENDIENTES POR MAG
# ============================================================

bacteriocin_loci = defaultdict(
    list
)

bacteriocin_candidates = defaultdict(
    list
)


for row in bact:

    mag = row[
        "source_MAG"
    ]

    if mag not in MAG_SET:
        raise SystemExit(
            "ERROR: MAG fuente de bacteriocina "
            f"no está entre los 18: {mag}"
        )

    bacteriocin_loci[
        mag
    ].append(
        row["locus_id"]
    )

    bacteriocin_candidates[
        mag
    ].append(
        row["candidate"]
    )


# ============================================================
# TABLA MAESTRA
# ============================================================

FEATURE_TYPES = [
    "KEGG_KO",
    "KEGG_Module",
    "KEGG_Pathway",
    "KEGG_Reaction",
    "EC",
    "CAZy",
]


master_out = []

for row in master:

    mag = row[
        "representative_MAG"
    ]

    rec = dict(
        row
    )

    rec[
        "MAG"
    ] = mag

    for ft in FEATURE_TYPES:

        rec[
            f"{ft}_unique"
        ] = functional[
            mag
        ].get(
            f"{ft}_unique",
            "0"
        )

        rec[
            f"{ft}_genes"
        ] = functional[
            mag
        ].get(
            f"{ft}_genes",
            "0"
        )

    loci = sorted(
        bacteriocin_loci[
            mag
        ]
    )

    candidates = sorted(
        bacteriocin_candidates[
            mag
        ]
    )

    rec[
        "n_independent_bacteriocin_loci"
    ] = len(
        loci
    )

    rec[
        "independent_bacteriocin_loci"
    ] = ";".join(
        loci
    )

    rec[
        "bacteriocin_candidates"
    ] = ";".join(
        candidates
    )

    master_out.append(
        rec
    )


extra_fields = [
    "MAG"
]

for ft in FEATURE_TYPES:
    extra_fields.extend([
        f"{ft}_unique",
        f"{ft}_genes"
    ])

extra_fields.extend([
    "n_independent_bacteriocin_loci",
    "independent_bacteriocin_loci",
    "bacteriocin_candidates"
])


write_tsv(
    os.path.join(
        OUTDIR,
        "MAG_functional_master.tsv"
    ),
    master_out,
    list(master[0].keys())
    + extra_fields
)


# ============================================================
# PANELES DE SCREENING
#
# IMPORTANTE:
# Son búsquedas sensibles por Preferred_name + Description.
# Un hit = candidato funcional.
# NO implica automáticamente una vía completa.
# ============================================================

PANELS = {

    "lactose_galactose": [
        r"\blacz\b",
        r"\blacg\b",
        r"\blacy\b",
        r"\blacs\b",
        r"\blace\b",
        r"\blacf\b",
        r"\bgalk\b",
        r"\bgalt\b",
        r"\bgale\b",
        r"\bgalm\b",
        r"beta[- ]galactosidase",
        r"phospho[- ]?beta[- ]galactosidase",
        r"lactose[- ]specific",
        r"lactose permease",
        r"galactokinase",
        r"galactose[- ]1[- ]phosphate",
        r"tagatose"
    ],

    "lactate_pyruvate_fermentation": [
        r"\bldh[a-z0-9]*\b",
        r"lactate dehydrogenase",
        r"pyruvate oxidase",
        r"pyruvate formate[- ]lyase",
        r"\bpox[a-z0-9]*\b",
        r"\bpfl[a-z0-9]*\b",
        r"phosphotransacetylase",
        r"\bpta\b",
        r"acetate kinase",
        r"\backa\b"
    ],

    "citrate_acetoin_diacetyl": [
        r"\bcitp\b",
        r"\bcit[a-z0-9]*\b",
        r"citrate transporter",
        r"citrate lyase",
        r"oxaloacetate decarboxylase",
        r"acetolactate synthase",
        r"alpha[- ]acetolactate",
        r"acetolactate decarboxylase",
        r"acetoin reductase",
        r"butanediol dehydrogenase",
        r"2,3[- ]butanediol"
    ],

    "proteolysis_peptide_utilization": [
        r"\bprtp\b",
        r"\bpepn\b",
        r"\bpepc\b",
        r"\bpepx\b",
        r"\bpepo\b",
        r"\bpepf\b",
        r"\bpepq\b",
        r"\bopp[a-z0-9]*\b",
        r"\bdpp[a-z0-9]*\b",
        r"aminopeptidase",
        r"endopeptidase",
        r"oligopeptide",
        r"dipeptide",
        r"tripeptide",
        r"peptidase"
    ],

    "amino_acid_flavor": [
        r"branched[- ]chain amino acid aminotransferase",
        r"aromatic amino acid aminotransferase",
        r"aminotransferase",
        r"methionine gamma[- ]lyase",
        r"methionine.*lyase",
        r"cysteine.*lyase"
    ],

    "lipolysis_esterases": [
        r"\blipase\b",
        r"\besterase\b",
        r"lipolytic",
        r"phospholipase",
        r"carboxylesterase"
    ],

    "exopolysaccharide_EPS": [
        r"exopolysaccharide",
        r"extracellular polysaccharide",
        r"\beps[a-z0-9]*\b",
        r"polysaccharide polymerase",
        r"polysaccharide export",
        r"\bwzx\b",
        r"\bwzy\b",
        r"glycosyltransferase"
    ],

    "acid_stress": [
        r"glutamate decarboxylase",
        r"\bgad[a-z0-9]*\b",
        r"arginine deiminase",
        r"\barc[a-z0-9]*\b",
        r"urease",
        r"proton[- ]translocating.*atp",
        r"f0f1.*atp",
        r"f-type.*atp"
    ],

    "osmotic_salt_stress": [
        r"glycine betaine",
        r"betaine transporter",
        r"compatible solute",
        r"\bopu[a-z0-9]*\b",
        r"\bbett\b",
        r"\bprop\b",
        r"sodium.*proton antiporter",
        r"Na\+/H\+.*antiporter"
    ],

    "oxidative_stress": [
        r"\bkat[a-z0-9]*\b",
        r"catalase",
        r"superoxide dismutase",
        r"\bsod[a-z0-9]*\b",
        r"peroxidase",
        r"thioredoxin",
        r"glutaredoxin",
        r"peroxiredoxin"
    ],

    "biogenic_amine_potential": [
        r"histidine decarboxylase",
        r"\bhdc[a-z0-9]*\b",
        r"tyrosine decarboxylase",
        r"\btdc[a-z0-9]*\b",
        r"ornithine decarboxylase",
        r"\bodc[a-z0-9]*\b",
        r"lysine decarboxylase",
        r"\bcada\b",
        r"arginine decarboxylase",
        r"agmatine deiminase"
    ]
}


COMPILED = {
    panel: [
        re.compile(
            pattern,
            re.IGNORECASE
        )
        for pattern in patterns
    ]
    for panel, patterns in PANELS.items()
}


# ============================================================
# SCREENING
# ============================================================

candidate_rows = []

counts = defaultdict(
    lambda: defaultdict(set)
)


for row in genes:

    if row[
        "has_eggnog_record"
    ] != "1":
        continue

    mag = row[
        "MAG"
    ]

    if mag not in MAG_SET:
        raise SystemExit(
            f"ERROR: MAG desconocido en inventario: {mag}"
        )

    preferred = clean(
        row.get(
            "Preferred_name"
        )
    )

    description = clean(
        row.get(
            "Description"
        )
    )

    text = (
        preferred
        + " "
        + description
    )

    for panel, regexes in COMPILED.items():

        matched = []

        for rgx in regexes:

            if rgx.search(
                text
            ):
                matched.append(
                    rgx.pattern
                )

        if not matched:
            continue

        counts[
            mag
        ][
            panel
        ].add(
            row["gene_id"]
        )

        candidate_rows.append({
            "MAG":
                mag,

            "panel":
                panel,

            "gene_id":
                row[
                    "gene_id"
                ],

            "Preferred_name":
                preferred,

            "Description":
                description,

            "KEGG_ko":
                clean(
                    row.get(
                        "KEGG_ko"
                    )
                ),

            "EC":
                clean(
                    row.get(
                        "EC"
                    )
                ),

            "KEGG_Pathway":
                clean(
                    row.get(
                        "KEGG_Pathway"
                    )
                ),

            "KEGG_Module":
                clean(
                    row.get(
                        "KEGG_Module"
                    )
                ),

            "CAZy":
                clean(
                    row.get(
                        "CAZy"
                    )
                ),

            "matched_patterns":
                ";".join(
                    matched
                )
        })


candidate_rows.sort(
    key=lambda x: (
        MAG_ORDER.index(
            x["MAG"]
        ),
        x["panel"],
        x["gene_id"]
    )
)


write_tsv(
    os.path.join(
        OUTDIR,
        "cheese_metabolism_candidate_genes.tsv"
    ),
    candidate_rows,
    [
        "MAG",
        "panel",
        "gene_id",
        "Preferred_name",
        "Description",
        "KEGG_ko",
        "EC",
        "KEGG_Pathway",
        "KEGG_Module",
        "CAZy",
        "matched_patterns"
    ]
)


# ============================================================
# MATRIZ DE NÚMERO DE CANDIDATOS
# ============================================================

PANEL_ORDER = list(
    PANELS.keys()
)


matrix_rows = []

for mag in MAG_ORDER:

    rec = {
        "MAG":
            mag
    }

    for panel in PANEL_ORDER:

        rec[
            panel
        ] = len(
            counts[
                mag
            ][
                panel
            ]
        )

    matrix_rows.append(
        rec
    )


write_tsv(
    os.path.join(
        OUTDIR,
        "cheese_metabolism_candidate_counts_18x11.tsv"
    ),
    matrix_rows,
    [
        "MAG",
        *PANEL_ORDER
    ]
)


# ============================================================
# LONG SUMMARY MAG × PANEL
# ============================================================

long_rows = []

for mag in MAG_ORDER:

    for panel in PANEL_ORDER:

        n = len(
            counts[
                mag
            ][
                panel
            ]
        )

        long_rows.append({
            "MAG":
                mag,

            "panel":
                panel,

            "n_candidate_genes":
                n,

            "candidate_signal":
                int(
                    n > 0
                )
        })


write_tsv(
    os.path.join(
        OUTDIR,
        "cheese_metabolism_candidate_summary_long.tsv"
    ),
    long_rows,
    [
        "MAG",
        "panel",
        "n_candidate_genes",
        "candidate_signal"
    ]
)


# ============================================================
# PREVALENCIA DE SCREENING
# ============================================================

prevalence_rows = []

for panel in PANEL_ORDER:

    mags_positive = [
        mag
        for mag in MAG_ORDER
        if len(
            counts[
                mag
            ][
                panel
            ]
        ) > 0
    ]

    n_genes = sum(
        len(
            counts[
                mag
            ][
                panel
            ]
        )
        for mag in MAG_ORDER
    )

    prevalence_rows.append({
        "panel":
            panel,

        "n_MAGs_with_candidate":
            len(
                mags_positive
            ),

        "percent_MAGs":
            f"{100*len(mags_positive)/18:.2f}",

        "n_candidate_genes":
            n_genes,

        "MAGs":
            ";".join(
                mags_positive
            )
    })


write_tsv(
    os.path.join(
        OUTDIR,
        "cheese_metabolism_screening_prevalence.tsv"
    ),
    prevalence_rows,
    [
        "panel",
        "n_MAGs_with_candidate",
        "percent_MAGs",
        "n_candidate_genes",
        "MAGs"
    ]
)


# ============================================================
# FOUR BACTERIOCIN-CONTEXT MAGs
# ============================================================

BACT_SOURCE_MAGS = {
    "L2__L2_maxbin2.004_sub",
    "L3__concoct_29",
    "M2__M2_maxbin2.004"
}

# L2__L2_maxbin2.011_sub is retained separately as
# competitive/shared context control.
COMPETITIVE_CONTROL_MAG = (
    "L2__L2_maxbin2.011_sub"
)


focused = [
    row
    for row in candidate_rows
    if (
        row["MAG"]
        in BACT_SOURCE_MAGS
        or row["MAG"]
        == COMPETITIVE_CONTROL_MAG
    )
]


write_tsv(
    os.path.join(
        OUTDIR,
        "bacteriocin_context_MAGs_metabolic_candidates.tsv"
    ),
    focused,
    [
        "MAG",
        "panel",
        "gene_id",
        "Preferred_name",
        "Description",
        "KEGG_ko",
        "EC",
        "KEGG_Pathway",
        "KEGG_Module",
        "CAZy",
        "matched_patterns"
    ]
)


# ============================================================
# README
# ============================================================

with open(
    os.path.join(
        OUTDIR,
        "README_step76.txt"
    ),
    "w"
) as f:

    f.write(
        "Paso 76 - Integracion funcional y screening "
        "de genes relevantes para queso\n"
    )

    f.write(
        "====================================================\n\n"
    )

    f.write(
        "La tabla MAG_functional_master.tsv integra "
        "calidad, taxonomia, anotacion funcional y "
        "loci bacteriocinicos independientes.\n\n"
    )

    f.write(
        "El screening metabolico usa Preferred_name "
        "y Description de eggNOG para recuperar "
        "genes candidatos.\n\n"
    )

    f.write(
        "IMPORTANTE: un candidate_signal NO significa "
        "que la ruta metabolica este completa.\n"
    )

    f.write(
        "Los candidatos deben ser revisados y "
        "clasificados por funcion/ruta en el paso 77.\n"
    )


# ============================================================
# CONSOLA
# ============================================================

completenesses = [
    float(
        row["completeness"]
    )
    for row in master
]

contaminations = [
    float(
        row["contamination"]
    )
    for row in master
]


print(
    "============================================================"
)

print(
    "PASO 76 COMPLETADO"
)

print(
    "============================================================"
)

print(
    f"MAGs integrados:              {len(master)}"
)

print(
    f"Completitud mínima:           {min(completenesses):.2f}%"
)

print(
    f"Completitud máxima:           {max(completenesses):.2f}%"
)

print(
    f"Contaminación máxima:         {max(contaminations):.2f}%"
)

print(
    f"Genes candidatos recuperados: {len(candidate_rows)}"
)

print()

print(
    "SCREENING POR PANEL"
)

for row in prevalence_rows:

    print(
        f"{row['panel']:34s}"
        f"{row['n_MAGs_with_candidate']:2d}/18 MAGs  "
        f"{row['n_candidate_genes']:4d} genes"
    )

print()

print(
    f"Salida: {OUTDIR}"
)
