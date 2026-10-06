#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import csv
import math
import os
import sys

# ============================================================
# PASO 80b
# Proyección del potencial funcional recuperado de los MAGs
# sobre las 18 muestras shotgun.
#
# PRINCIPIO:
#   Una capacidad funcional solo contribuye a una muestra si:
#
#   1. el MAG posee evidencia funcional curada para esa capacidad;
#   2. el MAG está detectado en esa muestra según Step 64/65.
#
# Se generan:
#   - proyección PRINCIPAL
#   - proyección ESTRICTA
#   - análisis de SENSIBILIDAD funcional
#
# IMPORTANTE:
#   * RA = abundancia relativa DENTRO de los MAGs recuperados/
#     aceptados, NO de toda la microbiota.
#   * No representa expresión génica.
#   * No representa flujo metabólico.
#   * No representa abundancia directa de genes/pathways.
#   * Las bacteriocinas NO se infieren desde el MAG fuente:
#     la evidencia locus-específica se mantiene separada.
# ============================================================

USER = os.environ["USER"]

ROOT = Path(
    f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"
)

FUNC_FILE = (
    ROOT
    / "61_MAG_functional_master"
    / "MAG_functional_master_18MAGs.tsv"
)

ABUND_FILE = (
    ROOT
    / "43_MAG_abundance_taxonomy"
    / "MAG_abundance_taxonomy_long.tsv"
)

OUT = (
    ROOT
    / "62_sample_functional_projection"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# HELPERS
# ============================================================

def die(msg):
    print(
        f"ERROR: {msg}",
        file=sys.stderr
    )
    sys.exit(1)


def read_tsv(path):
    if not path.exists():
        die(
            f"No existe archivo requerido: {path}"
        )

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        if reader.fieldnames is None:
            die(
                f"Archivo sin encabezado: {path}"
            )

        rows = list(reader)

    return reader.fieldnames, rows


def require_columns(
    fieldnames,
    required,
    source_name
):
    missing = [
        x
        for x in required
        if x not in fieldnames
    ]

    if missing:
        die(
            f"{source_name}: faltan columnas: "
            + ", ".join(missing)
        )


def num(x):
    x = str(x).strip()

    if (
        x == ""
        or x.upper() == "NA"
        or x.lower() == "nan"
    ):
        return 0.0

    try:
        return float(x)
    except ValueError:
        die(
            f"No se pudo convertir a número: {x}"
        )


def integer(x):
    return int(
        round(
            num(x)
        )
    )


def fmt(x):
    return f"{float(x):.10g}"


def write_tsv(
    path,
    rows,
    fields
):
    with path.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as fh:

        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=fields,
            extrasaction="ignore",
            lineterminator="\n"
        )

        writer.writeheader()
        writer.writerows(rows)


def split_states(x):
    return {
        value
        for value in x
        if value is not None
    }


# ============================================================
# REGLAS FUNCIONALES
#
# primary_states:
#   evidencia que aceptaremos como soporte principal.
#
# sensitivity_states:
#   incluye estados parciales/adicionales cuando tienen
#   interpretación biológica defendible.
#
# No se utiliza un score global.
# ============================================================

RULES = [

    # --------------------------------------------------------
    # CARBOHIDRATOS / FERMENTACIÓN
    # --------------------------------------------------------

    {
        "feature":
            "galactose_Leloir_core",
        "domain":
            "carbohydrate_fermentation",
        "field":
            "galactose_Leloir",
        "primary_states":
            split_states([
                "complete_core"
            ]),
        "sensitivity_states":
            split_states([
                "complete_core",
                "partial"
            ]),
    },

    {
        "feature":
            "lactose_transport_betaGal",
        "domain":
            "carbohydrate_fermentation",
        "field":
            "lactose_transport_betaGal",
        "primary_states":
            split_states([
                "supported"
            ]),
        "sensitivity_states":
            split_states([
                "supported",
                "partial"
            ]),
    },

    {
        "feature":
            "lactose_PTS_LacEFG",
        "domain":
            "carbohydrate_fermentation",
        "field":
            "lactose_PTS_LacEFG",
        "primary_states":
            split_states([
                "complete_marker_set"
            ]),
        "sensitivity_states":
            split_states([
                "complete_marker_set",
                "partial"
            ]),
    },

    {
        "feature":
            "lactate_LDH",
        "domain":
            "carbohydrate_fermentation",
        "field":
            "lactate_LDH",
        "primary_states":
            split_states([
                "supported"
            ]),
        "sensitivity_states":
            split_states([
                "supported"
            ]),
    },

    {
        "feature":
            "acetate_PtaAckA",
        "domain":
            "fermentation_end_products",
        "field":
            "acetate_PtaAckA",
        "primary_states":
            split_states([
                "complete_marker_set"
            ]),
        "sensitivity_states":
            split_states([
                "complete_marker_set",
                "partial"
            ]),
    },

    {
        "feature":
            "formate_PFL",
        "domain":
            "fermentation_end_products",
        "field":
            "formate_PFL",
        "primary_states":
            split_states([
                "complete_marker_set"
            ]),
        "sensitivity_states":
            split_states([
                "complete_marker_set"
            ]),
    },

    {
        "feature":
            "citrate_fermentation",
        "domain":
            "citrate_aroma",
        "field":
            "citrate_fermentation",
        "primary_states":
            split_states([
                "strong_system"
            ]),
        "sensitivity_states":
            split_states([
                "strong_system",
                "partial"
            ]),
    },

    {
        "feature":
            "acetoin_branch",
        "domain":
            "citrate_aroma",
        "field":
            "acetoin_branch",
        "primary_states":
            split_states([
                "supported_marker_pair"
            ]),
        "sensitivity_states":
            split_states([
                "supported_marker_pair",
                "partial"
            ]),
    },

    {
        "feature":
            "butanediol_branch",
        "domain":
            "citrate_aroma",
        "field":
            "butanediol_branch",
        "primary_states":
            split_states([
                "supported"
            ]),
        "sensitivity_states":
            split_states([
                "supported"
            ]),
    },


    # --------------------------------------------------------
    # AMINAS BIÓGENAS
    # --------------------------------------------------------

    {
        "feature":
            "histamine_potential",
        "domain":
            "biogenic_amines",
        "field":
            "histamine",
        "primary_states":
            split_states([
                "enzyme_plus_transport_context"
            ]),
        "sensitivity_states":
            split_states([
                "enzyme_plus_transport_context"
            ]),
    },

    {
        "feature":
            "tyramine_potential",
        "domain":
            "biogenic_amines",
        "field":
            "tyramine",
        "primary_states":
            split_states([
                "enzyme_plus_transport_context"
            ]),
        "sensitivity_states":
            split_states([
                "enzyme_plus_transport_context"
            ]),
    },

    {
        "feature":
            "cadaverine_potential",
        "domain":
            "biogenic_amines",
        "field":
            "cadaverine",
        "primary_states":
            split_states([
                "system_supported"
            ]),
        "sensitivity_states":
            split_states([
                "system_supported"
            ]),
    },

    {
        "feature":
            "putrescine_potential",
        "domain":
            "biogenic_amines",
        "field":
            "putrescine",
        "primary_states":
            split_states([
                "system_supported"
            ]),
        "sensitivity_states":
            split_states([
                "system_supported",
                "enzyme_only"
            ]),
    },


    # --------------------------------------------------------
    # PROTEÓLISIS
    # --------------------------------------------------------

    {
        "feature":
            "peptide_utilization_system",
        "domain":
            "proteolysis",
        "field":
            "peptide_utilization_evidence",
        "primary_states":
            split_states([
                "strong_complete_transport_plus_peptidases"
            ]),
        "sensitivity_states":
            split_states([
                "strong_complete_transport_plus_peptidases",
                "supported_transport_plus_peptidases"
            ]),
    },

    {
        "feature":
            "CEP_like_surface_proteinase",
        "domain":
            "proteolysis",
        "field":
            "best_surface_proteinase_evidence",
        "primary_states":
            split_states([
                "high_confidence_CEP_like_recovered"
            ]),
        "sensitivity_states":
            split_states([
                "high_confidence_CEP_like_recovered",
                "probable_CEP_like_recovered"
            ]),
    },


    # --------------------------------------------------------
    # LIPÓLISIS
    # --------------------------------------------------------

    {
        "feature":
            "lipolytic_enzyme_candidate",
        "domain":
            "lipolysis",
        "field":
            "lipolysis_evidence_curated",
        "primary_states":
            split_states([
                "specific_lipase_candidate_recovered"
            ]),
        "sensitivity_states":
            split_states([
                "specific_lipase_candidate_recovered",
                "esterase_lipase_candidate_recovered"
            ]),
    },


    # --------------------------------------------------------
    # AROMA DE AMINOÁCIDOS
    # --------------------------------------------------------

    {
        "feature":
            "amino_acid_transamination_aroma_potential",
        "domain":
            "amino_acid_aroma",
        "field":
            "amino_acid_aroma_evidence_curated",
        "primary_states":
            split_states([
                "sulfur_aroma_candidate_plus_amino_acid_transamination",
                "aromatic_amino_acid_transamination_potential",
                "branched_chain_transamination_only"
            ]),
        "sensitivity_states":
            split_states([
                "sulfur_aroma_candidate_plus_amino_acid_transamination",
                "aromatic_amino_acid_transamination_potential",
                "branched_chain_transamination_only"
            ]),
    },

    {
        "feature":
            "sulfur_aroma_candidate",
        "domain":
            "amino_acid_aroma",
        "field":
            "amino_acid_aroma_evidence_curated",
        "primary_states":
            split_states([
                "sulfur_aroma_candidate_plus_amino_acid_transamination"
            ]),
        "sensitivity_states":
            split_states([
                "sulfur_aroma_candidate_plus_amino_acid_transamination"
            ]),
    },


    # --------------------------------------------------------
    # EPS / CÁPSULA
    #
    # Se mantienen separadas porque cápsula != EPS tecnológico.
    # --------------------------------------------------------

    {
        "feature":
            "EPS_like_biosynthesis_context",
        "domain":
            "surface_polysaccharides",
        "field":
            "EPS_capsule_evidence_curated",
        "primary_states":
            split_states([
                "EPS_biosynthesis_like_locus_recovered",
                "mixed_EPS_capsule_like_locus_recovered",
                "separate_EPS_like_and_capsule_like_loci"
            ]),
        "sensitivity_states":
            split_states([
                "EPS_biosynthesis_like_locus_recovered",
                "mixed_EPS_capsule_like_locus_recovered",
                "separate_EPS_like_and_capsule_like_loci"
            ]),
    },

    {
        "feature":
            "capsule_like_context",
        "domain":
            "surface_polysaccharides",
        "field":
            "EPS_capsule_evidence_curated",
        "primary_states":
            split_states([
                "capsule_biosynthesis_like_locus_recovered",
                "mixed_EPS_capsule_like_locus_recovered",
                "separate_EPS_like_and_capsule_like_loci"
            ]),
        "sensitivity_states":
            split_states([
                "capsule_biosynthesis_like_locus_recovered",
                "mixed_EPS_capsule_like_locus_recovered",
                "separate_EPS_like_and_capsule_like_loci"
            ]),
    },


    # --------------------------------------------------------
    # ESTRÉS
    # --------------------------------------------------------

    {
        "feature":
            "acid_stress_system",
        "domain":
            "stress_adaptation",
        "field":
            "acid_stress_evidence",
        "primary_states":
            split_states([
                "system_supported"
            ]),
        "sensitivity_states":
            split_states([
                "system_supported",
                "marker_only"
            ]),
    },

    {
        "feature":
            "osmoadaptation_multigene_system",
        "domain":
            "stress_adaptation",
        "field":
            "osmotic_stress_evidence_curated",
        "primary_states":
            split_states([
                "strong_multigene_osmoadaptation_system"
            ]),
        "sensitivity_states":
            split_states([
                "strong_multigene_osmoadaptation_system"
            ]),
    },

    {
        "feature":
            "oxidative_stress_repertoire",
        "domain":
            "stress_adaptation",
        "field":
            "oxidative_stress_evidence",
        "primary_states":
            split_states([
                "broad_antioxidant_repertoire"
            ]),
        "sensitivity_states":
            split_states([
                "broad_antioxidant_repertoire",
                "supported"
            ]),
    },
]


# ============================================================
# LEER ARCHIVOS
# ============================================================

func_fields, func_rows = read_tsv(
    FUNC_FILE
)

ab_fields, ab_rows = read_tsv(
    ABUND_FILE
)


# ============================================================
# VALIDACIÓN DE COLUMNAS
# ============================================================

require_columns(
    func_fields,
    [
        "MAG",
        "genus",
        "species",
        "completeness",
        "contamination",
        "bacteriocin_analysis_role",
        "independent_bacteriocin_loci_n",
        "independent_bacteriocin_loci",
    ]
    + sorted({
        r["field"]
        for r in RULES
    }),
    "MAG_functional_master_18MAGs.tsv"
)

require_columns(
    ab_fields,
    [
        "sample",
        "producer",
        "biological_unit",
        "week",
        "MAG",
        "detection_class",
        "main_detection",
        "strict_detection",
        "relative_abundance_main_pct",
        "relative_abundance_strict_pct",
    ],
    "MAG_abundance_taxonomy_long.tsv"
)


# ============================================================
# ÍNDICE FUNCIONAL
# ============================================================

func = {}

for r in func_rows:

    mag = r["MAG"].strip()

    if mag in func:
        die(
            f"MAG duplicado en master funcional: {mag}"
        )

    func[mag] = r


if len(func) != 18:
    die(
        f"Se esperaban 18 MAGs funcionales; "
        f"se encontraron {len(func)}"
    )


# ============================================================
# VALIDAR TABLA DE ABUNDANCIA
# ============================================================

if len(ab_rows) != 324:
    die(
        f"Se esperaban 324 combinaciones "
        f"(18 muestras x 18 MAGs); "
        f"se encontraron {len(ab_rows)}"
    )


ab_mags = {
    r["MAG"].strip()
    for r in ab_rows
}

if ab_mags != set(func):
    die(
        "Los MAGs de abundancia no coinciden "
        "con los 18 MAGs funcionales."
    )


samples = {}

for r in ab_rows:

    sample = r["sample"].strip()

    meta = (
        r["producer"].strip(),
        r["biological_unit"].strip(),
        r["week"].strip()
    )

    if sample not in samples:
        samples[sample] = meta

    elif samples[sample] != meta:
        die(
            f"Metadata inconsistente para {sample}"
        )


if len(samples) != 18:
    die(
        f"Se esperaban 18 muestras; "
        f"se encontraron {len(samples)}"
    )


def numeric_sort_value(x):
    try:
        return int(float(x))
    except Exception:
        return 999999


sample_order = sorted(
    samples,
    key=lambda s: (
        samples[s][0],
        numeric_sort_value(
            samples[s][1]
        ),
        numeric_sort_value(
            samples[s][2]
        ),
        s
    )
)


# ============================================================
# ÍNDICE SAMPLE × MAG
# ============================================================

ab_lookup = {}

for r in ab_rows:

    key = (
        r["sample"].strip(),
        r["MAG"].strip()
    )

    if key in ab_lookup:
        die(
            f"Combinación duplicada: {key}"
        )

    ab_lookup[key] = r


for sample in sample_order:

    n = sum(
        1
        for key in ab_lookup
        if key[0] == sample
    )

    if n != 18:
        die(
            f"{sample}: se esperaban 18 MAGs; "
            f"se encontraron {n}"
        )


# ============================================================
# VALIDAR NORMALIZACIONES DE ABUNDANCIA
# ============================================================

for sample in sample_order:

    rows = [
        ab_lookup[
            (sample, mag)
        ]
        for mag in sorted(func)
    ]

    main_detected = sum(
        integer(
            r["main_detection"]
        )
        for r in rows
    )

    strict_detected = sum(
        integer(
            r["strict_detection"]
        )
        for r in rows
    )

    main_sum = sum(
        num(
            r["relative_abundance_main_pct"]
        )
        for r in rows
    )

    strict_sum = sum(
        num(
            r["relative_abundance_strict_pct"]
        )
        for r in rows
    )

    if main_detected > 0:
        if not math.isclose(
            main_sum,
            100.0,
            abs_tol=1e-5
        ):
            die(
                f"{sample}: RA main suma "
                f"{main_sum}, no ~100"
            )

    if strict_detected > 0:
        if not math.isclose(
            strict_sum,
            100.0,
            abs_tol=1e-5
        ):
            die(
                f"{sample}: RA strict suma "
                f"{strict_sum}, no ~100"
            )


# ============================================================
# EXPORTAR REGLAS
# ============================================================

rule_rows = []

for rule in RULES:

    rule_rows.append({
        "feature":
            rule["feature"],
        "domain":
            rule["domain"],
        "source_field":
            rule["field"],
        "primary_states":
            ";".join(
                sorted(
                    rule["primary_states"]
                )
            ),
        "sensitivity_states":
            ";".join(
                sorted(
                    rule["sensitivity_states"]
                )
            ),
    })


write_tsv(
    OUT / "functional_projection_rules.tsv",
    rule_rows,
    [
        "feature",
        "domain",
        "source_field",
        "primary_states",
        "sensitivity_states"
    ]
)


# ============================================================
# PROYECCIÓN SAMPLE × FEATURE
# ============================================================

projection = []

for sample in sample_order:

    producer, biological_unit, week = (
        samples[sample]
    )

    for rule in RULES:

        primary_main = []
        primary_strict = []

        sensitivity_main = []
        sensitivity_strict = []

        primary_RA_main = 0.0
        primary_RA_strict = 0.0

        sensitivity_RA_main = 0.0
        sensitivity_RA_strict = 0.0

        for mag in sorted(func):

            f = func[mag]
            a = ab_lookup[
                (sample, mag)
            ]

            state = (
                f[
                    rule["field"]
                ].strip()
            )

            main_detected = (
                integer(
                    a["main_detection"]
                ) == 1
            )

            strict_detected = (
                integer(
                    a["strict_detection"]
                ) == 1
            )

            ra_main = num(
                a[
                    "relative_abundance_main_pct"
                ]
            )

            ra_strict = num(
                a[
                    "relative_abundance_strict_pct"
                ]
            )

            # --------------------------------------------
            # EVIDENCIA PRINCIPAL
            # --------------------------------------------

            if (
                state
                in rule["primary_states"]
            ):

                if main_detected:
                    primary_main.append(
                        mag
                    )
                    primary_RA_main += (
                        ra_main
                    )

                if strict_detected:
                    primary_strict.append(
                        mag
                    )
                    primary_RA_strict += (
                        ra_strict
                    )

            # --------------------------------------------
            # SENSIBILIDAD FUNCIONAL
            # --------------------------------------------

            if (
                state
                in rule["sensitivity_states"]
            ):

                if main_detected:
                    sensitivity_main.append(
                        mag
                    )
                    sensitivity_RA_main += (
                        ra_main
                    )

                if strict_detected:
                    sensitivity_strict.append(
                        mag
                    )
                    sensitivity_RA_strict += (
                        ra_strict
                    )

        projection.append({

            "sample":
                sample,

            "producer":
                producer,

            "biological_unit":
                biological_unit,

            "week":
                week,

            "functional_domain":
                rule["domain"],

            "feature":
                rule["feature"],

            "source_field":
                rule["field"],

            "primary_supported_main":
                int(
                    len(primary_main) > 0
                ),

            "n_primary_carriers_main":
                len(primary_main),

            "primary_carriers_main":
                ";".join(
                    primary_main
                ),

            "primary_carrier_RA_main_pct":
                fmt(
                    primary_RA_main
                ),

            "primary_supported_strict":
                int(
                    len(primary_strict) > 0
                ),

            "n_primary_carriers_strict":
                len(primary_strict),

            "primary_carriers_strict":
                ";".join(
                    primary_strict
                ),

            "primary_carrier_RA_strict_pct":
                fmt(
                    primary_RA_strict
                ),

            "sensitivity_supported_main":
                int(
                    len(
                        sensitivity_main
                    ) > 0
                ),

            "n_sensitivity_carriers_main":
                len(
                    sensitivity_main
                ),

            "sensitivity_carriers_main":
                ";".join(
                    sensitivity_main
                ),

            "sensitivity_carrier_RA_main_pct":
                fmt(
                    sensitivity_RA_main
                ),

            "sensitivity_supported_strict":
                int(
                    len(
                        sensitivity_strict
                    ) > 0
                ),

            "n_sensitivity_carriers_strict":
                len(
                    sensitivity_strict
                ),

            "sensitivity_carriers_strict":
                ";".join(
                    sensitivity_strict
                ),

            "sensitivity_carrier_RA_strict_pct":
                fmt(
                    sensitivity_RA_strict
                ),
        })


LONG_FIELDS = [

    "sample",
    "producer",
    "biological_unit",
    "week",

    "functional_domain",
    "feature",
    "source_field",

    "primary_supported_main",
    "n_primary_carriers_main",
    "primary_carriers_main",
    "primary_carrier_RA_main_pct",

    "primary_supported_strict",
    "n_primary_carriers_strict",
    "primary_carriers_strict",
    "primary_carrier_RA_strict_pct",

    "sensitivity_supported_main",
    "n_sensitivity_carriers_main",
    "sensitivity_carriers_main",
    "sensitivity_carrier_RA_main_pct",

    "sensitivity_supported_strict",
    "n_sensitivity_carriers_strict",
    "sensitivity_carriers_strict",
    "sensitivity_carrier_RA_strict_pct",
]


write_tsv(
    OUT
    / "sample_functional_projection_long.tsv",
    projection,
    LONG_FIELDS
)


# ============================================================
# MATRICES SAMPLE × FEATURE
# ============================================================

feature_order = [
    r["feature"]
    for r in RULES
]

proj_lookup = {
    (
        r["sample"],
        r["feature"]
    ): r
    for r in projection
}


def write_matrix(
    field,
    filename
):

    path = (
        OUT
        / filename
    )

    with path.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as fh:

        writer = csv.writer(
            fh,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow(
            ["sample"]
            + feature_order
        )

        for sample in sample_order:

            writer.writerow(
                [sample]
                + [
                    proj_lookup[
                        (
                            sample,
                            feature
                        )
                    ][field]
                    for feature
                    in feature_order
                ]
            )


write_matrix(
    "primary_supported_main",
    "functional_primary_presence_main_18x23.tsv"
)

write_matrix(
    "primary_supported_strict",
    "functional_primary_presence_strict_18x23.tsv"
)

write_matrix(
    "primary_carrier_RA_main_pct",
    "functional_primary_carrier_RA_main_18x23.tsv"
)

write_matrix(
    "primary_carrier_RA_strict_pct",
    "functional_primary_carrier_RA_strict_18x23.tsv"
)

write_matrix(
    "sensitivity_supported_main",
    "functional_sensitivity_presence_main_18x23.tsv"
)

write_matrix(
    "sensitivity_carrier_RA_main_pct",
    "functional_sensitivity_carrier_RA_main_18x23.tsv"
)


# ============================================================
# PREVALENCIA POR FEATURE ENTRE LAS 18 MUESTRAS
# ============================================================

prev_rows = []

for rule in RULES:

    feature = rule["feature"]

    rr = [
        r
        for r in projection
        if r["feature"] == feature
    ]

    n_primary_main = sum(
        int(
            r[
                "primary_supported_main"
            ]
        )
        for r in rr
    )

    n_primary_strict = sum(
        int(
            r[
                "primary_supported_strict"
            ]
        )
        for r in rr
    )

    n_sens_main = sum(
        int(
            r[
                "sensitivity_supported_main"
            ]
        )
        for r in rr
    )

    mean_ra_main = sum(
        num(
            r[
                "primary_carrier_RA_main_pct"
            ]
        )
        for r in rr
    ) / 18.0

    mean_ra_strict = sum(
        num(
            r[
                "primary_carrier_RA_strict_pct"
            ]
        )
        for r in rr
    ) / 18.0

    prev_rows.append({

        "functional_domain":
            rule["domain"],

        "feature":
            feature,

        "samples_primary_main":
            n_primary_main,

        "percent_samples_primary_main":
            fmt(
                100.0
                * n_primary_main
                / 18.0
            ),

        "samples_primary_strict":
            n_primary_strict,

        "percent_samples_primary_strict":
            fmt(
                100.0
                * n_primary_strict
                / 18.0
            ),

        "samples_sensitivity_main":
            n_sens_main,

        "percent_samples_sensitivity_main":
            fmt(
                100.0
                * n_sens_main
                / 18.0
            ),

        "mean_primary_carrier_RA_main_pct":
            fmt(
                mean_ra_main
            ),

        "mean_primary_carrier_RA_strict_pct":
            fmt(
                mean_ra_strict
            ),
    })


write_tsv(
    OUT
    / "functional_projection_prevalence.tsv",
    prev_rows,
    [
        "functional_domain",
        "feature",
        "samples_primary_main",
        "percent_samples_primary_main",
        "samples_primary_strict",
        "percent_samples_primary_strict",
        "samples_sensitivity_main",
        "percent_samples_sensitivity_main",
        "mean_primary_carrier_RA_main_pct",
        "mean_primary_carrier_RA_strict_pct",
    ]
)


# ============================================================
# PROYECCIÓN DE LOS MAGs DEL ANÁLISIS DE BACTERIOCINAS
#
# ATENCIÓN:
# Esto NO significa que el locus bacteriocina esté presente.
# Es únicamente la detección/RA del MAG de referencia.
# ============================================================

bact_rows = []

focus_mags = [
    mag
    for mag, row in func.items()
    if (
        row[
            "bacteriocin_analysis_role"
        ].strip()
        != "none"
    )
]

for sample in sample_order:

    producer, biological_unit, week = (
        samples[sample]
    )

    for mag in sorted(
        focus_mags
    ):

        f = func[mag]
        a = ab_lookup[
            (
                sample,
                mag
            )
        ]

        bact_rows.append({

            "sample":
                sample,

            "producer":
                producer,

            "biological_unit":
                biological_unit,

            "week":
                week,

            "MAG":
                mag,

            "genus":
                f["genus"],

            "species":
                f["species"],

            "bacteriocin_analysis_role":
                f[
                    "bacteriocin_analysis_role"
                ],

            "independent_bacteriocin_loci_n":
                f[
                    "independent_bacteriocin_loci_n"
                ],

            "independent_bacteriocin_loci":
                f[
                    "independent_bacteriocin_loci"
                ],

            "detection_class":
                a[
                    "detection_class"
                ],

            "main_detection":
                a[
                    "main_detection"
                ],

            "strict_detection":
                a[
                    "strict_detection"
                ],

            "relative_abundance_main_pct":
                a[
                    "relative_abundance_main_pct"
                ],

            "relative_abundance_strict_pct":
                a[
                    "relative_abundance_strict_pct"
                ],

            "interpretation":
                (
                    "source_MAG_signal_only_"
                    "NOT_direct_locus_detection"
                ),
        })


write_tsv(
    OUT
    / "bacteriocin_focus_MAG_sample_projection.tsv",
    bact_rows,
    [
        "sample",
        "producer",
        "biological_unit",
        "week",
        "MAG",
        "genus",
        "species",
        "bacteriocin_analysis_role",
        "independent_bacteriocin_loci_n",
        "independent_bacteriocin_loci",
        "detection_class",
        "main_detection",
        "strict_detection",
        "relative_abundance_main_pct",
        "relative_abundance_strict_pct",
        "interpretation",
    ]
)


# ============================================================
# QC PRINCIPAL VS ESTRICTO
# ============================================================

qc_rows = []

for rule in RULES:

    feature = rule["feature"]

    rr = [
        r
        for r in projection
        if r["feature"] == feature
    ]

    discordant = [
        r["sample"]
        for r in rr
        if (
            r[
                "primary_supported_main"
            ]
            !=
            r[
                "primary_supported_strict"
            ]
        )
    ]

    qc_rows.append({

        "feature":
            feature,

        "n_samples_primary_main":
            sum(
                int(
                    r[
                        "primary_supported_main"
                    ]
                )
                for r in rr
            ),

        "n_samples_primary_strict":
            sum(
                int(
                    r[
                        "primary_supported_strict"
                    ]
                )
                for r in rr
            ),

        "n_main_strict_discordant":
            len(
                discordant
            ),

        "discordant_samples":
            ";".join(
                discordant
            ),
    })


write_tsv(
    OUT
    / "QC_function_projection_main_vs_strict.tsv",
    qc_rows,
    [
        "feature",
        "n_samples_primary_main",
        "n_samples_primary_strict",
        "n_main_strict_discordant",
        "discordant_samples",
    ]
)


# ============================================================
# README
# ============================================================

README = (
    OUT
    / "README_step80b.txt"
)

with README.open(
    "w",
    encoding="utf-8"
) as fh:

    fh.write(
        "PASO 80b - PROYECCION FUNCIONAL MAG -> MUESTRA\n"
        "================================================\n\n"

        "OBJETIVO\n"
        "--------\n"
        "Proyectar las capacidades funcionales curadas de los "
        "18 MAGs sobre las 18 muestras shotgun utilizando "
        "unicamente MAGs con evidencia de deteccion aceptada.\n\n"

        "INTERPRETACION DE ABUNDANCIA\n"
        "----------------------------\n"
        "relative_abundance_main_pct y relative_abundance_strict_pct "
        "son abundancias relativas dentro del conjunto de MAGs "
        "recuperados y aceptados de cada muestra. NO representan "
        "porcentaje de toda la microbiota.\n\n"

        "PRIMARY vs SENSITIVITY\n"
        "----------------------\n"
        "PRIMARY utiliza solamente los estados funcionales curados "
        "considerados suficientemente fuertes para interpretacion "
        "principal.\n\n"

        "SENSITIVITY incorpora determinados estados parciales cuando "
        "son biologicamente informativos, pero estos no sustituyen "
        "la interpretacion principal.\n\n"

        "ARCHIVOS PRINCIPALES\n"
        "--------------------\n"
        "sample_functional_projection_long.tsv\n"
        "  Tabla central muestra x capacidad funcional.\n\n"

        "functional_primary_presence_main_18x23.tsv\n"
        "  Presencia funcional principal basada en deteccion main.\n\n"

        "functional_primary_carrier_RA_main_18x23.tsv\n"
        "  Suma de abundancia relativa de MAGs portadores de cada "
        "capacidad en cada muestra.\n\n"

        "functional_primary_presence_strict_18x23.tsv\n"
        "  Analisis con deteccion estricta de MAGs.\n\n"

        "functional_projection_prevalence.tsv\n"
        "  Resumen descriptivo de prevalencia por capacidad.\n\n"

        "functional_projection_rules.tsv\n"
        "  Reglas exactas utilizadas para convertir estados "
        "funcionales de Step 80a en soporte primary/sensitivity.\n\n"

        "bacteriocin_focus_MAG_sample_projection.tsv\n"
        "  Presencia y abundancia de los MAGs vinculados al analisis "
        "de bacteriocinas. NO equivale a deteccion del locus.\n\n"

        "LIMITACIONES\n"
        "------------\n"
        "1. Esta es una proyeccion de POTENCIAL GENOMICO recuperado.\n"
        "2. No demuestra expresion o actividad metabolica.\n"
        "3. La suma de RA de MAGs portadores no es abundancia de "
        "genes ni de metabolitos.\n"
        "4. Una misma MAG puede aportar a multiples capacidades; "
        "por tanto las columnas funcionales no son composicionales "
        "y no deben sumarse entre si.\n"
        "5. Las bacteriocinas requieren su evidencia locus-especifica "
        "del analisis de reclutamiento competitivo.\n"
        "6. La ausencia de una capacidad en la proyeccion puede "
        "resultar de ausencia/no recuperacion del MAG o del marcador "
        "funcional y no demuestra ausencia absoluta en la microbiota.\n"
    )


# ============================================================
# RESUMEN FINAL
# ============================================================

print(
    "=" * 60
)

print(
    "PASO 80b COMPLETADO"
)

print(
    "=" * 60
)

print(
    f"Muestras integradas:              "
    f"{len(sample_order)}"
)

print(
    f"MAGs integrados:                  "
    f"{len(func)}"
)

print(
    f"Capacidades funcionales:          "
    f"{len(RULES)}"
)

print(
    f"Combinaciones muestra-función:    "
    f"{len(projection)}"
)

print(
    f"MAGs foco bacteriocina:           "
    f"{len(focus_mags)}"
)

print(
    f"Filas foco bacteriocina:          "
    f"{len(bact_rows)}"
)

print()

print(
    f"Salida: {OUT}"
)

print(
    "PASO 80b FINALIZÓ CORRECTAMENTE; "
    "ES SEGURO SALIR."
)
