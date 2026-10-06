#!/usr/bin/env python3

# ============================================================
# PASO 82a
# Síntesis final de evidencia para los 6 loci candidatos
# bacteriocina/RiPP y sus 3 MAGs fuente independientes.
#
# NO:
# - recalcula mapping
# - recalcula inStrain
# - recalcula anotación funcional
# - genera scores arbitrarios
#
# INTEGRA:
# 1. Evidencia locus-específica Step73
# 2. Integración locus/MAG/función Step81b
# 3. Perfil funcional MAG Step80a
# 4. Trayectorias longitudinales ya curadas
#
# PRODUCTOS:
# - tabla tesis por locus
# - tabla tesis por MAG fuente
# - trayectorias locus por unidad biológica
# - contribución de funciones clave
# - resumen Markdown listo para revisión
# ============================================================

from pathlib import Path
from collections import defaultdict
import csv
import os
import sys
import math


# ============================================================
# RUTAS
# ============================================================

USER = os.environ["USER"]

ROOT = Path(
    f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"
)

STEP73 = (
    ROOT
    / "51_bacteriocin_final_evidence"
    / "final_independent_locus_evidence_18x6.tsv"
)

STEP81B_LOCUS = (
    ROOT
    / "64_bacteriocin_functional_integration"
    / "integrated_summary_by_locus.tsv"
)

STEP81B_PROFILE = (
    ROOT
    / "64_bacteriocin_functional_integration"
    / "source_MAG_primary_functional_profile_summary.tsv"
)

STEP81B_PW = (
    ROOT
    / "64_bacteriocin_functional_integration"
    / "source_MAG_locus_summary_by_producer_week.tsv"
)

STEP81B_DOM = (
    ROOT
    / "64_bacteriocin_functional_integration"
    / "source_MAG_functional_dominance_summary.tsv"
)

MASTER = (
    ROOT
    / "61_MAG_functional_master"
    / "MAG_functional_master_18MAGs.tsv"
)

OUT = (
    ROOT
    / "65_thesis_candidate_synthesis"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# CONSTANTES
# ============================================================

EXPECTED_LOCI = [
    "ATTRLOC001",
    "ATTRLOC003",
    "ATTRLOC004",
    "ATTRLOC005",
    "ATTRLOC006",
    "ATTRLOC007",
]

EXPECTED_SOURCE_MAGS = [
    "L2__L2_maxbin2.004_sub",
    "L3__concoct_29",
    "M2__M2_maxbin2.004",
]

KEY_TECH_FUNCTIONS = [
    "galactose_Leloir_core",
    "lactose_transport_betaGal",
    "lactose_PTS_LacEFG",
    "lactate_LDH",
    "acetate_PtaAckA",
    "formate_PFL",
    "citrate_fermentation",
    "acetoin_branch",
    "butanediol_branch",
    "peptide_utilization_system",
    "CEP_like_surface_proteinase",
    "lipolytic_enzyme_candidate",
    "sulfur_aroma_candidate",
    "EPS_like_biosynthesis_context",
    "acid_stress_system",
    "osmoadaptation_multigene_system",
    "oxidative_stress_repertoire",
]

FUNCTION_LABELS = {
    "galactose_Leloir_core":
        "Leloir galactose core",
    "lactose_transport_betaGal":
        "Lactose transport + beta-gal",
    "lactose_PTS_LacEFG":
        "Lactose PTS LacEFG",
    "lactate_LDH":
        "Lactate dehydrogenase",
    "acetate_PtaAckA":
        "Acetate Pta-AckA",
    "formate_PFL":
        "Pyruvate formate-lyase",
    "citrate_fermentation":
        "Citrate fermentation",
    "acetoin_branch":
        "Acetoin branch",
    "butanediol_branch":
        "2,3-butanediol branch",
    "peptide_utilization_system":
        "Peptide utilization",
    "CEP_like_surface_proteinase":
        "CEP-like surface proteinase",
    "lipolytic_enzyme_candidate":
        "Lipolytic enzyme candidate",
    "sulfur_aroma_candidate":
        "Sulfur-aroma candidate",
    "EPS_like_biosynthesis_context":
        "EPS-like biosynthesis context",
    "acid_stress_system":
        "Acid-stress system",
    "osmoadaptation_multigene_system":
        "Multigene osmoadaptation",
    "oxidative_stress_repertoire":
        "Oxidative-stress repertoire",
}


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
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n"
        )

        writer.writeheader()
        writer.writerows(rows)


def require_columns(
    fields,
    required,
    filename
):

    missing = [
        c
        for c in required
        if c not in fields
    ]

    if missing:
        die(
            f"{filename}: faltan columnas: "
            + ", ".join(missing)
        )


def fnum(x):

    x = str(x).strip()

    if (
        x == ""
        or x.upper() == "NA"
        or x.lower() == "nan"
    ):
        return math.nan

    return float(x)


def fint(x):

    x = str(x).strip()

    if x == "":
        return 0

    return int(
        round(
            float(x)
        )
    )


def fmt(x, digits=6):

    if x is None:
        return ""

    try:
        if math.isnan(float(x)):
            return ""
    except Exception:
        return str(x)

    return (
        f"{float(x):.{digits}f}"
        .rstrip("0")
        .rstrip(".")
    )


def unique_nonempty(
    rows,
    field,
    expected_one=True
):

    values = sorted({
        str(r.get(field, "")).strip()
        for r in rows
        if str(
            r.get(field, "")
        ).strip() != ""
    })

    if (
        expected_one
        and len(values) != 1
    ):
        die(
            f"Campo {field} no es invariante: "
            + ";".join(values)
        )

    if expected_one:
        return values[0]

    return values


def md_escape(x):

    return (
        str(x)
        .replace("|", "\\|")
        .replace("\n", " ")
    )


# ============================================================
# CARGA
# ============================================================

f73, r73 = read_tsv(
    STEP73
)

f_locus, r_locus = read_tsv(
    STEP81B_LOCUS
)

f_profile, r_profile = read_tsv(
    STEP81B_PROFILE
)

f_pw, r_pw = read_tsv(
    STEP81B_PW
)

f_dom, r_dom = read_tsv(
    STEP81B_DOM
)

f_master, r_master = read_tsv(
    MASTER
)


# ============================================================
# VALIDACIÓN DE COLUMNAS
# ============================================================

require_columns(
    f73,
    [
        "sample",
        "producer",
        "producer_name",
        "biological_unit",
        "week",
        "subject",
        "locus_id",
        "source_MAG",
        "locus_context_length",
        "locus_contig_length",
        "locus_left_edge_truncated",
        "locus_right_edge_truncated",
        "locus_n_candidate_hits",
        "locus_n_unique_candidates",
        "locus_candidates",
        "locus_contains_interrupted_candidate",
        "primary_source_linked_detection",
        "strict_source_linked_detection",
        "sensitivity_locus_signal",
        "source_MAG_main_detection",
    ],
    STEP73.name
)

require_columns(
    f_locus,
    [
        "locus_id",
        "source_MAG",
        "source_MAG_genus",
        "source_MAG_species",
        "n_samples",
        "primary_source_linked",
        "strict_source_linked",
        "sensitivity_locus_signal",
        "A_source_linked",
        "B_locus_signal_source_unresolved",
        "C_source_MAG_without_locus_support",
        "D_no_support",
        "median_source_MAG_RA_when_primary",
    ],
    STEP81B_LOCUS.name
)

require_columns(
    f_profile,
    [
        "source_MAG",
        "genus",
        "species",
        "n_primary_functional_features",
        "primary_functional_features",
        "n_independent_loci",
        "independent_loci",
    ],
    STEP81B_PROFILE.name
)

require_columns(
    f_pw,
    [
        "producer",
        "producer_name",
        "week",
        "source_MAG",
        "source_MAG_genus",
        "source_MAG_species",
        "n_independent_loci",
        "n_biological_units",
        "n_units_source_MAG_main_detected",
        "prevalence_source_MAG_main",
        "median_source_MAG_RA_main_pct",
        "total_primary_locus_detections",
        "possible_locus_detections",
        "fraction_possible_locus_detections",
        "n_units_any_locus_primary",
        "n_units_all_loci_primary",
    ],
    STEP81B_PW.name
)

require_columns(
    f_dom,
    [
        "source_MAG",
        "genus",
        "species",
        "functional_domain",
        "feature",
        "n_samples_source_MAG_detected",
        "median_source_MAG_RA_main_pct",
        "median_community_functional_RA_main_pct",
        "median_source_MAG_share_of_functional_RA_pct",
        "min_source_MAG_share_of_functional_RA_pct",
        "max_source_MAG_share_of_functional_RA_pct",
    ],
    STEP81B_DOM.name
)

require_columns(
    f_master,
    [
        "MAG",
        "genus",
        "species",
        "completeness",
        "contamination",
        "galactose_Leloir",
        "lactose_transport_betaGal",
        "lactose_PTS_LacEFG",
        "lactate_LDH",
        "citrate_fermentation",
        "acetoin_branch",
        "butanediol_branch",
        "peptide_utilization_evidence",
        "best_surface_proteinase_evidence",
        "lipolysis_evidence_curated",
        "amino_acid_aroma_evidence_curated",
        "EPS_capsule_evidence_curated",
        "acid_stress_evidence",
        "acid_stress_systems",
        "osmotic_stress_evidence_curated",
        "oxidative_stress_evidence",
        "biogenic_amine_evidence_summary",
        "histamine",
        "tyramine",
        "cadaverine",
        "putrescine",
    ],
    MASTER.name
)


# ============================================================
# QC ESTRUCTURAL
# ============================================================

if len(r73) != 108:
    die(
        f"Step73: se esperaban 108 filas; "
        f"hay {len(r73)}"
    )

loci_set = sorted({
    r["locus_id"]
    for r in r73
})

if loci_set != sorted(
    EXPECTED_LOCI
):
    die(
        "Los loci de Step73 no coinciden con "
        "los seis loci independientes esperados."
    )

source_set = sorted({
    r["source_MAG"]
    for r in r73
})

if source_set != sorted(
    EXPECTED_SOURCE_MAGS
):
    die(
        "Los MAGs fuente no coinciden con "
        "los tres esperados."
    )

if len(r_locus) != 6:
    die(
        "integrated_summary_by_locus.tsv "
        "debe contener 6 filas."
    )

if len(r_profile) != 3:
    die(
        "source_MAG_primary_functional_profile_summary.tsv "
        "debe contener 3 filas."
    )

if len(r_master) != 18:
    die(
        "MAG_functional_master_18MAGs.tsv "
        "debe contener 18 MAGs."
    )


# ============================================================
# ÍNDICES
# ============================================================

locus_summary = {
    r["locus_id"]: r
    for r in r_locus
}

profile_by_mag = {
    r["source_MAG"]: r
    for r in r_profile
}

master_by_mag = {
    r["MAG"]: r
    for r in r_master
}

pw_by_key = {
    (
        r["source_MAG"],
        r["producer"],
        fint(
            r["week"]
        )
    ): r
    for r in r_pw
}


# ============================================================
# 1. TRAYECTORIAS DIRECTAS POR LOCUS
# ============================================================

trajectory_rows = []

for locus in EXPECTED_LOCI:

    rows_locus = [
        r
        for r in r73
        if r["locus_id"] == locus
    ]

    for producer in [
        "L",
        "M"
    ]:

        rows_prod = [
            r
            for r in rows_locus
            if r["producer"] == producer
        ]

        units = sorted({
            fint(
                r["biological_unit"]
            )
            for r in rows_prod
        })

        if units != [
            1,
            2,
            3
        ]:
            die(
                f"{locus}/{producer}: "
                "unidades biológicas inesperadas."
            )

        for unit in units:

            rr = sorted(
                [
                    r
                    for r in rows_prod
                    if fint(
                        r[
                            "biological_unit"
                        ]
                    ) == unit
                ],
                key=lambda x:
                    fint(
                        x["week"]
                    )
            )

            if [
                fint(
                    x["week"]
                )
                for x in rr
            ] != [
                0,
                1,
                2
            ]:
                die(
                    f"{locus}/{producer}{unit}: "
                    "faltan semanas 0,1,2."
                )

            primary = [
                fint(
                    x[
                        "primary_source_linked_detection"
                    ]
                )
                for x in rr
            ]

            sensitivity = [
                fint(
                    x[
                        "sensitivity_locus_signal"
                    ]
                )
                for x in rr
            ]

            trajectory_rows.append({
                "locus_id":
                    locus,
                "source_MAG":
                    rr[0]["source_MAG"],
                "producer":
                    producer,
                "producer_name":
                    rr[0]["producer_name"],
                "biological_unit":
                    unit,
                "subject":
                    rr[0]["subject"],
                "primary_week0":
                    primary[0],
                "primary_week1":
                    primary[1],
                "primary_week2":
                    primary[2],
                "primary_trajectory":
                    "".join(
                        str(x)
                        for x in primary
                    ),
                "sensitivity_week0":
                    sensitivity[0],
                "sensitivity_week1":
                    sensitivity[1],
                "sensitivity_week2":
                    sensitivity[2],
                "sensitivity_trajectory":
                    "".join(
                        str(x)
                        for x in sensitivity
                    ),
            })


write_tsv(
    OUT
    / "thesis_locus_longitudinal_trajectories.tsv",
    trajectory_rows,
    [
        "locus_id",
        "source_MAG",
        "producer",
        "producer_name",
        "biological_unit",
        "subject",
        "primary_week0",
        "primary_week1",
        "primary_week2",
        "primary_trajectory",
        "sensitivity_week0",
        "sensitivity_week1",
        "sensitivity_week2",
        "sensitivity_trajectory",
    ]
)


# ============================================================
# 2. SÍNTESIS POR LOCUS
# ============================================================

locus_synthesis = []

for locus in EXPECTED_LOCI:

    rr73 = [
        r
        for r in r73
        if r["locus_id"] == locus
    ]

    s = locus_summary[
        locus
    ]

    source_mag = s[
        "source_MAG"
    ]

    left_trunc = fint(
        unique_nonempty(
            rr73,
            "locus_left_edge_truncated"
        )
    )

    right_trunc = fint(
        unique_nonempty(
            rr73,
            "locus_right_edge_truncated"
        )
    )

    interrupted = fint(
        unique_nonempty(
            rr73,
            "locus_contains_interrupted_candidate"
        )
    )

    if (
        left_trunc == 1
        and right_trunc == 1
    ):
        context_edge_status = (
            "both_edges_truncated"
        )

    elif (
        left_trunc == 1
        or right_trunc == 1
    ):
        context_edge_status = (
            "one_edge_truncated"
        )

    else:
        context_edge_status = (
            "no_edge_truncation_detected"
        )

    A = fint(
        s[
            "A_source_linked"
        ]
    )

    B = fint(
        s[
            "B_locus_signal_source_unresolved"
        ]
    )

    C = fint(
        s[
            "C_source_MAG_without_locus_support"
        ]
    )

    D = fint(
        s[
            "D_no_support"
        ]
    )

    source_detected_samples = (
        A + C
    )

    primary_locus_signal_samples = (
        A + B
    )

    if source_detected_samples > 0:
        frac_source_with_locus = (
            A
            / source_detected_samples
        )
    else:
        frac_source_with_locus = math.nan

    if primary_locus_signal_samples > 0:
        frac_locus_source_linked = (
            A
            / primary_locus_signal_samples
        )
    else:
        frac_locus_source_linked = math.nan

    # --------------------------------------------
    # Patrón objetivo locus ↔ MAG
    # --------------------------------------------

    if (
        B == 0
        and C == 0
    ):
        linkage_pattern = (
            "exclusive_source_linked_concordance"
        )

    elif (
        C == 0
        and B > 0
    ):
        linkage_pattern = (
            "all_source_MAG_detections_supported_"
            "plus_unresolved_locus_signals"
        )

    elif (
        B == 0
        and C > 0
    ):
        linkage_pattern = (
            "incomplete_locus_support_among_"
            "source_MAG_detections"
        )

    else:
        linkage_pattern = (
            "mixed_bidirectional_discordance"
        )

    # --------------------------------------------
    # Trayectorias por productor
    # --------------------------------------------

    tr = [
        r
        for r in trajectory_rows
        if r["locus_id"] == locus
    ]

    prod_trajectories = {}

    for producer in [
        "L",
        "M"
    ]:

        x = sorted(
            [
                r
                for r in tr
                if r["producer"] ==
                producer
            ],
            key=lambda z:
                int(
                    z[
                        "biological_unit"
                    ]
                )
        )

        prod_trajectories[
            producer
        ] = ";".join(
            f"{producer}"
            f"{r['biological_unit']}:"
            f"{r['primary_trajectory']}"
            for r in x
        )

    # --------------------------------------------
    # Prevalencia locus por productor/semana
    # --------------------------------------------

    prevalence = {}

    for producer in [
        "L",
        "M"
    ]:

        for week in [
            0,
            1,
            2
        ]:

            x = [
                r
                for r in rr73
                if (
                    r["producer"]
                    == producer
                    and fint(
                        r["week"]
                    ) == week
                )
            ]

            detected = sum(
                fint(
                    r[
                        "primary_source_linked_detection"
                    ]
                )
                for r in x
            )

            prevalence[
                (
                    producer,
                    week
                )
            ] = (
                detected,
                detected / 3.0
            )

    # --------------------------------------------
    # Source MAG producer/week
    # --------------------------------------------

    source_long = {}

    for producer in [
        "L",
        "M"
    ]:

        for week in [
            0,
            1,
            2
        ]:

            key = (
                source_mag,
                producer,
                week
            )

            if key not in pw_by_key:
                die(
                    f"Falta producer-week para {key}"
                )

            x = pw_by_key[key]

            source_long[
                (
                    producer,
                    week
                )
            ] = x

    locus_synthesis.append({

        "locus_id":
            locus,

        "source_MAG":
            source_mag,

        "source_MAG_genus":
            s[
                "source_MAG_genus"
            ],

        "source_MAG_species":
            s[
                "source_MAG_species"
            ],

        "candidate_names":
            unique_nonempty(
                rr73,
                "locus_candidates"
            ),

        "locus_context_length_bp":
            unique_nonempty(
                rr73,
                "locus_context_length"
            ),

        "locus_contig_length_bp":
            unique_nonempty(
                rr73,
                "locus_contig_length"
            ),

        "left_edge_truncated":
            left_trunc,

        "right_edge_truncated":
            right_trunc,

        "context_edge_status":
            context_edge_status,

        "contains_interrupted_candidate":
            interrupted,

        "n_candidate_hits":
            unique_nonempty(
                rr73,
                "locus_n_candidate_hits"
            ),

        "n_unique_candidates":
            unique_nonempty(
                rr73,
                "locus_n_unique_candidates"
            ),

        "primary_source_linked":
            fint(
                s[
                    "primary_source_linked"
                ]
            ),

        "strict_source_linked":
            fint(
                s[
                    "strict_source_linked"
                ]
            ),

        "sensitivity_locus_signal":
            fint(
                s[
                    "sensitivity_locus_signal"
                ]
            ),

        "A_source_linked":
            A,

        "B_locus_signal_source_unresolved":
            B,

        "C_source_MAG_without_locus_support":
            C,

        "D_no_support":
            D,

        "source_MAG_detected_samples":
            source_detected_samples,

        "primary_locus_signal_samples":
            primary_locus_signal_samples,

        "fraction_source_MAG_detections_with_locus_support":
            fmt(
                frac_source_with_locus
            ),

        "fraction_primary_locus_signals_source_linked":
            fmt(
                frac_locus_source_linked
            ),

        "locus_MAG_linkage_pattern":
            linkage_pattern,

        "median_source_MAG_RA_when_primary_pct":
            fmt(
                fnum(
                    s[
                        "median_source_MAG_RA_when_primary"
                    ]
                )
            ),

        "L_primary_trajectories":
            prod_trajectories[
                "L"
            ],

        "M_primary_trajectories":
            prod_trajectories[
                "M"
            ],

        "L_locus_detected_w0":
            prevalence[
                (
                    "L",
                    0
                )
            ][0],

        "L_locus_detected_w1":
            prevalence[
                (
                    "L",
                    1
                )
            ][0],

        "L_locus_detected_w2":
            prevalence[
                (
                    "L",
                    2
                )
            ][0],

        "M_locus_detected_w0":
            prevalence[
                (
                    "M",
                    0
                )
            ][0],

        "M_locus_detected_w1":
            prevalence[
                (
                    "M",
                    1
                )
            ][0],

        "M_locus_detected_w2":
            prevalence[
                (
                    "M",
                    2
                )
            ][0],

        "L_source_MAG_detected_units_w0":
            source_long[
                (
                    "L",
                    0
                )
            ][
                "n_units_source_MAG_main_detected"
            ],

        "L_source_MAG_detected_units_w1":
            source_long[
                (
                    "L",
                    1
                )
            ][
                "n_units_source_MAG_main_detected"
            ],

        "L_source_MAG_detected_units_w2":
            source_long[
                (
                    "L",
                    2
                )
            ][
                "n_units_source_MAG_main_detected"
            ],

        "M_source_MAG_detected_units_w0":
            source_long[
                (
                    "M",
                    0
                )
            ][
                "n_units_source_MAG_main_detected"
            ],

        "M_source_MAG_detected_units_w1":
            source_long[
                (
                    "M",
                    1
                )
            ][
                "n_units_source_MAG_main_detected"
            ],

        "M_source_MAG_detected_units_w2":
            source_long[
                (
                    "M",
                    2
                )
            ][
                "n_units_source_MAG_main_detected"
            ],

        "L_source_MAG_median_RA_w0_pct":
            source_long[
                (
                    "L",
                    0
                )
            ][
                "median_source_MAG_RA_main_pct"
            ],

        "L_source_MAG_median_RA_w1_pct":
            source_long[
                (
                    "L",
                    1
                )
            ][
                "median_source_MAG_RA_main_pct"
            ],

        "L_source_MAG_median_RA_w2_pct":
            source_long[
                (
                    "L",
                    2
                )
            ][
                "median_source_MAG_RA_main_pct"
            ],

        "M_source_MAG_median_RA_w0_pct":
            source_long[
                (
                    "M",
                    0
                )
            ][
                "median_source_MAG_RA_main_pct"
            ],

        "M_source_MAG_median_RA_w1_pct":
            source_long[
                (
                    "M",
                    1
                )
            ][
                "median_source_MAG_RA_main_pct"
            ],

        "M_source_MAG_median_RA_w2_pct":
            source_long[
                (
                    "M",
                    2
                )
            ][
                "median_source_MAG_RA_main_pct"
            ],
    })


LOCUS_FIELDS = list(
    locus_synthesis[0].keys()
)

write_tsv(
    OUT
    / "thesis_candidate_locus_synthesis_6.tsv",
    locus_synthesis,
    LOCUS_FIELDS
)


# ============================================================
# 3. PERFIL DE LOS 3 MAGs FUENTE
# ============================================================

source_mag_synthesis = []

for mag in EXPECTED_SOURCE_MAGS:

    p = profile_by_mag[
        mag
    ]

    m = master_by_mag[
        mag
    ]

    loci_for_mag = [
        r
        for r in locus_synthesis
        if r["source_MAG"] == mag
    ]

    A_total = sum(
        int(
            r[
                "A_source_linked"
            ]
        )
        for r in loci_for_mag
    )

    B_total = sum(
        int(
            r[
                "B_locus_signal_source_unresolved"
            ]
        )
        for r in loci_for_mag
    )

    C_total = sum(
        int(
            r[
                "C_source_MAG_without_locus_support"
            ]
        )
        for r in loci_for_mag
    )

    D_total = sum(
        int(
            r[
                "D_no_support"
            ]
        )
        for r in loci_for_mag
    )

    source_detected_counts = {
        int(
            r[
                "source_MAG_detected_samples"
            ]
        )
        for r in loci_for_mag
    }

    if len(
        source_detected_counts
    ) != 1:
        die(
            f"{mag}: número de muestras con "
            "MAG fuente detectado difiere entre loci."
        )

    source_detected_samples = (
        list(
            source_detected_counts
        )[0]
    )

    ba_types = []

    for field in [
        "histamine",
        "tyramine",
        "cadaverine",
        "putrescine",
    ]:

        state = str(
            m[field]
        ).strip()

        if state != "not_detected":
            ba_types.append(
                f"{field}:{state}"
            )

    if ba_types:
        ba_detail = ";".join(
            ba_types
        )
    else:
        ba_detail = (
            "no_specific_BA_evidence"
        )

    # Producer-week details

    longitudinal = {}

    for producer in [
        "L",
        "M"
    ]:

        for week in [
            0,
            1,
            2
        ]:

            key = (
                mag,
                producer,
                week
            )

            x = pw_by_key[
                key
            ]

            longitudinal[
                (
                    producer,
                    week
                )
            ] = x

    source_mag_synthesis.append({

        "source_MAG":
            mag,

        "genus":
            p[
                "genus"
            ],

        "species":
            p[
                "species"
            ],

        "completeness":
            m[
                "completeness"
            ],

        "contamination":
            m[
                "contamination"
            ],

        "n_independent_loci":
            p[
                "n_independent_loci"
            ],

        "independent_loci":
            p[
                "independent_loci"
            ],

        "source_MAG_detected_samples":
            source_detected_samples,

        "total_A_source_linked_across_loci":
            A_total,

        "total_B_unresolved_across_loci":
            B_total,

        "total_C_MAG_without_locus_across_loci":
            C_total,

        "total_D_no_support_across_loci":
            D_total,

        "possible_sample_locus_pairs":
            18 * int(
                p[
                    "n_independent_loci"
                ]
            ),

        "fraction_all_possible_pairs_A_source_linked":
            fmt(
                A_total
                /
                (
                    18
                    *
                    int(
                        p[
                            "n_independent_loci"
                        ]
                    )
                )
            ),

        "n_primary_functional_features":
            p[
                "n_primary_functional_features"
            ],

        "primary_functional_features":
            p[
                "primary_functional_features"
            ],

        "galactose_Leloir":
            m[
                "galactose_Leloir"
            ],

        "lactose_transport_betaGal":
            m[
                "lactose_transport_betaGal"
            ],

        "lactose_PTS_LacEFG":
            m[
                "lactose_PTS_LacEFG"
            ],

        "lactate_LDH":
            m[
                "lactate_LDH"
            ],

        "citrate_fermentation":
            m[
                "citrate_fermentation"
            ],

        "acetoin_branch":
            m[
                "acetoin_branch"
            ],

        "butanediol_branch":
            m[
                "butanediol_branch"
            ],

        "peptide_utilization_evidence":
            m[
                "peptide_utilization_evidence"
            ],

        "best_surface_proteinase_evidence":
            m[
                "best_surface_proteinase_evidence"
            ],

        "lipolysis_evidence_curated":
            m[
                "lipolysis_evidence_curated"
            ],

        "amino_acid_aroma_evidence_curated":
            m[
                "amino_acid_aroma_evidence_curated"
            ],

        "EPS_capsule_evidence_curated":
            m[
                "EPS_capsule_evidence_curated"
            ],

        "acid_stress_evidence":
            m[
                "acid_stress_evidence"
            ],

        "acid_stress_systems":
            m[
                "acid_stress_systems"
            ],

        "osmotic_stress_evidence_curated":
            m[
                "osmotic_stress_evidence_curated"
            ],

        "oxidative_stress_evidence":
            m[
                "oxidative_stress_evidence"
            ],

        "biogenic_amine_evidence_summary":
            m[
                "biogenic_amine_evidence_summary"
            ],

        "biogenic_amine_detail":
            ba_detail,

        "L_source_MAG_detected_units_w0":
            longitudinal[
                (
                    "L",
                    0
                )
            ][
                "n_units_source_MAG_main_detected"
            ],

        "L_source_MAG_detected_units_w1":
            longitudinal[
                (
                    "L",
                    1
                )
            ][
                "n_units_source_MAG_main_detected"
            ],

        "L_source_MAG_detected_units_w2":
            longitudinal[
                (
                    "L",
                    2
                )
            ][
                "n_units_source_MAG_main_detected"
            ],

        "L_source_MAG_median_RA_w0_pct":
            longitudinal[
                (
                    "L",
                    0
                )
            ][
                "median_source_MAG_RA_main_pct"
            ],

        "L_source_MAG_median_RA_w1_pct":
            longitudinal[
                (
                    "L",
                    1
                )
            ][
                "median_source_MAG_RA_main_pct"
            ],

        "L_source_MAG_median_RA_w2_pct":
            longitudinal[
                (
                    "L",
                    2
                )
            ][
                "median_source_MAG_RA_main_pct"
            ],

        "M_source_MAG_detected_units_w0":
            longitudinal[
                (
                    "M",
                    0
                )
            ][
                "n_units_source_MAG_main_detected"
            ],

        "M_source_MAG_detected_units_w1":
            longitudinal[
                (
                    "M",
                    1
                )
            ][
                "n_units_source_MAG_main_detected"
            ],

        "M_source_MAG_detected_units_w2":
            longitudinal[
                (
                    "M",
                    2
                )
            ][
                "n_units_source_MAG_main_detected"
            ],

        "M_source_MAG_median_RA_w0_pct":
            longitudinal[
                (
                    "M",
                    0
                )
            ][
                "median_source_MAG_RA_main_pct"
            ],

        "M_source_MAG_median_RA_w1_pct":
            longitudinal[
                (
                    "M",
                    1
                )
            ][
                "median_source_MAG_RA_main_pct"
            ],

        "M_source_MAG_median_RA_w2_pct":
            longitudinal[
                (
                    "M",
                    2
                )
            ][
                "median_source_MAG_RA_main_pct"
            ],
    })


SOURCE_FIELDS = list(
    source_mag_synthesis[0].keys()
)

write_tsv(
    OUT
    / "thesis_source_MAG_synthesis_3.tsv",
    source_mag_synthesis,
    SOURCE_FIELDS
)


# ============================================================
# 4. CONTRIBUCIÓN DE FUNCIONES CLAVE
# ============================================================

key_function_rows = []

for r in r_dom:

    if (
        r["source_MAG"]
        not in EXPECTED_SOURCE_MAGS
    ):
        continue

    if (
        r["feature"]
        not in KEY_TECH_FUNCTIONS
    ):
        continue

    key_function_rows.append({

        "source_MAG":
            r[
                "source_MAG"
            ],

        "genus":
            r[
                "genus"
            ],

        "species":
            r[
                "species"
            ],

        "functional_domain":
            r[
                "functional_domain"
            ],

        "feature":
            r[
                "feature"
            ],

        "feature_label":
            FUNCTION_LABELS.get(
                r[
                    "feature"
                ],
                r[
                    "feature"
                ]
            ),

        "n_samples_source_MAG_detected":
            r[
                "n_samples_source_MAG_detected"
            ],

        "median_source_MAG_RA_main_pct":
            r[
                "median_source_MAG_RA_main_pct"
            ],

        "median_community_functional_RA_main_pct":
            r[
                "median_community_functional_RA_main_pct"
            ],

        "median_source_MAG_share_of_functional_RA_pct":
            r[
                "median_source_MAG_share_of_functional_RA_pct"
            ],

        "min_source_MAG_share_of_functional_RA_pct":
            r[
                "min_source_MAG_share_of_functional_RA_pct"
            ],

        "max_source_MAG_share_of_functional_RA_pct":
            r[
                "max_source_MAG_share_of_functional_RA_pct"
            ],
    })


key_function_rows.sort(
    key=lambda x: (
        EXPECTED_SOURCE_MAGS.index(
            x[
                "source_MAG"
            ]
        ),
        KEY_TECH_FUNCTIONS.index(
            x[
                "feature"
            ]
        )
    )
)

write_tsv(
    OUT
    / "thesis_source_MAG_key_function_contributions.tsv",
    key_function_rows,
    [
        "source_MAG",
        "genus",
        "species",
        "functional_domain",
        "feature",
        "feature_label",
        "n_samples_source_MAG_detected",
        "median_source_MAG_RA_main_pct",
        "median_community_functional_RA_main_pct",
        "median_source_MAG_share_of_functional_RA_pct",
        "min_source_MAG_share_of_functional_RA_pct",
        "max_source_MAG_share_of_functional_RA_pct",
    ]
)


# ============================================================
# 5. EJES DE EVIDENCIA, SIN SCORE
# ============================================================

evidence_axes = []

for r in locus_synthesis:

    source_mag = r[
        "source_MAG"
    ]

    m = master_by_mag[
        source_mag
    ]

    evidence_axes.append({

        "locus_id":
            r[
                "locus_id"
            ],

        "source_MAG":
            source_mag,

        "species":
            r[
                "source_MAG_species"
            ],

        "axis_locus_MAG_linkage":
            r[
                "locus_MAG_linkage_pattern"
            ],

        "axis_source_detection_with_locus_fraction":
            r[
                "fraction_source_MAG_detections_with_locus_support"
            ],

        "axis_locus_signal_source_linked_fraction":
            r[
                "fraction_primary_locus_signals_source_linked"
            ],

        "axis_context_edges":
            r[
                "context_edge_status"
            ],

        "axis_interrupted_candidate":
            r[
                "contains_interrupted_candidate"
            ],

        "axis_primary_source_linked_n18":
            r[
                "primary_source_linked"
            ],

        "axis_sensitivity_signal_n18":
            r[
                "sensitivity_locus_signal"
            ],

        "axis_MAG_completeness":
            m[
                "completeness"
            ],

        "axis_MAG_contamination":
            m[
                "contamination"
            ],

        "axis_primary_functional_features":
            profile_by_mag[
                source_mag
            ][
                "n_primary_functional_features"
            ],

        "axis_CEP_like":
            m[
                "best_surface_proteinase_evidence"
            ],

        "axis_peptide_utilization":
            m[
                "peptide_utilization_evidence"
            ],

        "axis_lactose_direct":
            m[
                "lactose_transport_betaGal"
            ],

        "axis_lactose_PTS":
            m[
                "lactose_PTS_LacEFG"
            ],

        "axis_EPS_capsule":
            m[
                "EPS_capsule_evidence_curated"
            ],

        "axis_biogenic_amine_summary":
            m[
                "biogenic_amine_evidence_summary"
            ],
    })


write_tsv(
    OUT
    / "thesis_candidate_evidence_axes_no_score.tsv",
    evidence_axes,
    list(
        evidence_axes[0].keys()
    )
)


# ============================================================
# 6. MARKDOWN PARA REVISIÓN DE TESIS
# ============================================================

md = []

md.append(
    "# Síntesis de candidatos bacteriocina/RiPP y MAGs fuente"
)

md.append("")

md.append(
    "Esta tabla resume evidencia genómica y ecológica. "
    "No demuestra expresión, producción de bacteriocina, "
    "actividad antimicrobiana ni actividad metabólica."
)

md.append("")

md.append(
    "Las abundancias relativas de MAG corresponden al conjunto "
    "de MAGs recuperados y aceptados, no a toda la microbiota."
)

md.append("")

md.append(
    "## Seis loci independientes"
)

md.append("")

md.append(
    "| Locus | MAG fuente | Especie | A/B/C/D | "
    "Asociación locus–MAG | Contexto | L trayectorias | "
    "M trayectorias |"
)

md.append(
    "|---|---|---|---|---|---|---|---|"
)

for r in locus_synthesis:

    abcd = (
        f"{r['A_source_linked']}/"
        f"{r['B_locus_signal_source_unresolved']}/"
        f"{r['C_source_MAG_without_locus_support']}/"
        f"{r['D_no_support']}"
    )

    md.append(
        "| "
        + md_escape(
            r[
                "locus_id"
            ]
        )
        + " | "
        + md_escape(
            r[
                "source_MAG"
            ]
        )
        + " | "
        + md_escape(
            r[
                "source_MAG_species"
            ]
        )
        + " | "
        + abcd
        + " | "
        + md_escape(
            r[
                "locus_MAG_linkage_pattern"
            ]
        )
        + " | "
        + md_escape(
            r[
                "context_edge_status"
            ]
        )
        + " | "
        + md_escape(
            r[
                "L_primary_trajectories"
            ]
        )
        + " | "
        + md_escape(
            r[
                "M_primary_trajectories"
            ]
        )
        + " |"
    )

md.append("")

md.append(
    "## Tres MAGs fuente independientes"
)

md.append("")

md.append(
    "| MAG | Especie | Comp. | Cont. | Loci | "
    "Funciones primary | BA específicas | "
    "L RA mediana w0→w1→w2 | "
    "M RA mediana w0→w1→w2 |"
)

md.append(
    "|---|---|---:|---:|---|---:|---|---|---|"
)

for r in source_mag_synthesis:

    L_ra = (
        f"{r['L_source_MAG_median_RA_w0_pct']}→"
        f"{r['L_source_MAG_median_RA_w1_pct']}→"
        f"{r['L_source_MAG_median_RA_w2_pct']}"
    )

    M_ra = (
        f"{r['M_source_MAG_median_RA_w0_pct']}→"
        f"{r['M_source_MAG_median_RA_w1_pct']}→"
        f"{r['M_source_MAG_median_RA_w2_pct']}"
    )

    md.append(
        "| "
        + md_escape(
            r[
                "source_MAG"
            ]
        )
        + " | "
        + md_escape(
            r[
                "species"
            ]
        )
        + " | "
        + md_escape(
            r[
                "completeness"
            ]
        )
        + " | "
        + md_escape(
            r[
                "contamination"
            ]
        )
        + " | "
        + md_escape(
            r[
                "independent_loci"
            ]
        )
        + " | "
        + md_escape(
            r[
                "n_primary_functional_features"
            ]
        )
        + " | "
        + md_escape(
            r[
                "biogenic_amine_detail"
            ]
        )
        + " | "
        + md_escape(
            L_ra
        )
        + " | "
        + md_escape(
            M_ra
        )
        + " |"
    )

md.append("")

md.append(
    "## Reglas de interpretación"
)

md.append("")

md.append(
    "- `A_source_linked`: señal del locus y MAG fuente "
    "coherentemente soportados."
)

md.append(
    "- `B_locus_signal_source_unresolved`: señal del locus sin "
    "soporte suficiente para atribuirla al MAG fuente."
)

md.append(
    "- `C_source_MAG_without_locus_support`: MAG fuente detectado "
    "sin soporte principal para el locus."
)

md.append(
    "- `D_no_support`: sin soporte principal del locus ni del "
    "MAG fuente."
)

md.append("")

md.append(
    "No se aplica un ranking numérico porque los ejes de evidencia "
    "no son equivalentes ni aditivos."
)

with (
    OUT
    / "thesis_candidate_synthesis.md"
).open(
    "w",
    encoding="utf-8"
) as fh:

    fh.write(
        "\n".join(
            md
        )
        + "\n"
    )


# ============================================================
# 7. README
# ============================================================

README = """
PASO 82a - SINTESIS FINAL DE CANDIDATOS PARA TESIS
==================================================

OBJETIVO
--------
Integrar la evidencia existente de los seis loci candidatos
bacteriocina/RiPP y sus tres MAGs fuente independientes.

NO SE CALCULA NINGUN SCORE GLOBAL.

RAZON
-----
Los siguientes ejes no son equivalentes ni deben sumarse:

- asociación locus-MAG
- prevalencia
- truncamiento del contexto
- calidad del MAG
- repertorio funcional
- contribución ecológica del MAG
- evidencia de aminas biogénicas

ARCHIVOS PRINCIPALES
--------------------

thesis_candidate_locus_synthesis_6.tsv
    Tabla final por cada uno de los seis loci.

thesis_source_MAG_synthesis_3.tsv
    Tabla final por cada uno de los tres MAGs fuente.

thesis_locus_longitudinal_trajectories.tsv
    36 trayectorias:
    6 loci x 6 unidades biológicas.

thesis_source_MAG_key_function_contributions.tsv
    Contribución relativa de cada MAG fuente a funciones
    tecnológicamente relevantes.

thesis_candidate_evidence_axes_no_score.tsv
    Ejes independientes de evidencia, sin ranking numérico.

thesis_candidate_synthesis.md
    Resumen compacto para revisión y posterior escritura
    del capítulo de resultados.

INTERPRETACION
--------------
La evidencia de locus indica recuperación de un contexto genómico
candidato. No demuestra expresión, producción ni actividad.

La evidencia funcional indica potencial genómico recuperado.
No demuestra actividad metabólica.

Las abundancias MAG corresponden únicamente al conjunto de MAGs
recuperados y aceptados, no a toda la microbiota.

ATTRLOC002 / Atlantibacter queda fuera de los seis loci biológicos
independientes y no se reincorpora como candidato independiente.
"""

with (
    OUT
    / "README_step82a.txt"
).open(
    "w",
    encoding="utf-8"
) as fh:

    fh.write(
        README.strip()
        + "\n"
    )


# ============================================================
# RESUMEN FINAL
# ============================================================

print(
    "=" * 60
)

print(
    "PASO 82a COMPLETADO"
)

print(
    "=" * 60
)

print(
    f"Loci independientes integrados:      "
    f"{len(locus_synthesis)}"
)

print(
    f"MAGs fuente independientes:          "
    f"{len(source_mag_synthesis)}"
)

print(
    f"Trayectorias locus × unidad:          "
    f"{len(trajectory_rows)}"
)

print(
    f"Filas funciones clave MAG-fuente:     "
    f"{len(key_function_rows)}"
)

print(
    f"Total A_source_linked:                "
    f"{sum(int(r['A_source_linked']) for r in locus_synthesis)}"
)

print(
    f"Total B_source_unresolved:            "
    f"{sum(int(r['B_locus_signal_source_unresolved']) for r in locus_synthesis)}"
)

print(
    f"Total C_MAG_without_locus:            "
    f"{sum(int(r['C_source_MAG_without_locus_support']) for r in locus_synthesis)}"
)

print(
    f"Total D_no_support:                   "
    f"{sum(int(r['D_no_support']) for r in locus_synthesis)}"
)

print()

print(
    f"Salida: {OUT}"
)

print(
    "PASO 82a FINALIZÓ CORRECTAMENTE; "
    "ES SEGURO SALIR."
)
