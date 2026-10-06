#!/usr/bin/env python3

import csv
import math
import os
from collections import Counter

USER = os.environ["USER"]
ROOT = f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"

LOCUS_FILE = os.path.join(
    ROOT,
    "49_bacteriocin_locus_matrices",
    "bacteriocin_locus_detection_long_18x7.tsv"
)

MAG_FILE = os.path.join(
    ROOT,
    "42_MAG_detection_classification",
    "MAG_detection_classification_18x18.tsv"
)

OUTDIR = os.path.join(
    ROOT,
    "50_bacteriocin_locus_MAG_concordance"
)

os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# CONFIGURACIÓN
# ============================================================

SAMPLES = [
    "L1_1", "L1_2", "L1_3",
    "L2_1", "L2_2", "L2_3",
    "L3_1", "L3_2", "L3_3",
    "M1_1", "M1_2", "M1_3",
    "M2_1", "M2_2", "M2_3",
    "M3_1", "M3_2", "M3_3",
]

LOCI = [
    f"ATTRLOC{i:03d}"
    for i in range(1, 8)
]

SOURCE_MAGS = [
    "L2__L2_maxbin2.004_sub",
    "L2__L2_maxbin2.011_sub",
    "L3__concoct_29",
    "M2__M2_maxbin2.004",
]


# ATTRLOC002 es una referencia competitiva adicional producida
# por la coincidencia exacta del precursor de ATTRLOC006
# en otro contexto genómico.
REFERENCE_ROLE = {
    "ATTRLOC001": "candidate_source_context",
    "ATTRLOC002": "competitive_shared_context",
    "ATTRLOC003": "candidate_source_context",
    "ATTRLOC004": "candidate_source_context",
    "ATTRLOC005": "candidate_source_context",
    "ATTRLOC006": "candidate_source_context",
    "ATTRLOC007": "candidate_source_context",
}

COUNT_AS_INDEPENDENT_LOCUS = {
    locus: int(locus != "ATTRLOC002")
    for locus in LOCI
}


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


def num(value):
    if value is None:
        return None

    value = str(value).strip()

    if value == "" or value.upper() in {
        "NA", "N/A", "NAN", "NONE"
    }:
        return None

    try:
        return float(value)

    except ValueError:
        return None


def integer(value, default=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def fmt(value):
    if value is None:
        return "NA"

    return f"{value:.10g}"


# ============================================================
# LEER
# ============================================================

locus_rows = read_tsv(
    LOCUS_FILE
)

mag_rows = read_tsv(
    MAG_FILE
)


if len(locus_rows) != 126:
    raise SystemExit(
        "ERROR: tabla ATTRLOC tiene "
        f"{len(locus_rows)} filas; "
        "se esperaban 126."
    )


if len(mag_rows) != 324:
    raise SystemExit(
        "ERROR: tabla MAG tiene "
        f"{len(mag_rows)} filas; "
        "se esperaban 324."
    )


# ============================================================
# VALIDAR KEYS
# ============================================================

locus_lookup = {}

for row in locus_rows:

    key = (
        row["sample"],
        row["locus_id"]
    )

    if key in locus_lookup:
        raise SystemExit(
            "ERROR: ATTRLOC duplicado: "
            f"{key}"
        )

    locus_lookup[key] = row


if len(locus_lookup) != 126:
    raise SystemExit(
        "ERROR: no existen 126 combinaciones "
        "sample × ATTRLOC únicas."
    )


mag_lookup = {}

for row in mag_rows:

    key = (
        row["sample"],
        row["MAG"]
    )

    if key in mag_lookup:
        raise SystemExit(
            "ERROR: MAG duplicado: "
            f"{key}"
        )

    mag_lookup[key] = row


if len(mag_lookup) != 324:
    raise SystemExit(
        "ERROR: no existen 324 combinaciones "
        "sample × MAG únicas."
    )


# ============================================================
# COMPROBAR LOS 4 MAGs FUENTE
# ============================================================

for sample in SAMPLES:

    for mag in SOURCE_MAGS:

        if (
            sample,
            mag
        ) not in mag_lookup:

            raise SystemExit(
                "ERROR: falta combinación "
                f"{sample} × {mag}"
            )


# ============================================================
# CAMPOS MAG QUE CONSERVAREMOS
# ============================================================

MAG_FIELDS = [
    "reported_by_instrain",
    "coverage",
    "breadth",
    "breadth_minCov",
    "coverage_raw",
    "breadth_raw",
    "breadth_expected",
    "conANI_reference",
    "popANI_reference",
    "filtered_read_pair_count",
    "reads_unfiltered_pairs",
    "reads_mean_PID",
    "length",
    "breadth_to_expected_ratio",
    "detection_class",
    "main_detection",
    "strict_detection",
]


# ============================================================
# INTEGRAR LOCUS + MAG FUENTE
# ============================================================

joint_rows = []


for locus_row in locus_rows:

    sample = locus_row["sample"]
    locus = locus_row["locus_id"]
    source_mag = locus_row["locus_MAG"]

    key = (
        sample,
        source_mag
    )

    if key not in mag_lookup:

        raise SystemExit(
            "ERROR: no existe MAG fuente "
            f"{key}"
        )

    mag_row = mag_lookup[key]


    # --------------------------------------------------------
    # Comprobar metadata experimental
    # --------------------------------------------------------

    for field in [
        "producer",
        "biological_unit",
        "week"
    ]:

        if (
            str(locus_row[field])
            != str(mag_row[field])
        ):

            raise SystemExit(
                "ERROR: metadata no coincide "
                f"para {sample}, {locus}, "
                f"campo {field}."
            )


    locus_main = integer(
        locus_row[
            "context_evidence_main"
        ]
    )

    locus_strict = integer(
        locus_row[
            "context_evidence_strict"
        ]
    )

    mag_main = integer(
        mag_row[
            "main_detection"
        ]
    )

    mag_strict = integer(
        mag_row[
            "strict_detection"
        ]
    )


    # ========================================================
    # CONCORDANCIA PRINCIPAL
    # ========================================================

    if locus_main == 1 and mag_main == 1:

        joint_main = (
            "locus_and_MAG_supported"
        )

    elif locus_main == 1 and mag_main == 0:

        joint_main = (
            "locus_supported_MAG_not_supported"
        )

    elif locus_main == 0 and mag_main == 1:

        joint_main = (
            "MAG_supported_locus_not_supported"
        )

    else:

        joint_main = (
            "neither_supported"
        )


    # ========================================================
    # CONCORDANCIA ESTRICTA
    # ========================================================

    if (
        locus_strict == 1
        and mag_strict == 1
    ):

        joint_strict = (
            "locus_and_MAG_supported"
        )

    elif (
        locus_strict == 1
        and mag_strict == 0
    ):

        joint_strict = (
            "locus_supported_MAG_not_supported"
        )

    elif (
        locus_strict == 0
        and mag_strict == 1
    ):

        joint_strict = (
            "MAG_supported_locus_not_supported"
        )

    else:

        joint_strict = (
            "neither_supported"
        )


    # ========================================================
    # COVERAGE LOCUS / MAG
    #
    # SOLO interpretable como concordancia cuando
    # ambos pasan el criterio principal.
    #
    # NO es estimación de copy number.
    # ========================================================

    locus_cov = num(
        locus_row.get(
            "coverage"
        )
    )

    mag_cov = num(
        mag_row.get(
            "coverage"
        )
    )

    coverage_ratio = None
    log2_coverage_ratio = None
    ratio_interpretable = 0


    if (
        locus_main == 1
        and mag_main == 1
        and locus_cov is not None
        and mag_cov is not None
        and locus_cov >= 0
        and mag_cov > 0
    ):

        coverage_ratio = (
            locus_cov / mag_cov
        )

        ratio_interpretable = 1

        if coverage_ratio > 0:

            log2_coverage_ratio = (
                math.log2(
                    coverage_ratio
                )
            )


    # ========================================================
    # CONSTRUIR FILA
    # ========================================================

    rec = dict(
        locus_row
    )

    rec[
        "reference_role"
    ] = REFERENCE_ROLE[locus]

    rec[
        "count_as_independent_locus"
    ] = COUNT_AS_INDEPENDENT_LOCUS[
        locus
    ]

    rec[
        "source_MAG"
    ] = source_mag


    for field in MAG_FIELDS:

        rec[
            f"source_MAG_{field}"
        ] = mag_row.get(
            field,
            "NA"
        )


    rec[
        "joint_support_main"
    ] = joint_main

    rec[
        "joint_support_strict"
    ] = joint_strict

    rec[
        "coverage_locus_to_MAG_ratio"
    ] = fmt(
        coverage_ratio
    )

    rec[
        "log2_coverage_locus_to_MAG_ratio"
    ] = fmt(
        log2_coverage_ratio
    )

    rec[
        "coverage_ratio_interpretable"
    ] = ratio_interpretable


    joint_rows.append(
        rec
    )


# ============================================================
# VALIDAR RESULTADO
# ============================================================

if len(joint_rows) != 126:

    raise SystemExit(
        "ERROR: resultado integrado "
        f"tiene {len(joint_rows)} filas."
    )


# ============================================================
# TABLA LARGA COMPLETA
# ============================================================

new_fields = [
    "reference_role",
    "count_as_independent_locus",
    "source_MAG",

    *[
        f"source_MAG_{field}"
        for field in MAG_FIELDS
    ],

    "joint_support_main",
    "joint_support_strict",

    "coverage_locus_to_MAG_ratio",
    "log2_coverage_locus_to_MAG_ratio",
    "coverage_ratio_interpretable",
]


fields = (
    list(
        locus_rows[0].keys()
    )
    + new_fields
)


write_tsv(
    os.path.join(
        OUTDIR,
        "locus_MAG_concordance_long_18x7.tsv"
    ),
    joint_rows,
    fields
)


# ============================================================
# TABLA ÚNICA DE LOS 4 MAGs FUENTE
# 18 × 4 = 72 filas
# ============================================================

source_mag_rows = []


for sample in SAMPLES:

    for mag in SOURCE_MAGS:

        row = dict(
            mag_lookup[
                (
                    sample,
                    mag
                )
            ]
        )

        source_mag_rows.append(
            row
        )


if len(source_mag_rows) != 72:

    raise SystemExit(
        "ERROR: tabla de MAGs fuente "
        "no tiene 72 filas."
    )


write_tsv(
    os.path.join(
        OUTDIR,
        "source_MAG_status_18x4.tsv"
    ),
    source_mag_rows,
    list(
        mag_rows[0].keys()
    )
)


# ============================================================
# RESUMEN POR ATTRLOC
# ============================================================

JOINT_CLASSES = [
    "locus_and_MAG_supported",
    "locus_supported_MAG_not_supported",
    "MAG_supported_locus_not_supported",
    "neither_supported",
]


locus_summary = []


for locus in LOCI:

    rows = [
        row
        for row in joint_rows
        if row["locus_id"] == locus
    ]

    rec = {
        "locus_id":
            locus,

        "reference_role":
            REFERENCE_ROLE[
                locus
            ],

        "count_as_independent_locus":
            COUNT_AS_INDEPENDENT_LOCUS[
                locus
            ],

        "source_MAG":
            rows[0][
                "source_MAG"
            ],

        "n_samples":
            len(rows),

        "locus_main_evidence":
            sum(
                integer(
                    row[
                        "context_evidence_main"
                    ]
                )
                for row in rows
            ),

        "MAG_main_evidence":
            sum(
                integer(
                    row[
                        "source_MAG_main_detection"
                    ]
                )
                for row in rows
            ),

        "locus_strict_evidence":
            sum(
                integer(
                    row[
                        "context_evidence_strict"
                    ]
                )
                for row in rows
            ),

        "MAG_strict_evidence":
            sum(
                integer(
                    row[
                        "source_MAG_strict_detection"
                    ]
                )
                for row in rows
            ),
    }


    for state in JOINT_CLASSES:

        rec[
            f"main_{state}"
        ] = sum(
            row[
                "joint_support_main"
            ] == state
            for row in rows
        )


    for state in JOINT_CLASSES:

        rec[
            f"strict_{state}"
        ] = sum(
            row[
                "joint_support_strict"
            ] == state
            for row in rows
        )


    ratios = [
        num(
            row[
                "coverage_locus_to_MAG_ratio"
            ]
        )
        for row in rows
        if integer(
            row[
                "coverage_ratio_interpretable"
            ]
        ) == 1
    ]

    ratios = [
        x
        for x in ratios
        if x is not None
    ]


    if ratios:

        ratios_sorted = sorted(
            ratios
        )

        n = len(
            ratios_sorted
        )

        if n % 2 == 1:

            median = (
                ratios_sorted[
                    n // 2
                ]
            )

        else:

            median = (
                ratios_sorted[
                    n // 2 - 1
                ]
                +
                ratios_sorted[
                    n // 2
                ]
            ) / 2

    else:

        median = None


    rec[
        "n_ratio_interpretable"
    ] = len(
        ratios
    )

    rec[
        "median_locus_MAG_coverage_ratio"
    ] = fmt(
        median
    )


    locus_summary.append(
        rec
    )


summary_fields = [
    "locus_id",
    "reference_role",
    "count_as_independent_locus",
    "source_MAG",
    "n_samples",
    "locus_main_evidence",
    "MAG_main_evidence",
    "locus_strict_evidence",
    "MAG_strict_evidence",

    *[
        f"main_{state}"
        for state in JOINT_CLASSES
    ],

    *[
        f"strict_{state}"
        for state in JOINT_CLASSES
    ],

    "n_ratio_interpretable",
    "median_locus_MAG_coverage_ratio",
]


write_tsv(
    os.path.join(
        OUTDIR,
        "locus_MAG_concordance_summary_by_ATTRLOC.tsv"
    ),
    locus_summary,
    summary_fields
)


# ============================================================
# RESUMEN POR MUESTRA
#
# Excluye ATTRLOC002 del conteo de loci biológicos
# ============================================================

sample_summary = []


for sample in SAMPLES:

    rows = [
        row
        for row in joint_rows
        if (
            row["sample"] == sample
            and integer(
                row[
                    "count_as_independent_locus"
                ]
            ) == 1
        )
    ]


    rec = {
        "sample":
            sample,

        "producer":
            rows[0][
                "producer"
            ],

        "producer_name":
            rows[0][
                "producer_name"
            ],

        "biological_unit":
            rows[0][
                "biological_unit"
            ],

        "week":
            rows[0][
                "week"
            ],

        "n_independent_loci":
            len(rows),

        "locus_main_evidence":
            sum(
                integer(
                    row[
                        "context_evidence_main"
                    ]
                )
                for row in rows
            ),

        "locus_strict_evidence":
            sum(
                integer(
                    row[
                        "context_evidence_strict"
                    ]
                )
                for row in rows
            ),

        "locus_and_MAG_supported_main":
            sum(
                row[
                    "joint_support_main"
                ]
                ==
                "locus_and_MAG_supported"
                for row in rows
            ),

        "locus_supported_MAG_not_supported_main":
            sum(
                row[
                    "joint_support_main"
                ]
                ==
                "locus_supported_MAG_not_supported"
                for row in rows
            ),

        "MAG_supported_locus_not_supported_main":
            sum(
                row[
                    "joint_support_main"
                ]
                ==
                "MAG_supported_locus_not_supported"
                for row in rows
            ),

        "neither_supported_main":
            sum(
                row[
                    "joint_support_main"
                ]
                ==
                "neither_supported"
                for row in rows
            ),
    }


    sample_summary.append(
        rec
    )


write_tsv(
    os.path.join(
        OUTDIR,
        "locus_MAG_concordance_summary_by_sample.tsv"
    ),
    sample_summary,
    list(
        sample_summary[0].keys()
    )
)


# ============================================================
# CASOS DISCORDANTES
# ============================================================

discordant = [
    row
    for row in joint_rows
    if (
        row[
            "joint_support_main"
        ]
        in {
            "locus_supported_MAG_not_supported",
            "MAG_supported_locus_not_supported",
        }
    )
]


discordant.sort(
    key=lambda row: (
        row[
            "joint_support_main"
        ],
        row[
            "locus_id"
        ],
        row[
            "sample"
        ]
    )
)


discordant_fields = [
    "sample",
    "producer",
    "biological_unit",
    "week",

    "locus_id",
    "reference_role",
    "source_MAG",

    "detection_class",
    "context_evidence_main",
    "coverage",
    "breadth",
    "breadth_to_expected_ratio",

    "source_MAG_detection_class",
    "source_MAG_main_detection",
    "source_MAG_coverage",
    "source_MAG_breadth",
    "source_MAG_breadth_to_expected_ratio",

    "joint_support_main",

    "coverage_locus_to_MAG_ratio",
    "coverage_ratio_interpretable",
]


write_tsv(
    os.path.join(
        OUTDIR,
        "discordant_locus_MAG_support.tsv"
    ),
    discordant,
    discordant_fields
)


# ============================================================
# CONCORDANTES CON RATIO INTERPRETABLE
# ============================================================

ratio_rows = [
    row
    for row in joint_rows
    if integer(
        row[
            "coverage_ratio_interpretable"
        ]
    ) == 1
]


write_tsv(
    os.path.join(
        OUTDIR,
        "concordant_main_coverage_ratios.tsv"
    ),
    ratio_rows,
    [
        "sample",
        "producer",
        "biological_unit",
        "week",
        "locus_id",
        "reference_role",
        "source_MAG",

        "coverage",
        "breadth",
        "detection_class",

        "source_MAG_coverage",
        "source_MAG_breadth",
        "source_MAG_detection_class",

        "coverage_locus_to_MAG_ratio",
        "log2_coverage_locus_to_MAG_ratio",
    ]
)


# ============================================================
# RESUMEN GLOBAL DE CONCORDANCIA
# ============================================================

main_counts = Counter(
    row[
        "joint_support_main"
    ]
    for row in joint_rows
)

strict_counts = Counter(
    row[
        "joint_support_strict"
    ]
    for row in joint_rows
)


global_rows = []


for state in JOINT_CLASSES:

    global_rows.append({
        "criterion":
            "main",

        "joint_support":
            state,

        "count":
            main_counts[
                state
            ],

        "percent_of_126":
            f"{100*main_counts[state]/126:.3f}"
    })


for state in JOINT_CLASSES:

    global_rows.append({
        "criterion":
            "strict",

        "joint_support":
            state,

        "count":
            strict_counts[
                state
            ],

        "percent_of_126":
            f"{100*strict_counts[state]/126:.3f}"
    })


write_tsv(
    os.path.join(
        OUTDIR,
        "joint_support_global_summary.tsv"
    ),
    global_rows,
    [
        "criterion",
        "joint_support",
        "count",
        "percent_of_126"
    ]
)


# ============================================================
# README
# ============================================================

with open(
    os.path.join(
        OUTDIR,
        "README_step72.txt"
    ),
    "w"
) as f:

    f.write(
        "Paso 72 - Concordancia ATTRLOC / MAG fuente\n"
    )

    f.write(
        "=========================================\n\n"
    )

    f.write(
        "El objetivo es comprobar si la evidencia "
        "de un contexto bacteriocinico/RiPP es "
        "concordante con la deteccion del MAG que "
        "origino ese contexto.\n\n"
    )

    f.write(
        "ATTRLOC002 se conserva como referencia "
        "competitiva compartida y NO se cuenta "
        "como un locus biologico independiente.\n\n"
    )

    f.write(
        "coverage_locus_to_MAG_ratio solo se "
        "calcula cuando locus y MAG pasan el "
        "criterio principal. Es un indice "
        "diagnostico de concordancia y NO una "
        "estimacion de copy number.\n"
    )


# ============================================================
# MOSTRAR RESULTADO
# ============================================================

print(
    "============================================================"
)

print(
    "PASO 72 COMPLETADO"
)

print(
    "============================================================"
)

print(
    f"Filas locus x muestra:      {len(joint_rows)}"
)

print(
    f"MAGs fuente únicos:         {len(SOURCE_MAGS)}"
)

print(
    f"Filas MAG fuente:           {len(source_mag_rows)}"
)

print()

print(
    "CONCORDANCIA PRINCIPAL"
)

for state in JOINT_CLASSES:

    print(
        f"{state:38s}"
        f"{main_counts[state]:4d}"
    )

print()

print(
    "CONCORDANCIA ESTRICTA"
)

for state in JOINT_CLASSES:

    print(
        f"{state:38s}"
        f"{strict_counts[state]:4d}"
    )

print()

print(
    f"Salida: {OUTDIR}"
)
