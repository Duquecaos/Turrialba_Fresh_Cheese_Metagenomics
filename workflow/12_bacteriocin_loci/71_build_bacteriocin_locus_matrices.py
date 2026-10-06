#!/usr/bin/env python3

import csv
import os
import re

USER = os.environ["USER"]

ROOT = f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"

REFMETA = os.path.join(
    ROOT,
    "47_bacteriocin_locus_reference",
    "all_match_locus_metadata.tsv"
)

MAPBASE = os.path.join(
    ROOT,
    "48_bacteriocin_locus_mapping"
)

OUTDIR = os.path.join(
    ROOT,
    "49_bacteriocin_locus_matrices"
)

os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# DISEÑO EXPERIMENTAL
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

CATEGORIES = [
    "robust_context_wide",
    "coherent_limited_breadth",
    "intermediate",
    "very_low_breadth",
    "partial_inconsistent",
    "insufficient_expected",
    "not_reported",
]


# ============================================================
# CAMPOS QUE CONSERVAREMOS
# ============================================================

META_FIELDS = [
    "MAG",
    "contig",
    "context_start",
    "context_end",
    "context_length",
    "contig_length",
    "left_edge_truncated",
    "right_edge_truncated",
    "n_candidate_hits",
    "n_unique_candidates",
    "candidates",
    "contains_interrupted_candidate",
    "context_sha256",
]

METRICS = [
    "coverage",
    "breadth",
    "breadth_expected",
    "coverage_median",
    "coverage_std",
    "breadth_minCov",
    "nucl_diversity",
    "conANI_reference",
    "popANI_reference",
    "filtered_read_pair_count",
    "reads_unfiltered_pairs",
    "reads_mean_PID",
    "reads_unfiltered_reads",
    "divergent_site_count",
]


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
        "NA",
        "N/A",
        "NAN",
        "NONE"
    }:
        return None

    try:
        return float(value)

    except ValueError:
        return None


def fmt(value):

    if value is None:
        return "NA"

    return f"{value:.10g}"


def sample_meta(sample):

    match = re.fullmatch(
        r"([LM])([1-3])_([1-3])",
        sample
    )

    if not match:

        raise SystemExit(
            f"ERROR: sample ID inválido: {sample}"
        )

    producer = match.group(1)

    biological_unit = int(
        match.group(2)
    )

    week = (
        int(match.group(3))
        - 1
    )

    producer_name = (
        "Lidieth"
        if producer == "L"
        else "Marino"
    )

    return {
        "producer": producer,
        "producer_name": producer_name,
        "biological_unit": biological_unit,
        "week": week,
        "subject": (
            f"{producer}"
            f"{biological_unit}"
        )
    }


# ============================================================
# METADATA DE LOS 7 ATTRLOC
# ============================================================

meta_rows = read_tsv(
    REFMETA
)

if len(meta_rows) != 7:

    raise SystemExit(
        "ERROR: se esperaban 7 ATTRLOC "
        f"en metadata; hay {len(meta_rows)}."
    )

meta = {
    row["locus_id"]: row
    for row in meta_rows
}

if sorted(meta) != LOCI:

    raise SystemExit(
        "ERROR: los IDs de metadata no son "
        "ATTRLOC001-ATTRLOC007."
    )


# ============================================================
# CONSTRUIR TABLA COMPLETA
# 18 muestras × 7 ATTRLOC = 126 combinaciones
# ============================================================

long_rows = []

reported_total = 0


for sample in SAMPLES:

    marker = os.path.join(
        MAPBASE,
        sample,
        "LOCUS_MAPPING_COMPLETED"
    )

    info = os.path.join(
        MAPBASE,
        sample,
        f"{sample}.locus_genome_info.tsv"
    )

    if not os.path.isfile(marker):

        raise SystemExit(
            "ERROR: falta "
            "LOCUS_MAPPING_COMPLETED "
            f"para {sample}."
        )

    if (
        not os.path.isfile(info)
        or os.path.getsize(info) == 0
    ):

        raise SystemExit(
            "ERROR: falta genome_info "
            f"para {sample}."
        )

    observed_rows = read_tsv(
        info
    )

    observed = {}

    for row in observed_rows:

        locus = row["genome"]

        if locus not in LOCI:

            raise SystemExit(
                f"ERROR: locus inesperado "
                f"en {sample}: {locus}"
            )

        if locus in observed:

            raise SystemExit(
                f"ERROR: locus duplicado "
                f"en {sample}: {locus}"
            )

        observed[locus] = row

    reported_total += len(
        observed
    )

    smeta = sample_meta(
        sample
    )


    for locus in LOCI:

        reported = int(
            locus in observed
        )

        raw = observed.get(
            locus,
            {}
        )

        coverage = (
            num(
                raw.get("coverage")
            )
            if reported
            else 0.0
        )

        breadth = (
            num(
                raw.get("breadth")
            )
            if reported
            else 0.0
        )

        expected = (
            num(
                raw.get(
                    "breadth_expected"
                )
            )
            if reported
            else None
        )

        if coverage is None:
            coverage = 0.0

        if breadth is None:
            breadth = 0.0

        ratio = (
            breadth / expected
            if (
                expected is not None
                and expected > 0
            )
            else None
        )


        # ====================================================
        # CLASIFICACIÓN OPERACIONAL
        #
        # IMPORTANTE:
        # describe coherencia del reclutamiento del contexto.
        # NO demuestra producción de bacteriocina.
        # ====================================================

        if (
            reported == 0
            or breadth == 0
        ):

            category = (
                "not_reported"
            )

        elif (
            ratio is not None
            and ratio < 0.50
        ):

            category = (
                "partial_inconsistent"
            )

        elif breadth < 0.10:

            category = (
                "very_low_breadth"
            )

        elif ratio is None:

            category = (
                "insufficient_expected"
            )

        elif ratio < 0.80:

            category = (
                "intermediate"
            )

        elif breadth < 0.50:

            category = (
                "coherent_limited_breadth"
            )

        else:

            category = (
                "robust_context_wide"
            )


        main_evidence = int(
            category in {
                "coherent_limited_breadth",
                "robust_context_wide"
            }
        )

        strict_evidence = int(
            category
            == "robust_context_wide"
        )


        rec = {
            "sample": sample,
            **smeta,

            "locus_id": locus,

            "reported_by_instrain":
                reported,

            "breadth_to_expected_ratio":
                fmt(ratio),

            "detection_class":
                category,

            "context_evidence_main":
                main_evidence,

            "context_evidence_strict":
                strict_evidence
        }


        # ====================================================
        # METADATA DEL CONTEXTO ATTRLOC
        # ====================================================

        for field in META_FIELDS:

            rec[
                f"locus_{field}"
            ] = (
                meta[locus].get(
                    field,
                    "NA"
                )
                or "NA"
            )


        # ====================================================
        # MÉTRICAS inStrain
        # ====================================================

        for field in METRICS:

            if reported:

                value = raw.get(
                    field,
                    ""
                )

                rec[field] = (
                    value
                    if str(value).strip()
                    else "NA"
                )

            else:

                rec[field] = "NA"


        # ====================================================
        # AUSENCIAS EXPLÍCITAS
        #
        # Para matrices:
        # coverage y breadth = 0
        #
        # Otros parámetros permanecen NA.
        # ====================================================

        if not reported:

            rec["coverage"] = "0"
            rec["breadth"] = "0"

            rec[
                "filtered_read_pair_count"
            ] = "0"

            rec[
                "reads_unfiltered_pairs"
            ] = "0"

            rec[
                "reads_unfiltered_reads"
            ] = "0"


        long_rows.append(
            rec
        )


# ============================================================
# QC GLOBAL
# ============================================================

if len(long_rows) != 126:

    raise SystemExit(
        "ERROR: se esperaban 126 "
        f"combinaciones; hay {len(long_rows)}."
    )


# El preflight inmediatamente anterior dio 96.
if reported_total != 96:

    raise SystemExit(
        "ERROR: el preflight mostró "
        "96 combinaciones reportadas, "
        f"pero ahora hay {reported_total}."
    )


lookup = {
    (
        row["sample"],
        row["locus_id"]
    ): row

    for row in long_rows
}


if len(lookup) != 126:

    raise SystemExit(
        "ERROR: existen combinaciones "
        "sample × ATTRLOC duplicadas."
    )


# ============================================================
# TABLA LARGA
# ============================================================

BASE_FIELDS = [
    "sample",
    "producer",
    "producer_name",
    "biological_unit",
    "week",
    "subject",
    "locus_id",
    "reported_by_instrain",
    "breadth_to_expected_ratio",
    "detection_class",
    "context_evidence_main",
    "context_evidence_strict",
]


LONG_FIELDS = (
    BASE_FIELDS
    + [
        f"locus_{field}"
        for field in META_FIELDS
    ]
    + METRICS
)


write_tsv(
    os.path.join(
        OUTDIR,
        "bacteriocin_locus_detection_long_18x7.tsv"
    ),
    long_rows,
    LONG_FIELDS
)


# ============================================================
# MATRICES
# ============================================================

def write_matrix(
    field,
    filename
):

    with open(
        os.path.join(
            OUTDIR,
            filename
        ),
        "w",
        newline=""
    ) as f:

        writer = csv.writer(
            f,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow(
            ["sample"] + LOCI
        )

        for sample in SAMPLES:

            writer.writerow(
                [sample]
                + [
                    lookup[
                        (
                            sample,
                            locus
                        )
                    ][field]

                    for locus in LOCI
                ]
            )


write_matrix(
    "coverage",
    "coverage_matrix_18x7.tsv"
)

write_matrix(
    "breadth",
    "breadth_matrix_18x7.tsv"
)

write_matrix(
    "breadth_to_expected_ratio",
    "breadth_consistency_matrix_18x7.tsv"
)

write_matrix(
    "detection_class",
    "detection_class_matrix_18x7.tsv"
)

write_matrix(
    "context_evidence_main",
    "context_evidence_main_matrix_18x7.tsv"
)

write_matrix(
    "context_evidence_strict",
    "context_evidence_strict_matrix_18x7.tsv"
)


# ============================================================
# RESUMEN GLOBAL
# ============================================================

global_summary = []


for category in CATEGORIES:

    n = sum(
        row["detection_class"]
        == category

        for row in long_rows
    )

    global_summary.append({
        "detection_class":
            category,

        "count":
            n,

        "percent_of_126":
            f"{100*n/126:.3f}"
    })


write_tsv(
    os.path.join(
        OUTDIR,
        "detection_class_summary.tsv"
    ),
    global_summary,
    [
        "detection_class",
        "count",
        "percent_of_126"
    ]
)


# ============================================================
# RESUMEN POR MUESTRA
# ============================================================

sample_summary = []


for sample in SAMPLES:

    rows = [
        lookup[
            (
                sample,
                locus
            )
        ]

        for locus in LOCI
    ]

    rec = {
        "sample":
            sample,

        **sample_meta(
            sample
        ),

        "reported_by_instrain":
            sum(
                int(
                    row[
                        "reported_by_instrain"
                    ]
                )
                for row in rows
            )
    }


    for category in CATEGORIES:

        rec[category] = sum(
            row[
                "detection_class"
            ]
            == category

            for row in rows
        )


    rec[
        "context_evidence_main"
    ] = sum(
        int(
            row[
                "context_evidence_main"
            ]
        )
        for row in rows
    )


    rec[
        "context_evidence_strict"
    ] = sum(
        int(
            row[
                "context_evidence_strict"
            ]
        )
        for row in rows
    )


    sample_summary.append(
        rec
    )


write_tsv(
    os.path.join(
        OUTDIR,
        "locus_summary_by_sample.tsv"
    ),
    sample_summary,
    [
        "sample",
        "producer",
        "producer_name",
        "biological_unit",
        "week",
        "subject",
        "reported_by_instrain",
        *CATEGORIES,
        "context_evidence_main",
        "context_evidence_strict"
    ]
)


# ============================================================
# RESUMEN POR ATTRLOC
# ============================================================

locus_summary = []


for locus in LOCI:

    rows = [
        lookup[
            (
                sample,
                locus
            )
        ]

        for sample in SAMPLES
    ]

    m = meta[locus]

    rec = {
        "locus_id":
            locus,

        "MAG":
            m.get(
                "MAG",
                "NA"
            ),

        "contig":
            m.get(
                "contig",
                "NA"
            ),

        "context_length":
            m.get(
                "context_length",
                "NA"
            ),

        "candidates":
            m.get(
                "candidates",
                "NA"
            ),

        "contains_interrupted_candidate":
            m.get(
                "contains_interrupted_candidate",
                "NA"
            ),

        "left_edge_truncated":
            m.get(
                "left_edge_truncated",
                "NA"
            ),

        "right_edge_truncated":
            m.get(
                "right_edge_truncated",
                "NA"
            ),

        "reported_by_instrain":
            sum(
                int(
                    row[
                        "reported_by_instrain"
                    ]
                )
                for row in rows
            )
    }


    for category in CATEGORIES:

        rec[category] = sum(
            row[
                "detection_class"
            ]
            == category

            for row in rows
        )


    rec[
        "context_evidence_main"
    ] = sum(
        int(
            row[
                "context_evidence_main"
            ]
        )
        for row in rows
    )


    rec[
        "context_evidence_strict"
    ] = sum(
        int(
            row[
                "context_evidence_strict"
            ]
        )
        for row in rows
    )


    locus_summary.append(
        rec
    )


write_tsv(
    os.path.join(
        OUTDIR,
        "locus_summary_by_ATTRLOC.tsv"
    ),
    locus_summary,
    [
        "locus_id",
        "MAG",
        "contig",
        "context_length",
        "candidates",
        "contains_interrupted_candidate",
        "left_edge_truncated",
        "right_edge_truncated",
        "reported_by_instrain",
        *CATEGORIES,
        "context_evidence_main",
        "context_evidence_strict"
    ]
)


# ============================================================
# COBERTURA ALTA PERO CONTEXTO INCONSISTENTE
# ============================================================

discordant = []


for row in long_rows:

    coverage = (
        num(
            row["coverage"]
        )
        or 0.0
    )

    if (
        coverage >= 1.0
        and row[
            "detection_class"
        ] in {
            "partial_inconsistent",
            "intermediate"
        }
    ):

        discordant.append(
            row
        )


discordant.sort(
    key=lambda row:
        -(
            num(
                row["coverage"]
            )
            or 0.0
        )
)


write_tsv(
    os.path.join(
        OUTDIR,
        "discordant_high_coverage_contexts.tsv"
    ),
    discordant,
    [
        "sample",
        "producer",
        "biological_unit",
        "week",
        "locus_id",
        "locus_MAG",
        "locus_candidates",
        "coverage",
        "breadth",
        "breadth_expected",
        "breadth_to_expected_ratio",
        "filtered_read_pair_count",
        "reads_mean_PID",
        "detection_class"
    ]
)


# ============================================================
# README
# ============================================================

with open(
    os.path.join(
        OUTDIR,
        "README_step71.txt"
    ),
    "w"
) as f:

    f.write(
        "Paso 71 - Reclutamiento competitivo "
        "de loci bacteriocinicos/RiPP\n"
    )

    f.write(
        "===================================="
        "=========================\n\n"
    )

    f.write(
        "Muestras: 18\n"
    )

    f.write(
        "ATTRLOC competitivos: 7\n"
    )

    f.write(
        "Combinaciones: 126\n"
    )

    f.write(
        "Reportadas por inStrain: "
        f"{reported_total}\n"
    )

    f.write(
        "No reportadas: "
        f"{126-reported_total}\n\n"
    )

    f.write(
        "Las categorias son criterios "
        "operacionales de coherencia de "
        "reclutamiento del contexto genomico. "
        "No demuestran expresion, produccion "
        "de bacteriocina ni funcionalidad "
        "del BGC.\n"
    )


# ============================================================
# RESULTADO
# ============================================================

print(
    "============================================================"
)

print(
    "PASO 71 COMPLETADO"
)

print(
    "============================================================"
)

print(
    f"Muestras:                {len(SAMPLES)}"
)

print(
    f"ATTRLOC:                 {len(LOCI)}"
)

print(
    f"Combinaciones:           {len(long_rows)}"
)

print(
    f"Reportadas por inStrain: {reported_total}"
)

print(
    f"No reportadas:           {len(long_rows)-reported_total}"
)

print()


for row in global_summary:

    print(
        f"{row['detection_class']:28s}"
        f"{int(row['count']):4d}  "
        f"({row['percent_of_126']}%)"
    )


print()

print(
    f"Salida: {OUTDIR}"
)
