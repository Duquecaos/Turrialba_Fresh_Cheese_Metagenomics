#!/usr/bin/env python3

import csv
import os
import re
from collections import Counter, defaultdict

USER = os.environ["USER"]
ROOT = f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"

INFILE = os.path.join(
    ROOT,
    "50_bacteriocin_locus_MAG_concordance",
    "locus_MAG_concordance_long_18x7.tsv"
)

MAPBASE = os.path.join(
    ROOT,
    "48_bacteriocin_locus_mapping"
)

OUTDIR = os.path.join(
    ROOT,
    "51_bacteriocin_final_evidence"
)

os.makedirs(OUTDIR, exist_ok=True)

SAMPLES = [
    "L1_1", "L1_2", "L1_3",
    "L2_1", "L2_2", "L2_3",
    "L3_1", "L3_2", "L3_3",
    "M1_1", "M1_2", "M1_3",
    "M2_1", "M2_2", "M2_3",
    "M3_1", "M3_2", "M3_3",
]

LOCI = [
    "ATTRLOC001",
    "ATTRLOC003",
    "ATTRLOC004",
    "ATTRLOC005",
    "ATTRLOC006",
    "ATTRLOC007",
]


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


def num(x):
    if x is None:
        return None
    x = str(x).strip()
    if x == "" or x.upper() in {"NA", "N/A", "NAN", "NONE"}:
        return None
    try:
        return float(x)
    except ValueError:
        return None


def integer(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return 0


def fmt(x):
    return "NA" if x is None else f"{x:.10g}"


# ============================================================
# LEER PASO 72
# ============================================================

all_rows = read_tsv(INFILE)

if len(all_rows) != 126:
    raise SystemExit(
        f"ERROR: se esperaban 126 filas; hay {len(all_rows)}"
    )

# Excluir ATTRLOC002 de análisis biológico independiente
rows = [
    r for r in all_rows
    if integer(r["count_as_independent_locus"]) == 1
]

if len(rows) != 108:
    raise SystemExit(
        f"ERROR: se esperaban 108 combinaciones independientes; "
        f"hay {len(rows)}"
    )


# ============================================================
# PROFUNDIDAD DE SECUENCIACIÓN
# Bowtie2: primera línea "N reads; of these:"
#
# En nuestros datos todos fueron paired.
# Se conserva el nombre input_paired_fragments para no
# confundirlo con número de mates individuales.
# ============================================================

input_pairs = {}

for sample in SAMPLES:

    log = os.path.join(
        MAPBASE,
        sample,
        f"{sample}.bowtie2.log"
    )

    if not os.path.isfile(log):
        raise SystemExit(
            f"ERROR: falta log Bowtie2 para {sample}: {log}"
        )

    n = None
    paired_confirmed = False

    with open(log) as f:
        for line in f:

            m = re.match(
                r"^\s*([0-9]+)\s+reads; of these:",
                line
            )

            if m and n is None:
                n = int(m.group(1))

            if re.search(
                r"\(100\.00%\)\s+were paired",
                line
            ):
                paired_confirmed = True

    if n is None:
        raise SystemExit(
            f"ERROR: no pude leer profundidad para {sample}"
        )

    if not paired_confirmed:
        raise SystemExit(
            f"ERROR: Bowtie2 no confirma 100% paired en {sample}"
        )

    input_pairs[sample] = n


# ============================================================
# CLASIFICACIÓN FINAL
# ============================================================

final_rows = []

for r in rows:

    sample = r["sample"]

    locus_main = integer(
        r["context_evidence_main"]
    )

    locus_strict = integer(
        r["context_evidence_strict"]
    )

    mag_main = integer(
        r["source_MAG_main_detection"]
    )

    mag_strict = integer(
        r["source_MAG_strict_detection"]
    )

    joint_main = r["joint_support_main"]
    joint_strict = r["joint_support_strict"]

    # --------------------------------------------------------
    # TIER PRINCIPAL
    # --------------------------------------------------------

    if joint_main == "locus_and_MAG_supported":
        tier_main = "A_source_linked"

    elif joint_main == "locus_supported_MAG_not_supported":
        tier_main = "B_locus_signal_source_unresolved"

    elif joint_main == "MAG_supported_locus_not_supported":
        tier_main = "C_source_MAG_without_locus_support"

    else:
        tier_main = "D_no_support"

    # --------------------------------------------------------
    # TIER ESTRICTO
    # --------------------------------------------------------

    if joint_strict == "locus_and_MAG_supported":
        tier_strict = "A_source_linked"

    elif joint_strict == "locus_supported_MAG_not_supported":
        tier_strict = "B_locus_signal_source_unresolved"

    elif joint_strict == "MAG_supported_locus_not_supported":
        tier_strict = "C_source_MAG_without_locus_support"

    else:
        tier_strict = "D_no_support"

    # --------------------------------------------------------
    # TRES DEFINICIONES
    # --------------------------------------------------------

    primary_source_linked = int(
        joint_main == "locus_and_MAG_supported"
    )

    strict_source_linked = int(
        joint_strict == "locus_and_MAG_supported"
    )

    sensitivity_locus_signal = int(
        locus_main == 1
    )

    # --------------------------------------------------------
    # NORMALIZACIÓN POR PROFUNDIDAD DE MUESTRA
    #
    # coverage ya está normalizado por longitud del contexto.
    # Dividimos por millones de fragmentos paired de entrada.
    #
    # Esto NO es abundancia relativa de toda la microbiota.
    # --------------------------------------------------------

    cov = num(r["coverage"])

    filtered_pairs = num(
        r["filtered_read_pair_count"]
    )

    n_input = input_pairs[sample]

    coverage_per_million = (
        cov * 1_000_000 / n_input
        if cov is not None
        else None
    )

    filtered_pairs_cpm = (
        filtered_pairs * 1_000_000 / n_input
        if filtered_pairs is not None
        else None
    )

    primary_cov_pm = (
        coverage_per_million
        if primary_source_linked
        else None
    )

    primary_pairs_cpm = (
        filtered_pairs_cpm
        if primary_source_linked
        else None
    )

    rec = dict(r)

    rec["input_paired_fragments"] = n_input

    rec["evidence_tier_main"] = tier_main
    rec["evidence_tier_strict"] = tier_strict

    rec["primary_source_linked_detection"] = (
        primary_source_linked
    )

    rec["strict_source_linked_detection"] = (
        strict_source_linked
    )

    rec["sensitivity_locus_signal"] = (
        sensitivity_locus_signal
    )

    rec["coverage_per_million_input_pairs"] = (
        fmt(coverage_per_million)
    )

    rec["filtered_pairs_CPM"] = (
        fmt(filtered_pairs_cpm)
    )

    rec["primary_coverage_per_million"] = (
        fmt(primary_cov_pm)
    )

    rec["primary_filtered_pairs_CPM"] = (
        fmt(primary_pairs_cpm)
    )

    final_rows.append(rec)


# ============================================================
# VALIDACIONES
# ============================================================

lookup = {
    (r["sample"], r["locus_id"]): r
    for r in final_rows
}

if len(lookup) != 108:
    raise SystemExit(
        "ERROR: combinaciones sample × locus duplicadas."
    )


# ============================================================
# TABLA LARGA FINAL
# ============================================================

NEW_FIELDS = [
    "input_paired_fragments",
    "evidence_tier_main",
    "evidence_tier_strict",
    "primary_source_linked_detection",
    "strict_source_linked_detection",
    "sensitivity_locus_signal",
    "coverage_per_million_input_pairs",
    "filtered_pairs_CPM",
    "primary_coverage_per_million",
    "primary_filtered_pairs_CPM",
]

FIELDS = list(rows[0].keys()) + NEW_FIELDS

write_tsv(
    os.path.join(
        OUTDIR,
        "final_independent_locus_evidence_18x6.tsv"
    ),
    final_rows,
    FIELDS
)


# ============================================================
# MATRICES
# ============================================================

def matrix(field, filename):

    with open(
        os.path.join(OUTDIR, filename),
        "w",
        newline=""
    ) as f:

        w = csv.writer(
            f,
            delimiter="\t",
            lineterminator="\n"
        )

        w.writerow(
            ["sample"] + LOCI
        )

        for sample in SAMPLES:

            w.writerow(
                [sample] +
                [
                    lookup[(sample, locus)][field]
                    for locus in LOCI
                ]
            )


matrix(
    "evidence_tier_main",
    "evidence_tier_main_matrix_18x6.tsv"
)

matrix(
    "primary_source_linked_detection",
    "primary_source_linked_detection_matrix_18x6.tsv"
)

matrix(
    "strict_source_linked_detection",
    "strict_source_linked_detection_matrix_18x6.tsv"
)

matrix(
    "sensitivity_locus_signal",
    "sensitivity_locus_signal_matrix_18x6.tsv"
)

matrix(
    "primary_coverage_per_million",
    "primary_normalized_coverage_matrix_18x6.tsv"
)

matrix(
    "primary_filtered_pairs_CPM",
    "primary_filtered_pairs_CPM_matrix_18x6.tsv"
)


# ============================================================
# RESUMEN POR LOCUS
# ============================================================

TIER_ORDER = [
    "A_source_linked",
    "B_locus_signal_source_unresolved",
    "C_source_MAG_without_locus_support",
    "D_no_support",
]

locus_summary = []

for locus in LOCI:

    rr = [
        r for r in final_rows
        if r["locus_id"] == locus
    ]

    counts = Counter(
        r["evidence_tier_main"]
        for r in rr
    )

    strict_counts = Counter(
        r["evidence_tier_strict"]
        for r in rr
    )

    rec = {
        "locus_id": locus,
        "source_MAG": rr[0]["source_MAG"],
        "candidate": rr[0]["locus_candidates"],
        "n_samples": 18,
        "primary_source_linked":
            sum(integer(r["primary_source_linked_detection"])
                for r in rr),
        "strict_source_linked":
            sum(integer(r["strict_source_linked_detection"])
                for r in rr),
        "sensitivity_locus_signal":
            sum(integer(r["sensitivity_locus_signal"])
                for r in rr),
    }

    for tier in TIER_ORDER:
        rec[f"main_{tier}"] = counts[tier]

    for tier in TIER_ORDER:
        rec[f"strict_{tier}"] = strict_counts[tier]

    locus_summary.append(rec)


write_tsv(
    os.path.join(
        OUTDIR,
        "final_summary_by_locus.tsv"
    ),
    locus_summary,
    list(locus_summary[0].keys())
)


# ============================================================
# RESUMEN POR MUESTRA
# ============================================================

sample_summary = []

for sample in SAMPLES:

    rr = [
        r for r in final_rows
        if r["sample"] == sample
    ]

    rec = {
        "sample": sample,
        "producer": rr[0]["producer"],
        "producer_name": rr[0]["producer_name"],
        "biological_unit": rr[0]["biological_unit"],
        "week": rr[0]["week"],
        "primary_source_linked_loci":
            sum(integer(r["primary_source_linked_detection"])
                for r in rr),
        "strict_source_linked_loci":
            sum(integer(r["strict_source_linked_detection"])
                for r in rr),
        "sensitivity_locus_signals":
            sum(integer(r["sensitivity_locus_signal"])
                for r in rr),
        "unresolved_locus_signals":
            sum(
                r["evidence_tier_main"]
                == "B_locus_signal_source_unresolved"
                for r in rr
            ),
    }

    sample_summary.append(rec)


write_tsv(
    os.path.join(
        OUTDIR,
        "final_summary_by_sample.tsv"
    ),
    sample_summary,
    list(sample_summary[0].keys())
)


# ============================================================
# PREVALENCIA DESCRIPTIVA
# productor × semana × locus
#
# Cada grupo contiene 3 unidades biológicas.
# ============================================================

prev_rows = []

for producer in ["L", "M"]:

    for week in [0, 1, 2]:

        for locus in LOCI:

            rr = [
                r for r in final_rows
                if (
                    r["producer"] == producer
                    and integer(r["week"]) == week
                    and r["locus_id"] == locus
                )
            ]

            if len(rr) != 3:
                raise SystemExit(
                    f"ERROR: {producer}, semana {week}, "
                    f"{locus}: n={len(rr)}, esperado 3."
                )

            primary = sum(
                integer(
                    r["primary_source_linked_detection"]
                )
                for r in rr
            )

            strict = sum(
                integer(
                    r["strict_source_linked_detection"]
                )
                for r in rr
            )

            sensitivity = sum(
                integer(
                    r["sensitivity_locus_signal"]
                )
                for r in rr
            )

            prev_rows.append({
                "producer": producer,
                "producer_name": rr[0]["producer_name"],
                "week": week,
                "locus_id": locus,
                "n_biological_units": 3,
                "primary_detected": primary,
                "primary_prevalence": f"{primary/3:.6f}",
                "strict_detected": strict,
                "strict_prevalence": f"{strict/3:.6f}",
                "sensitivity_detected": sensitivity,
                "sensitivity_prevalence": f"{sensitivity/3:.6f}",
            })


write_tsv(
    os.path.join(
        OUTDIR,
        "prevalence_by_producer_week_locus.tsv"
    ),
    prev_rows,
    list(prev_rows[0].keys())
)


# ============================================================
# TIER GLOBAL
# ============================================================

tier_counts = Counter(
    r["evidence_tier_main"]
    for r in final_rows
)

tier_summary = []

for tier in TIER_ORDER:

    n = tier_counts[tier]

    tier_summary.append({
        "evidence_tier_main": tier,
        "count": n,
        "percent_of_108": f"{100*n/108:.3f}",
    })


write_tsv(
    os.path.join(
        OUTDIR,
        "final_tier_summary.tsv"
    ),
    tier_summary,
    [
        "evidence_tier_main",
        "count",
        "percent_of_108"
    ]
)


# ============================================================
# COMPARAR PRINCIPAL VS ESTRICTO
# ============================================================

disagreement = [
    r for r in final_rows
    if integer(r["primary_source_linked_detection"])
       != integer(r["strict_source_linked_detection"])
]

write_tsv(
    os.path.join(
        OUTDIR,
        "primary_vs_strict_disagreement.tsv"
    ),
    disagreement,
    [
        "sample",
        "locus_id",
        "source_MAG",
        "detection_class",
        "source_MAG_detection_class",
        "primary_source_linked_detection",
        "strict_source_linked_detection",
    ]
)


# ============================================================
# README
# ============================================================

with open(
    os.path.join(
        OUTDIR,
        "README_step73.txt"
    ),
    "w"
) as f:

    f.write(
        "Final evidence framework for six independent "
        "bacteriocin/RiPP candidate loci.\n\n"
    )

    f.write(
        "ATTRLOC002 is excluded from independent-locus "
        "analyses because it is a competitive shared context.\n\n"
    )

    f.write(
        "Primary evidence = locus context and source MAG "
        "both pass the main recruitment criterion.\n"
    )

    f.write(
        "Strict evidence = locus context and source MAG "
        "both pass strict criteria.\n"
    )

    f.write(
        "Sensitivity evidence = locus passes main criterion "
        "even when source MAG attribution is unresolved.\n\n"
    )

    f.write(
        "Coverage is normalized per million paired input "
        "fragments. This is NOT relative abundance of the "
        "whole microbiome and NOT locus copy number.\n"
    )


# ============================================================
# RESULTADO
# ============================================================

primary_n = sum(
    integer(r["primary_source_linked_detection"])
    for r in final_rows
)

strict_n = sum(
    integer(r["strict_source_linked_detection"])
    for r in final_rows
)

sensitivity_n = sum(
    integer(r["sensitivity_locus_signal"])
    for r in final_rows
)

print("============================================================")
print("PASO 73 COMPLETADO")
print("============================================================")
print(f"Combinaciones independientes: {len(final_rows)}")
print(f"Evidencia primaria:           {primary_n}")
print(f"Evidencia estricta:           {strict_n}")
print(f"Señales sensibilidad:         {sensitivity_n}")
print(f"Discordancias main/strict:    {len(disagreement)}")
print()
for r in tier_summary:
    print(
        f'{r["evidence_tier_main"]:38s}'
        f'{int(r["count"]):4d} '
        f'({r["percent_of_108"]}%)'
    )
print()
print(f"Salida: {OUTDIR}")
