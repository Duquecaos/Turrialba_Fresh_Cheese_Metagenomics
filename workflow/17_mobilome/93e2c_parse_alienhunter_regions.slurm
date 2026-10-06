#!/bin/bash

#SBATCH --job-name=93e2c_AH_parse
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=93e2c_AH_parse_%j.out
#SBATCH --error=93e2c_AH_parse_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

BASE="${ROOT}/95_mobilome/93E2_alienhunter_final18"
PRED="${BASE}/predictions"

AUDIT="${ROOT}/95_mobilome/93E2C0_alienhunter_audit"

OUT="${ROOT}/95_mobilome/93E2C_alienhunter_regions"

rm -rf "$OUT"
mkdir -p "$OUT"

echo "======================================================================"
echo "93E2C - PARSEO DEFINITIVO DE ALIENHUNTER"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "======================================================================"

python3 - \
    "$BASE" \
    "$PRED" \
    "$AUDIT" \
    "$OUT" <<'PY'

import csv
import math
import re
import statistics
import sys

from collections import Counter, defaultdict
from pathlib import Path


BASE = Path(sys.argv[1])
PRED = Path(sys.argv[2])
AUDIT = Path(sys.argv[3])
OUT = Path(sys.argv[4])

MANIFEST = BASE / "93E2A_manifest.tsv"


# ======================================================================
# Funciones
# ======================================================================

def clean_value(value):
    value = value.strip()

    if (
        len(value) >= 2
        and value[0] == '"'
        and value[-1] == '"'
    ):
        value = value[1:-1]

    return value


def parse_float(value):
    value = clean_value(value)

    try:
        return float(value)
    except Exception:
        return None


def merge_intervals(intervals):
    """
    intervals: [(start, end), ...], coordenadas 1-based inclusivas
    """

    if not intervals:
        return []

    intervals = sorted(intervals)

    merged = [
        list(intervals[0])
    ]

    for start, end in intervals[1:]:

        last = merged[-1]

        if start <= last[1] + 1:
            last[1] = max(last[1], end)

        else:
            merged.append(
                [start, end]
            )

    return [
        tuple(x)
        for x in merged
    ]


def overlap_pair_count(intervals):

    intervals = sorted(intervals)

    n = 0

    for i in range(len(intervals)):

        a_start, a_end = intervals[i]

        for j in range(i + 1, len(intervals)):

            b_start, b_end = intervals[j]

            if b_start > a_end:
                break

            if (
                a_start <= b_end
                and b_start <= a_end
            ):
                n += 1

    return n


def median_or_blank(values):

    values = [
        x
        for x in values
        if x is not None
    ]

    if not values:
        return ""

    return statistics.median(values)


# ======================================================================
# Leer manifest
# ======================================================================

with MANIFEST.open() as fh:

    manifest = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )


if len(manifest) != 271:

    raise RuntimeError(
        f"Manifest esperado=271 observado={len(manifest)}"
    )


# ======================================================================
# Parsear cada EMBL
# ======================================================================

all_regions = []

parse_issues = []

qualifier_presence = Counter()

algorithms = Counter()

notes = Counter()

colours = Counter()

region_counter = 0


for row in manifest:

    iid = row["input_id"]

    embl = (
        PRED
        / iid
        / f"{iid}.alienhunter.embl"
    )

    if not embl.is_file():

        raise RuntimeError(
            f"Falta output AlienHunter: {embl}"
        )


    text = embl.read_text(
        errors="replace"
    )

    current = None
    current_qualifier = None

    features = []




    for raw in text.splitlines():

        # --------------------------------------------------------------
        # Feature principal
        # --------------------------------------------------------------

        m_feature = re.match(
            r"^FT\s{3}(\S+)\s+(.+?)\s*$",
            raw
        )

        if m_feature:

            if current is not None:
                features.append(current)

            current = {
                "feature_key": m_feature.group(1),
                "location": m_feature.group(2).strip(),
                "qualifiers": {},
            }

            current_qualifier = None

            continue


        # --------------------------------------------------------------
        # Qualifier
        # --------------------------------------------------------------

        m_qual = re.match(
            r"^FT\s+/([^=\s]+)(?:=(.*))?\s*$",
            raw
        )

        if m_qual and current is not None:

            key = m_qual.group(1)

            value = (
                m_qual.group(2)
                if m_qual.group(2) is not None
                else ""
            )

            current["qualifiers"][key] = (
                value.strip()
            )

            current_qualifier = key

            continue


        # --------------------------------------------------------------
        # Continuación potencial de qualifier
        # --------------------------------------------------------------

        m_cont = re.match(
            r"^FT\s+(.+?)\s*$",
            raw
        )

        if (
            m_cont
            and current is not None
            and current_qualifier is not None
        ):

            piece = m_cont.group(1).strip()

            # No incorporar nuevas líneas de feature por error
            if not piece.startswith("/"):

                old = current["qualifiers"].get(
                    current_qualifier,
                    ""
                )

                if old:
                    current["qualifiers"][
                        current_qualifier
                    ] = old + " " + piece
                else:
                    current["qualifiers"][
                        current_qualifier
                    ] = piece


    if current is not None:
        features.append(current)


    # ==================================================================
    # Convertir features
    # ==================================================================

    for feature in features:

        if feature["feature_key"] != "misc_feature":

            parse_issues.append({
                "input_id": iid,
                "issue": "unexpected_feature_key",
                "value": feature["feature_key"],
            })

            continue


        location = feature["location"]

        m_loc = re.fullmatch(
            r"(\d+)\.\.(\d+)",
            location
        )

        if not m_loc:

            parse_issues.append({
                "input_id": iid,
                "issue": "unexpected_location",
                "value": location,
            })

            continue


        start = int(
            m_loc.group(1)
        )

        end = int(
            m_loc.group(2)
        )

        contig_length = int(
            row["contig_length_bp"]
        )


        if start < 1 or end < start:

            parse_issues.append({
                "input_id": iid,
                "issue": "invalid_coordinates",
                "value": location,
            })

            continue


        if end > contig_length:

            parse_issues.append({
                "input_id": iid,
                "issue": "coordinate_exceeds_contig",
                "value": (
                    f"{location};"
                    f"contig_length={contig_length}"
                ),
            })

            continue


        q = feature["qualifiers"]


        for expected in [
            "algorithm",
            "note",
            "score",
            "colour",
        ]:

            if expected in q:
                qualifier_presence[
                    expected
                ] += 1


        algorithm = clean_value(
            q.get("algorithm", "")
        )

        note = clean_value(
            q.get("note", "")
        )

        colour = clean_value(
            q.get("colour", "")
        )

        score_raw = clean_value(
            q.get("score", "")
        )

        score = parse_float(
            score_raw
        )


        if algorithm:
            algorithms[algorithm] += 1

        if note:
            notes[note] += 1

        if colour:
            colours[colour] += 1


        region_counter += 1

        region_id = (
            f"GICAND{region_counter:05d}"
        )

        length_bp = (
            end - start + 1
        )

        touches_left = (
            start == 1
        )

        touches_right = (
            end == contig_length
        )

        edge_touching = (
            touches_left
            or touches_right
        )


        all_regions.append({
            "GI_candidate_id": region_id,
            "array_index": row["array_index"],
            "input_id": iid,
            "coassembly": row["coassembly"],
            "MAG": row["MAG"],
            "contig": row["contig"],
            "contig_length_bp": contig_length,
            "start": start,
            "end": end,
            "length_bp": length_bp,
            "fraction_of_contig": (
                length_bp / contig_length
            ),
            "touches_left_edge": (
                "YES"
                if touches_left
                else "NO"
            ),
            "touches_right_edge": (
                "YES"
                if touches_right
                else "NO"
            ),
            "edge_truncation_possible": (
                "YES"
                if edge_touching
                else "NO"
            ),
            "score_raw": score_raw,
            "score_numeric": (
                ""
                if score is None
                else score
            ),
            "algorithm": algorithm,
            "note": note,
            "colour": colour,
            "prediction_source": (
                "AlienHunter_1.7_default_IVOM"
            ),
            "interpretation": (
                "composition_based_candidate_GI_region"
            ),
        })


# ======================================================================
# Validación global contra 93E2C0
# ======================================================================

if len(all_regions) != 1251:

    raise RuntimeError(
        "Número de regiones parseadas distinto de auditoría: "
        f"esperado=1251 observado={len(all_regions)}"
    )


if parse_issues:

    raise RuntimeError(
        f"Se detectaron {len(parse_issues)} problemas de parseo"
    )


for q in [
    "algorithm",
    "note",
    "score",
    "colour",
]:

    if qualifier_presence[q] != 1251:

        raise RuntimeError(
            f"Qualifier {q}: esperado=1251 "
            f"observado={qualifier_presence[q]}"
        )


# ======================================================================
# Tabla completa de candidatos
# ======================================================================

region_fields = list(
    all_regions[0].keys()
)

with (
    OUT
    / "93E2C_candidate_GI_regions.tsv"
).open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=region_fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(
        all_regions
    )


# ======================================================================
# Resumen por contig
# ======================================================================

regions_by_contig = defaultdict(
    list
)

for r in all_regions:

    key = (
        r["coassembly"],
        r["MAG"],
        r["contig"],
    )

    regions_by_contig[
        key
    ].append(r)


contig_summary = []


for row in manifest:

    key = (
        row["coassembly"],
        row["MAG"],
        row["contig"],
    )

    rr = regions_by_contig.get(
        key,
        []
    )

    intervals = [
        (
            int(r["start"]),
            int(r["end"])
        )
        for r in rr
    ]

    merged = merge_intervals(
        intervals
    )

    raw_bp = sum(
        end - start + 1
        for start, end in intervals
    )

    merged_bp = sum(
        end - start + 1
        for start, end in merged
    )

    contig_len = int(
        row["contig_length_bp"]
    )

    scores = []

    for r in rr:

        value = r["score_numeric"]

        if value != "":
            scores.append(
                float(value)
            )


    contig_summary.append({
        "coassembly": row["coassembly"],
        "MAG": row["MAG"],
        "contig": row["contig"],
        "contig_length_bp": contig_len,
        "n_candidate_regions": len(rr),
        "raw_candidate_bp": raw_bp,
        "merged_candidate_bp": merged_bp,
        "fraction_contig_candidate": (
            merged_bp / contig_len
        ),
        "candidate_regions_per_Mb": (
            len(rr)
            / (contig_len / 1_000_000)
        ),
        "n_merged_candidate_blocks": len(merged),
        "n_overlapping_region_pairs": (
            overlap_pair_count(intervals)
        ),
        "n_edge_touching_regions": sum(
            r["edge_truncation_possible"]
            == "YES"
            for r in rr
        ),
        "score_min": (
            min(scores)
            if scores
            else ""
        ),
        "score_median": (
            median_or_blank(scores)
        ),
        "score_max": (
            max(scores)
            if scores
            else ""
        ),
    })


contig_fields = list(
    contig_summary[0].keys()
)

with (
    OUT
    / "93E2C_contig_summary.tsv"
).open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=contig_fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(
        contig_summary
    )


# ======================================================================
# Resumen por MAG
# ======================================================================

def summarize_group(
    key_name,
    key_function,
    filename
):

    groups = defaultdict(list)

    for row in contig_summary:

        groups[
            key_function(row)
        ].append(row)


    output = []

    for key in sorted(groups):

        rows = groups[key]

        screened_bp = sum(
            int(x["contig_length_bp"])
            for x in rows
        )

        n_regions = sum(
            int(x["n_candidate_regions"])
            for x in rows
        )

        merged_bp = sum(
            int(x["merged_candidate_bp"])
            for x in rows
        )

        edge_n = sum(
            int(x["n_edge_touching_regions"])
            for x in rows
        )

        positive_contigs = sum(
            int(x["n_candidate_regions"]) > 0
            for x in rows
        )


        output.append({
            key_name: key,
            "screened_contigs_ge50kb": len(rows),
            "screened_bp": screened_bp,
            "contigs_with_candidate_regions": (
                positive_contigs
            ),
            "candidate_regions": n_regions,
            "candidate_regions_per_Mb": (
                n_regions
                / (screened_bp / 1_000_000)
            ),
            "merged_candidate_bp": merged_bp,
            "fraction_screened_bp_candidate": (
                merged_bp / screened_bp
            ),
            "edge_touching_candidate_regions": (
                edge_n
            ),
        })


    fields = list(
        output[0].keys()
    )

    with (
        OUT
        / filename
    ).open(
        "w",
        newline=""
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(output)


    return output


mag_summary = summarize_group(
    "MAG",
    lambda x: x["MAG"],
    "93E2C_MAG_summary.tsv"
)


coassembly_summary = summarize_group(
    "coassembly",
    lambda x: x["coassembly"],
    "93E2C_coassembly_summary.tsv"
)


# ======================================================================
# Score summary
# ======================================================================

scores = []

for r in all_regions:

    if r["score_numeric"] != "":
        scores.append(
            float(
                r["score_numeric"]
            )
        )


edge_count = sum(
    r["edge_truncation_possible"]
    == "YES"
    for r in all_regions
)


total_screened_bp = sum(
    int(x["contig_length_bp"])
    for x in contig_summary
)


total_merged_candidate_bp = sum(
    int(x["merged_candidate_bp"])
    for x in contig_summary
)


positive_contigs = sum(
    int(x["n_candidate_regions"]) > 0
    for x in contig_summary
)


overlap_pairs = sum(
    int(x["n_overlapping_region_pairs"])
    for x in contig_summary
)


# ======================================================================
# Resumen global
# ======================================================================

global_metrics = [
    ("screened_contigs_ge50kb", len(contig_summary)),
    ("screened_MAGs", len(set(x["MAG"] for x in contig_summary))),
    ("screened_bp", total_screened_bp),
    ("candidate_regions", len(all_regions)),
    ("contigs_with_candidate_regions", positive_contigs),
    (
        "candidate_regions_per_Mb",
        len(all_regions)
        / (total_screened_bp / 1_000_000)
    ),
    (
        "merged_candidate_bp",
        total_merged_candidate_bp
    ),
    (
        "fraction_screened_bp_candidate",
        total_merged_candidate_bp
        / total_screened_bp
    ),
    (
        "edge_touching_candidate_regions",
        edge_count
    ),
    (
        "overlapping_candidate_region_pairs",
        overlap_pairs
    ),
    (
        "score_numeric_count",
        len(scores)
    ),
    (
        "score_min",
        min(scores)
        if scores
        else ""
    ),
    (
        "score_median",
        statistics.median(scores)
        if scores
        else ""
    ),
    (
        "score_mean",
        statistics.mean(scores)
        if scores
        else ""
    ),
    (
        "score_max",
        max(scores)
        if scores
        else ""
    ),
    (
        "distinct_algorithm_values",
        len(algorithms)
    ),
    (
        "distinct_note_values",
        len(notes)
    ),
    (
        "distinct_colour_values",
        len(colours)
    ),
    (
        "parse_issues",
        len(parse_issues)
    ),
]


with (
    OUT
    / "93E2C_global_summary.tsv"
).open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    for key, value in global_metrics:
        fh.write(
            f"{key}\t{value}\n"
        )


# ======================================================================
# Valores de qualifiers
# ======================================================================

def write_counter(
    counter,
    filename,
    fieldname
):

    with (
        OUT
        / filename
    ).open("w") as fh:

        fh.write(
            f"{fieldname}\tcount\n"
        )

        for value, count in counter.most_common():

            safe = value.replace(
                "\t",
                " "
            )

            fh.write(
                f"{safe}\t{count}\n"
            )


write_counter(
    algorithms,
    "93E2C_algorithm_values.tsv",
    "algorithm"
)

write_counter(
    notes,
    "93E2C_note_values.tsv",
    "note"
)

write_counter(
    colours,
    "93E2C_colour_values.tsv",
    "colour"
)


# ======================================================================
# Edge-touching candidates
# ======================================================================

edge_regions = [
    r
    for r in all_regions
    if r["edge_truncation_possible"]
    == "YES"
]

with (
    OUT
    / "93E2C_edge_touching_regions.tsv"
).open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=region_fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(
        edge_regions
    )


# ======================================================================
# Parse issues table
# ======================================================================

with (
    OUT
    / "93E2C_parse_issues.tsv"
).open(
    "w",
    newline=""
) as fh:

    fields = [
        "input_id",
        "issue",
        "value",
    ]

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()

    for issue in parse_issues:
        writer.writerow(issue)


# ======================================================================
# Methodological interpretation
# ======================================================================

with (
    OUT
    / "93E2C_methodological_scope.tsv"
).open("w") as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        "prediction_tool\tAlienHunter_1.7\n"
    )

    fh.write(
        "prediction_mode\tdefault_IVOM\n"
    )

    fh.write(
        "input_scope\tfinal18_MAG_contigs_ge_50000bp\n"
    )

    fh.write(
        "prediction_term\t"
        "composition_based_candidate_GI_region\n"
    )

    fh.write(
        "candidate_region_equals_validated_genomic_island\tNO\n"
    )

    fh.write(
        "candidate_region_equals_demonstrated_HGT\tNO\n"
    )

    fh.write(
        "negative_result_equals_GI_absence\tNO\n"
    )

    fh.write(
        "edge_touching_region_complete_boundary_guaranteed\tNO\n"
    )

    fh.write(
        "raw_candidate_bp_used_for_fraction\tNO\n"
    )

    fh.write(
        "merged_nonredundant_candidate_bp_used_for_fraction\tYES\n"
    )


# ======================================================================
# Final validation
# ======================================================================

required = [
    "93E2C_candidate_GI_regions.tsv",
    "93E2C_contig_summary.tsv",
    "93E2C_MAG_summary.tsv",
    "93E2C_coassembly_summary.tsv",
    "93E2C_global_summary.tsv",
    "93E2C_algorithm_values.tsv",
    "93E2C_note_values.tsv",
    "93E2C_colour_values.tsv",
    "93E2C_edge_touching_regions.tsv",
    "93E2C_parse_issues.tsv",
    "93E2C_methodological_scope.tsv",
]


for filename in required:

    path = OUT / filename

    if not path.is_file():

        raise RuntimeError(
            f"Salida faltante: {path}"
        )


print(
    f"CANDIDATE_REGIONS={len(all_regions)}"
)

print(
    f"SCREENED_CONTIGS={len(contig_summary)}"
)

print(
    f"SCREENED_MAGS={len(set(x['MAG'] for x in contig_summary))}"
)

print(
    f"EDGE_TOUCHING={edge_count}"
)

print(
    f"PARSE_ISSUES={len(parse_issues)}"
)

print(
    "93E2C_PARSE=PASS"
)

PY


# ======================================================================
# Marker
# ======================================================================

printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"${SLURM_JOB_ID}" \
> "${OUT}/93E2C_COMPLETE.ok"


echo
echo "======================================================================"
echo "GLOBAL SUMMARY"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93E2C_global_summary.tsv"


echo
echo "======================================================================"
echo "MAG SUMMARY"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93E2C_MAG_summary.tsv"


echo
echo "======================================================================"
echo "COASSEMBLY SUMMARY"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93E2C_coassembly_summary.tsv"


echo
echo "======================================================================"
echo "ALGORITHM VALUES"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93E2C_algorithm_values.tsv"


echo
echo "======================================================================"
echo "93E2C FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "======================================================================"

