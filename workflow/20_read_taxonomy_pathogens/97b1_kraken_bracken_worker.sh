#!/bin/bash

set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

BASE="${ROOT}/98_read_taxonomy"

ENVFILE="${BASE}/kraken_bracken_environment.sh"

OUT="${BASE}/97B1_all_samples"

MANIFEST="${OUT}/97B1_manifest.tsv"


if [[ ! -s "$ENVFILE" ]]
then
    echo "ERROR: environment missing: $ENVFILE" >&2
    exit 10
fi

source "$ENVFILE"


if [[ ! -s "$MANIFEST" ]]
then
    echo "ERROR: manifest missing" >&2
    exit 11
fi


# ============================================================
# Recover array row
# ============================================================

LINE=$(
    sed -n \
    "$((SLURM_ARRAY_TASK_ID + 2))p" \
    "$MANIFEST"
)

LINE="${LINE//$'\r'/}"


if [[ -z "$LINE" ]]
then
    echo "ERROR: empty manifest row" >&2
    exit 12
fi


IFS=$'\t' read -r \
    SAMPLE \
    R1 \
    R2 \
    <<< "$LINE"


if [[ ! -s "$R1" || ! -s "$R2" ]]
then
    echo "ERROR: FASTQ missing for $SAMPLE" >&2
    exit 13
fi


SOUT="${OUT}/${SAMPLE}"

mkdir -p "$SOUT"


# ============================================================
# Paths
# ============================================================

KREPORT="${SOUT}/${SAMPLE}.kraken2.report"

KLOG="${SOUT}/${SAMPLE}.kraken2.log"

BG="${SOUT}/${SAMPLE}.bracken.genus.tsv"

BG_REPORT="${SOUT}/${SAMPLE}.bracken.genus.report"

BS="${SOUT}/${SAMPLE}.bracken.species.tsv"

BS_REPORT="${SOUT}/${SAMPLE}.bracken.species.report"


echo "============================================================"
echo "97B1 - KRAKEN2 + BRACKEN"
echo "Sample: ${SAMPLE}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Start: $(date -Iseconds)"
echo "============================================================"


# ============================================================
# Metadata
# ============================================================

KRAKEN_VERSION_OUTPUT="$(
    "$KRAKEN2" --version 2>&1 |
    head -n 1
)"


BRACKEN_VERSION_OUTPUT="$(
    "$BRACKEN" -v 2>&1 |
    head -n 1
)"


{
    printf "field\tvalue\n"

    printf "sample\t%s\n" \
        "$SAMPLE"

    printf "R1\t%s\n" \
        "$R1"

    printf "R2\t%s\n" \
        "$R2"

    printf "R1_size_bytes\t%s\n" \
        "$(stat -c '%s' "$R1")"

    printf "R2_size_bytes\t%s\n" \
        "$(stat -c '%s' "$R2")"

    printf "kraken2_reported_version\t%s\n" \
        "$KRAKEN_VERSION_OUTPUT"

    printf "kraken2_installation_source_tag\t2.17.2\n"

    printf "bracken_reported_version\t%s\n" \
        "$BRACKEN_VERSION_OUTPUT"

    printf "bracken_installation_source_tag\t3.1\n"

    printf "database\tPlusPF\n"

    printf "database_release\t2026-06-26\n"

    printf "read_length_actual_bp\t151\n"

    printf "bracken_distribution_bp\t%s\n" \
        "$BRACKEN_READ_LEN"

    printf "bracken_threshold_reads\t10\n"

    printf "paired_end\tYES\n"

    printf "host_removed_input\tYES\n"

    printf "raw_per_read_kraken_output_saved\tNO\n"

} > "${SOUT}/${SAMPLE}.metadata.tsv"


# ============================================================
# Kraken2
# ============================================================

"$KRAKEN2" \
    --db "$KRAKEN_DB" \
    --threads "${SLURM_CPUS_PER_TASK:-8}" \
    --paired \
    --gzip-compressed \
    --report "$KREPORT" \
    --output /dev/null \
    "$R1" \
    "$R2" \
    > /dev/null \
    2> "$KLOG"


if [[ ! -s "$KREPORT" ]]
then
    echo "ERROR: Kraken report missing for ${SAMPLE}" >&2
    exit 20
fi


# ============================================================
# Bracken genus
# ============================================================

"$BRACKEN" \
    -d "$KRAKEN_DB" \
    -i "$KREPORT" \
    -o "$BG" \
    -w "$BG_REPORT" \
    -r "$BRACKEN_READ_LEN" \
    -l G \
    -t 10 \
    > "${SOUT}/${SAMPLE}.bracken.genus.log" \
    2>&1


if [[ ! -s "$BG" ]]
then
    echo "ERROR: Bracken genus missing for ${SAMPLE}" >&2
    exit 21
fi


# ============================================================
# Bracken species
# ============================================================

"$BRACKEN" \
    -d "$KRAKEN_DB" \
    -i "$KREPORT" \
    -o "$BS" \
    -w "$BS_REPORT" \
    -r "$BRACKEN_READ_LEN" \
    -l S \
    -t 10 \
    > "${SOUT}/${SAMPLE}.bracken.species.log" \
    2>&1


if [[ ! -s "$BS" ]]
then
    echo "ERROR: Bracken species missing for ${SAMPLE}" >&2
    exit 22
fi


# ============================================================
# Parse outputs
# ============================================================

python3 - \
"$KREPORT" \
"$BG" \
"$BS" \
"$SOUT" \
"$SAMPLE" <<'PY'
import csv
import sys
from pathlib import Path


kreport = Path(sys.argv[1])
genus_file = Path(sys.argv[2])
species_file = Path(sys.argv[3])
out = Path(sys.argv[4])
sample = sys.argv[5]


classified = None
unclassified = None


with kreport.open(
    errors="replace"
) as fh:

    for line in fh:

        parts = line.rstrip(
            "\n"
        ).split("\t")

        if len(parts) < 6:
            continue

        pct, clade, direct, rank, taxid, name = parts[:6]

        if (
            rank == "U"
            and
            taxid == "0"
        ):
            unclassified = int(
                clade
            )

        if (
            rank == "R"
            and
            taxid == "1"
        ):
            classified = int(
                clade
            )


if classified is None:
    raise RuntimeError(
        "Kraken classified count not recovered"
    )

if unclassified is None:
    raise RuntimeError(
        "Kraken unclassified count not recovered"
    )


total = (
    classified
    + unclassified
)


if total <= 0:
    raise RuntimeError(
        "Total sequence pairs <= 0"
    )


def load_bracken(path):

    with path.open(
        newline="",
        errors="replace"
    ) as fh:

        rows = list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )

    if not rows:
        raise RuntimeError(
            f"Empty Bracken output: {path}"
        )

    return rows


genus = load_bracken(
    genus_file
)

species = load_bracken(
    species_file
)


def process_rank(
    rows,
    rank_name
):

    processed = []

    for row in rows:

        est = int(
            float(
                row[
                    "new_est_reads"
                ]
            )
        )

        rec = dict(
            row
        )

        rec[
            "fraction_all_input_pairs"
        ] = (
            f"{est / total:.10f}"
        )

        processed.append(
            rec
        )


    processed.sort(
        key=lambda r:
            int(
                float(
                    r[
                        "new_est_reads"
                    ]
                )
            ),
        reverse=True
    )


    est_sum = sum(
        int(
            float(
                r[
                    "new_est_reads"
                ]
            )
        )
        for r in processed
    )


    outfile = (
        out
        / f"{sample}.bracken.{rank_name}.augmented.tsv"
    )


    fields = (
        list(
            processed[0].keys()
        )
    )


    with outfile.open(
        "w",
        newline=""
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()

        writer.writerows(
            processed
        )


    topfile = (
        out
        / f"{sample}.top20_{rank_name}.tsv"
    )


    topfields = [
        "rank",
        "name",
        "taxonomy_id",
        "taxonomy_lvl",
        "new_est_reads",
        "fraction_total_reads",
        "fraction_all_input_pairs",
    ]


    with topfile.open(
        "w",
        newline=""
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=topfields,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()

        for i, row in enumerate(
            processed[:20],
            1
        ):

            writer.writerow({
                "rank":
                    i,

                "name":
                    row[
                        "name"
                    ],

                "taxonomy_id":
                    row[
                        "taxonomy_id"
                    ],

                "taxonomy_lvl":
                    row[
                        "taxonomy_lvl"
                    ],

                "new_est_reads":
                    row[
                        "new_est_reads"
                    ],

                "fraction_total_reads":
                    row[
                        "fraction_total_reads"
                    ],

                "fraction_all_input_pairs":
                    row[
                        "fraction_all_input_pairs"
                    ],
            })


    return {
        "rows":
            len(processed),

        "estimated_reads_sum":
            est_sum,

        "resolved_fraction_all_pairs":
            (
                est_sum
                / total
            ),
    }


g = process_rank(
    genus,
    "genus"
)

s = process_rank(
    species,
    "species"
)


# Explicit host-like Homo genus signal
homo_est = 0

for row in genus:

    if (
        row["name"].strip()
        == "Homo"
    ):

        homo_est += int(
            float(
                row[
                    "new_est_reads"
                ]
            )
        )


summary = [
    (
        "sample",
        sample
    ),
    (
        "total_sequence_pairs",
        total
    ),
    (
        "kraken_classified_pairs",
        classified
    ),
    (
        "kraken_unclassified_pairs",
        unclassified
    ),
    (
        "kraken_classified_pct",
        f"{100 * classified / total:.6f}"
    ),
    (
        "kraken_unclassified_pct",
        f"{100 * unclassified / total:.6f}"
    ),
    (
        "bracken_genus_rows",
        g["rows"]
    ),
    (
        "bracken_genus_estimated_reads_sum",
        g[
            "estimated_reads_sum"
        ]
    ),
    (
        "bracken_genus_resolved_pct_all_pairs",
        f"{100 * g['resolved_fraction_all_pairs']:.6f}"
    ),
    (
        "bracken_species_rows",
        s["rows"]
    ),
    (
        "bracken_species_estimated_reads_sum",
        s[
            "estimated_reads_sum"
        ]
    ),
    (
        "bracken_species_resolved_pct_all_pairs",
        f"{100 * s['resolved_fraction_all_pairs']:.6f}"
    ),
    (
        "Homo_genus_estimated_reads",
        homo_est
    ),
    (
        "Homo_genus_pct_all_pairs",
        f"{100 * homo_est / total:.8f}"
    ),
    (
        "bracken_distribution_bp",
        150
    ),
    (
        "actual_read_length_bp",
        151
    ),
    (
        "raw_kraken_per_read_output_saved",
        "NO"
    ),
]


with (
    out
    / f"{sample}.classification_summary.tsv"
).open(
    "w",
    newline=""
) as fh:

    writer = csv.writer(
        fh,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writerow([
        "metric",
        "value",
    ])

    writer.writerows(
        summary
    )


print(
    f"SAMPLE={sample}"
)

print(
    f"TOTAL_PAIRS={total}"
)

print(
    f"CLASSIFIED_PCT={100 * classified / total:.6f}"
)

print(
    f"GENUS_RESOLVED_PCT_ALL_PAIRS="
    f"{100 * g['resolved_fraction_all_pairs']:.6f}"
)

print(
    f"SPECIES_RESOLVED_PCT_ALL_PAIRS="
    f"{100 * s['resolved_fraction_all_pairs']:.6f}"
)

print(
    f"HOMO_PCT_ALL_PAIRS="
    f"{100 * homo_est / total:.8f}"
)

print(
    "97B1_PARSE=PASS"
)
PY


# ============================================================
# Validation
# ============================================================

for F in \
    "${SOUT}/${SAMPLE}.classification_summary.tsv" \
    "${SOUT}/${SAMPLE}.bracken.genus.augmented.tsv" \
    "${SOUT}/${SAMPLE}.bracken.species.augmented.tsv" \
    "${SOUT}/${SAMPLE}.top20_genus.tsv" \
    "${SOUT}/${SAMPLE}.top20_species.tsv"
do

    if [[ ! -s "$F" ]]
    then
        echo "ERROR: missing output: $F" >&2
        exit 30
    fi

done


printf \
"COMPLETED\t%s\tjob=%s\tarray_task=%s\n" \
"$(date -Iseconds)" \
"${SLURM_ARRAY_JOB_ID:-NA}" \
"${SLURM_ARRAY_TASK_ID:-NA}" \
> "${SOUT}/97B1_${SAMPLE}.ok"


echo
echo "============================================================"
echo "97B1 ${SAMPLE} FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"

