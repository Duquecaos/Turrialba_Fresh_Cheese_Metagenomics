#!/bin/bash

#SBATCH --job-name=92g2_votucov
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --array=0-17%4

#SBATCH --mail-type=END,FAIL

#SBATCH --output=92g2_votucov_%A_%a.out
#SBATCH --error=92g2_votucov_%A_%a.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"
SAMTOOLS_SIF="${SAMTOOLS_SIF:-/opt/ohpc/pub/containers/BIO/samtools-1.21.sif}"

MAPROOT="${ROOT}/92_votu_abundance/mapping"
REF="${ROOT}/92_votu_abundance/reference/votu_4491.fna"

OUTROOT="${ROOT}/93_votu_coverage"

mkdir -p "$OUTROOT"

SAMPLES=(
    L1_1
    L1_2
    L1_3
    L2_1
    L2_2
    L2_3
    L3_1
    L3_2
    L3_3
    M1_1
    M1_2
    M1_3
    M2_1
    M2_2
    M2_3
    M3_1
    M3_2
    M3_3
)

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

BAM="${MAPROOT}/${SAMPLE}/${SAMPLE}.votu.sorted.bam"
BAI="${BAM}.bai"

BT2LOG="${MAPROOT}/${SAMPLE}/${SAMPLE}.bowtie2.log"
OLDMET="${MAPROOT}/${SAMPLE}/${SAMPLE}.mapping_metrics.tsv"

SOUT="${OUTROOT}/${SAMPLE}"

mkdir -p "$SOUT"

DONE="${SOUT}/${SAMPLE}.92G2_COMPLETE.ok"

echo "======================================================================"
echo "92G2 - vOTU COVERAGE"
echo "Muestra: ${SAMPLE}"
echo "Job: ${SLURM_JOB_ID}"
echo "Inicio: $(date)"
echo "======================================================================"

for F in \
    "$SAMTOOLS_SIF" \
    "$REF" \
    "$BAM" \
    "$BAI" \
    "$BT2LOG" \
    "$OLDMET"
do
    [[ -s "$F" ]] || {
        echo "ERROR: falta $F" >&2
        exit 2
    }
done

[[ -x "$APPTAINER" ]] || {
    echo "ERROR: falta Apptainer." >&2
    exit 3
}


# ======================================================================
# 1. DENOMINADORES
# ======================================================================

INPUT_PAIRS="$(
    awk '
        NR==1 {
            print $1
            exit
        }
    ' "$BT2LOG"
)"

[[ "$INPUT_PAIRS" =~ ^[0-9]+$ ]] || {
    echo "ERROR: no pude leer pares de Bowtie2." >&2
    exit 4
}

INPUT_READ_RECORDS=$(( INPUT_PAIRS * 2 ))


EXPECTED_PRIMARY="$(
    awk -F '\t' \
        '$1=="primary_mapped_alignment_records"{print $2}' \
        "$OLDMET"
)"

EXPECTED_MAPQ20="$(
    awk -F '\t' \
        '$1=="mapq20_primary_mapped_alignment_records"{print $2}' \
        "$OLDMET"
)"

for X in \
    "$EXPECTED_PRIMARY" \
    "$EXPECTED_MAPQ20"
do
    [[ "$X" =~ ^[0-9]+$ ]] || {
        echo "ERROR: denominador inválido: $X" >&2
        exit 5
    }
done


# ======================================================================
# 2. FUNCIÓN: CONTEO DE REGISTROS POR vOTU
# ======================================================================

count_records () {

    local MODE="$1"
    local MQ="$2"
    local OUT="$3"

    {
        echo -e "votu_id\tmapped_records"

        "$APPTAINER" exec \
            --bind "${ROOT}:${ROOT}" \
            "$SAMTOOLS_SIF" \
            samtools view \
                -F 3844 \
                -q "$MQ" \
                "$BAM" \
        | awk '
            BEGIN {
                OFS="\t"
            }
            {
                n[$3]++
            }
            END {
                for (r in n)
                    print r, n[r]
            }
        ' \
        | sort -k1,1

    } > "$OUT"
}


# ======================================================================
# 3. FUNCIÓN: DEPTH / BREADTH
#
# -aa = incluye todos los vOTUs, aun sin cobertura
# -s  = evita doble conteo en solapamiento de mates
# -q0 = sin filtro adicional de calidad de base
# -Q  = filtro de mapping quality
# ======================================================================

depth_metrics () {

    local MODE="$1"
    local MQ="$2"
    local OUT="$3"

    "$APPTAINER" exec \
        --bind "${ROOT}:${ROOT}" \
        "$SAMTOOLS_SIF" \
        samtools depth \
            -aa \
            -s \
            -q 0 \
            -Q "$MQ" \
            -G UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY \
            "$BAM" \
    | awk \
        -v SAMPLE="$SAMPLE" \
        -v MODE="$MODE" '
        BEGIN {
            OFS="\t"

            print \
                "sample", \
                "mode", \
                "votu_id", \
                "length_bp", \
                "depth_sum", \
                "mean_depth", \
                "covered_bp_1x", \
                "breadth_1x", \
                "breadth_1x_percent", \
                "covered_bp_5x", \
                "breadth_5x", \
                "breadth_5x_percent", \
                "max_depth"
        }

        function emit() {

            if (ref == "")
                return

            mean = sumdepth / npos

            b1 = cov1 / npos
            b5 = cov5 / npos

            print \
                SAMPLE, \
                MODE, \
                ref, \
                npos, \
                sumdepth, \
                mean, \
                cov1, \
                b1, \
                b1 * 100, \
                cov5, \
                b5, \
                b5 * 100, \
                maxdepth
        }

        {
            if (ref != "" && $1 != ref) {

                emit()

                npos=0
                sumdepth=0
                cov1=0
                cov5=0
                maxdepth=0
            }

            if ($1 != ref)
                ref=$1

            d=$3

            npos++
            sumdepth += d

            if (d >= 1)
                cov1++

            if (d >= 5)
                cov5++

            if (d > maxdepth)
                maxdepth=d
        }

        END {
            emit()
        }
    ' > "$OUT"
}


# ======================================================================
# 4. EJECUTAR DOS VISTAS
# ======================================================================

ALL_COUNTS="${SOUT}/${SAMPLE}.ALL_PRIMARY.counts.tsv"
Q20_COUNTS="${SOUT}/${SAMPLE}.MAPQ20.counts.tsv"

ALL_DEPTH="${SOUT}/${SAMPLE}.ALL_PRIMARY.depth.tsv"
Q20_DEPTH="${SOUT}/${SAMPLE}.MAPQ20.depth.tsv"


echo
echo "===== ALL_PRIMARY COUNTS ====="

count_records \
    "ALL_PRIMARY" \
    0 \
    "$ALL_COUNTS"


echo
echo "===== MAPQ20 COUNTS ====="

count_records \
    "MAPQ20" \
    20 \
    "$Q20_COUNTS"


echo
echo "===== ALL_PRIMARY DEPTH ====="

depth_metrics \
    "ALL_PRIMARY" \
    0 \
    "$ALL_DEPTH"


echo
echo "===== MAPQ20 DEPTH ====="

depth_metrics \
    "MAPQ20" \
    20 \
    "$Q20_DEPTH"


# ======================================================================
# 5. VALIDACIONES DE FILAS Y LONGITUD TOTAL
# ======================================================================

for F in \
    "$ALL_DEPTH" \
    "$Q20_DEPTH"
do

    N="$(
        awk 'END{print NR-1}' "$F"
    )"

    BP="$(
        awk -F '\t' \
            'NR>1{s+=$4} END{printf "%.0f\n",s}' \
            "$F"
    )"

    echo "$(basename "$F"): rows=${N} bp=${BP}"

    [[ "$N" -eq 4491 ]] || {
        echo "ERROR: $F no tiene 4491 vOTUs." >&2
        exit 6
    }

    [[ "$BP" -eq 26784056 ]] || {
        echo "ERROR: longitud total inesperada en $F: $BP" >&2
        exit 7
    }

done


# ======================================================================
# 6. INTEGRACIÓN
# ======================================================================

FINAL="${SOUT}/${SAMPLE}.votu_metrics.tsv"
SUMMARY="${SOUT}/${SAMPLE}.92G2_summary.tsv"

python3 - \
    "$SAMPLE" \
    "$INPUT_PAIRS" \
    "$INPUT_READ_RECORDS" \
    "$EXPECTED_PRIMARY" \
    "$EXPECTED_MAPQ20" \
    "$ALL_COUNTS" \
    "$Q20_COUNTS" \
    "$ALL_DEPTH" \
    "$Q20_DEPTH" \
    "$FINAL" \
    "$SUMMARY" <<'PY'

import csv
import sys
from collections import Counter

(
    sample,
    input_pairs,
    input_read_records,
    expected_primary,
    expected_mapq20,
    all_counts_file,
    q20_counts_file,
    all_depth_file,
    q20_depth_file,
    final_file,
    summary_file,
) = sys.argv[1:]

input_pairs = int(input_pairs)
input_read_records = int(input_read_records)

expected = {
    "ALL_PRIMARY": int(expected_primary),
    "MAPQ20": int(expected_mapq20),
}


def read_counts(path):
    d = {}

    with open(path) as fh:
        r = csv.DictReader(
            fh,
            delimiter="\t"
        )

        for row in r:
            d[row["votu_id"]] = int(
                row["mapped_records"]
            )

    return d


def read_depth(path):
    with open(path) as fh:
        return list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


counts_by_mode = {
    "ALL_PRIMARY":
        read_counts(all_counts_file),

    "MAPQ20":
        read_counts(q20_counts_file),
}


rows_by_mode = {
    "ALL_PRIMARY":
        read_depth(all_depth_file),

    "MAPQ20":
        read_depth(q20_depth_file),
}


final_rows = []
summary_rows = []


for mode in [
    "ALL_PRIMARY",
    "MAPQ20",
]:

    depth_rows = rows_by_mode[mode]
    counts = counts_by_mode[mode]

    if len(depth_rows) != 4491:
        raise RuntimeError(
            f"{mode}: depth rows={len(depth_rows)}"
        )

    mapped_total = sum(
        counts.values()
    )

    if mapped_total != expected[mode]:
        raise RuntimeError(
            f"{mode}: mapped records "
            f"{mapped_total} != "
            f"92G1 {expected[mode]}"
        )


    for row in depth_rows:

        votu = row["votu_id"]

        length = int(
            row["length_bp"]
        )

        mapped = counts.get(
            votu,
            0
        )

        mean_depth = float(
            row["mean_depth"]
        )

        b1 = float(
            row["breadth_1x"]
        )

        b5 = float(
            row["breadth_5x"]
        )


        final_rows.append({

            "sample":
                sample,

            "mode":
                mode,

            "votu_id":
                votu,

            "length_bp":
                length,

            "mapped_records":
                mapped,

            "mapped_records_per_million_input_read_records":
                (
                    mapped
                    * 1_000_000
                    / input_read_records
                ),

            "mapped_records_per_million_mode_mapped_records":
                (
                    mapped
                    * 1_000_000
                    / mapped_total
                    if mapped_total
                    else 0
                ),

            "covered_bp_1x":
                int(
                    row["covered_bp_1x"]
                ),

            "breadth_1x":
                b1,

            "breadth_1x_percent":
                float(
                    row[
                        "breadth_1x_percent"
                    ]
                ),

            "covered_bp_5x":
                int(
                    row["covered_bp_5x"]
                ),

            "breadth_5x":
                b5,

            "breadth_5x_percent":
                float(
                    row[
                        "breadth_5x_percent"
                    ]
                ),

            "depth_sum":
                int(
                    float(
                        row["depth_sum"]
                    )
                ),

            "mean_depth":
                mean_depth,

            "max_depth":
                int(
                    float(
                        row["max_depth"]
                    )
                ),

            "mean_depth_per_million_input_read_records":
                (
                    mean_depth
                    * 1_000_000
                    / input_read_records
                ),

            "mean_depth_per_million_mode_mapped_records":
                (
                    mean_depth
                    * 1_000_000
                    / mapped_total
                    if mapped_total
                    else 0
                ),

            "input_read_pairs":
                input_pairs,

            "input_read_records":
                input_read_records,

            "mode_total_mapped_records":
                mapped_total,
        })


    # ------------------------------------------------------
    # Descriptivo, NO presencia biológica
    # ------------------------------------------------------

    summary_rows.extend([

        {
            "sample": sample,
            "mode": mode,
            "metric": "mode_total_mapped_records",
            "value": mapped_total,
        },

        {
            "sample": sample,
            "mode": mode,
            "metric": "votu_mapped_records_gt0",
            "value": sum(
                counts.get(
                    x["votu_id"],
                    0
                ) > 0
                for x in depth_rows
            ),
        },

        {
            "sample": sample,
            "mode": mode,
            "metric": "breadth1_ge_0.10",
            "value": sum(
                float(
                    x["breadth_1x"]
                ) >= 0.10
                for x in depth_rows
            ),
        },

        {
            "sample": sample,
            "mode": mode,
            "metric": "breadth1_ge_0.50",
            "value": sum(
                float(
                    x["breadth_1x"]
                ) >= 0.50
                for x in depth_rows
            ),
        },

        {
            "sample": sample,
            "mode": mode,
            "metric": "breadth1_ge_0.80",
            "value": sum(
                float(
                    x["breadth_1x"]
                ) >= 0.80
                for x in depth_rows
            ),
        },

        {
            "sample": sample,
            "mode": mode,
            "metric": "breadth1_ge_0.95",
            "value": sum(
                float(
                    x["breadth_1x"]
                ) >= 0.95
                for x in depth_rows
            ),
        },

        {
            "sample": sample,
            "mode": mode,
            "metric": "breadth5_ge_0.10",
            "value": sum(
                float(
                    x["breadth_5x"]
                ) >= 0.10
                for x in depth_rows
            ),
        },

        {
            "sample": sample,
            "mode": mode,
            "metric": "breadth5_ge_0.50",
            "value": sum(
                float(
                    x["breadth_5x"]
                ) >= 0.50
                for x in depth_rows
            ),
        },

        {
            "sample": sample,
            "mode": mode,
            "metric": "breadth5_ge_0.80",
            "value": sum(
                float(
                    x["breadth_5x"]
                ) >= 0.80
                for x in depth_rows
            ),
        },
    ])


with open(
    final_file,
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=list(
            final_rows[0].keys()
        ),
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(
        final_rows
    )


with open(
    summary_file,
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=[
            "sample",
            "mode",
            "metric",
            "value",
        ],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(
        summary_rows
    )


print(
    "sample",
    sample,
    sep="\t"
)

print(
    "input_read_pairs",
    input_pairs,
    sep="\t"
)

print(
    "input_read_records",
    input_read_records,
    sep="\t"
)

for mode in [
    "ALL_PRIMARY",
    "MAPQ20",
]:

    print(
        mode,
        expected[mode],
        sep="\t"
    )

PY


# ======================================================================
# 7. VALIDACIÓN FINAL
# ======================================================================

NFINAL="$(
    awk 'END{print NR-1}' "$FINAL"
)"

[[ "$NFINAL" -eq 8982 ]] || {
    echo "ERROR: final tiene $NFINAL filas; esperaba 8982." >&2
    exit 8
}


{
    echo -e "field\tvalue"
    echo -e "analysis\tvOTU_breadth_depth"
    echo -e "sample\t${SAMPLE}"
    echo -e "reference_vOTUs\t4491"
    echo -e "modes\tALL_PRIMARY;MAPQ20"
    echo -e "samtools_version\t1.21"
    echo -e "depth_overlap_handling\t-s_first_read_only_in_overlap"
    echo -e "base_quality_min\t0"
    echo -e "ALL_PRIMARY_min_MAPQ\t0"
    echo -e "MAPQ20_min_MAPQ\t20"
    echo -e "excluded_flags\tUNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY"
    echo -e "breadth_thresholds_descriptive_only\tYES"
    echo -e "biological_presence_threshold\tNOT_DEFINED"
    echo -e "date\t$(date -Iseconds)"
} > "${SOUT}/${SAMPLE}.92G2_metadata.tsv"


printf \
"COMPLETED\t%s\tjob=%s\tsample=%s\n" \
"$(date -Iseconds)" \
"$SLURM_JOB_ID" \
"$SAMPLE" \
> "$DONE"


echo
echo "======================================================================"
echo "RESUMEN ${SAMPLE}"
echo "======================================================================"

column -t -s $'\t' \
    "$SUMMARY"


echo
echo "FINAL_ROWS=${NFINAL}"

echo
echo "======================================================================"
echo "92G2 ${SAMPLE} FINALIZÓ CORRECTAMENTE."
echo "NO se definió presencia/ausencia biológica."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "======================================================================"

