#!/bin/bash

#SBATCH --job-name=92g3_id95
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --array=0-17%4

#SBATCH --mail-type=END,FAIL

#SBATCH --output=92g3_id95_%A_%a.out
#SBATCH --error=92g3_id95_%A_%a.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"
SAMTOOLS_SIF="${SAMTOOLS_SIF:-/opt/ohpc/pub/containers/BIO/samtools-1.21.sif}"

MAPROOT="${ROOT}/92_votu_abundance/mapping"
PREVROOT="${ROOT}/93_votu_coverage"
OUTROOT="${ROOT}/93_votu_coverage/identity95"

SAMPLES=(
    L1_1 L1_2 L1_3
    L2_1 L2_2 L2_3
    L3_1 L3_2 L3_3
    M1_1 M1_2 M1_3
    M2_1 M2_2 M2_3
    M3_1 M3_2 M3_3
)

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

BAM="${MAPROOT}/${SAMPLE}/${SAMPLE}.votu.sorted.bam"
PREV="${PREVROOT}/${SAMPLE}/${SAMPLE}.votu_metrics.tsv"

SOUT="${OUTROOT}/${SAMPLE}"
mkdir -p "$SOUT"

FILTER="${SOUT}/filter_identity95.py"
AUDIT="${SOUT}/${SAMPLE}.identity95_filter_audit.tsv"

FBAM="${SOUT}/${SAMPLE}.MAPQ20_ID95.bam"
IDXSTATS="${SOUT}/${SAMPLE}.MAPQ20_ID95.idxstats.tsv"
DEPTH="${SOUT}/${SAMPLE}.MAPQ20_ID95.depth.tsv"

FINAL="${SOUT}/${SAMPLE}.MAPQ20_ID95.votu_metrics.tsv"
SUMMARY="${SOUT}/${SAMPLE}.MAPQ20_ID95.summary.tsv"
DONE="${SOUT}/${SAMPLE}.92G3_COMPLETE.ok"

echo "======================================================================"
echo "92G3 MAPQ20 + EDIT IDENTITY >=95%"
echo "Muestra: ${SAMPLE}"
echo "Inicio: $(date)"
echo "======================================================================"

for F in \
    "$BAM" \
    "${BAM}.bai" \
    "$PREV" \
    "$SAMTOOLS_SIF"
do
    [[ -s "$F" ]] || {
        echo "ERROR: falta $F" >&2
        exit 2
    }
done


# ======================================================================
# 1. FILTRO SAM POR NM/CIGAR
# ======================================================================

cat > "$FILTER" <<'PY'
#!/usr/bin/env python3

import re
import sys

audit_file = sys.argv[1]

cigar_re = re.compile(
    r"(\d+)([MIDNSHP=X])"
)

seen = 0
retained = 0
dropped = 0
missing_nm = 0

bins = {
    "lt_90": 0,
    "90_to_lt95": 0,
    "95_to_lt97": 0,
    "97_to_lt99": 0,
    "99_to_lt100": 0,
    "eq100": 0,
}

for line in sys.stdin:

    if line.startswith("@"):
        sys.stdout.write(line)
        continue

    seen += 1

    fields = line.rstrip("\n").split("\t")

    if len(fields) < 12:
        raise RuntimeError(
            "Registro SAM incompleto"
        )

    cigar = fields[5]

    nm = None

    for tag in fields[11:]:
        if tag.startswith("NM:i:"):
            nm = int(
                tag.split(":")[-1]
            )
            break

    if nm is None:
        missing_nm += 1
        continue

    alignment_columns = 0

    for n, op in cigar_re.findall(cigar):

        n = int(n)

        if op in {
            "M",
            "I",
            "D",
            "=",
            "X",
        }:
            alignment_columns += n

    if alignment_columns <= 0:
        raise RuntimeError(
            f"CIGAR sin columnas alineadas: {cigar}"
        )

    identity = (
        1.0
        - nm / alignment_columns
    )

    if identity < 0.90:
        bins["lt_90"] += 1

    elif identity < 0.95:
        bins["90_to_lt95"] += 1

    elif identity < 0.97:
        bins["95_to_lt97"] += 1

    elif identity < 0.99:
        bins["97_to_lt99"] += 1

    elif identity < 1.0:
        bins["99_to_lt100"] += 1

    else:
        bins["eq100"] += 1


    if identity >= 0.95:

        retained += 1
        sys.stdout.write(line)

    else:

        dropped += 1


with open(
    audit_file,
    "w"
) as fh:

    fh.write(
        "metric\tvalue\n"
    )

    fh.write(
        f"mapq20_records_seen\t{seen}\n"
    )

    fh.write(
        f"retained_identity_ge95\t{retained}\n"
    )

    fh.write(
        f"dropped_identity_lt95\t{dropped}\n"
    )

    fh.write(
        f"missing_NM\t{missing_nm}\n"
    )

    for key, value in bins.items():

        fh.write(
            f"identity_bin_{key}\t{value}\n"
        )


if missing_nm != 0:

    raise RuntimeError(
        f"{missing_nm} registros sin NM"
    )

if retained + dropped != seen:

    raise RuntimeError(
        "Balance de registros inconsistente"
    )

PY


# ======================================================================
# 2. FILTRAR BAM
#
# -F 3844:
# UNMAP + SECONDARY + QCFAIL + DUP + SUPPLEMENTARY
#
# -q 20:
# MAPQ >=20
# ======================================================================

rm -f \
    "$FBAM" \
    "${FBAM}.bai" \
    "$AUDIT"


export \
    APPTAINER \
    ROOT \
    SAMTOOLS_SIF \
    BAM \
    FILTER \
    AUDIT \
    FBAM


srun --exact \
    -n 1 \
    -c 2 \
    bash -c '
        set -euo pipefail

        "$APPTAINER" exec \
            --bind "${ROOT}:${ROOT}" \
            "$SAMTOOLS_SIF" \
            samtools view \
                -h \
                -F 3844 \
                -q 20 \
                "$BAM" \
        | \
        python3 \
            "$FILTER" \
            "$AUDIT" \
        | \
        "$APPTAINER" exec \
            --bind "${ROOT}:${ROOT}" \
            "$SAMTOOLS_SIF" \
            samtools view \
                -b \
                -o "$FBAM" \
                -
    '


[[ -s "$FBAM" ]] || {
    echo "ERROR: BAM ID95 vacío" >&2
    exit 3
}

[[ -s "$AUDIT" ]] || {
    echo "ERROR: falta auditoría ID95" >&2
    exit 4
}


# ======================================================================
# 3. QUICKCHECK + INDEX
# ======================================================================

srun --exact \
    -n 1 \
    -c 1 \
    "$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$SAMTOOLS_SIF" \
    samtools quickcheck \
        -v \
        "$FBAM"


srun --exact \
    -n 1 \
    -c 2 \
    "$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$SAMTOOLS_SIF" \
    samtools index \
        -@ 2 \
        "$FBAM"


[[ -s "${FBAM}.bai" ]] || {
    echo "ERROR: falta índice BAM" >&2
    exit 5
}


# ======================================================================
# 4. IDXSTATS
# ======================================================================

srun --exact \
    -n 1 \
    -c 1 \
    "$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$SAMTOOLS_SIF" \
    samtools idxstats \
        "$FBAM" \
    > "$IDXSTATS"


REFS="$(
    awk '$1!="*"{n++} END{print n+0}' \
        "$IDXSTATS"
)"

[[ "$REFS" -eq 4491 ]] || {
    echo "ERROR: referencias=${REFS}" >&2
    exit 6
}


# ======================================================================
# 5. DEPTH
#
# Ya está filtrado:
# primary + MAPQ20 + identity>=95%.
#
# Inicialización explícita corrige el bug de 92G2.
# ======================================================================

export FBAM DEPTH SAMPLE

srun --exact \
    -n 1 \
    -c 1 \
    bash -c '
        set -euo pipefail

        "$APPTAINER" exec \
            --bind "${ROOT}:${ROOT}" \
            "$SAMTOOLS_SIF" \
            samtools depth \
                -aa \
                -s \
                -q 0 \
                -Q 0 \
                "$FBAM" \
        | awk \
            -v SAMPLE="$SAMPLE" "
            BEGIN {
                OFS=\"\t\"

                print \
                    \"sample\", \
                    \"votu_id\", \
                    \"length_bp\", \
                    \"depth_sum\", \
                    \"covered_bp_1x\", \
                    \"covered_bp_5x\", \
                    \"max_depth\"
            }

            function emit() {

                if (ref == \"\")
                    return

                print \
                    SAMPLE, \
                    ref, \
                    npos, \
                    sumdepth, \
                    cov1, \
                    cov5, \
                    maxdepth
            }

            {
                if (\$1 != ref) {

                    if (ref != \"\")
                        emit()

                    ref = \$1

                    npos = 0
                    sumdepth = 0
                    cov1 = 0
                    cov5 = 0
                    maxdepth = 0
                }

                d = \$3

                npos++
                sumdepth += d

                if (d >= 1)
                    cov1++

                if (d >= 5)
                    cov5++

                if (d > maxdepth)
                    maxdepth = d
            }

            END {
                emit()
            }
        " > "$DEPTH"
    '


NDEPTH="$(
    awk 'END{print NR-1}' "$DEPTH"
)"

[[ "$NDEPTH" -eq 4491 ]] || {
    echo "ERROR: depth rows=${NDEPTH}" >&2
    exit 7
}


# ======================================================================
# 6. INTEGRACIÓN + COMPARACIÓN CON MAPQ20
# ======================================================================

python3 - \
    "$SAMPLE" \
    "$AUDIT" \
    "$IDXSTATS" \
    "$DEPTH" \
    "$PREV" \
    "$FINAL" \
    "$SUMMARY" <<'PY'

import csv
import sys

(
    sample,
    audit_file,
    idxstats_file,
    depth_file,
    previous_file,
    final_file,
    summary_file,
) = sys.argv[1:]


# ------------------------------------------------------------
# Auditoría filtro
# ------------------------------------------------------------

audit = {}

with open(audit_file) as fh:

    r = csv.reader(
        fh,
        delimiter="\t"
    )

    next(r)

    for key, value in r:
        audit[key] = int(value)


# ------------------------------------------------------------
# Counts ID95
# ------------------------------------------------------------

counts = {}

with open(idxstats_file) as fh:

    for line in fh:

        ref, length, mapped, unmapped = (
            line.rstrip("\n").split("\t")
        )

        if ref == "*":
            continue

        counts[ref] = int(mapped)


if len(counts) != 4491:
    raise RuntimeError(
        f"idxstats refs={len(counts)}"
    )


retained = sum(
    counts.values()
)

if retained != audit[
    "retained_identity_ge95"
]:

    raise RuntimeError(
        "idxstats total != filtro retained"
    )


# ------------------------------------------------------------
# Previous MAPQ20 metrics
# ------------------------------------------------------------

previous = {}

input_read_records = None
mapq20_total = None

with open(previous_file) as fh:

    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in r:

        if row["mode"] != "MAPQ20":
            continue

        previous[
            row["votu_id"]
        ] = row

        input_read_records = int(
            row[
                "input_read_records"
            ]
        )

        mapq20_total = int(
            row[
                "mode_total_mapped_records"
            ]
        )


if len(previous) != 4491:
    raise RuntimeError(
        f"previous MAPQ20 rows={len(previous)}"
    )

if mapq20_total != audit[
    "mapq20_records_seen"
]:

    raise RuntimeError(
        "MAPQ20 total != filter records seen"
    )


# ------------------------------------------------------------
# Depth ID95
# ------------------------------------------------------------

rows = []

with open(depth_file) as fh:

    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in r:

        votu = row["votu_id"]

        length = int(
            row["length_bp"]
        )

        depth_sum = int(
            row["depth_sum"]
        )

        cov1 = int(
            row["covered_bp_1x"]
        )

        cov5 = int(
            row["covered_bp_5x"]
        )

        maxdepth = int(
            row["max_depth"]
        )

        if not (
            0 <= cov5 <= cov1 <= length
        ):
            raise RuntimeError(
                f"coverage inconsistente: {votu}"
            )

        mean_depth = (
            depth_sum / length
        )

        b1 = (
            cov1 / length
        )

        b5 = (
            cov5 / length
        )

        prev = previous[votu]

        # ID95 debe ser subconjunto de MAPQ20.
        if counts[votu] > int(
            prev["mapped_records"]
        ):
            raise RuntimeError(
                f"{votu}: ID95 records > MAPQ20"
            )

        if cov1 > int(
            prev["covered_bp_1x"]
        ):
            raise RuntimeError(
                f"{votu}: ID95 cov1 > MAPQ20"
            )

        if cov5 > int(
            prev["covered_bp_5x"]
        ):
            raise RuntimeError(
                f"{votu}: ID95 cov5 > MAPQ20"
            )

        if depth_sum > int(
            prev["depth_sum"]
        ):
            raise RuntimeError(
                f"{votu}: ID95 depth > MAPQ20"
            )


        rows.append({

            "sample":
                sample,

            "mode":
                "MAPQ20_ID95",

            "votu_id":
                votu,

            "length_bp":
                length,

            "mapped_records":
                counts[votu],

            "mapped_records_per_million_input_read_records":
                (
                    counts[votu]
                    * 1_000_000
                    / input_read_records
                ),

            "mapped_records_per_million_mode_mapped_records":
                (
                    counts[votu]
                    * 1_000_000
                    / retained
                    if retained
                    else 0
                ),

            "covered_bp_1x":
                cov1,

            "breadth_1x":
                b1,

            "breadth_1x_percent":
                b1 * 100,

            "covered_bp_5x":
                cov5,

            "breadth_5x":
                b5,

            "breadth_5x_percent":
                b5 * 100,

            "depth_sum":
                depth_sum,

            "mean_depth":
                mean_depth,

            "max_depth":
                maxdepth,

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
                    / retained
                    if retained
                    else 0
                ),

            "input_read_records":
                input_read_records,

            "mode_total_mapped_records":
                retained,
        })


if len(rows) != 4491:
    raise RuntimeError(
        f"final rows={len(rows)}"
    )


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

summary = {

    "mapq20_records":
        mapq20_total,

    "identity95_records":
        retained,

    "identity95_retention_fraction":
        (
            retained / mapq20_total
            if mapq20_total
            else 0
        ),

    "votu_mapped_records_gt0":
        sum(
            x["mapped_records"] > 0
            for x in rows
        ),

    "breadth1_ge_0.10":
        sum(
            x["breadth_1x"] >= 0.10
            for x in rows
        ),

    "breadth1_ge_0.50":
        sum(
            x["breadth_1x"] >= 0.50
            for x in rows
        ),

    "breadth1_ge_0.80":
        sum(
            x["breadth_1x"] >= 0.80
            for x in rows
        ),

    "breadth1_ge_0.95":
        sum(
            x["breadth_1x"] >= 0.95
            for x in rows
        ),

    "breadth5_ge_0.10":
        sum(
            x["breadth_5x"] >= 0.10
            for x in rows
        ),

    "breadth5_ge_0.50":
        sum(
            x["breadth_5x"] >= 0.50
            for x in rows
        ),

    "breadth5_ge_0.80":
        sum(
            x["breadth_5x"] >= 0.80
            for x in rows
        ),
}


with open(
    final_file,
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            rows[0].keys()
        ),
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(rows)


with open(
    summary_file,
    "w",
    newline=""
) as fh:

    w = csv.writer(
        fh,
        delimiter="\t"
    )

    w.writerow([
        "sample",
        "mode",
        "metric",
        "value",
    ])

    for key, value in summary.items():

        w.writerow([
            sample,
            "MAPQ20_ID95",
            key,
            value,
        ])


print(
    "sample",
    sample,
    sep="\t"
)

for key, value in summary.items():

    print(
        key,
        value,
        sep="\t"
    )

PY


# ======================================================================
# 7. METADATA
# ======================================================================

{
    echo -e "field\tvalue"
    echo -e "analysis\tvOTU_MAPQ20_edit_identity95_sensitivity"
    echo -e "source_mapping\t92G1"
    echo -e "source_MAPQ20_coverage\t92G2b"
    echo -e "mapping_recomputed\tNO"
    echo -e "minimum_MAPQ\t20"
    echo -e "minimum_edit_identity\t0.95"
    echo -e "identity_definition\t1-NM/(M+I+D+EQ+X)"
    echo -e "identity_is_ANI\tNO"
    echo -e "primary_only\tYES"
    echo -e "presence_threshold\tNOT_DEFINED"
    echo -e "date\t$(date -Iseconds)"
} > "${SOUT}/${SAMPLE}.92G3_metadata.tsv"


printf \
"COMPLETED\t%s\tjob=%s\tsample=%s\n" \
"$(date -Iseconds)" \
"$SLURM_JOB_ID" \
"$SAMPLE" \
> "$DONE"


echo
echo "===== FILTER AUDIT ====="

column -t -s $'\t' "$AUDIT"

echo
echo "===== SUMMARY ====="

column -t -s $'\t' "$SUMMARY"

echo
echo "======================================================================"
echo "92G3 ${SAMPLE} FINALIZO CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "======================================================================"

