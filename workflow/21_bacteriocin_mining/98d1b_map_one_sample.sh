#!/bin/bash

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C


ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

BASE="${ROOT}/98_bacteriocin_mining/98D1_competitive_mapping_v2"

REFDIR="${BASE}/reference"

REF="${REFDIR}/98D1AR_competitive_locus_windows_30.fna"

META="${REFDIR}/98D1AR_competitive_locus_windows_30_metadata.tsv"

SAMPLE_MANIFEST="${REFDIR}/98D1AR_sample_manifest.tsv"

INDEX="${REFDIR}/bowtie2_index/bacteriocin_locuswindows30"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"

BOWTIE2="${BOWTIE2:-/opt/ohpc/pub/containers/BIO/bowtie2-2.5.4.sif}"

SAMTOOLS="${SAMTOOLS:-/opt/ohpc/pub/containers/BIO/samtools-1.21.sif}"

PY="${SCRIPT_DIR}/98d1b_summarize_mapping.py"


ARRAY_INDEX="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID not set}"


SAMPLE=$(
    awk -F $'\t' \
        -v IDX="$ARRAY_INDEX" '
        NR > 1 && $1 == IDX {
            print $2
        }
    ' "$SAMPLE_MANIFEST"
)


if [[ -z "$SAMPLE" ]]
then

    echo \
    "ERROR: no sample for array index ${ARRAY_INDEX}" \
    >&2

    exit 20

fi


R1="${ROOT}/04_host_removed/${SAMPLE}/${SAMPLE}_1.hostfree.fastq.gz"

R2="${ROOT}/04_host_removed/${SAMPLE}/${SAMPLE}_2.hostfree.fastq.gz"


SAMPLEBASE="${BASE}/samples"

OUT="${SAMPLEBASE}/${SAMPLE}"

WORK="${BASE}/.working_${SAMPLE}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"


mkdir -p "$SAMPLEBASE"


# Safe restart behavior
if [[ -s "${OUT}/98D1B_COMPLETE.ok" ]]
then

    echo "SAMPLE_ALREADY_COMPLETE=${SAMPLE}"

    exit 0

fi


rm -rf \
    "$WORK" \
    "$OUT"

mkdir -p "$WORK"


BAM="${WORK}/${SAMPLE}.contexts30.sorted.bam"

BTLOG="${WORK}/${SAMPLE}.bowtie2.log"


BOWTIE_THREADS=3
SORT_THREADS=1


echo "============================================================"
echo "98D1B - COMPETITIVE MAPPING"
echo "SAMPLE=${SAMPLE}"
echo "ARRAY_INDEX=${ARRAY_INDEX}"
echo "START=$(date -Iseconds)"
echo "NODE=$(hostname)"
echo "============================================================"


# ============================================================
# 1. Preflight
# ============================================================

for F in \
    "$REF" \
    "$META" \
    "$SAMPLE_MANIFEST" \
    "$R1" \
    "$R2" \
    "$BOWTIE2" \
    "$SAMTOOLS" \
    "$PY"
do

    [[ -s "$F" ]] || {
        echo "ERROR: missing $F" >&2
        exit 21
    }

done


NREF=$(
    grep -c '^>' "$REF"
)


[[ "$NREF" -eq 30 ]] || {
    echo "ERROR: expected 30 references; found $NREF" >&2
    exit 22
}


NINDEX=$(
    find "${REFDIR}/bowtie2_index" \
        -maxdepth 1 \
        -type f \
        \( \
            -name 'bacteriocin_locuswindows30*.bt2' \
            -o \
            -name 'bacteriocin_locuswindows30*.bt2l' \
        \) \
    | wc -l
)


[[ "$NINDEX" -eq 6 ]] || {
    echo "ERROR: Bowtie2 index incomplete: $NINDEX/6" >&2
    exit 23
}


echo "REFERENCE_PREFLIGHT=30/30"
echo "INDEX_PREFLIGHT=6/6"


# ============================================================
# 2. Competitive Bowtie2 mapping
#
# No -a / -k:
# one primary/best alignment representation per read.
# Ambiguous recruitment is assessed via MAPQ.
# ============================================================

"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$BOWTIE2" \
    bowtie2 \
    --very-sensitive \
    --threads "$BOWTIE_THREADS" \
    --rg-id "$SAMPLE" \
    --rg "SM:${SAMPLE}" \
    --no-unal \
    -x "$INDEX" \
    -1 "$R1" \
    -2 "$R2" \
    2> "$BTLOG" \
| \
"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$SAMTOOLS" \
    samtools sort \
    -@ "$SORT_THREADS" \
    -m 2G \
    -o "$BAM" \
    -


[[ -s "$BAM" ]] || {
    echo "ERROR: BAM empty" >&2
    exit 30
}


# ============================================================
# 3. BAM QC / index
# ============================================================

"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$SAMTOOLS" \
    samtools quickcheck \
    -v \
    "$BAM"


"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$SAMTOOLS" \
    samtools index \
    -@ 2 \
    "$BAM"


[[ -s "${BAM}.bai" ]] || {
    echo "ERROR: BAM index absent" >&2
    exit 31
}


# ============================================================
# 4. IDXSTATS
# ============================================================

"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$SAMTOOLS" \
    samtools idxstats \
    "$BAM" \
> "${WORK}/${SAMPLE}.idxstats.tsv"


NIDX=$(
    awk -F $'\t' '
        $1 != "*" {
            n++
        }
        END {
            print n+0
        }
    ' "${WORK}/${SAMPLE}.idxstats.tsv"
)


[[ "$NIDX" -eq 30 ]] || {
    echo "ERROR: idxstats references=$NIDX" >&2
    exit 32
}


# ============================================================
# 5. Mapping counts
# ============================================================

MAPPED_ALL=$(
    "$APPTAINER" exec \
        --bind "${ROOT}:${ROOT}" \
        "$SAMTOOLS" \
        samtools view \
        -c \
        -F 4 \
        "$BAM"
)


MAPPED_Q10=$(
    "$APPTAINER" exec \
        --bind "${ROOT}:${ROOT}" \
        "$SAMTOOLS" \
        samtools view \
        -c \
        -F 4 \
        -q 10 \
        "$BAM"
)


{
    printf "metric\tvalue\n"

    printf "mapped_alignment_records_all\t%s\n" \
        "$MAPPED_ALL"

    printf "mapped_alignment_records_MAPQ10\t%s\n" \
        "$MAPPED_Q10"

} > "${WORK}/${SAMPLE}.mapping_counts.tsv"


echo "MAPPED_ALIGNMENT_RECORDS_ALL=$MAPPED_ALL"
echo "MAPPED_ALIGNMENT_RECORDS_MAPQ10=$MAPPED_Q10"


# ============================================================
# 6. MAPQ distribution
# ============================================================

{
    printf "MAPQ\talignment_records\n"

    "$APPTAINER" exec \
        --bind "${ROOT}:${ROOT}" \
        "$SAMTOOLS" \
        samtools view \
        "$BAM" \
    | awk '
        {
            n[$5]++
        }
        END {
            for (q in n) {
                print q "\t" n[q]
            }
        }
    ' \
    | sort -n -k1,1

} > "${WORK}/${SAMPLE}.MAPQ_distribution.tsv"


# ============================================================
# 7. Depth — all primary mappings
#
# -aa ensures zero-covered positions are retained.
# -q = minimum base quality
# -Q = minimum mapping quality
# ============================================================

"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$SAMTOOLS" \
    samtools depth \
    -aa \
    -q 0 \
    -Q 0 \
    "$BAM" \
| gzip -c \
> "${WORK}/${SAMPLE}.depth_ALL.tsv.gz"


# ============================================================
# 8. Depth — conservative MAPQ >=10
# ============================================================

"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$SAMTOOLS" \
    samtools depth \
    -aa \
    -q 0 \
    -Q 10 \
    "$BAM" \
| gzip -c \
> "${WORK}/${SAMPLE}.depth_MAPQ10.tsv.gz"


[[ -s "${WORK}/${SAMPLE}.depth_ALL.tsv.gz" ]] || {
    echo "ERROR: ALL depth absent" >&2
    exit 33
}


[[ -s "${WORK}/${SAMPLE}.depth_MAPQ10.tsv.gz" ]] || {
    echo "ERROR: MAPQ10 depth absent" >&2
    exit 34
}


# ============================================================
# 9. Summarize all 30 contexts
# ============================================================

python3 \
    "$PY" \
    "$META" \
    "$SAMPLE" \
    "${WORK}/${SAMPLE}.depth_ALL.tsv.gz" \
    "${WORK}/${SAMPLE}.depth_MAPQ10.tsv.gz" \
    "${WORK}/${SAMPLE}.idxstats.tsv" \
    "$BTLOG" \
    "${WORK}/${SAMPLE}.mapping_counts.tsv" \
    "$WORK"


[[ -s "${WORK}/${SAMPLE}.98D1B_locus_metrics.tsv" ]] || {
    echo "ERROR: locus metrics absent" >&2
    exit 35
}


ROWS=$(
    awk '
        NR > 1 {
            n++
        }
        END {
            print n+0
        }
    ' "${WORK}/${SAMPLE}.98D1B_locus_metrics.tsv"
)


[[ "$ROWS" -eq 30 ]] || {
    echo "ERROR: locus metric rows=$ROWS" >&2
    exit 36
}


# ============================================================
# 10. Provenance
# ============================================================

cat > "${WORK}/${SAMPLE}.98D1B_provenance.tsv" <<EOF
field	value
sample	${SAMPLE}
reference	${REF}
competitive_contexts	30
existing_ATTRLOC	7
additional_specialized	20
additional_exploratory	3
mapper	Bowtie2_2.5.4
bowtie2_preset	very-sensitive
bowtie2_a_or_k	NO
unaligned_records_written	NO
samtools	1.21
primary_depth_MAPQ_threshold	10
all_depth_MAPQ_threshold	0
base_quality_threshold	0
operational_robust_rule	breadth_ALL_ge0.90_AND_breadth_MAPQ10_ge0.75
depth_normalization	per_million_input_read_pairs
normalized_depth_equivalent_to_cell_abundance	NO
DNA_recruitment_equivalent_to_expression	NO
functional_bacteriocin_confirmed	NO
amplicon_data_used	NO
job_id	${SLURM_ARRAY_JOB_ID}
array_task	${SLURM_ARRAY_TASK_ID}
EOF


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${WORK}/98D1B_COMPLETE.ok"


# ============================================================
# 11. Atomic-ish finalization
# ============================================================

mv \
    "$WORK" \
    "$OUT"


echo
echo "============================================================"
echo "SAMPLE SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/${SAMPLE}.98D1B_sample_summary.tsv"


echo
echo "============================================================"
echo "98D1B SAMPLE ${SAMPLE} FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
