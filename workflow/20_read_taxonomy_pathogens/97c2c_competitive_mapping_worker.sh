#!/bin/bash

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

TASK_ID="${1:?TASK_ID required}"

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

REFROOT="${ROOT}/98_read_taxonomy/97C2B_reference_panels"
OUTROOT="${ROOT}/98_read_taxonomy/97C2C_competitive_mapping"

MANIFEST="${OUTROOT}/97C2C_manifest.tsv"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"
BOWTIE2_SIF="${BOWTIE2_SIF:-/opt/ohpc/pub/containers/BIO/bowtie2-2.5.4.sif}"
SAMTOOLS_SIF="${SAMTOOLS_SIF:-/opt/ohpc/pub/containers/BIO/samtools-1.21.sif}"

DEPTH_PARSER="${SCRIPT_DIR}/97c2c_depth_metrics.py"
SPECIES_PARSER="${SCRIPT_DIR}/97c2c_species_counts.py"
SUMMARY_PARSER="${SCRIPT_DIR}/97c2c_build_task_summary.py"

CPUS="${SLURM_CPUS_PER_TASK:-6}"

mkdir -p "$OUTROOT"


if [[ ! -s "$MANIFEST" ]]
then

    cat > "$MANIFEST" <<'EOF'
task_id	sample	panel	target_slug	target_species	priority
0	M2_3	salmonella	Salmonella_enterica	Salmonella enterica	VERY_HIGH
1	M3_3	salmonella	Salmonella_enterica	Salmonella enterica	HIGH
2	L2_3	listeria	Listeria_monocytogenes	Listeria monocytogenes	VERY_HIGH
3	L3_1	listeria	Listeria_monocytogenes	Listeria monocytogenes	HIGH
4	M1_1	staph_aureus	Staphylococcus_aureus	Staphylococcus aureus	VERY_HIGH
5	M2_1	ecoli	Escherichia_coli	Escherichia coli	HIGH
6	M2_3	yersinia	Yersinia_enterocolitica	Yersinia enterocolitica	VERY_HIGH
7	L2_2	cronobacter	Cronobacter_sakazakii	Cronobacter sakazakii	HIGH
8	L2_3	bcereus	Bacillus_cereus	Bacillus cereus	MEDIUM
9	L1_3	cperfringens	Clostridium_perfringens	Clostridium perfringens	MEDIUM
10	L2_2	campylobacter	Campylobacter_jejuni	Campylobacter jejuni	LOW
EOF

fi


ROW=$(
    awk -F'\t' \
        -v id="$TASK_ID" \
        'NR > 1 && $1 == id {print; exit}' \
        "$MANIFEST"
)

ROW="${ROW//$'\r'/}"


if [[ -z "$ROW" ]]
then
    echo "ERROR: task not found: $TASK_ID" >&2
    exit 10
fi


IFS=$'\t' read -r \
    TASK \
    SAMPLE \
    PANEL \
    TARGET_SLUG \
    TARGET_SPECIES \
    PRIORITY \
<<< "$ROW"


TASKDIR="${OUTROOT}/${TASK}_${SAMPLE}_${PANEL}"

rm -rf "$TASKDIR"

mkdir -p "$TASKDIR"


R1="${ROOT}/04_host_removed/${SAMPLE}/${SAMPLE}_1.hostfree.fastq.gz"
R2="${ROOT}/04_host_removed/${SAMPLE}/${SAMPLE}_2.hostfree.fastq.gz"

INDEX="${REFROOT}/bowtie2_indexes/${PANEL}"
PANEL_FASTA="${REFROOT}/panels/${PANEL}.competitive.fna"

CLASSIFICATION="${ROOT}/98_read_taxonomy/97B1_all_samples/${SAMPLE}/${SAMPLE}.classification_summary.tsv"

BAM="${TASKDIR}/${SAMPLE}.${PANEL}.competitive.sorted.bam"
BT2LOG="${TASKDIR}/${SAMPLE}.${PANEL}.bowtie2.log"


for F in \
    "$R1" \
    "$R2" \
    "$PANEL_FASTA" \
    "$CLASSIFICATION" \
    "$DEPTH_PARSER" \
    "$SPECIES_PARSER" \
    "$SUMMARY_PARSER"

do

    if [[ ! -s "$F" ]]
    then
        echo "ERROR: missing input $F" >&2
        exit 11
    fi

done


if [[ ! -s "${INDEX}.1.bt2" && ! -s "${INDEX}.1.bt2l" ]]
then
    echo "ERROR: Bowtie2 index absent: $INDEX" >&2
    exit 12
fi


echo "============================================================"
echo "97C2C COMPETITIVE MAPPING"
echo "task=$TASK"
echo "sample=$SAMPLE"
echo "panel=$PANEL"
echo "target=$TARGET_SPECIES"
echo "priority=$PRIORITY"
echo "start=$(date -Iseconds)"
echo "============================================================"


# ============================================================
# Competitive mapping
#
# Default Bowtie2: one best alignment reported.
# MAPQ carries ambiguity information relative to the panel.
# ============================================================

"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$BOWTIE2_SIF" \
    bowtie2 \
    --very-sensitive \
    --no-unal \
    -p 4 \
    -x "$INDEX" \
    -1 "$R1" \
    -2 "$R2" \
    2> "$BT2LOG" \
| "$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$SAMTOOLS_SIF" \
    samtools view \
    -b \
    -F 4 \
    - \
| "$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$SAMTOOLS_SIF" \
    samtools sort \
    -@ 2 \
    -o "$BAM" \
    -


"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$SAMTOOLS_SIF" \
    samtools quickcheck \
    -v \
    "$BAM"


"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$SAMTOOLS_SIF" \
    samtools index \
    -@ 2 \
    "$BAM"


# ============================================================
# Resolve target chromosome
# ============================================================

TARGET_INFO=$(
    "$APPTAINER" exec \
        --bind "${ROOT}:${ROOT}" \
        "$SAMTOOLS_SIF" \
        samtools idxstats \
        "$BAM" \
    | awk -F'\t' \
        -v prefix="${TARGET_SLUG}|" \
        'index($1,prefix)==1 {
            print $1 "\t" $2
        }'
)


TARGET_N=$(
    printf '%s\n' "$TARGET_INFO" \
    | awk 'NF {n++} END {print n+0}'
)


if [[ "$TARGET_N" -ne 1 ]]
then
    echo "ERROR: expected exactly one target chromosome, found $TARGET_N" >&2
    printf '%s\n' "$TARGET_INFO" >&2
    exit 20
fi


TARGET_REF=$(
    printf '%s\n' "$TARGET_INFO" \
    | cut -f1
)

TARGET_LENGTH=$(
    printf '%s\n' "$TARGET_INFO" \
    | cut -f2
)


echo "target_ref=$TARGET_REF"
echo "target_length=$TARGET_LENGTH"


# ============================================================
# Mapping counts
# ============================================================

sam_count () {

    "$APPTAINER" exec \
        --bind "${ROOT}:${ROOT}" \
        "$SAMTOOLS_SIF" \
        samtools view \
        "$@"
}


PANEL_ALL=$(
    sam_count \
        -c \
        -F 2304 \
        "$BAM"
)

PANEL_Q10=$(
    sam_count \
        -c \
        -F 2304 \
        -q 10 \
        "$BAM"
)

PANEL_Q20=$(
    sam_count \
        -c \
        -F 2304 \
        -q 20 \
        "$BAM"
)


TARGET_ALL=$(
    sam_count \
        -c \
        -F 2304 \
        "$BAM" \
        "$TARGET_REF"
)

TARGET_Q10=$(
    sam_count \
        -c \
        -F 2304 \
        -q 10 \
        "$BAM" \
        "$TARGET_REF"
)

TARGET_Q20=$(
    sam_count \
        -c \
        -F 2304 \
        -q 20 \
        "$BAM" \
        "$TARGET_REF"
)


TARGET_PP_ALL=$(
    sam_count \
        -c \
        -F 2304 \
        -f 2 \
        "$BAM" \
        "$TARGET_REF"
)

TARGET_PP_Q10=$(
    sam_count \
        -c \
        -F 2304 \
        -f 2 \
        -q 10 \
        "$BAM" \
        "$TARGET_REF"
)

TARGET_PP_Q20=$(
    sam_count \
        -c \
        -F 2304 \
        -f 2 \
        -q 20 \
        "$BAM" \
        "$TARGET_REF"
)


{
    printf "metric\tvalue\n"

    printf "panel_mapped_records_all\t%s\n" "$PANEL_ALL"
    printf "panel_mapped_records_q10\t%s\n" "$PANEL_Q10"
    printf "panel_mapped_records_q20\t%s\n" "$PANEL_Q20"

    printf "target_mapped_records_all\t%s\n" "$TARGET_ALL"
    printf "target_mapped_records_q10\t%s\n" "$TARGET_Q10"
    printf "target_mapped_records_q20\t%s\n" "$TARGET_Q20"

    printf "target_proper_pair_records_all\t%s\n" "$TARGET_PP_ALL"
    printf "target_proper_pair_records_q10\t%s\n" "$TARGET_PP_Q10"
    printf "target_proper_pair_records_q20\t%s\n" "$TARGET_PP_Q20"

} > "${TASKDIR}/mapping_counts.tsv"


# ============================================================
# Species-level competitive counts
# ============================================================

for Q in 0 10 20
do

    sam_count \
        -F 2304 \
        -q "$Q" \
        "$BAM" \
    | python3 \
        "$SPECIES_PARSER" \
        "$Q" \
        "${TASKDIR}/species_counts_q${Q}.tsv"

done


# ============================================================
# Target coverage
#
# -q = minimum base quality
# -Q = minimum mapping quality
# ============================================================

for Q in 0 10 20
do

    "$APPTAINER" exec \
        --bind "${ROOT}:${ROOT}" \
        "$SAMTOOLS_SIF" \
        samtools depth \
        -aa \
        -q 0 \
        -Q "$Q" \
        -r "$TARGET_REF" \
        "$BAM" \
    | python3 \
        "$DEPTH_PARSER" \
        "$SAMPLE" \
        "$PANEL" \
        "$TARGET_SPECIES" \
        "$TARGET_REF" \
        "$TARGET_LENGTH" \
        "$Q" \
        10000 \
        "${TASKDIR}/depth_summary_q${Q}.tsv" \
        "${TASKDIR}/depth_windows_q${Q}.tsv"

done


# ============================================================
# Final task summary
# ============================================================

python3 \
    "$SUMMARY_PARSER" \
    "$TASK" \
    "$SAMPLE" \
    "$PANEL" \
    "$TARGET_SLUG" \
    "$TARGET_SPECIES" \
    "$PRIORITY" \
    "$CLASSIFICATION" \
    "${TASKDIR}/mapping_counts.tsv" \
    "${TASKDIR}/depth_summary_q0.tsv" \
    "${TASKDIR}/depth_summary_q10.tsv" \
    "${TASKDIR}/depth_summary_q20.tsv" \
    "${TASKDIR}/species_counts_q10.tsv" \
    "${TASKDIR}/task_summary.tsv"


# ============================================================
# Validation
# ============================================================

for F in \
    "$BAM" \
    "${BAM}.bai" \
    "$BT2LOG" \
    "${TASKDIR}/mapping_counts.tsv" \
    "${TASKDIR}/species_counts_q0.tsv" \
    "${TASKDIR}/species_counts_q10.tsv" \
    "${TASKDIR}/species_counts_q20.tsv" \
    "${TASKDIR}/depth_summary_q0.tsv" \
    "${TASKDIR}/depth_summary_q10.tsv" \
    "${TASKDIR}/depth_summary_q20.tsv" \
    "${TASKDIR}/task_summary.tsv"

do

    if [[ ! -s "$F" ]]
    then
        echo "ERROR: missing output $F" >&2
        exit 30
    fi

done


printf \
"COMPLETED\t%s\ttask=%s\tsample=%s\tpanel=%s\n" \
"$(date -Iseconds)" \
"$TASK" \
"$SAMPLE" \
"$PANEL" \
> "${TASKDIR}/97C2C_COMPLETE.ok"


echo
echo "============================================================"
echo "TASK SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${TASKDIR}/task_summary.tsv"

echo
echo "============================================================"
echo "97C2C TASK $TASK COMPLETE"
echo "ES SEGURO SALIR"
echo "============================================================"
