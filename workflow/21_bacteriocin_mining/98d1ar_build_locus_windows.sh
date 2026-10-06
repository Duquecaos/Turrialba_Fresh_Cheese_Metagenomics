#!/bin/bash

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

D0="${ROOT}/98_bacteriocin_mining/98D0_deduplication"

C2="${ROOT}/98_bacteriocin_mining/98C2_positive_context"

OLDREF="${ROOT}/47_bacteriocin_locus_reference/all_match_locus_contexts_5kb.fna"

CONTIGS="${ROOT}/69_incomplete_unique_contigs/unique_contigs_lt90_all.fna"

OUT="${ROOT}/98_bacteriocin_mining/98D1_competitive_mapping_v2"

REFDIR="${OUT}/reference"

INDEXDIR="${REFDIR}/bowtie2_index"

PY="${SCRIPT_DIR}/98d1ar_build_locus_windows.py"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"

BOWTIE2="${BOWTIE2:-/opt/ohpc/pub/containers/BIO/bowtie2-2.5.4.sif}"

THREADS="${SLURM_CPUS_PER_TASK:-4}"


rm -rf "$OUT"

mkdir -p \
    "$REFDIR" \
    "$INDEXDIR"


echo "============================================================"
echo "98D1AR - CORRECTED COMPETITIVE LOCUS-WINDOW REFERENCE"
echo "START=$(date -Iseconds)"
echo "============================================================"


# ============================================================
# Inputs
# ============================================================

for F in \
    "${D0}/98D0_additional_context_mapping_manifest.tsv" \
    "${C2}/98C2_positive_context_summary.tsv" \
    "$CONTIGS" \
    "$OLDREF" \
    "$PY" \
    "$BOWTIE2"
do

    [[ -s "$F" ]] || {
        echo "ERROR: missing $F" >&2
        exit 20
    }

done


# ============================================================
# Build corrected reference
# ============================================================

python3 \
    "$PY" \
    "${D0}/98D0_additional_context_mapping_manifest.tsv" \
    "${C2}/98C2_positive_context_summary.tsv" \
    "$CONTIGS" \
    "$OLDREF" \
    "$REFDIR"


REF="${REFDIR}/98D1AR_competitive_locus_windows_30.fna"

META="${REFDIR}/98D1AR_competitive_locus_windows_30_metadata.tsv"

SAMPLES="${REFDIR}/98D1AR_sample_manifest.tsv"


NREF=$(
    grep -c '^>' "$REF"
)


NSAMP=$(
    awk -F $'\t' '
        NR > 1 {n++}
        END {print n+0}
    ' "$SAMPLES"
)


[[ "$NREF" -eq 30 ]] || {
    echo "ERROR: reference count=$NREF" >&2
    exit 21
}


[[ "$NSAMP" -eq 18 ]] || {
    echo "ERROR: sample count=$NSAMP" >&2
    exit 22
}


echo "REFERENCE_SEQUENCES=$NREF"
echo "SAMPLES=$NSAMP"


# ============================================================
# FASTQ preflight
# ============================================================

MISSING=0


while IFS=$'\t' read -r IDX SAMPLE PRODUCER UNIT TIMECODE
do

    [[ "$IDX" == "array_index" ]] && continue

    R1="${ROOT}/04_host_removed/${SAMPLE}/${SAMPLE}_1.hostfree.fastq.gz"
    R2="${ROOT}/04_host_removed/${SAMPLE}/${SAMPLE}_2.hostfree.fastq.gz"

    for F in "$R1" "$R2"
    do

        if [[ ! -s "$F" ]]
        then
            echo "MISSING_FASTQ=$F" >&2
            MISSING=$((MISSING + 1))
        fi

    done

done < "$SAMPLES"


[[ "$MISSING" -eq 0 ]] || {
    echo "ERROR: missing FASTQ=$MISSING" >&2
    exit 23
}


echo "HOSTFREE_FASTQ_PREFLIGHT=36/36_PASS"


# ============================================================
# Build Bowtie2 index
# ============================================================

INDEX="${INDEXDIR}/bacteriocin_locuswindows30"


"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$BOWTIE2" \
    bowtie2-build \
    --threads "$THREADS" \
    "$REF" \
    "$INDEX" \
    > "${REFDIR}/98D1AR_bowtie2_build.stdout" \
    2> "${REFDIR}/98D1AR_bowtie2_build.stderr"


NINDEX=$(
    find "$INDEXDIR" \
        -maxdepth 1 \
        -type f \
        \( \
            -name 'bacteriocin_locuswindows30*.bt2' \
            -o \
            -name 'bacteriocin_locuswindows30*.bt2l' \
        \) \
    | wc -l
)


echo "BOWTIE2_INDEX_FILES=$NINDEX"


[[ "$NINDEX" -eq 6 ]] || {
    echo "ERROR: Bowtie2 index files=$NINDEX" >&2
    exit 24
}


# ============================================================
# Checksums / scope
# ============================================================

sha256sum \
    "$REF" \
    "$META" \
    "$SAMPLES" \
> "${REFDIR}/98D1AR_reference_checksums.sha256"


cat > "${REFDIR}/98D1AR_methodological_scope.tsv" <<'EOF'
field	value
analysis	corrected_competitive_bacteriocin_RiPP_reference
known_context_reference	original_7_ATTRLOC_local_contexts_from_step47
additional_context_reference	candidate_gene_span_plus_5000bp_each_side_clipped_to_contig
additional_deduplication	representatives_frozen_from_98D0
known_ATTRLOC_count	7
additional_specialized_count	20
additional_exploratory_count	3
total_reference_count	30
full_long_additional_contigs_used_directly	NO
reason_for_windowing	reduce_context_length_imbalance_and_focus_recruitment_around_candidate_locus
ATTRLOC002_and_ATTRLOC006_retained	YES
shared_sequence_handling	competitive_mapping_plus_MAPQ
job_133413_reference	technically_valid_but_superseded_before_read_mapping
reads_mapped	NO
novelty_inferred	NO
functional_activity_inferred	NO
amplicon_data_used	NO
EOF


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${REFDIR}/98D1AR_COMPLETE.ok"


echo
echo "============================================================"
echo "98D1AR SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${REFDIR}/98D1AR_reference_summary.tsv"


echo
echo "============================================================"
echo "98D1AR FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
