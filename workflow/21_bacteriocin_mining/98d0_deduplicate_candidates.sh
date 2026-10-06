#!/bin/bash

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

OUT="${ROOT}/98_bacteriocin_mining/98D0_deduplication"

C="${ROOT}/98_bacteriocin_mining/98C_FINAL_synthesis"

C2="${ROOT}/98_bacteriocin_mining/98C2_positive_context"

C1="${ROOT}/98_bacteriocin_mining/98C1R_positive_screen_corrected"

B2="${ROOT}/98_bacteriocin_mining/98B2R_strict_comparippson"

ASROOT="${ROOT}/26_antismash_rep18/results"

UNIQUE="${ROOT}/69_incomplete_unique_contigs/unique_contigs_lt90_all.fna"

MASTER="${ROOT}/38_bacteriocin_master/structural_candidates_protein.faa"

PREP="${SCRIPT_DIR}/98d0_prepare_dedup.py"

PARSE="${SCRIPT_DIR}/98d0_parse_dedup.py"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"

BLAST="${BLAST:-/opt/ohpc/pub/containers/BIO/blast-2.16.0.sif}"

THREADS="${SLURM_CPUS_PER_TASK:-4}"


rm -rf "$OUT"

mkdir -p \
    "$OUT/protein_db" \
    "$OUT/context_db"


echo "============================================================"
echo "98D0 - DEDUPLICATION OF ADDITIONAL BACTERIOCIN/RiPP CONTEXTS"
echo "START=$(date -Iseconds)"
echo "============================================================"


# ============================================================
# Inputs
# ============================================================

INPUTS=(
    "${C}/98C_FINAL_primary_protein_synthesis.tsv"
    "${C}/98C_FINAL_primary_contig_synthesis.tsv"
    "${C2}/98C2_positive_context_summary.tsv"
    "${C1}/98C1_positive_proteins.faa"
    "$MASTER"
    "$UNIQUE"
    "${B2}/98B2R_locus_strict_summary.tsv"
    "$PREP"
    "$PARSE"
    "$BLAST"
)


for F in "${INPUTS[@]}"
do

    [[ -s "$F" ]] || {
        echo "ERROR: missing input $F" >&2
        exit 20
    }

done


[[ -d "$ASROOT" ]] || {
    echo "ERROR: antiSMASH results root missing: $ASROOT" >&2
    exit 21
}


# ============================================================
# 1. Prepare panels
# ============================================================

python3 \
    "$PREP" \
    "${C}/98C_FINAL_primary_protein_synthesis.tsv" \
    "${C}/98C_FINAL_primary_contig_synthesis.tsv" \
    "${C2}/98C2_positive_context_summary.tsv" \
    "${C1}/98C1_positive_proteins.faa" \
    "$MASTER" \
    "$UNIQUE" \
    "${B2}/98B2R_locus_strict_summary.tsv" \
    "$ASROOT" \
    "$OUT"


# ============================================================
# 2. Protein all-vs-all
# ============================================================

"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$BLAST" \
    makeblastdb \
    -in "${OUT}/98D0_protein_redundancy_panel.faa" \
    -dbtype prot \
    -parse_seqids \
    -out "${OUT}/protein_db/panel" \
    > "${OUT}/98D0_makeblastdb_protein.log" \
    2>&1


"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$BLAST" \
    blastp \
    -task blastp-short \
    -query "${OUT}/98D0_protein_redundancy_panel.faa" \
    -db "${OUT}/protein_db/panel" \
    -seg no \
    -comp_based_stats 0 \
    -evalue 10 \
    -max_target_seqs 1000 \
    -max_hsps 1 \
    -num_threads "$THREADS" \
    -outfmt \
    "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
    -out "${OUT}/98D0_protein_all_vs_all.tsv"


# ============================================================
# 3. Context all-vs-all
# ============================================================

"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$BLAST" \
    makeblastdb \
    -in "${OUT}/98D0_context_redundancy_panel.fna" \
    -dbtype nucl \
    -parse_seqids \
    -out "${OUT}/context_db/panel" \
    > "${OUT}/98D0_makeblastdb_context.log" \
    2>&1


"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$BLAST" \
    blastn \
    -task megablast \
    -query "${OUT}/98D0_context_redundancy_panel.fna" \
    -db "${OUT}/context_db/panel" \
    -dust no \
    -perc_identity 80 \
    -evalue 1e-20 \
    -max_target_seqs 1000 \
    -max_hsps 1 \
    -num_threads "$THREADS" \
    -outfmt \
    "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
    -out "${OUT}/98D0_context_all_vs_all.tsv"


# ============================================================
# 4. Parse redundancy
# ============================================================

python3 \
    "$PARSE" \
    "$OUT"


# ============================================================
# 5. Scope
# ============================================================

cat > "${OUT}/98D0_methodological_scope.tsv" <<'EOF'
field	value
analysis	operational_redundancy_analysis_before_optional_read_mapping
existing_reference_proteins	8_master_candidate_sequences
existing_reference_contexts	7_ATTRLOC_contigs
additional_context_definition	no_direct_homology_link_to_existing_candidate_in_98C
protein_near_redundancy_threshold	PID_ge95_and_bidirectional_coverage_ge90
context_near_redundancy_threshold	PID_ge95_and_alignment_ge85_percent_of_shorter_sequence
exact_context_hash	orientation_insensitive
near_redundancy_group	biological_identity_NO
context_homology_group	same_MGE_or_same_strain_NO
mapping_representative	computational_device_not_biological_representative
profile_only_contexts	retained_as_exploratory
novelty_assessed	NO
functional_activity_assessed	NO
reads_remapped	NO
amplicon_data_used	NO
EOF


# ============================================================
# 6. Validate
# ============================================================

OUTPUTS=(
    "${OUT}/98D0_preparation_summary.tsv"
    "${OUT}/98D0_context_classification.tsv"
    "${OUT}/98D0_additional_candidate_proteins.tsv"
    "${OUT}/98D0_additional_candidate_CDS.fna"
    "${OUT}/98D0_additional_CDS_metadata.tsv"
    "${OUT}/98D0_protein_similarity_pairs.tsv"
    "${OUT}/98D0_protein_redundancy_groups.tsv"
    "${OUT}/98D0_CDS_exact_clusters.tsv"
    "${OUT}/98D0_context_similarity_pairs.tsv"
    "${OUT}/98D0_context_redundancy_groups.tsv"
    "${OUT}/98D0_additional_context_mapping_manifest.tsv"
    "${OUT}/98D0_global_summary.tsv"
)


for F in "${OUTPUTS[@]}"
do

    [[ -s "$F" ]] || {
        echo "ERROR: missing output $F" >&2
        exit 30
    }

done


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/98D0_COMPLETE.ok"


echo
echo "============================================================"
echo "98D0 GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/98D0_global_summary.tsv"


echo
echo "============================================================"
echo "98D0 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
