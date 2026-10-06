#!/bin/bash

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

OUT="${ROOT}/98_bacteriocin_mining/98C1_positive_screen"

CATALOG="${ROOT}/70_incomplete_gene_catalog/incomplete_lt90_new.prodigal.faa"

EGGNOG="${ROOT}/72_incomplete_eggnog_consolidated/incomplete_lt90_eggnog_annotations_698054.tsv"

MASTER="${ROOT}/38_bacteriocin_master/structural_candidates_protein.faa"

HMMROOT="${ROOT}/databases/antismash_7.1.0_package/antismash/detection/hmm_detection/data"

PY="${SCRIPT_DIR}/98c1_integrate_positive_screen.py"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"

HMMER="${HMMER:-/opt/ohpc/pub/containers/BIO/hmmer-3.4.sif}"

BLAST="${BLAST:-/opt/ohpc/pub/containers/BIO/blast-2.16.0.sif}"

THREADS="${SLURM_CPUS_PER_TASK:-16}"


rm -rf "$OUT"

mkdir -p \
    "$OUT/hmm" \
    "$OUT/blastdb"


echo "============================================================"
echo "98C1 - POSITIVE BACTERIOCIN/RiPP SCREEN"
echo "START=$(date -Iseconds)"
echo "THREADS=$THREADS"
echo "============================================================"


# ============================================================
# 1. Required inputs
# ============================================================

for F in \
    "$CATALOG" \
    "$EGGNOG" \
    "$MASTER" \
    "$PY" \
    "$HMMER" \
    "$BLAST"
do

    if [[ ! -s "$F" ]]
    then
        echo "ERROR: missing $F" >&2
        exit 20
    fi

done


# ============================================================
# 2. Build HMM databases
# ============================================================

GA_DB="${OUT}/hmm/bacteriocin_GA.hmm"
NOGA_DB="${OUT}/hmm/bacteriocin_noGA.hmm"


GA_MODELS=(
    "BacteriocIIc_cy.hmm"
    "Bacteriocin_II.hmm"
    "Bacteriocin_IIc.hmm"
    "Bacteriocin_IId.hmm"
    "Bacteriocin_IIi.hmm"
    "Lactococcin.hmm"
    "Lactococcin_972.hmm"
    "LcnG-beta.hmm"
)


NOGA_MODELS=(
    "glycocin.hmm"
    "lasso.hmm"
    "thio_amide.hmm"
    "thiostrepton.hmm"
)


: > "$GA_DB"
: > "$NOGA_DB"


for H in "${GA_MODELS[@]}"
do

    F="${HMMROOT}/${H}"

    [[ -s "$F" ]] || {
        echo "ERROR: missing GA HMM $F" >&2
        exit 21
    }

    cat "$F" >> "$GA_DB"

done


for H in "${NOGA_MODELS[@]}"
do

    F="${HMMROOT}/${H}"

    [[ -s "$F" ]] || {
        echo "ERROR: missing noGA HMM $F" >&2
        exit 22
    }

    cat "$F" >> "$NOGA_DB"

done


"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$HMMER" \
    hmmpress \
    -f \
    "$GA_DB"


"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$HMMER" \
    hmmpress \
    -f \
    "$NOGA_DB"


# ============================================================
# 3. Curated GA HMM scan
# ============================================================

echo
echo "[98C1] GA-profile hmmscan: $(date -Iseconds)"


"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$HMMER" \
    hmmscan \
    --cpu "$THREADS" \
    --cut_ga \
    --tblout "${OUT}/98C1_GA_hmmscan.tbl" \
    --domtblout "${OUT}/98C1_GA_hmmscan.domtbl" \
    "$GA_DB" \
    "$CATALOG" \
    > "${OUT}/98C1_GA_hmmscan.stdout"


# ============================================================
# 4. Supplementary HMM scan without curated GA
#
# These are screening hits only.
# ============================================================

echo
echo "[98C1] supplementary hmmscan: $(date -Iseconds)"


"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$HMMER" \
    hmmscan \
    --cpu "$THREADS" \
    -E 1e-5 \
    --domE 1e-5 \
    --tblout "${OUT}/98C1_noGA_hmmscan.tbl" \
    --domtblout "${OUT}/98C1_noGA_hmmscan.domtbl" \
    "$NOGA_DB" \
    "$CATALOG" \
    > "${OUT}/98C1_noGA_hmmscan.stdout"


# ============================================================
# 5. BLASTP-short against 8 unique master proteins
# ============================================================

echo
echo "[98C1] makeblastdb: $(date -Iseconds)"


"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$BLAST" \
    makeblastdb \
    -in "$CATALOG" \
    -dbtype prot \
    -parse_seqids \
    -out "${OUT}/blastdb/incomplete_catalog" \
    > "${OUT}/98C1_makeblastdb.log" \
    2>&1


echo
echo "[98C1] blastp-short: $(date -Iseconds)"


"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$BLAST" \
    blastp \
    -task blastp-short \
    -query "$MASTER" \
    -db "${OUT}/blastdb/incomplete_catalog" \
    -seg no \
    -comp_based_stats 0 \
    -evalue 10 \
    -max_target_seqs 100000 \
    -max_hsps 1 \
    -num_threads "$THREADS" \
    -outfmt \
    "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
    -out "${OUT}/98C1_blastp_short_raw.tsv"


# ============================================================
# 6. Integrate
# ============================================================

echo
echo "[98C1] integration: $(date -Iseconds)"


python3 \
    "$PY" \
    "$CATALOG" \
    "$EGGNOG" \
    "${OUT}/98C1_GA_hmmscan.tbl" \
    "${OUT}/98C1_noGA_hmmscan.tbl" \
    "${OUT}/98C1_blastp_short_raw.tsv" \
    "$MASTER" \
    "$OUT"


# ============================================================
# 7. Scope
# ============================================================

cat > "${OUT}/98C1_methodological_scope.tsv" <<'EOF'
field	value
analysis	positive_screen_of_protein_catalog_from_incomplete_bins
catalog	incomplete_lt90_new.prodigal.faa
GA_HMM_models	8
GA_HMM_threshold	model_specific_GA_cutoff
supplementary_HMM_models	4
supplementary_HMM_threshold	full_and_domain_Evalue_1e-5
short_protein_similarity	BLASTP_short
BLAST_query_set	8_unique_existing_structural_candidate_sequences
BLAST_raw_evalue_cutoff	10_for_sensitive_short_peptide_search
promoted_exact_full_length	100pct_identity_100pct_query_coverage_same_length
promoted_exact_contained	100pct_identity_100pct_query_coverage
promoted_near_exact	identity_ge90_query_coverage_ge90_bitscore_ge25
promoted_strong_homology	identity_ge70_query_coverage_ge80_bitscore_ge25_evalue_le1e-2
eggNOG	role_as_auxiliary_annotation_only_after_sequence_screen
HMM_or_BLAST_positive	equivalent_to_validated_bacteriocin	NO
incomplete_bin_negative_result	equivalent_to_biological_absence	NO
expression	NOT_ASSESSED
activity	NOT_ASSESSED
reads_remapped	NO
amplicon_data_used	NO
EOF


# ============================================================
# 8. Validate
# ============================================================

for F in \
    "${OUT}/98C1_GA_hmmscan.tbl" \
    "${OUT}/98C1_noGA_hmmscan.tbl" \
    "${OUT}/98C1_blastp_short_raw.tsv" \
    "${OUT}/98C1_positive_proteins.tsv" \
    "${OUT}/98C1_positive_proteins.faa" \
    "${OUT}/98C1_HMM_hits.tsv" \
    "${OUT}/98C1_promoted_BLAST_hits.tsv" \
    "${OUT}/98C1_global_summary.tsv" \
    "${OUT}/98C1_methodological_scope.tsv"
do

    [[ -e "$F" ]] || {
        echo "ERROR: missing output $F" >&2
        exit 30
    }

done


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/98C1_COMPLETE.ok"


echo
echo "============================================================"
echo "98C1 GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/98C1_global_summary.tsv"


echo
echo "============================================================"
echo "98C1 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
