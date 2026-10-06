#!/bin/bash

set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"
OUT="${ROOT}/98_bacteriocin_mining/98A_audit"

mkdir -p "$OUT"

STAGE="initialization"

trap '
RC=$?
echo "98A_ERROR_STAGE=${STAGE}" >&2
echo "98A_ERROR_RC=${RC}" >&2
exit "$RC"
' ERR


echo "============================================================"
echo "98A - BACTERIOCIN MINING AUDIT"
echo "START=$(date -Iseconds)"
echo "ROOT=$ROOT"
echo "============================================================"


# ============================================================
# 1. Known tool resources
# ============================================================

STAGE="known_tools"

cat > "${OUT}/98A_known_tool_resources.tsv" <<'EOF'
tool	path	expected_use
Apptainer	/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer	container_runtime
BLAST	/opt/ohpc/pub/containers/BIO/blast-2.16.0.sif	sequence_similarity
DIAMOND	/opt/ohpc/pub/containers/BIO/diamond-2.1.9.sif	protein_similarity
HMMER	/opt/ohpc/pub/containers/BIO/hmmer-3.4.sif	profile_HMM_search
Prodigal	/opt/ohpc/pub/containers/BIO/prodigal-2.6.3.sif	protein_prediction
eggNOG	/opt/ohpc/pub/containers/BIO/eggnog-mapper-2.1.12.sif	functional_annotation
EOF


printf \
"tool\tpath\texpected_use\tstatus\n" \
> "${OUT}/98A_known_tool_status.tsv"


tail -n +2 \
"${OUT}/98A_known_tool_resources.tsv" \
| while IFS=$'\t' read -r TOOL PATH USE
do

    if [[ -e "$PATH" ]]
    then
        STATUS="FOUND"
    else
        STATUS="NOT_FOUND"
    fi

    printf \
"%s\t%s\t%s\t%s\n" \
"$TOOL" \
"$PATH" \
"$USE" \
"$STATUS" \
>> "${OUT}/98A_known_tool_status.tsv"

done


# ============================================================
# 2. Specialized containers
#
# All discovery searches are NON-FATAL.
# ============================================================

STAGE="container_audit"

: > "${OUT}/98A_all_container_files.txt"


for BASE in \
    "${OHPC_PUBLIC_ROOT:-/opt/ohpc/pub}/containers/BIO" \
    "${OHPC_PUBLIC_ROOT:-/opt/ohpc/pub}/containers"
do

    if [[ -d "$BASE" ]]
    then

        {
            find "$BASE" \
                -maxdepth 4 \
                -type f \
                -print \
                2>/dev/null \
            || true
        } \
        >> "${OUT}/98A_all_container_files.txt"

    fi

done


sort -u \
"${OUT}/98A_all_container_files.txt" \
-o "${OUT}/98A_all_container_files.txt"


grep -Ei \
'antismash|bagel|bacterioc|ripp|lanthi|lasso|sacti|thiopep|hmmer|blast|diamond|prodigal|interpro|pfam' \
"${OUT}/98A_all_container_files.txt" \
> "${OUT}/98A_relevant_containers.txt" \
|| true


# ============================================================
# 3. Local databases
# ============================================================

STAGE="database_audit"

: > "${OUT}/98A_relevant_database_paths.txt"


for BASE in \
    "${ROOT}/databases" \
    "${OHPC_PUBLIC_ROOT:-/opt/ohpc/pub}/databases" \
    "${OHPC_PUBLIC_ROOT:-/opt/ohpc/pub}/data"
do

    [[ -d "$BASE" ]] || continue


    {
        find "$BASE" \
            -maxdepth 7 \
            \( -type f -o -type d \) \
            -print \
            2>/dev/null \
        || true
    } \
    | grep -Ei \
        'bacterioc|bagel|mibig|ripp|lant|nisin|pediocin|enterocin|lactococcin|plantaricin|sakacin|pfam|tigrfam' \
    >> "${OUT}/98A_relevant_database_paths.txt" \
    || true

done


sort -u \
"${OUT}/98A_relevant_database_paths.txt" \
-o "${OUT}/98A_relevant_database_paths.txt"


# ============================================================
# 4. Project top-level inventory
# ============================================================

STAGE="project_inventory"

{
    find "$ROOT" \
        -maxdepth 1 \
        -mindepth 1 \
        -type d \
        -printf '%f\n' \
        2>/dev/null \
    || true
} \
| sort \
> "${OUT}/98A_project_directories.txt"


# ============================================================
# 5. Existing bacteriocin evidence
# ============================================================

STAGE="existing_bacteriocin_search"

{
    find "$ROOT" \
        -maxdepth 6 \
        \( -type f -o -type d \) \
        -print \
        2>/dev/null \
    || true
} \
| grep -Ei \
    'bacteriocin|ATTRLOC|lactococcin|lantipeptide|lanthipeptide|pediocin|enterocin|nisin|plantaricin|sakacin' \
| sort -u \
> "${OUT}/98A_existing_bacteriocin_paths.txt" \
|| true


# ============================================================
# 6. Existing ATTRLOC workflow
# ============================================================

STAGE="ATTRLOC_snapshot"

ATTRDIR="${ROOT}/49_bacteriocin_locus_matrices"

SNAP="${OUT}/existing_ATTRLOC_snapshot"

mkdir -p "$SNAP"

ATTRLOC_DIR_PRESENT="NO"
ATTRLOC_COUNT="0"


if [[ -d "$ATTRDIR" ]]
then

    ATTRLOC_DIR_PRESENT="YES"


    for NAME in \
        detection_class_summary.tsv \
        locus_summary_by_ATTRLOC.tsv \
        locus_summary_by_sample.tsv \
        discordant_high_coverage_contexts.tsv
    do

        SRC="${ATTRDIR}/${NAME}"

        if [[ -s "$SRC" ]]
        then

            cp -p \
                "$SRC" \
                "${SNAP}/${NAME}"

        fi

    done


    LOCUS_FILE="${ATTRDIR}/locus_summary_by_ATTRLOC.tsv"

    if [[ -s "$LOCUS_FILE" ]]
    then

        ATTRLOC_COUNT=$(
            awk '
                NR > 1 && NF {
                    n++
                }
                END {
                    print n+0
                }
            ' "$LOCUS_FILE"
        )

    fi

fi


# ============================================================
# 7. eggNOG outputs
# ============================================================

STAGE="eggnog_search"

{
    find "$ROOT" \
        -maxdepth 8 \
        -type f \
        \( \
            -iname '*.emapper.annotations' \
            -o -iname '*eggnog*.tsv' \
            -o -iname '*emapper*.tsv' \
        \) \
        -print \
        2>/dev/null \
    || true
} \
| sort -u \
> "${OUT}/98A_eggnog_annotation_files.txt"


# ============================================================
# 8. Protein FASTA files
# ============================================================

STAGE="protein_fasta_search"

{
    find "$ROOT" \
        -maxdepth 8 \
        -type f \
        \( \
            -iname '*.faa' \
            -o -iname '*proteins*.fa' \
            -o -iname '*proteins*.fasta' \
            -o -iname '*protein*.faa' \
        \) \
        -print \
        2>/dev/null \
    || true
} \
| sort -u \
> "${OUT}/98A_protein_fasta_files.txt"


# ============================================================
# 9. Genomic-context files
# ============================================================

STAGE="context_file_search"

{
    find "$ROOT" \
        -maxdepth 8 \
        -type f \
        \( \
            -iname '*.gff' \
            -o -iname '*.gff3' \
            -o -iname '*genes*.tsv' \
            -o -iname '*orfs*.tsv' \
        \) \
        -print \
        2>/dev/null \
    || true
} \
| sort -u \
> "${OUT}/98A_gene_context_files.txt"


# ============================================================
# 10. Candidate files from steps 86-90
# ============================================================

STAGE="steps86_90_search"

printf \
"path\tbytes\n" \
> "${OUT}/98A_steps86_90_candidate_files.tsv"


{
    find "$ROOT" \
        -maxdepth 8 \
        -type f \
        \( \
            -iname '*.faa' \
            -o -iname '*.fna' \
            -o -iname '*.fa' \
            -o -iname '*.fasta' \
            -o -iname '*.gff' \
            -o -iname '*.gff3' \
            -o -iname '*.emapper.annotations' \
            -o -iname '*eggnog*.tsv' \
            -o -iname '*annotation*.tsv' \
            -o -iname '*protein*.tsv' \
            -o -iname '*gene*.tsv' \
        \) \
        -print \
        2>/dev/null \
    || true
} \
| grep -E \
    '/(86_|87_|88_|89_|90_)' \
| sort -u \
> "${OUT}/.98A_steps86_90_paths.tmp" \
|| true


while IFS= read -r F
do

    [[ -n "$F" ]] || continue
    [[ -e "$F" ]] || continue

    SIZE=$(
        stat -c '%s' "$F" \
        2>/dev/null \
        || echo 0
    )

    printf \
"%s\t%s\n" \
"$F" \
"$SIZE" \
>> "${OUT}/98A_steps86_90_candidate_files.tsv"

done < "${OUT}/.98A_steps86_90_paths.tmp"


rm -f \
"${OUT}/.98A_steps86_90_paths.tmp"


# ============================================================
# 11. Existing specialized outputs
# ============================================================

STAGE="specialized_output_search"

{
    find "$ROOT" \
        -maxdepth 8 \
        \( -type f -o -type d \) \
        -print \
        2>/dev/null \
    || true
} \
| grep -Ei \
    'antismash|bagel|ripp|lanthipeptide|lantipeptide|lassopeptide|sactipeptide' \
| sort -u \
> "${OUT}/98A_existing_specialized_tool_outputs.txt" \
|| true


# ============================================================
# 12. Build list of annotation-like files for keyword search
# ============================================================

STAGE="annotation_inventory"

{
    find "$ROOT" \
        -maxdepth 8 \
        -type f \
        \( \
            -iname '*.tsv' \
            -o -iname '*.csv' \
            -o -iname '*.txt' \
            -o -iname '*.gff' \
            -o -iname '*.gff3' \
            -o -iname '*.annotations' \
        \) \
        -print \
        2>/dev/null \
    || true
} \
| grep -E \
    '/(49_bacteriocin|86_|87_|88_|89_|90_)' \
| sort -u \
> "${OUT}/98A_annotation_candidate_files.txt" \
|| true


# ============================================================
# 13. Keyword search
# ============================================================

STAGE="keyword_search"

KEYWORDS='bacteriocin|bacteriocidal|lantipeptide|lanthipeptide|nisin|pediocin|lactococcin|enterocin|plantaricin|sakacin|carnobacteriocin|ribosomally synthesized|immunity protein'


printf \
"path\thit_lines\n" \
> "${OUT}/98A_annotation_keyword_file_hits.tsv"


printf \
"source\tline_number\tannotation_text\n" \
> "${OUT}/98A_annotation_keyword_examples.tsv"


while IFS= read -r F
do

    [[ -s "$F" ]] || continue


    SIZE=$(
        stat -c '%s' "$F" \
        2>/dev/null \
        || echo 0
    )


    # Do not grep annotation files >512 MiB during audit.
    if [[ "$SIZE" -gt 536870912 ]]
    then
        continue
    fi


    N=$(
        grep -Eic \
            "$KEYWORDS" \
            "$F" \
            2>/dev/null \
        || true
    )


    if [[ "$N" -gt 0 ]]
    then

        printf \
"%s\t%s\n" \
"$F" \
"$N" \
>> "${OUT}/98A_annotation_keyword_file_hits.tsv"


        grep -Ein \
            "$KEYWORDS" \
            "$F" \
            2>/dev/null \
        | head -n 100 \
        | awk -v src="$F" '
            {
                pos=index($0,":")
                if (pos > 0) {
                    ln=substr($0,1,pos-1)
                    txt=substr($0,pos+1)

                    gsub(/\t/," ",txt)

                    print src "\t" ln "\t" txt
                }
            }
        ' \
        >> "${OUT}/98A_annotation_keyword_examples.tsv" \
        || true

    fi

done < "${OUT}/98A_annotation_candidate_files.txt"


# ============================================================
# 14. Counts
# ============================================================

STAGE="summary_counts"

count_nonempty_lines () {

    local F="$1"

    awk '
        NF {
            n++
        }
        END {
            print n+0
        }
    ' "$F"
}


count_data_rows () {

    local F="$1"

    awk '
        NR > 1 && NF {
            n++
        }
        END {
            print n+0
        }
    ' "$F"
}


RELEVANT_CONTAINERS=$(
    count_nonempty_lines \
        "${OUT}/98A_relevant_containers.txt"
)

RELEVANT_DATABASE_PATHS=$(
    count_nonempty_lines \
        "${OUT}/98A_relevant_database_paths.txt"
)

EXISTING_BACT_PATHS=$(
    count_nonempty_lines \
        "${OUT}/98A_existing_bacteriocin_paths.txt"
)

EGGNOG_FILES=$(
    count_nonempty_lines \
        "${OUT}/98A_eggnog_annotation_files.txt"
)

PROTEIN_FASTAS=$(
    count_nonempty_lines \
        "${OUT}/98A_protein_fasta_files.txt"
)

CONTEXT_FILES=$(
    count_nonempty_lines \
        "${OUT}/98A_gene_context_files.txt"
)

SPECIALIZED_OUTPUTS=$(
    count_nonempty_lines \
        "${OUT}/98A_existing_specialized_tool_outputs.txt"
)

ANNOT_FILES_WITH_HITS=$(
    count_data_rows \
        "${OUT}/98A_annotation_keyword_file_hits.tsv"
)

ANNOT_EXAMPLE_HITS=$(
    count_data_rows \
        "${OUT}/98A_annotation_keyword_examples.tsv"
)


# ============================================================
# 15. Global summary
# ============================================================

STAGE="global_summary"

{
    printf "metric\tvalue\n"

    printf "existing_ATTRLOC_directory\t%s\n" \
        "$ATTRLOC_DIR_PRESENT"

    printf "existing_ATTRLOC_count\t%s\n" \
        "$ATTRLOC_COUNT"

    printf "relevant_container_paths\t%s\n" \
        "$RELEVANT_CONTAINERS"

    printf "relevant_database_paths\t%s\n" \
        "$RELEVANT_DATABASE_PATHS"

    printf "existing_bacteriocin_related_paths\t%s\n" \
        "$EXISTING_BACT_PATHS"

    printf "eggnog_annotation_files\t%s\n" \
        "$EGGNOG_FILES"

    printf "protein_fasta_files\t%s\n" \
        "$PROTEIN_FASTAS"

    printf "gene_context_files\t%s\n" \
        "$CONTEXT_FILES"

    printf "existing_specialized_tool_outputs\t%s\n" \
        "$SPECIALIZED_OUTPUTS"

    printf "annotation_files_with_bacteriocin_keywords\t%s\n" \
        "$ANNOT_FILES_WITH_HITS"

    printf "annotation_keyword_example_lines\t%s\n" \
        "$ANNOT_EXAMPLE_HITS"

    printf "new_bacteriocin_prediction_performed\tNO\n"

    printf "reads_remapped\tNO\n"

    printf "amplicon_data_used\tNO\n"

    printf "next_step\t98B_select_specific_mining_strategy\n"

} > "${OUT}/98A_global_summary.tsv"


# ============================================================
# 16. Methodological scope
# ============================================================

STAGE="methodological_scope"

cat > "${OUT}/98A_methodological_scope.tsv" <<'EOF'
field	value
workflow	98_specific_bacteriocin_mining
step	98A_resource_and_existing_evidence_audit
existing_ATTRLOC_loci	reused_not_replaced
primary_goal	find_additional_bacteriocin_or_RiPP_candidates_and_strengthen_existing_candidates
candidate_unit	gene_or_gene_cluster_with_sequence_and_genomic_context
annotation_keyword_hit	insufficient_alone_for_bacteriocin_claim
transporter_or_immunity_gene	insufficient_alone_for_bacteriocin_claim
high_confidence_goal	precursor_or_core_peptide_plus_specific_biosynthetic_or_transport_or_immunity_context
community_contigs	positive_evidence_allowed
incomplete_MAG_absence	not_interpreted_as_biological_absence
sample_recurrence_not_yet_assessed	YES
new_mapping_performed	NO
amplicon_data_used	NO
EOF


# ============================================================
# 17. Final validation
# ============================================================

STAGE="final_validation"

REQUIRED=(
    "98A_known_tool_status.tsv"
    "98A_relevant_containers.txt"
    "98A_relevant_database_paths.txt"
    "98A_existing_bacteriocin_paths.txt"
    "98A_steps86_90_candidate_files.tsv"
    "98A_eggnog_annotation_files.txt"
    "98A_protein_fasta_files.txt"
    "98A_gene_context_files.txt"
    "98A_existing_specialized_tool_outputs.txt"
    "98A_annotation_keyword_file_hits.tsv"
    "98A_annotation_keyword_examples.tsv"
    "98A_global_summary.tsv"
    "98A_methodological_scope.tsv"
)


for NAME in "${REQUIRED[@]}"
do

    if [[ ! -e "${OUT}/${NAME}" ]]
    then

        echo \
        "ERROR: required audit artifact missing: ${NAME}" \
        >&2

        exit 20

    fi

done


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/98A_COMPLETE.ok"


echo
echo "============================================================"
echo "98A GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/98A_global_summary.tsv"


echo
echo "============================================================"
echo "98A FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
