#!/bin/bash

set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

OUT="${ROOT}/98_read_taxonomy/97C2A_validation_audit"

mkdir -p "$OUT"

REPORT="${OUT}/97C2A_tool_audit.tsv"


echo "============================================================"
echo "97C2A - PATHOGEN VALIDATION AUDIT"
echo "START: $(date -Iseconds)"
echo "============================================================"


# ============================================================
# Helpers
# ============================================================

tool_row () {

    local NAME="$1"
    local CMD="$2"

    local PATH_FOUND="NOT_FOUND"
    local VERSION="NA"

    if command -v "$CMD" >/dev/null 2>&1
    then

        PATH_FOUND="$(command -v "$CMD")"

        case "$CMD" in

            datasets)
                VERSION="$("$CMD" version 2>&1 | head -n 1 || true)"
                ;;

            dataformat)
                VERSION="$("$CMD" --version 2>&1 | head -n 1 || true)"
                ;;

            bowtie2)
                VERSION="$("$CMD" --version 2>&1 | head -n 1 || true)"
                ;;

            samtools)
                VERSION="$("$CMD" --version 2>&1 | head -n 1 || true)"
                ;;

            blastn)
                VERSION="$("$CMD" -version 2>&1 | head -n 1 || true)"
                ;;

            seqkit)
                VERSION="$("$CMD" version 2>&1 | head -n 1 || true)"
                ;;

            curl)
                VERSION="$("$CMD" --version 2>&1 | head -n 1 || true)"
                ;;

            wget)
                VERSION="$("$CMD" --version 2>&1 | head -n 1 || true)"
                ;;

            *)
                VERSION="$("$CMD" --help 2>&1 | head -n 1 || true)"
                ;;

        esac

    fi


    printf "%s\t%s\t%s\n" \
        "$NAME" \
        "$PATH_FOUND" \
        "$VERSION"
}


# ============================================================
# Tool audit
# ============================================================

printf "tool\tpath\tversion\n" > "$REPORT"

tool_row "NCBI_Datasets" "datasets" >> "$REPORT"
tool_row "NCBI_dataformat" "dataformat" >> "$REPORT"

tool_row "Entrez_esearch" "esearch" >> "$REPORT"
tool_row "Entrez_efetch" "efetch" >> "$REPORT"

tool_row "curl" "curl" >> "$REPORT"
tool_row "wget" "wget" >> "$REPORT"

tool_row "bowtie2" "bowtie2" >> "$REPORT"
tool_row "bowtie2_build" "bowtie2-build" >> "$REPORT"

tool_row "samtools" "samtools" >> "$REPORT"

tool_row "blastn" "blastn" >> "$REPORT"
tool_row "makeblastdb" "makeblastdb" >> "$REPORT"

tool_row "seqkit" "seqkit" >> "$REPORT"


# ============================================================
# Known containers
# ============================================================

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"

BOWTIE2_SIF="${BOWTIE2_SIF:-/opt/ohpc/pub/containers/BIO/bowtie2-2.5.4.sif}"
SAMTOOLS_SIF="${SAMTOOLS_SIF:-/opt/ohpc/pub/containers/BIO/samtools-1.21.sif}"
BLAST_SIF="${BLAST_SIF:-/opt/ohpc/pub/containers/BIO/blast-2.16.0.sif}"


{
    printf "resource\tstatus\tpath\n"

    for P in \
        "$APPTAINER" \
        "$BOWTIE2_SIF" \
        "$SAMTOOLS_SIF" \
        "$BLAST_SIF"
    do

        if [[ -e "$P" ]]
        then
            printf "%s\tFOUND\t%s\n" "$(basename "$P")" "$P"
        else
            printf "%s\tNOT_FOUND\t%s\n" "$(basename "$P")" "$P"
        fi

    done

} > "${OUT}/97C2A_container_audit.tsv"


# ============================================================
# Search for Datasets / Entrez in software tree
# ============================================================

find "${OHPC_PUBLIC_ROOT:-/opt/ohpc/pub}" \
    -type f \
    \( \
        -name datasets \
        -o -name dataformat \
        -o -name esearch \
        -o -name efetch \
    \) \
    -print \
    2>/dev/null \
    | head -n 100 \
    > "${OUT}/97C2A_ncbi_tools_search.txt" \
    || true


# ============================================================
# Search relevant local databases
# ============================================================

find "${OHPC_PUBLIC_ROOT:-/opt/ohpc/pub}" \
    /scratch/global/"$USER"/Shotgun_MAGs_Turrialba/databases \
    -maxdepth 6 \
    \( \
        -iname '*refseq*' \
        -o -iname '*genbank*' \
        -o -iname '*kraken*' \
        -o -iname '*nt*.nal' \
        -o -iname '*nt*.ndb' \
        -o -iname '*bacteria*.fna*' \
        -o -iname '*assembly*summary*' \
    \) \
    -print \
    2>/dev/null \
    | head -n 300 \
    > "${OUT}/97C2A_reference_database_search.txt" \
    || true


# ============================================================
# Validate host-free reads for proposed samples
# ============================================================

SAMPLES=(
    M2_3
    M3_3
    L2_3
    L3_1
    M1_1
    M2_1
    L2_2
    L1_3
)

printf \
"sample\tR1\tR2\tR1_bytes\tR2_bytes\tstatus\n" \
> "${OUT}/97C2A_selected_sample_reads.tsv"


for SAMPLE in "${SAMPLES[@]}"
do

    R1="${ROOT}/04_host_removed/${SAMPLE}/${SAMPLE}_1.hostfree.fastq.gz"
    R2="${ROOT}/04_host_removed/${SAMPLE}/${SAMPLE}_2.hostfree.fastq.gz"

    if [[ -s "$R1" && -s "$R2" ]]
    then

        STATUS="PASS"

        B1=$(stat -c '%s' "$R1")
        B2=$(stat -c '%s' "$R2")

    else

        STATUS="FAIL"

        B1=0
        B2=0

    fi


    printf \
"%s\t%s\t%s\t%s\t%s\t%s\n" \
"$SAMPLE" \
"$R1" \
"$R2" \
"$B1" \
"$B2" \
"$STATUS" \
>> "${OUT}/97C2A_selected_sample_reads.tsv"

done


# ============================================================
# Candidate validation plan
# ============================================================

cat > "${OUT}/97C2A_candidate_validation_plan.tsv" <<'EOF'
priority	target_species	sample	validation_priority	rationale
tier1	Salmonella enterica	M2_3	VERY_HIGH	strongest_Salmonella_classifier_signal
tier1	Salmonella enterica	M3_3	HIGH	second_strongest_Salmonella_signal
tier1	Listeria monocytogenes	L2_3	VERY_HIGH	high_Kraken_to_Bracken_ratio_and_thesis_pathogen
tier1	Listeria monocytogenes	L3_1	HIGH	recurrent_Listeria_signal_and_thesis_pathogen
tier1	Staphylococcus aureus	M1_1	VERY_HIGH	strongest_S_aureus_signal_and_thesis_pathogen
tier1	Escherichia coli	M2_1	HIGH	strong_classifier_signal_but_close_taxonomic_neighbors
tier2	Yersinia enterocolitica	M2_3	VERY_HIGH	Kraken_and_Bracken_nearly_agree
tier2	Cronobacter sakazakii	L2_2	HIGH	strong_direct_classifier_support
tier2	Bacillus cereus	L2_3	MEDIUM	B_cereus_group_species_resolution_is_difficult
tier2	Clostridium perfringens	L1_3	MEDIUM	lower_but_direct_signal
tier2	Campylobacter jejuni	L2_2	LOW	very_low_absolute_classifier_signal
EOF


# ============================================================
# Filesystem
# ============================================================

df -h \
    "$ROOT" \
    > "${OUT}/97C2A_filesystem.txt"


# ============================================================
# Network probe: lightweight only
#
# Do not download genomes.
# ============================================================

NETWORK_STATUS="NOT_TESTED"

if command -v curl >/dev/null 2>&1
then

    if curl \
        -L \
        -sS \
        --connect-timeout 15 \
        --max-time 30 \
        -o /dev/null \
        -w '%{http_code}' \
        "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/Salmonella%20enterica/dataset_report" \
        > "${OUT}/97C2A_ncbi_http_code.txt" \
        2> "${OUT}/97C2A_ncbi_network.err"
    then

        CODE=$(cat "${OUT}/97C2A_ncbi_http_code.txt")

        if [[ "$CODE" =~ ^2|3 ]]
        then
            NETWORK_STATUS="PASS"
        else
            NETWORK_STATUS="HTTP_${CODE}"
        fi

    else

        NETWORK_STATUS="FAIL"

    fi

else

    printf "curl_not_in_PATH\n" \
        > "${OUT}/97C2A_ncbi_network.err"

    NETWORK_STATUS="NO_CURL"

fi


# ============================================================
# Global summary
# ============================================================

READ_FAILS=$(
    awk -F'\t' '
        NR > 1 && $6 != "PASS" {
            n++
        }
        END {
            print n+0
        }
    ' "${OUT}/97C2A_selected_sample_reads.tsv"
)


{
    printf "metric\tvalue\n"

    printf "selected_unique_samples\t8\n"
    printf "candidate_sample_target_pairs\t11\n"

    printf "tier1_pairs\t6\n"
    printf "tier2_pairs\t5\n"

    printf "hostfree_read_failures\t%s\n" "$READ_FAILS"

    printf "NCBI_network_probe\t%s\n" "$NETWORK_STATUS"

    printf "reference_genomes_downloaded\tNO\n"

    printf "mapping_performed\tNO\n"

    printf "preferred_validation_strategy\tcompetitive_reference_panel_mapping\n"

    printf "next_step\t97C2B_reference_panel_design_and_download\n"

} > "${OUT}/97C2A_global_summary.tsv"


# ============================================================
# Methodological scope
# ============================================================

cat > "${OUT}/97C2A_methodological_scope.tsv" <<'EOF'
field	value
analysis	resource_audit_for_high_specificity_pathogen_validation
upstream_screening	97C1_Kraken2_Bracken_classifier_support
validation_goal	discriminate_target_species_signal_from_close_taxonomic_neighbors
planned_primary_method	competitive_reference_panel_mapping
planned_secondary_metrics	MAPQ_filtered_depth;breadth;mean_depth;reference_specificity
viability_inference	NO
pathogenicity_inference	NO
food_safety_confirmation	NO
short_read_limitation	close_species_may_remain_unresolvable
Bacillus_cereus_note	B_cereus_group_requires_extra_caution
E_coli_note	E_coli_and_Shigella_can_be_difficult_to_separate_with_short_reads
amplicon_data_used	NO
EOF


# ============================================================
# Final validation
# ============================================================

for F in \
    "${OUT}/97C2A_tool_audit.tsv" \
    "${OUT}/97C2A_container_audit.tsv" \
    "${OUT}/97C2A_selected_sample_reads.tsv" \
    "${OUT}/97C2A_candidate_validation_plan.tsv" \
    "${OUT}/97C2A_global_summary.tsv" \
    "${OUT}/97C2A_methodological_scope.tsv"
do

    [[ -s "$F" ]] || {
        echo "ERROR: missing output $F" >&2
        exit 20
    }

done


if [[ "$READ_FAILS" -ne 0 ]]
then
    echo "ERROR: selected sample FASTQ audit failed" >&2
    exit 21
fi


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/97C2A_COMPLETE.ok"


echo
echo "============================================================"
echo "97C2A GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/97C2A_global_summary.tsv"


echo
echo "============================================================"
echo "97C2A FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
