#!/bin/bash

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

ASROOT="${ROOT}/26_antismash_rep18"
RESULTS="${ASROOT}/results"

ATTRLOC="${ROOT}/49_bacteriocin_locus_matrices/locus_summary_by_ATTRLOC.tsv"

OUT="${ROOT}/98_bacteriocin_mining/98B0_existing_antismash"

PARSER="${SCRIPT_DIR}/98b0_parse_existing_antismash.py"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"

ANTISMASH_SIF="${ANTISMASH_SIF:-/opt/ohpc/pub/containers/BIO/antismash-7.1.0.sif}"

ASDB="${ROOT}/databases/antismash_7.1.0"

BAGEL_ENV="${ROOT}/envs/bagel4-1.2"

mkdir -p "$OUT"


echo "============================================================"
echo "98B0 - EXISTING ANTISMASH + BAGEL AUDIT"
echo "START=$(date -Iseconds)"
echo "============================================================"


# ============================================================
# 1. Required inputs
# ============================================================

for F in \
    "$PARSER" \
    "$ATTRLOC"
do

    if [[ ! -s "$F" ]]
    then
        echo "ERROR: missing $F" >&2
        exit 20
    fi

done


if [[ ! -d "$RESULTS" ]]
then
    echo "ERROR: antiSMASH results directory absent: $RESULTS" >&2
    exit 21
fi


# ============================================================
# 2. antiSMASH executable/version
# ============================================================

{
    printf "field\tvalue\n"

    if [[ -s "$ANTISMASH_SIF" ]]
    then

        printf \
        "container\t%s\n" \
        "$ANTISMASH_SIF"

        VERSION=$(
            "$APPTAINER" exec \
                "$ANTISMASH_SIF" \
                antismash \
                --version \
                2>&1 \
            | head -n 1 \
            || true
        )

        printf \
        "version\t%s\n" \
        "${VERSION:-UNRESOLVED}"

    else

        printf "container\tNOT_FOUND\n"
        printf "version\tNOT_TESTED\n"

    fi

    printf \
    "database_root\t%s\n" \
    "$ASDB"

    printf \
    "database_root_present\t%s\n" \
    "$(
        [[ -d "$ASDB" ]] \
        && echo YES \
        || echo NO
    )"

} > "${OUT}/98B0_antismash_environment.tsv"


# ============================================================
# 3. Existing run summary
# ============================================================

RUNSUMMARY="${ASROOT}/antismash_rep18_run_summary.tsv"

if [[ -s "$RUNSUMMARY" ]]
then

    cp -p \
        "$RUNSUMMARY" \
        "${OUT}/98B0_existing_run_summary.tsv"

else

    printf \
    "status\nRUN_SUMMARY_NOT_FOUND\n" \
    > "${OUT}/98B0_existing_run_summary.tsv"

fi


# ============================================================
# 4. Parse existing antiSMASH output
# ============================================================

python3 \
    "$PARSER" \
    "$RESULTS" \
    "$ATTRLOC" \
    "$OUT"


# ============================================================
# 5. antiSMASH bacteriocin-related HMM inventory
# ============================================================

HMMROOT="${ROOT}/databases/antismash_7.1.0_package/antismash/detection/hmm_detection/data"

printf \
"hmm_file\n" \
> "${OUT}/98B0_bacteriocin_HMM_inventory.tsv"


if [[ -d "$HMMROOT" ]]
then

    find "$HMMROOT" \
        -maxdepth 1 \
        -type f \
        -name '*.hmm' \
        -printf '%f\n' \
    | grep -Ei \
        'bacterioc|lactococc|lcn|lanthi|lasso|sacti|thio|glycocin|ripp' \
    | sort \
    | awk '
        {
            print $0
        }
    ' \
    >> "${OUT}/98B0_bacteriocin_HMM_inventory.tsv" \
    || true

fi


# ============================================================
# 6. BAGEL environment audit
# ============================================================

printf \
"metric\tvalue\n" \
> "${OUT}/98B0_bagel_environment.tsv"


if [[ -d "$BAGEL_ENV" ]]
then

    printf \
    "bagel_environment_present\tYES\n" \
    >> "${OUT}/98B0_bagel_environment.tsv"

else

    printf \
    "bagel_environment_present\tNO\n" \
    >> "${OUT}/98B0_bagel_environment.tsv"

fi


BAGEL_CANDIDATES="${OUT}/98B0_bagel_candidate_executables.txt"

: > "$BAGEL_CANDIDATES"


if [[ -d "${BAGEL_ENV}/bin" ]]
then

    find "${BAGEL_ENV}/bin" \
        -maxdepth 1 \
        -type f \
        -printf '%f\n' \
        2>/dev/null \
    | grep -Ei \
        'bagel|bacterioc' \
    | sort \
    > "$BAGEL_CANDIDATES" \
    || true

fi


N_BAGEL_EXE=$(
    awk '
        NF {
            n++
        }
        END {
            print n+0
        }
    ' "$BAGEL_CANDIDATES"
)


printf \
"candidate_bagel_executables\t%s\n" \
"$N_BAGEL_EXE" \
>> "${OUT}/98B0_bagel_environment.tsv"


# Search for BAGEL-named source/scripts,
# but do not treat environment libraries as biological outputs.

{
    find "$BAGEL_ENV" \
        -maxdepth 3 \
        -type f \
        -print \
        2>/dev/null \
    || true
} \
| grep -Ei \
    '/[^/]*bagel[^/]*$' \
| sort \
> "${OUT}/98B0_bagel_named_files.txt" \
|| true


# ============================================================
# 7. Scope
# ============================================================

cat > "${OUT}/98B0_methodological_scope.tsv" <<'EOF'
field	value
analysis	reanalysis_of_existing_antismash_predictions
new_antismash_predictions	NO
primary_antismash_version	7.1.0
reason_for_version_choice	matches_existing_MAG_results_and_local_database_installation
bacteriocin_RiPP_filter	product_or_detection_rule_contains_bacteriocin_or_RiPP_related_terms
ATTRLOC_comparison	level_of_comparison_is_same_MAG_and_same_contig
same_contig_overlap	not_equivalent_to_exact_gene_overlap
BAGEL_status	audited_only_not_used_for_biological_calls
keyword_only_evidence	not_sufficient
absence_from_fragmented_MAG	not_biological_absence
amplicon_data_used	NO
EOF


# ============================================================
# 8. Required outputs
# ============================================================

for F in \
    "${OUT}/98B0_antismash_environment.tsv" \
    "${OUT}/98B0_existing_run_summary.tsv" \
    "${OUT}/98B0_antismash_MAG_inventory.tsv" \
    "${OUT}/98B0_all_antismash_regions.tsv" \
    "${OUT}/98B0_bacteriocin_RiPP_regions.tsv" \
    "${OUT}/98B0_product_summary.tsv" \
    "${OUT}/98B0_ATTRLOC_antismash_comparison.tsv" \
    "${OUT}/98B0_global_summary.tsv" \
    "${OUT}/98B0_bacteriocin_HMM_inventory.tsv" \
    "${OUT}/98B0_bagel_environment.tsv" \
    "${OUT}/98B0_bagel_candidate_executables.txt" \
    "${OUT}/98B0_methodological_scope.tsv"
do

    if [[ ! -e "$F" ]]
    then
        echo "ERROR: output missing: $F" >&2
        exit 30
    fi

done


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/98B0_COMPLETE.ok"


echo
echo "============================================================"
echo "GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/98B0_global_summary.tsv"


echo
echo "============================================================"
echo "98B0 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
