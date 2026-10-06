#!/bin/bash

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

BASE="${ROOT}/98_bacteriocin_mining/98C3BR2_targeted_BAGEL"

SESSION="${BASE}/session"
QUERYFOLDER="${SESSION}/queryfolder"
SUMMARY="${BASE}/summary"

B0="${ROOT}/98_bacteriocin_mining/98C3B0_bagel_audit"

TARGET_FASTA="${B0}/98C3B0_antismash_negative7.fna"
MANIFEST="${B0}/98C3B0_negative7_manifest.tsv"

CONTEXT="${ROOT}/98_bacteriocin_mining/98C2_positive_context/98C2_positive_context_summary.tsv"

BAGEL="${ROOT}/tools/BAGEL4_hpc"
ENV="${ROOT}/envs/bagel4-1.2"

PERL="${ENV}/bin/perl"

PY="${SCRIPT_DIR}/98c3b_bagel_parse.py"


rm -rf "$BASE"

mkdir -p \
    "$SESSION" \
    "$QUERYFOLDER" \
    "$SUMMARY"


echo "============================================================"
echo "98C3BR2 - TARGETED BAGEL4 WITH PREPOPULATED QUERYFOLDER"
echo "START=$(date -Iseconds)"
echo "============================================================"


# ============================================================
# 1. Required inputs
# ============================================================

for F in \
    "$TARGET_FASTA" \
    "$MANIFEST" \
    "$CONTEXT" \
    "$PERL" \
    "$PY" \
    "${BAGEL}/bagel4_wrapper.pl" \
    "${BAGEL}/bagel4.conf"
do

    [[ -s "$F" ]] || {
        echo "ERROR: missing $F" >&2
        exit 20
    }

done


[[ -d "${BAGEL}/lib" ]] || {
    echo "ERROR: BAGEL library directory missing" >&2
    exit 21
}


# ============================================================
# 2. Environment
# ============================================================

export PATH="${ENV}/bin:${PATH}"

export PERL5LIB="${BAGEL}/lib${PERL5LIB:+:${PERL5LIB}}"


echo "PERL=$PERL"
echo "PERL5LIB=$PERL5LIB"


"$PERL" \
    -I"${BAGEL}/lib" \
    -c \
    "${BAGEL}/bagel4_wrapper.pl"


echo "BAGEL_PERL_PREFLIGHT=PASS"


# ============================================================
# 3. IMPORTANT:
#    Put the seven FASTAs directly in:
#
#        session/queryfolder
#
#    because this patched BAGEL wrapper actually scans
#    $sessiondir/queryfolder.
# ============================================================

python3 \
    "$PY" \
    prepare \
    "$TARGET_FASTA" \
    "$MANIFEST" \
    "$QUERYFOLDER"


NQUERY=$(
    find "$QUERYFOLDER" \
        -maxdepth 1 \
        -type f \
        -name '*.fna' \
    | wc -l
)


echo "SESSION_QUERY_FASTAS_BEFORE_RUN=$NQUERY"


if [[ "$NQUERY" -ne 7 ]]
then

    echo \
    "ERROR: expected 7 FASTAs in session/queryfolder; found $NQUERY" \
    >&2

    exit 22

fi


echo
echo "Prepared BAGEL query files:"

find "$QUERYFOLDER" \
    -maxdepth 1 \
    -type f \
    -name '*.fna' \
    -printf '%f\n' \
| sort


# ============================================================
# 4. Session writability
# ============================================================

TEST="${SESSION}/.write_test"

printf "test\n" > "$TEST"

[[ -s "$TEST" ]] || {
    echo "ERROR: session is not writable" >&2
    exit 23
}

rm -f "$TEST"


echo "SESSION_PRECREATED_AND_WRITABLE=YES"


# ============================================================
# 5. Run BAGEL
#
# The wrapper's implementation scans session/queryfolder.
# We therefore also pass that same path via -query so
# provenance and actual scan location agree.
# ============================================================

echo
echo "============================================================"
echo "RUN BAGEL4"
echo "START=$(date -Iseconds)"
echo "============================================================"


set +e

"$PERL" \
    -I"${BAGEL}/lib" \
    "${BAGEL}/bagel4_wrapper.pl" \
    -s "$SESSION" \
    -query "$QUERYFOLDER" \
    -r '*.fna' \
    > "${BASE}/BAGEL_wrapper.stdout" \
    2> "${BASE}/BAGEL_wrapper.stderr"

RC=$?

set -e


echo "BAGEL_WRAPPER_EXIT_CODE=$RC"


if [[ "$RC" -ne 0 ]]
then

    echo "ERROR: BAGEL wrapper returned $RC" >&2

    cat \
        "${BASE}/BAGEL_wrapper.stderr" \
        >&2 \
        || true

    exit 30

fi


# ============================================================
# 6. Reject silent fatal errors
# ============================================================

if grep -Eqi \
'Could not write|Can.t stat|No such file or directory|Can.t locate|Permission denied|Segmentation fault|fatal error' \
"${BASE}/BAGEL_wrapper.stderr"
then

    echo "ERROR: fatal-looking BAGEL stderr detected" >&2

    cat \
        "${BASE}/BAGEL_wrapper.stderr" \
        >&2

    exit 31

fi


# ============================================================
# 7. Session initialization QC
# ============================================================

for F in \
    "${SESSION}/bagel4.conf.json" \
    "${SESSION}/BAGEL_wrapper.log" \
    "${SESSION}/00.all_contigs.table" \
    "${SESSION}/00.filenames_db.json" \
    "${SESSION}/system_command.log"
do

    [[ -s "$F" ]] || {
        echo "ERROR: BAGEL processing file missing/empty: $F" >&2
        exit 32
    }

done


# The query directory must still have exactly seven FASTAs.
NQUERY_AFTER=$(
    find "$QUERYFOLDER" \
        -maxdepth 1 \
        -type f \
        -name '*.fna' \
    | wc -l
)


echo "SESSION_QUERY_FASTAS_AFTER_RUN=$NQUERY_AFTER"


if [[ "$NQUERY_AFTER" -ne 7 ]]
then

    echo \
    "ERROR: session/queryfolder contains $NQUERY_AFTER FASTAs after run; expected 7" \
    >&2

    exit 33

fi


# ============================================================
# 8. Confirm BAGEL recognized all seven sequences
# ============================================================

INPUT_COUNT=$(
    grep -Eo \
        'InputfilesCount=[0-9]+' \
        "${BASE}/BAGEL_wrapper.stdout" \
        "${SESSION}/BAGEL_wrapper.log" \
        2>/dev/null \
    | sed 's/.*InputfilesCount=//' \
    | tail -n 1 \
    || true
)


echo "BAGEL_REPORTED_INPUTFILESCOUNT=${INPUT_COUNT:-UNRESOLVED}"


if [[ "${INPUT_COUNT:-}" != "7" ]]
then

    echo \
    "ERROR: BAGEL did not report InputfilesCount=7" \
    >&2

    echo "--- wrapper stdout ---" >&2
    cat "${BASE}/BAGEL_wrapper.stdout" >&2 || true

    echo "--- BAGEL_wrapper.log ---" >&2
    cat "${SESSION}/BAGEL_wrapper.log" >&2 || true

    exit 34

fi


# ============================================================
# 9. Require one AOI.table per target
#
# Empty AOI.table is allowed and means no AOI for that input.
# Missing AOI.table means that target was not processed fully.
# ============================================================

N_AOI_TABLES=0
MISSING_AOI=0


while IFS= read -r F
do

    CID=$(
        basename "$F" .fna
    )

    AOI="${SESSION}/${CID}.AOI.table"


    if [[ -e "$AOI" ]]
    then

        N_AOI_TABLES=$(
            (
                N_AOI_TABLES + 1
            )
        )

    else

        echo \
        "MISSING_AOI_TABLE=$CID" \
        >&2

        MISSING_AOI=$(
            (
                MISSING_AOI + 1
            )
        )

    fi

done < <(
    find "$QUERYFOLDER" \
        -maxdepth 1 \
        -type f \
        -name '*.fna' \
    | sort
)


echo "AOI_TABLES_PRESENT=$N_AOI_TABLES"


if [[ "$MISSING_AOI" -ne 0 ]]
then

    echo \
    "ERROR: $MISSING_AOI target(s) lack AOI.table" \
    >&2

    exit 35

fi


# ============================================================
# 10. HMM execution QC
# ============================================================

[[ -d "${SESSION}/hmm_rules_domtblout" ]] || {
    echo "ERROR: hmm_rules_domtblout absent" >&2
    exit 36
}


N_HMM_FILES=$(
    find "${SESSION}/hmm_rules_domtblout" \
        -type f \
    | wc -l
)


echo "HMM_OUTPUT_FILES=$N_HMM_FILES"


if [[ "$N_HMM_FILES" -lt 1 ]]
then

    echo "ERROR: no BAGEL HMM output files generated" >&2
    exit 37

fi


# ============================================================
# 11. Confirm commands were actually launched
# ============================================================

N_SYSTEM_COMMANDS=$(
    grep -cv \
        '^[[:space:]]*$' \
        "${SESSION}/system_command.log" \
        || true
)


echo "SYSTEM_COMMAND_LOG_LINES=$N_SYSTEM_COMMANDS"


if [[ "$N_SYSTEM_COMMANDS" -lt 7 ]]
then

    echo \
    "ERROR: suspiciously few BAGEL commands executed: $N_SYSTEM_COMMANDS" \
    >&2

    exit 38

fi


# ============================================================
# 12. Completion
# ============================================================

[[ -s "${SESSION}/sessionstop" ]] || {
    echo "ERROR: sessionstop absent" >&2
    exit 39
}


if ! grep -q \
    'Analysis done' \
    "${SESSION}/sessionstop"
then

    echo "ERROR: BAGEL did not report Analysis done" >&2
    exit 40

fi


echo "BAGEL_REAL_PROCESSING_VALIDATION=PASS"


# ============================================================
# 13. Parse biological results
# ============================================================

python3 \
    "$PY" \
    parse \
    "$MANIFEST" \
    "$CONTEXT" \
    "$SESSION" \
    "$SUMMARY"


# ============================================================
# 14. Execution QC
# ============================================================

{
    printf "metric\tvalue\n"

    printf "target_FASTAs_expected\t7\n"
    printf "target_FASTAs_before_run\t%s\n" "$NQUERY"
    printf "target_FASTAs_after_run\t%s\n" "$NQUERY_AFTER"
    printf "BAGEL_reported_InputfilesCount\t%s\n" "$INPUT_COUNT"
    printf "AOI_tables_present\t%s\n" "$N_AOI_TABLES"
    printf "HMM_output_files\t%s\n" "$N_HMM_FILES"
    printf "system_command_log_lines\t%s\n" "$N_SYSTEM_COMMANDS"
    printf "bagel4_conf_json_present\tYES\n"
    printf "filenames_db_present\tYES\n"
    printf "sessionstop_analysis_done\tYES\n"
    printf "wrapper_exit_code\t%s\n" "$RC"
    printf "real_processing_validation\tPASS\n"

} > "${SUMMARY}/98C3BR2_execution_QC.tsv"


# ============================================================
# 15. Methodological scope
# ============================================================

cat > "${SUMMARY}/98C3BR2_methodological_scope.tsv" <<'EOF'
field	value
analysis	validated_targeted_BAGEL4_rescue
previous_job_133394	scientifically_invalid_no_query_processing
previous_job_133399	failed_validation_before_query_processing
target_contigs	7
query_staging	session_queryfolder_prepopulated_before_wrapper
reason	patched_wrapper_scans_sessiondir_queryfolder_directly
BAGEL_reported_input_count_required	7
AOI_table_required_for_every_target	YES
empty_AOI_table	allowed_computational_negative
missing_AOI_table	invalid_processing
HMM_processing_required	YES
BAGEL_AOI_absence_after_valid_processing	biological_absence_NO
BAGEL_AOI_containing_target	functional_confirmation_NO
additional_incomplete_context	equivalent_to_novel_bacteriocin_NO
expression	NOT_ASSESSED
antimicrobial_activity	NOT_ASSESSED
reads_remapped	NO
amplicon_data_used	NO
EOF


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${SUMMARY}/98C3BR2_COMPLETE.ok"


echo
echo "============================================================"
echo "98C3BR2 EXECUTION QC"
echo "============================================================"

column -t -s $'\t' \
"${SUMMARY}/98C3BR2_execution_QC.tsv"


echo
echo "============================================================"
echo "98C3BR2 BIOLOGICAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${SUMMARY}/98C3B_global_summary.tsv"


echo
echo "============================================================"
echo "98C3BR2 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
