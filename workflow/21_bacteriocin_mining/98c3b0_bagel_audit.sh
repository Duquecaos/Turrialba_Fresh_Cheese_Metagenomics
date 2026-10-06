#!/bin/bash

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

OUT="${ROOT}/98_bacteriocin_mining/98C3B0_bagel_audit"

ANTI="${ROOT}/98_bacteriocin_mining/98C3A_targeted_antismash/summary/98C3A_primary_contig_antismash_summary.tsv"

CONTEXT="${ROOT}/98_bacteriocin_mining/98C2_positive_context/98C2_positive_context_summary.tsv"

FASTA="${ROOT}/69_incomplete_unique_contigs/unique_contigs_lt90_all.fna"

PY="${SCRIPT_DIR}/98c3b0_prepare_bagel_targets.py"

BAGEL="${ROOT}/tools/BAGEL4_hpc"

OLD="${ROOT}/36_bagel4_BAL4"

ENV="${ROOT}/envs/bagel4-1.2"


rm -rf "$OUT"
mkdir -p "$OUT"


for F in \
    "$ANTI" \
    "$CONTEXT" \
    "$FASTA" \
    "$PY" \
    "${BAGEL}/bagel4_wrapper.pl"
do

    [[ -s "$F" ]] || {
        echo "ERROR: missing $F" >&2
        exit 20
    }

done


# ============================================================
# 1. Prepare seven antiSMASH-negative targets
# ============================================================

python3 \
    "$PY" \
    "$ANTI" \
    "$CONTEXT" \
    "$FASTA" \
    "$OUT"


# ============================================================
# 2. BAGEL installation audit
# ============================================================

{
    printf "field\tvalue\n"

    printf "BAGEL_root\t%s\n" \
        "$BAGEL"

    printf "wrapper_present\t%s\n" \
        "$(
            [[ -s "${BAGEL}/bagel4_wrapper.pl" ]] \
            && echo YES \
            || echo NO
        )"

    printf "config_present\t%s\n" \
        "$(
            [[ -s "${BAGEL}/bagel4.conf" ]] \
            && echo YES \
            || echo NO
        )"

    printf "environment_present\t%s\n" \
        "$(
            [[ -d "$ENV" ]] \
            && echo YES \
            || echo NO
        )"

    printf "previous_completed_runs\t%s\n" \
        "$(
            find "$OLD" \
                -type f \
                -name BAGEL4_COMPLETED \
                2>/dev/null \
            | wc -l
        )"

} > "${OUT}/98C3B0_bagel_installation.tsv"


# ============================================================
# 3. Perl environment
# ============================================================

PERL=""

if [[ -x "${ENV}/bin/perl" ]]
then

    PERL="${ENV}/bin/perl"

elif command -v perl >/dev/null 2>&1
then

    PERL="$(
        command -v perl
    )"

fi


{
    echo "PERL=${PERL:-NOT_FOUND}"

    if [[ -n "$PERL" ]]
    then

        "$PERL" -v \
        | head -n 3 \
        || true

        echo

        "$PERL" -c \
            "${BAGEL}/bagel4_wrapper.pl" \
            2>&1 \
        || true

    fi

} > "${OUT}/98C3B0_perl_audit.txt"


# ============================================================
# 4. Recover previous BAGEL command evidence
# ============================================================

: > "${OUT}/98C3B0_previous_command_evidence.txt"


find "$OLD" \
    -type f \
    \( \
        -name 'system_command.log' \
        -o \
        -name 'BAGEL_wrapper.log' \
    \) \
    -print \
    2>/dev/null \
| sort \
| head -n 12 \
| while IFS= read -r F
do

    echo "============================================================"
    echo "FILE=$F"
    echo "============================================================"

    grep -Ein \
        'bagel4_wrapper|perl|queryfolder|bagel4\.conf|AOI_prerun|command' \
        "$F" \
        2>/dev/null \
    | head -n 80 \
    || true

    echo

done \
>> "${OUT}/98C3B0_previous_command_evidence.txt"


# ============================================================
# 5. Wrapper source — argument/usage audit
# ============================================================

grep -Ein \
    'Usage|GetOptions|ARGV|query|fasta|session|output|folder|conf|input' \
    "${BAGEL}/bagel4_wrapper.pl" \
    2>/dev/null \
| head -n 160 \
> "${OUT}/98C3B0_wrapper_interface.txt" \
|| true


# ============================================================
# 6. README/setup instructions
# ============================================================

: > "${OUT}/98C3B0_README_run_instructions.txt"


for F in \
    "${BAGEL}/README_BAGEL4_SETUP_and_RUN.txt" \
    "${BAGEL}/README.md" \
    "${BAGEL}/Pipeline_description.txt"
do

    [[ -s "$F" ]] || continue

    echo "============================================================" \
    >> "${OUT}/98C3B0_README_run_instructions.txt"

    echo "FILE=$F" \
    >> "${OUT}/98C3B0_README_run_instructions.txt"

    echo "============================================================" \
    >> "${OUT}/98C3B0_README_run_instructions.txt"

    grep -Ein \
        'wrapper|perl|run|query|fasta|input|output|session|command|usage' \
        "$F" \
        2>/dev/null \
    | head -n 160 \
    >> "${OUT}/98C3B0_README_run_instructions.txt" \
    || true

done


# ============================================================
# 7. Provenance
# ============================================================

if [[ -s "${OLD}/BAGEL4_HPC_provenance.txt" ]]
then

    cp -p \
        "${OLD}/BAGEL4_HPC_provenance.txt" \
        "${OUT}/98C3B0_previous_BAGEL4_provenance.txt"

else

    printf \
    "PROVENANCE_NOT_FOUND\n" \
    > "${OUT}/98C3B0_previous_BAGEL4_provenance.txt"

fi


# ============================================================
# 8. Scope
# ============================================================

cat > "${OUT}/98C3B0_methodological_scope.tsv" <<'EOF'
field	value
analysis	BAGEL_preflight_and_target_selection
target_set	7_primary_contigs_without_antiSMASH_RiPP_region
BAGEL_run_performed	NO
purpose	recover_exact_previously_working_HPC_invocation_before_new_BAGEL_run
antiSMASH_negative_on_fragmented_contig	biological_absence_NO
antiSMASH_and_98C1_GA_HMM	fully_independent_methods_NO
reason	98C1_profiles_were_taken_from_antiSMASH_profile_resources
BAGEL	role_complementary_specialized_bacteriocin_predictor
functional_activity	NOT_ASSESSED
amplicon_data_used	NO
EOF


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/98C3B0_COMPLETE.ok"


echo "============================================================"
echo "98C3B0 TARGET SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/98C3B0_target_summary.tsv"

echo
echo "============================================================"
echo "98C3B0 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
