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

OUT="${ROOT}/98_bacteriocin_mining/98C3BR2R_BAGEL_salvage"

MANIFEST="${ROOT}/98_bacteriocin_mining/98C3B0_bagel_audit/98C3B0_negative7_manifest.tsv"

CONTEXT="${ROOT}/98_bacteriocin_mining/98C2_positive_context/98C2_positive_context_summary.tsv"

PY="${SCRIPT_DIR}/98c3br2r_salvage_bagel.py"

STDOUT="${BASE}/BAGEL_wrapper.stdout"
STDERR="${BASE}/BAGEL_wrapper.stderr"
COMMANDLOG="${SESSION}/system_command.log"


rm -rf "$OUT"
mkdir -p "$OUT"


echo "============================================================"
echo "98C3BR2R - SALVAGE EXISTING BAGEL RUN"
echo "NO BAGEL REEXECUTION"
echo "START=$(date -Iseconds)"
echo "============================================================"


for F in \
    "$MANIFEST" \
    "$CONTEXT" \
    "$PY" \
    "$COMMANDLOG" \
    "${SESSION}/00.filenames_db.json" \
    "${SESSION}/00.all_contigs.table"
do

    [[ -s "$F" ]] || {
        echo "ERROR: missing $F" >&2
        exit 20
    }

done


[[ -d "${SESSION}/hmm_rules_domtblout" ]] || {
    echo "ERROR: HMM output directory absent" >&2
    exit 21
}


python3 \
    "$PY" \
    "$MANIFEST" \
    "$CONTEXT" \
    "$SESSION" \
    "$STDOUT" \
    "$STDERR" \
    "$COMMANDLOG" \
    "$OUT"


VALID=$(
    awk -F $'\t' '
        $1 == "real_BAGEL_processing_validation" {
            print $2
        }
    ' "${OUT}/98C3BR2R_global_summary.tsv"
)


if [[ "$VALID" != "PASS" ]]
then

    echo \
    "ERROR: BAGEL processing validation did not pass" \
    >&2

    column -t -s $'\t' \
        "${OUT}/98C3BR2R_target_execution_QC.tsv" \
        >&2 \
        || true

    exit 30

fi


cat > "${OUT}/98C3BR2R_methodological_scope.tsv" <<'EOF'
field	value
analysis	salvage_and_validation_of_existing_133401_BAGEL_outputs
BAGEL_reexecuted	NO
reason	job_133401_failed_after_BAGEL_due_to_shell_validation_bug
shell_failure	MISSING_AOI_counter_used_command_substitution_instead_of_arithmetic_expansion
BAGEL_target_processing_validation	per_target_HMM_plus_BLAST_plus_AOI_merge_commands
missing_AOI_table_after_complete_detection	interpreted_as_no_AOI_not_processing_failure
BAGEL_queryname_suffix_handling	prefix_mapping_to_original_ICONTIG
BAGEL_AOI_containing_target	functional_confirmation_NO
BAGEL_negative	biological_absence_NO
expression	NOT_ASSESSED
antimicrobial_activity	NOT_ASSESSED
reads_remapped	NO
amplicon_data_used	NO
EOF


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/98C3BR2R_COMPLETE.ok"


echo
echo "============================================================"
echo "GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/98C3BR2R_global_summary.tsv"


echo
echo "============================================================"
echo "98C3BR2R FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
