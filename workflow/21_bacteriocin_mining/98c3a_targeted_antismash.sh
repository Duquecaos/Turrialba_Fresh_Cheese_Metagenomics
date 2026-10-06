#!/bin/bash

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

BASE="${ROOT}/98_bacteriocin_mining/98C3A_targeted_antismash"

PREP="${BASE}/prepared"
ASOUT="${BASE}/antismash"
SUMMARY="${BASE}/summary"

RUNTIME_ROOT="${BASE}/runtime_antismash_package"
RUNTIME_PKG="${RUNTIME_ROOT}/antismash"

PRIMARY="${ROOT}/98_bacteriocin_mining/98C2_positive_context/98C2_primary_target_contigs.tsv"

PY="${SCRIPT_DIR}/98c3a_prepare_and_parse.py"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"

ANTISMASH="${ANTISMASH:-/opt/ohpc/pub/containers/BIO/antismash-7.1.0.sif}"

ASDB="${ROOT}/databases/antismash_7.1.0"

ASPKG_SOURCE="${ROOT}/databases/antismash_7.1.0_package/antismash"

THREADS="${SLURM_CPUS_PER_TASK:-16}"


echo "============================================================"
echo "98C3AR2 - TARGETED ANTISMASH WITH WRITABLE PACKAGE"
echo "START=$(date -Iseconds)"
echo "THREADS=$THREADS"
echo "============================================================"


# ============================================================
# 1. Preserve previous failed run
# ============================================================

if [[ -d "$BASE" ]]
then

    if [[ ! -s "${BASE}/summary/98C3A_COMPLETE.ok" ]]
    then

        ARCHIVE="${BASE}.FAILED_previous_$(date +%Y%m%dT%H%M%S)"

        mv \
            "$BASE" \
            "$ARCHIVE"

        echo "Archived previous failed run:"
        echo "$ARCHIVE"

    else

        echo "WARNING: previous completed 98C3A exists."
        echo "It will not be overwritten automatically." >&2
        exit 10

    fi

fi


mkdir -p \
    "$PREP" \
    "$ASOUT" \
    "$SUMMARY" \
    "$RUNTIME_ROOT"


# ============================================================
# 2. Required inputs
# ============================================================

for F in \
    "$PRIMARY" \
    "$PY" \
    "$ANTISMASH"
do

    [[ -s "$F" ]] || {
        echo "ERROR: missing $F" >&2
        exit 20
    }

done


[[ -d "$ASDB" ]] || {
    echo "ERROR: antiSMASH database root missing: $ASDB" >&2
    exit 21
}


[[ -d "$ASPKG_SOURCE" ]] || {
    echo "ERROR: prepared antiSMASH package missing: $ASPKG_SOURCE" >&2
    exit 22
}


# ============================================================
# 3. Validate important prepared files
# ============================================================

REQUIRED_PREPARED=(
    "${ASPKG_SOURCE}/modules/lanthipeptides/data/lanthipeptide.classifier.pkl"
    "${ASPKG_SOURCE}/modules/lanthipeptides/data/lanthipeptide.scaler.pkl"

    "${ASPKG_SOURCE}/modules/lassopeptides/data/lassopeptide.classifier.pkl"
    "${ASPKG_SOURCE}/modules/lassopeptides/data/lassopeptide.scaler.pkl"

    "${ASPKG_SOURCE}/modules/sactipeptides/data/sactipeptide.classifier.pkl"
    "${ASPKG_SOURCE}/modules/sactipeptides/data/sactipeptide.scaler.pkl"

    "${ASPKG_SOURCE}/modules/thiopeptides/data/thiopeptide.classifier.pkl"
    "${ASPKG_SOURCE}/modules/thiopeptides/data/thiopeptide.scaler.pkl"

    "${ASPKG_SOURCE}/modules/nrps_pks/data/cterm.fasta.brawn_cache"
    "${ASPKG_SOURCE}/modules/nrps_pks/data/nterm.fasta.brawn_cache"
)


for F in "${REQUIRED_PREPARED[@]}"
do

    [[ -s "$F" ]] || {
        echo "ERROR: prepared antiSMASH file missing: $F" >&2
        exit 23
    }

done


echo "Prepared antiSMASH files: PASS"


# ============================================================
# 4. Create job-private writable package copy
#
# Do not modify the master copy under databases/.
# ============================================================

echo
echo "============================================================"
echo "CREATE WRITABLE ANTISMASH PACKAGE"
echo "============================================================"


cp -a \
    "$ASPKG_SOURCE" \
    "$RUNTIME_PKG"


[[ -d "$RUNTIME_PKG" ]] || {
    echo "ERROR: runtime package copy failed" >&2
    exit 24
}


echo "SOURCE=$ASPKG_SOURCE"
echo "RUNTIME=$RUNTIME_PKG"

du -sh \
    "$RUNTIME_PKG" \
    || true


# ============================================================
# 5. Verify host-side writability
# ============================================================

TESTFILE="${RUNTIME_PKG}/.98C3AR2_write_test"

printf \
"98C3AR2_WRITE_TEST\n" \
> "$TESTFILE"


[[ -s "$TESTFILE" ]] || {
    echo "ERROR: runtime package is not writable" >&2
    exit 25
}


rm -f "$TESTFILE"

echo "HOST_RUNTIME_PACKAGE_WRITABLE=YES"


# ============================================================
# 6. Bind the ENTIRE antiSMASH package
# ============================================================

CONTAINER_PKG="/usr/local/lib/python3.10/site-packages/antismash"

BIND_ARGS=(
    --bind "${ROOT}:${ROOT}"
    --bind "${RUNTIME_PKG}:${CONTAINER_PKG}"
)


# ============================================================
# 7. Container preflight
# ============================================================

echo
echo "============================================================"
echo "CONTAINER PREFLIGHT"
echo "============================================================"


"$APPTAINER" exec \
    "${BIND_ARGS[@]}" \
    "$ANTISMASH" \
    python3 - <<'PY'
import os
import antismash
from pathlib import Path

pkg = Path(antismash.__file__).resolve().parent

print(f"ANTISMASH_IMPORTED_FROM={pkg}")

expected = Path(
    "/usr/local/lib/python3.10/site-packages/antismash"
)

if pkg != expected:
    raise SystemExit(
        f"ERROR: unexpected antiSMASH package path: {pkg}"
    )


required = [
    pkg / "modules/lanthipeptides/data/lanthipeptide.scaler.pkl",
    pkg / "modules/lassopeptides/data/lassopeptide.scaler.pkl",
    pkg / "modules/sactipeptides/data/sactipeptide.scaler.pkl",
    pkg / "modules/thiopeptides/data/thiopeptide.scaler.pkl",
    pkg / "modules/nrps_pks/data/cterm.fasta.brawn_cache",
    pkg / "modules/nrps_pks/data/nterm.fasta.brawn_cache",
]


for p in required:

    if not p.is_file():
        raise SystemExit(
            f"ERROR: missing inside container: {p}"
        )

    print(
        f"FOUND\t{p.stat().st_size}\t{p}"
    )


# Real write test from INSIDE the container.
test = pkg / ".98C3AR2_container_write_test"

test.write_text(
    "container_write_test\n"
)

if not test.is_file():
    raise SystemExit(
        "ERROR: package bind is not writable inside container"
    )

test.unlink()

print("ANTISMASH_FULL_PACKAGE_BIND_WRITABLE=YES")
print("98C3AR2_CONTAINER_PREFLIGHT=PASS")
PY


"$APPTAINER" exec \
    "${BIND_ARGS[@]}" \
    "$ANTISMASH" \
    antismash \
    --version


# ============================================================
# 8. Prepare exactly the 40 primary contigs
# ============================================================

echo
echo "============================================================"
echo "PREPARE PRIMARY CONTIGS"
echo "============================================================"


python3 \
    "$PY" \
    prepare \
    "$ROOT" \
    "$PRIMARY" \
    "$PREP"


FASTA="${PREP}/98C3A_primary40_contigs.fna"


[[ -s "$FASTA" ]] || {
    echo "ERROR: target FASTA missing" >&2
    exit 26
}


NSEQ=$(
    grep -c '^>' "$FASTA"
)


if [[ "$NSEQ" -ne 40 ]]
then

    echo "ERROR: expected 40 records; observed $NSEQ" >&2
    exit 27

fi


echo "TARGET_FASTA_RECORDS=$NSEQ"


# ============================================================
# 9. antiSMASH
# ============================================================

echo
echo "============================================================"
echo "RUN ANTISMASH"
echo "START=$(date -Iseconds)"
echo "============================================================"


"$APPTAINER" exec \
    "${BIND_ARGS[@]}" \
    "$ANTISMASH" \
    antismash \
    --taxon bacteria \
    --genefinding-tool prodigal \
    --cpus "$THREADS" \
    --databases "$ASDB" \
    --output-dir "$ASOUT" \
    --output-basename 98C3A_primary40 \
    --cb-knownclusters \
    "$FASTA"


echo
echo "ANTISMASH_END=$(date -Iseconds)"


# ============================================================
# 10. Validate main outputs
# ============================================================

JSON="${ASOUT}/98C3A_primary40.json"
GBK="${ASOUT}/98C3A_primary40.gbk"
HTML="${ASOUT}/index.html"


for F in \
    "$JSON" \
    "$GBK" \
    "$HTML"
do

    [[ -s "$F" ]] || {
        echo "ERROR: antiSMASH output missing: $F" >&2
        exit 30
    }

done


N_REGION=$(
    find "$ASOUT" \
        -maxdepth 1 \
        -type f \
        -name '*.region*.gbk' \
        | wc -l
)


echo "ANTISMASH_REGION_FILES=$N_REGION"


# ============================================================
# 11. Parse
# ============================================================

python3 \
    "$PY" \
    parse \
    "$PRIMARY" \
    "$ASOUT" \
    "$SUMMARY"


# ============================================================
# 12. Scope
# ============================================================

cat > "${SUMMARY}/98C3A_methodological_scope.tsv" <<'EOF'
field	value
analysis	targeted_antismash_on_primary_positive_incomplete_contigs
antiSMASH_version	7.1.0
input_contigs	40_primary_positive_unique_contigs
input_selection	homology_or_curated_GA_HMM_positive
supplementary_noGA_only_contigs_included	NO
gene_finding	Prodigal_by_antismash
runtime_package	job_private_writable_copy_of_prepared_antismash_7.1.0_package
container_package_bind	full_antismash_python_package
reason_for_full_package_bind	antiSMASH_prerequisite_steps_generate_or_update_internal_model_and_alignment_cache_files
master_prepared_package_modified	NO
same_contig_RiPP_prediction_equivalent_to_functional_confirmation	NO
antiSMASH_negative_on_fragmented_contig_equivalent_to_biological_absence	NO
contig_edge_truncation_requires_cautious_interpretation	YES
additional_to_final18_equivalent_to_novel_bacteriocin	NO
taxonomic_bin_attribution_from_incomplete_bins_requires_caution	YES
expression	NOT_ASSESSED
antimicrobial_activity	NOT_ASSESSED
reads_remapped	NO
amplicon_data_used	NO
EOF


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${SUMMARY}/98C3A_COMPLETE.ok"


echo
echo "============================================================"
echo "98C3A GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
"${SUMMARY}/98C3A_global_summary.tsv"


echo
echo "============================================================"
echo "98C3AR2 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
