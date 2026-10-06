#!/bin/bash

set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

BASE="${ROOT}/98_read_taxonomy"
SOFT="${BASE}/software"

KRAKEN_VERSION="2.17.2"
BRACKEN_VERSION="3.1"

KRAKEN_INSTALL="${SOFT}/kraken2-${KRAKEN_VERSION}"
BRACKEN_INSTALL="${SOFT}/Bracken-${BRACKEN_VERSION}"

DBROOT="${ROOT}/databases"
DB="${DBROOT}/kraken2_pluspf_20260626"

DOWNLOADS="${DBROOT}/downloads"
DB_ARCHIVE="${DOWNLOADS}/k2_pluspf_20260626.tar.gz"
DB_MD5="${DOWNLOADS}/pluspf_20260626.md5"

OLDOUT="${BASE}/97B0_prepare"
OUT="${BASE}/97B0R_finalize"

mkdir -p "$OUT"

echo "============================================================"
echo "97B0R - FINALIZE EXISTING KRAKEN2/BRACKEN/PLUSPF"
echo "Inicio: $(date -Iseconds)"
echo "============================================================"


# ============================================================
# 1. Validate executables
# ============================================================

KRAKEN="${KRAKEN_INSTALL}/kraken2"
INSPECT="${KRAKEN_INSTALL}/kraken2-inspect"
BRACKEN="${BRACKEN_INSTALL}/bracken"

for F in \
    "$KRAKEN" \
    "$INSPECT" \
    "$BRACKEN"
do

    if [[ ! -x "$F" ]]
    then
        echo "ERROR: executable missing: $F" >&2
        exit 20
    fi

done

echo "EXECUTABLES=PASS"


# ============================================================
# 2. Validate database core
# ============================================================

for F in \
    hash.k2d \
    opts.k2d \
    taxo.k2d
do

    if [[ ! -s "${DB}/${F}" ]]
    then
        echo "ERROR: database core missing: ${DB}/${F}" >&2
        exit 21
    fi

done

echo "DATABASE_CORE=PASS"


# ============================================================
# 3. Validate Bracken 150 bp distribution
# ============================================================

BRACKEN150="${DB}/database150mers.kmer_distrib"

if [[ ! -s "$BRACKEN150" ]]
then
    echo "ERROR: Bracken 150-mer distribution missing" >&2
    exit 22
fi

echo "BRACKEN_150MER=PASS"


# ============================================================
# 4. Validate already-generated inspect output
#
# IMPORTANT:
# Do NOT rerun kraken2-inspect.
# The original job already loaded/inspected the ~110 GB index.
# ============================================================

INSPECT_HEAD="${OLDOUT}/97B0_kraken2_inspect_head.txt"

if [[ ! -s "$INSPECT_HEAD" ]]
then
    echo "ERROR: previous inspect output missing" >&2
    exit 23
fi

if ! grep -q "Database options" "$INSPECT_HEAD"
then
    echo "ERROR: inspect output lacks database header" >&2
    exit 24
fi

if ! grep -q $'\troot$' "$INSPECT_HEAD"
then
    echo "WARNING: root line not matched exactly; checking text form"

    if ! grep -q "root" "$INSPECT_HEAD"
    then
        echo "ERROR: inspect output lacks root taxonomy line" >&2
        exit 25
    fi
fi

echo "PREVIOUS_KRAKEN_INSPECT=PASS"


# ============================================================
# 5. Validate downloaded archive
# ============================================================

if [[ ! -s "$DB_ARCHIVE" ]]
then
    echo "ERROR: downloaded PlusPF archive missing" >&2
    exit 26
fi

echo "DATABASE_ARCHIVE_PRESENT=PASS"


# ============================================================
# 6. Versions
# ============================================================

KRAKEN_VERSION_OUTPUT="$(
    "$KRAKEN" --version 2>&1 |
    head -n 1
)"

BRACKEN_HELP="$(
    "$BRACKEN" -h 2>&1 |
    head -n 1 ||
    true
)"


{
    printf "field\tvalue\n"

    printf "kraken2_version\t%s\n" \
        "$KRAKEN_VERSION_OUTPUT"

    printf "bracken_requested_version\t%s\n" \
        "$BRACKEN_VERSION"

    printf "bracken_executable\t%s\n" \
        "$BRACKEN"

    printf "database\tPlusPF\n"

    printf "database_release\t2026-06-26\n"

    printf "database_path\t%s\n" \
        "$DB"

    printf "database_archive\t%s\n" \
        "$DB_ARCHIVE"

    printf "bracken_distribution_bp\t150\n"

    printf "shotgun_read_length_bp\t151\n"

    printf "bracken_read_length_strategy\tclosest_prebuilt_150bp_distribution\n"

    printf "previous_prepare_job\t132825\n"

    printf "previous_prepare_job_status\tFAILED_after_successful_install_download_extract_and_inspect\n"

    printf "failure_interpretation\tbookkeeping_pipeline_failure_after_kraken_inspect_output\n"

} > "${OUT}/97B0R_versions_database.tsv"


# ============================================================
# 7. Database inventory
# ============================================================

{
    printf "file\tsize_bytes\n"

    for F in \
        hash.k2d \
        opts.k2d \
        taxo.k2d
    do
        printf "%s\t%s\n" \
            "$F" \
            "$(stat -c '%s' "${DB}/${F}")"
    done

    while IFS= read -r F
    do
        printf "%s\t%s\n" \
            "$(basename "$F")" \
            "$(stat -c '%s' "$F")"

    done < <(
        find "$DB" \
            -maxdepth 1 \
            -type f \
            -name '*mers.kmer_distrib' \
            -print |
        sort
    )

} > "${OUT}/97B0R_database_inventory.tsv"


# ============================================================
# 8. Disk usage
# ============================================================

{
    printf "item\tsize\n"

    printf "database_archive\t%s\n" \
        "$(du -sh "$DB_ARCHIVE" | awk '{print $1}')"

    printf "database_extracted\t%s\n" \
        "$(du -sh "$DB" | awk '{print $1}')"

    printf "kraken2_install\t%s\n" \
        "$(du -sh "$KRAKEN_INSTALL" | awk '{print $1}')"

    printf "bracken_install\t%s\n" \
        "$(du -sh "$BRACKEN_INSTALL" | awk '{print $1}')"

} > "${OUT}/97B0R_disk_usage.tsv"


# ============================================================
# 9. Preserve inspect evidence
# ============================================================

cp \
    "$INSPECT_HEAD" \
    "${OUT}/97B0R_kraken2_inspect_head.txt"


# ============================================================
# 10. Environment for 97B1+
# ============================================================

cat > "${OUT}/97B0R_environment.sh" <<EOF
export KRAKEN2="${KRAKEN}"
export KRAKEN2_INSPECT="${INSPECT}"
export BRACKEN="${BRACKEN}"
export KRAKEN_DB="${DB}"
export BRACKEN_READ_LEN="150"
EOF


# Also place canonical environment file where later scripts
# can find it without caring whether preparation needed recovery.

cp \
    "${OUT}/97B0R_environment.sh" \
    "${BASE}/kraken_bracken_environment.sh"


# ============================================================
# 11. Recovery summary
# ============================================================

{
    printf "metric\tvalue\n"
    printf "kraken2_executable\tPASS\n"
    printf "bracken_executable\tPASS\n"
    printf "database_core\tPASS\n"
    printf "bracken_150mer\tPASS\n"
    printf "previous_inspect_output\tPASS\n"
    printf "database_redownloaded\tNO\n"
    printf "database_reextracted\tNO\n"
    printf "kraken_inspect_rerun\tNO\n"
    printf "original_job\t132825\n"
    printf "original_job_state\tFAILED\n"
    printf "recovery_status\tVALIDATED_EXISTING_ARTIFACTS\n"

} > "${OUT}/97B0R_recovery_summary.tsv"


# ============================================================
# 12. Marker
# ============================================================

printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/97B0R_COMPLETE.ok"


echo
echo "============================================================"
echo "RECOVERY SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/97B0R_recovery_summary.tsv"

echo
echo "============================================================"
echo "97B0R COMPLETE"
echo "NO DATABASE REDOWNLOAD"
echo "NO DATABASE REEXTRACTION"
echo "NO KRAKEN INSPECT RERUN"
echo "ES SEGURO SALIR"
echo "============================================================"

