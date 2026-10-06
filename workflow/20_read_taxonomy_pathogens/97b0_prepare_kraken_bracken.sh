#!/bin/bash

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

BASE="${ROOT}/98_read_taxonomy"
SOFT="${BASE}/software"
SRC="${SOFT}/src"
DBROOT="${ROOT}/databases"
DB="${DBROOT}/kraken2_pluspf_20260626"
DOWNLOADS="${DBROOT}/downloads"

KRAKEN_VERSION="2.17.2"
BRACKEN_VERSION="3.1"

KRAKEN_INSTALL="${SOFT}/kraken2-${KRAKEN_VERSION}"
BRACKEN_INSTALL="${SOFT}/Bracken-${BRACKEN_VERSION}"

KRAKEN_SRC_ARCHIVE="${DOWNLOADS}/kraken2-v${KRAKEN_VERSION}.tar.gz"
BRACKEN_SRC_ARCHIVE="${DOWNLOADS}/bracken-v${BRACKEN_VERSION}.tar.gz"

DB_ARCHIVE="${DOWNLOADS}/k2_pluspf_20260626.tar.gz"
DB_MD5="${DOWNLOADS}/pluspf_20260626.md5"

KRAKEN_URL="https://github.com/DerrickWood/kraken2/archive/refs/tags/${KRAKEN_VERSION}.tar.gz"
BRACKEN_URL="https://github.com/jenniferlu717/Bracken/archive/refs/tags/v${BRACKEN_VERSION}.tar.gz"

DB_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20260626.tar.gz"
MD5_URL="https://genome-idx.s3.amazonaws.com/kraken/pluspf_20260626/pluspf.md5"

OUT="${BASE}/97B0_prepare"

mkdir -p \
    "$OUT" \
    "$SOFT" \
    "$SRC" \
    "$DBROOT" \
    "$DOWNLOADS"

echo "============================================================"
echo "97B0 - PREPARE KRAKEN2 + BRACKEN + PLUSPF"
echo "Inicio: $(date -Iseconds)"
echo "============================================================"

# ============================================================
# 1. Basic requirements
# ============================================================

for TOOL in \
    curl \
    tar \
    make \
    g++ \
    perl \
    python3
do

    if ! command -v "$TOOL" >/dev/null 2>&1
    then
        echo "ERROR: required tool missing: $TOOL" >&2
        exit 10
    fi

done

echo "Basic build tools: PASS"

# ============================================================
# 2. External connectivity
# ============================================================

echo
echo "Testing external HTTPS..."

if ! curl \
    -L \
    --connect-timeout 20 \
    --max-time 60 \
    --retry 2 \
    --retry-delay 5 \
    -fsSI \
    "$DB_URL" \
    >/dev/null
then
    echo "ERROR: compute node cannot reach database URL" >&2
    echo "NETWORK_ACCESS_FROM_COMPUTE=FAIL" >&2
    exit 11
fi

echo "NETWORK_ACCESS_FROM_COMPUTE=PASS"

# ============================================================
# 3. Kraken2 source + installation
# ============================================================

if [[ ! -x "${KRAKEN_INSTALL}/kraken2" ]]
then

    echo
    echo "Downloading Kraken2 v${KRAKEN_VERSION}..."

    curl \
        -L \
        --fail \
        --retry 5 \
        --retry-delay 10 \
        -C - \
        -o "$KRAKEN_SRC_ARCHIVE" \
        "$KRAKEN_URL"

    rm -rf \
        "${SRC}/kraken2-${KRAKEN_VERSION}"

    tar -xzf \
        "$KRAKEN_SRC_ARCHIVE" \
        -C "$SRC"

    rm -rf \
        "$KRAKEN_INSTALL"

    mkdir -p \
        "$KRAKEN_INSTALL"

    cd \
        "${SRC}/kraken2-${KRAKEN_VERSION}"

    ./install_kraken2.sh \
        "$KRAKEN_INSTALL"

fi


if [[ ! -x "${KRAKEN_INSTALL}/kraken2" ]]
then
    echo "ERROR: Kraken2 installation failed" >&2
    exit 20
fi


echo
echo "Kraken2 installed:"
"${KRAKEN_INSTALL}/kraken2" \
    --version


# ============================================================
# 4. Bracken source + installation
# ============================================================

if [[ ! -x "${BRACKEN_INSTALL}/bracken" ]]
then

    echo
    echo "Downloading Bracken v${BRACKEN_VERSION}..."

    curl \
        -L \
        --fail \
        --retry 5 \
        --retry-delay 10 \
        -C - \
        -o "$BRACKEN_SRC_ARCHIVE" \
        "$BRACKEN_URL"

    rm -rf \
        "${SRC}/Bracken-${BRACKEN_VERSION}"

    tar -xzf \
        "$BRACKEN_SRC_ARCHIVE" \
        -C "$SRC"

    rm -rf \
        "$BRACKEN_INSTALL"

    cp -a \
        "${SRC}/Bracken-${BRACKEN_VERSION}" \
        "$BRACKEN_INSTALL"

    export PATH="${KRAKEN_INSTALL}:$PATH"

    cd \
        "$BRACKEN_INSTALL"

    bash \
        install_bracken.sh

fi


if [[ ! -x "${BRACKEN_INSTALL}/bracken" ]]
then
    echo "ERROR: Bracken installation failed" >&2
    exit 21
fi


echo
echo "Bracken installed:"

"${BRACKEN_INSTALL}/bracken" \
    -h \
    2>&1 \
| head -n 5 \
|| true


# ============================================================
# 5. Download PlusPF June 2026
# ============================================================

echo
echo "============================================================"
echo "DATABASE"
echo "============================================================"

if [[ ! -s "$DB_ARCHIVE" ]]
then

    echo "Downloading PlusPF archive..."
    echo "Expected compressed size: approximately 84.8 GB"

    curl \
        -L \
        --fail \
        --retry 20 \
        --retry-delay 20 \
        -C - \
        -o "$DB_ARCHIVE" \
        "$DB_URL"

else

    echo "Database archive already present:"
    ls -lh \
        "$DB_ARCHIVE"

fi


# ============================================================
# 6. Download official checksum manifest
# ============================================================

curl \
    -L \
    --fail \
    --retry 5 \
    --retry-delay 5 \
    -o "$DB_MD5" \
    "$MD5_URL"


# ============================================================
# 7. Extract database
# ============================================================

if [[ ! -s "${DB}/hash.k2d" ]]
then

    echo
    echo "Extracting database..."

    rm -rf "$DB"

    mkdir -p "$DB"

    tar -xzf \
        "$DB_ARCHIVE" \
        -C "$DB"

fi


# ============================================================
# 8. Validate core Kraken2 files
# ============================================================

for F in \
    hash.k2d \
    opts.k2d \
    taxo.k2d
do

    if [[ ! -s "${DB}/${F}" ]]
    then
        echo "ERROR: Kraken DB core file missing: ${F}" >&2
        exit 30
    fi

done


# ============================================================
# 9. Validate Bracken 150-mer distribution
# ============================================================

BRACKEN150=$(
    find "$DB" \
        -maxdepth 1 \
        -type f \
        -name '*150mers.kmer_distrib' \
        -print \
    | head -n 1
)


if [[ -z "$BRACKEN150" ]]
then

    echo "ERROR: 150-mer Bracken distribution missing" >&2

    echo "Available Bracken distributions:" >&2

    find "$DB" \
        -maxdepth 1 \
        -type f \
        -name '*mers.kmer_distrib' \
        -print \
        >&2

    exit 31

fi


echo
echo "Bracken 150-mer DB:"
ls -lh \
    "$BRACKEN150"


# ============================================================
# 10. Kraken inspect sanity
# ============================================================

echo
echo "Running kraken2-inspect sanity check..."

"${KRAKEN_INSTALL}/kraken2-inspect" \
    --db "$DB" \
| head -n 25 \
> "${OUT}/97B0_kraken2_inspect_head.txt"


if [[ ! -s "${OUT}/97B0_kraken2_inspect_head.txt" ]]
then
    echo "ERROR: kraken2-inspect produced no output" >&2
    exit 32
fi


# ============================================================
# 11. Versions and database metadata
# ============================================================

{
    echo -e "field\tvalue"

    echo -e \
        "kraken2_version\t$(
            "${KRAKEN_INSTALL}/kraken2" \
                --version \
            2>&1 \
            | head -n 1
        )"

    echo -e \
        "bracken_requested_version\t${BRACKEN_VERSION}"

    echo -e \
        "database\tPlusPF"

    echo -e \
        "database_release\t2026-06-26"

    echo -e \
        "database_path\t${DB}"

    echo -e \
        "database_archive\t${DB_ARCHIVE}"

    echo -e \
        "bracken_read_distribution\t150"

    echo -e \
        "shotgun_read_length_bp\t151"

    echo -e \
        "bracken_length_interpretation\tclosest_prebuilt_distribution_150bp_for_151bp_reads"

} > "${OUT}/97B0_versions_database.tsv"


# ============================================================
# 12. Disk usage
# ============================================================

{
    echo -e "item\tsize"

    echo -e \
        "database_archive\t$(
            du -sh "$DB_ARCHIVE" \
            | awk '{print $1}'
        )"

    echo -e \
        "database_extracted\t$(
            du -sh "$DB" \
            | awk '{print $1}'
        )"

    echo -e \
        "kraken2_install\t$(
            du -sh "$KRAKEN_INSTALL" \
            | awk '{print $1}'
        )"

    echo -e \
        "bracken_install\t$(
            du -sh "$BRACKEN_INSTALL" \
            | awk '{print $1}'
        )"

} > "${OUT}/97B0_disk_usage.tsv"


# ============================================================
# 13. Environment file for next steps
# ============================================================

cat > "${OUT}/97B0_environment.sh" <<EOF
export KRAKEN2="${KRAKEN_INSTALL}/kraken2"
export KRAKEN2_INSPECT="${KRAKEN_INSTALL}/kraken2-inspect"
export BRACKEN="${BRACKEN_INSTALL}/bracken"
export KRAKEN_DB="${DB}"
export BRACKEN_READ_LEN="150"
EOF


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/97B0_COMPLETE.ok"


echo
echo "============================================================"
echo "97B0 COMPLETE"
echo "============================================================"

cat \
    "${OUT}/97B0_versions_database.tsv"

echo

cat \
    "${OUT}/97B0_disk_usage.tsv"

echo

echo "KRAKEN_DATABASE=PASS"
echo "BRACKEN_150MER=PASS"
echo "ES SEGURO SALIR"

