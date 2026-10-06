\
#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

PATCHER="${SCRIPT_DIR}/patch_read_gtdb_taxonomy.py"
OUTPUT="${MDMCLEANER_PATCH:-${SCRIPT_DIR}/read_gtdb_taxonomy.py}"

# Optional repository configuration.
CONFIG="${PROJECT_CONFIG:-${REPO_ROOT}/config/config.sh}"

if [[ -r "${CONFIG}" ]]; then
    # shellcheck disable=SC1090
    source "${CONFIG}"
fi

# Container may be supplied either through:
#
#   ./prepare_mdmcleaner_patch.sh /path/to/mdmcleaner.sif
#
# or:
#
#   MDMCLEANER=/path/to/mdmcleaner.sif \
#       ./prepare_mdmcleaner_patch.sh

if [[ $# -ge 1 ]]; then
    MDMCLEANER="$1"
fi

: "${MDMCLEANER:?Set MDMCLEANER or provide the MDMcleaner 0.8.7 container as argument}"

UPSTREAM_PATH="${MDMCLEANER_UPSTREAM_SOURCE:-/usr/local/lib/python3.11/site-packages/mdmcleaner/read_gtdb_taxonomy.py}"

if [[ ! -s "${PATCHER}" ]]; then
    echo "ERROR: patcher not found: ${PATCHER}" >&2
    exit 2
fi

if [[ ! -s "${MDMCLEANER}" ]]; then
    echo "ERROR: MDMcleaner container not found: ${MDMCLEANER}" >&2
    exit 2
fi

if command -v apptainer >/dev/null 2>&1; then
    CONTAINER_RUNTIME="$(command -v apptainer)"
elif command -v singularity >/dev/null 2>&1; then
    CONTAINER_RUNTIME="$(command -v singularity)"
else
    echo "ERROR: neither apptainer nor singularity was found." >&2
    exit 2
fi

TMP_SOURCE="$(mktemp)"

cleanup() {
    rm -f "${TMP_SOURCE}"
}
trap cleanup EXIT

echo "============================================================"
echo "Preparing MDMcleaner 0.8.7 compatibility patch"
echo "============================================================"
echo "Container:      ${MDMCLEANER}"
echo "Upstream file:  ${UPSTREAM_PATH}"
echo "Output:         ${OUTPUT}"
echo

"${CONTAINER_RUNTIME}" exec \
    "${MDMCLEANER}" \
    cat "${UPSTREAM_PATH}" \
    > "${TMP_SOURCE}"

if [[ ! -s "${TMP_SOURCE}" ]]; then
    echo "ERROR: extracted upstream source is empty." >&2
    exit 2
fi

echo "Upstream SHA-256:"
sha256sum "${TMP_SOURCE}"

echo

python3 \
    "${PATCHER}" \
    "${TMP_SOURCE}" \
    "${OUTPUT}"

python3 -m py_compile \
    "${OUTPUT}"

echo
echo "Patched SHA-256:"
sha256sum "${OUTPUT}"

echo
echo "Expected thesis patched SHA-256:"
echo "e9b683ce506941e89638691ac57bd66da6620d7e2a0defb8fa528df54cae69f3"

ACTUAL_SHA="$(
    sha256sum "${OUTPUT}" |
    awk '{print $1}'
)"

EXPECTED_SHA="e9b683ce506941e89638691ac57bd66da6620d7e2a0defb8fa528df54cae69f3"

echo

if [[ "${ACTUAL_SHA}" == "${EXPECTED_SHA}" ]]; then
    echo "MDMCLEANER_PATCH_HISTORICAL_MATCH=YES"
else
    echo "MDMCLEANER_PATCH_HISTORICAL_MATCH=NO"
    echo
    echo "WARNING:"
    echo "The transformation completed, but the resulting file is not"
    echo "byte-identical to the MDMcleaner source used in the thesis."
    echo "This may indicate a different upstream MDMcleaner build."
fi

echo
echo "MDMCLEANER_PATCH_PREPARATION=PASS"
