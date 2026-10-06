#!/bin/bash

set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"
OUT="${ROOT}/98_read_taxonomy/97B2B0R_R_audit"

mkdir -p "$OUT"

REPORT="${OUT}/97B2B0R_R_environment.tsv"

echo "============================================================"
echo "97B2B0R - R ENVIRONMENT AUDIT"
echo "Inicio: $(date -Iseconds)"
echo "============================================================"

# ============================================================
# 1. Load R module
# ============================================================

module purge

module load gnu15/15.1.0
module load R/4.5.0

echo "MODULE_LOAD=PASS"

echo
echo "R:"
command -v R
R --version | head -n 1

echo
echo "Rscript:"
command -v Rscript
Rscript --version


# ============================================================
# 2. Check required packages
# ============================================================

Rscript - "$REPORT" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
outfile <- args[1]

required <- c(
    "vegan",
    "permute",
    "ggplot2"
)

rows <- data.frame(
    package = required,
    installed = FALSE,
    version = NA_character_,
    stringsAsFactors = FALSE
)

for (i in seq_along(required)) {

    pkg <- required[i]

    ok <- requireNamespace(
        pkg,
        quietly = TRUE
    )

    rows$installed[i] <- ok

    if (ok) {
        rows$version[i] <- as.character(
            packageVersion(pkg)
        )
    }
}

write.table(
    rows,
    file = outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)

cat(
    paste0(
        rows$package,
        "=",
        rows$installed,
        " version=",
        rows$version,
        collapse = "\n"
    ),
    "\n"
)
RS


# ============================================================
# 3. Overall result
# ============================================================

MISSING=$(
    awk -F'\t' '
        NR > 1 && $2 != "TRUE" {
            print $1
        }
    ' "$REPORT"
)

if [[ -z "$MISSING" ]]
then

    STATUS="MODULE_R_READY"

else

    STATUS="MODULE_R_MISSING_PACKAGES"

fi


{
    printf "field\tvalue\n"
    printf "status\t%s\n" "$STATUS"
    printf "R_path\t%s\n" "$(command -v R)"
    printf "Rscript_path\t%s\n" "$(command -v Rscript)"
    printf "R_version\t%s\n" "$(R --version | head -n 1)"
    printf "compiler_module\tgnu15/15.1.0\n"
    printf "R_module\tR/4.5.0\n"

} > "${OUT}/97B2B0R_summary.tsv"


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/97B2B0R_COMPLETE.ok"


echo
echo "============================================================"
echo "SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/97B2B0R_summary.tsv"

echo
echo "============================================================"
echo "PACKAGES"
echo "============================================================"

column -t -s $'\t' \
    "$REPORT"

echo
echo "============================================================"
echo "97B2B0R COMPLETE"
echo "============================================================"

