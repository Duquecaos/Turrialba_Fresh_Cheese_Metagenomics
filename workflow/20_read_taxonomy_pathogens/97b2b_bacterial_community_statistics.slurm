#!/bin/bash

#SBATCH --job-name=97b2b_bact_stats
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=97b2b_bact_stats_%j.out
#SBATCH --error=97b2b_bact_stats_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

SCRIPT="${SCRIPT_DIR}/97b2b_bacterial_community_statistics.R"

OUT="${ROOT}/98_read_taxonomy/97B2B_bacterial_community"

RENV="${ROOT}/98_read_taxonomy/R_environment.sh"

# ============================================================
# 1. Load validated R environment
# ============================================================

if [[ ! -s "$RENV" ]]
then
    echo "ERROR: R environment file missing: $RENV" >&2
    exit 10
fi

source "$RENV"

echo "============================================================"
echo "R ENVIRONMENT"
echo "============================================================"

echo "R=$(command -v R)"
echo "Rscript=$(command -v Rscript)"
echo "R_LIBS_USER=${R_LIBS_USER}"

R --version | head -n 1

echo
echo "Required packages:"

Rscript - <<'RS'
pkgs <- c(
    "vegan",
    "permute",
    "ggplot2"
)

for (p in pkgs) {

    if (!requireNamespace(
        p,
        quietly = TRUE
    )) {
        stop(
            paste(
                "Missing package:",
                p
            )
        )
    }

    cat(
        p,
        "\t",
        as.character(
            packageVersion(p)
        ),
        "\n",
        sep = ""
    )
}
RS

# ============================================================
# 2. Clean previous failed output
# ============================================================

rm -rf "$OUT"

# ============================================================
# 3. Run analysis
# ============================================================

Rscript \
    "$SCRIPT" \
    "$ROOT"

# ============================================================
# 4. Validate major outputs
# ============================================================

[[ -s "${OUT}/97B2B_global_summary.tsv" ]] || {
    echo "ERROR: global summary missing" >&2
    exit 20
}

[[ -s "${OUT}/97B2B_PERMANOVA_time_primary.tsv" ]] || {
    echo "ERROR: time PERMANOVA missing" >&2
    exit 21
}

[[ -s "${OUT}/97B2B_PERMANOVA_producer_subject_level.tsv" ]] || {
    echo "ERROR: producer PERMANOVA missing" >&2
    exit 22
}

[[ -s "${OUT}/97B2B_alpha_diversity_primary.tsv" ]] || {
    echo "ERROR: alpha diversity output missing" >&2
    exit 23
}

[[ -s "${OUT}/97B2B_PCoA_primary_coordinates.tsv" ]] || {
    echo "ERROR: PCoA coordinates missing" >&2
    exit 24
}

# ============================================================
# 5. Completion marker
# ============================================================

printf \
"COMPLETED\t%s\tjob=%s\tR=4.5.0\tvegan=2.7.6\tpermute=0.9.10\tggplot2=4.0.3\n" \
"$(date -Iseconds)" \
"$SLURM_JOB_ID" \
> "${OUT}/97B2B_COMPLETE.ok"

echo
echo "============================================================"
echo "97B2B GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/97B2B_global_summary.tsv"

echo
echo "============================================================"
echo "97B2B FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"

