#!/bin/bash

set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

BASE="${ROOT}/98_read_taxonomy"
SOFT="${BASE}/software"

RLIB="${SOFT}/R_libs/4.5"

OUT="${BASE}/97B2B0I_R_packages"

mkdir -p \
    "$RLIB" \
    "$OUT"


echo "============================================================"
echo "97B2B0I - INSTALL PRIVATE R PACKAGES"
echo "Inicio: $(date -Iseconds)"
echo "============================================================"


# ============================================================
# 1. Load institutional R
# ============================================================

module purge
module load gnu15/15.1.0
module load R/4.5.0

export R_LIBS_USER="$RLIB"


echo "R=$(command -v R)"
echo "Rscript=$(command -v Rscript)"
echo "R_LIBS_USER=$R_LIBS_USER"

R --version | head -n 1


# ============================================================
# 2. Install required packages
# ============================================================

Rscript - \
    "$RLIB" \
    "${SLURM_CPUS_PER_TASK:-4}" <<'RS'
args <- commandArgs(trailingOnly = TRUE)

lib <- args[1]
ncpus <- as.integer(args[2])

dir.create(
    lib,
    recursive = TRUE,
    showWarnings = FALSE
)

.libPaths(
    c(
        lib,
        .libPaths()
    )
)

options(
    repos = c(
        CRAN = "https://cloud.r-project.org"
    )
)

options(
    timeout = 1200
)

required <- c(
    "vegan",
    "permute",
    "ggplot2"
)

cat(
    "Library paths:\n"
)

print(
    .libPaths()
)


for (pkg in required) {

    if (
        requireNamespace(
            pkg,
            quietly = TRUE
        )
    ) {

        cat(
            pkg,
            "already installed:",
            as.character(
                packageVersion(pkg)
            ),
            "\n"
        )

        next
    }


    cat(
        "\nInstalling:",
        pkg,
        "\n"
    )


    install.packages(
        pkg,
        lib = lib,
        dependencies = c(
            "Depends",
            "Imports",
            "LinkingTo"
        ),
        Ncpus = ncpus
    )


    if (
        !requireNamespace(
            pkg,
            quietly = TRUE
        )
    ) {

        stop(
            paste(
                "Installation failed:",
                pkg
            )
        )
    }


    cat(
        pkg,
        "installed:",
        as.character(
            packageVersion(pkg)
        ),
        "\n"
    )
}
RS


# ============================================================
# 3. Validate installed packages
# ============================================================

Rscript - \
    "$OUT/97B2B0I_package_versions.tsv" <<'RS'
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
    library_path = NA_character_,
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

        rows$library_path[i] <- dirname(
            find.package(pkg)
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


if (
    !all(
        rows$installed
    )
) {

    stop(
        "One or more required packages are still missing"
    )
}
RS


# ============================================================
# 4. Environment file for later R analyses
# ============================================================

cat > "${OUT}/97B2B0I_R_environment.sh" <<EOF
module purge
module load gnu15/15.1.0
module load R/4.5.0
export R_LIBS_USER="${RLIB}"
EOF


cp \
    "${OUT}/97B2B0I_R_environment.sh" \
    "${BASE}/R_environment.sh"


# ============================================================
# 5. Reproducibility snapshot
# ============================================================

Rscript - \
    "$OUT/97B2B0I_sessionInfo.txt" <<'RS'
args <- commandArgs(trailingOnly = TRUE)

outfile <- args[1]

sink(
    outfile
)

cat(
    "===== R sessionInfo =====\n\n"
)

print(
    sessionInfo()
)

cat(
    "\n===== Required packages =====\n\n"
)

for (
    pkg
    in c(
        "vegan",
        "permute",
        "ggplot2"
    )
) {

    cat(
        pkg,
        "\t",
        as.character(
            packageVersion(pkg)
        ),
        "\t",
        find.package(pkg),
        "\n",
        sep = ""
    )
}

sink()
RS


# ============================================================
# 6. Final checks
# ============================================================

for F in \
    "${OUT}/97B2B0I_package_versions.tsv" \
    "${OUT}/97B2B0I_R_environment.sh" \
    "${OUT}/97B2B0I_sessionInfo.txt"
do

    if [[ ! -s "$F" ]]
    then
        echo "ERROR: missing output: $F" >&2
        exit 20
    fi

done


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/97B2B0I_COMPLETE.ok"


echo
echo "============================================================"
echo "PACKAGE VERSIONS"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/97B2B0I_package_versions.tsv"


echo
echo "============================================================"
echo "97B2B0I FINALIZO CORRECTAMENTE"
echo "PRIVATE R LIBRARY READY"
echo "ES SEGURO SALIR"
echo "============================================================"

