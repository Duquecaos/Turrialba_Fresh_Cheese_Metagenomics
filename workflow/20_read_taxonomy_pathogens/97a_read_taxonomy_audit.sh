#!/bin/bash

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"
OUT="${ROOT}/98_read_taxonomy/97A_audit"

mkdir -p "$OUT"

REPORT="${OUT}/97A_audit.txt"

: > "$REPORT"

{
    echo "============================================================"
    echo "97A - READ-BASED TAXONOMY AUDIT"
    echo "Date: $(date -Iseconds)"
    echo "User: $USER"
    echo "============================================================"

    echo
    echo "============================================================"
    echo "1. HOST TOOLS"
    echo "============================================================"

    for TOOL in \
        kraken2 \
        bracken \
        metaphlan \
        kaiju \
        centrifuge
    do
        printf "%-15s\t" "$TOOL"

        if command -v "$TOOL" >/dev/null 2>&1
        then
            command -v "$TOOL"
        else
            echo "NOT_FOUND"
        fi
    done

    echo
    echo "============================================================"
    echo "2. BIO CONTAINERS"
    echo "============================================================"

    CONTAINER_ROOT="${CONTAINER_ROOT:-/opt/ohpc/pub/containers/BIO}"

    if [[ -d "$CONTAINER_ROOT" ]]
    then
        find "$CONTAINER_ROOT" \
            -maxdepth 2 \
            -type f \
            \( \
                -iname '*kraken*' \
                -o \
                -iname '*bracken*' \
                -o \
                -iname '*metaphlan*' \
                -o \
                -iname '*kaiju*' \
                -o \
                -iname '*centrifuge*' \
            \) \
            -print \
        | sort
    else
        echo "CONTAINER_ROOT_NOT_FOUND"
    fi

    echo
    echo "============================================================"
    echo "3. PROJECT DATABASE CANDIDATES"
    echo "============================================================"

    PROJECT_DB="${ROOT}/databases"

    if [[ -d "$PROJECT_DB" ]]
    then
        find "$PROJECT_DB" \
            -maxdepth 4 \
            \( \
                -type d \
                -o \
                -type f \
            \) \
            \( \
                -iname '*kraken*' \
                -o \
                -iname '*bracken*' \
                -o \
                -iname '*pluspf*' \
                -o \
                -iname '*standard*' \
                -o \
                -iname '*metaphlan*' \
                -o \
                -iname '*chocophlan*' \
                -o \
                -iname '*kaiju*' \
                -o \
                -iname '*centrifuge*' \
            \) \
            -print \
        | sort
    else
        echo "PROJECT_DATABASE_DIR_NOT_FOUND"
    fi

    echo
    echo "============================================================"
    echo "4. SHARED DATABASE ROOTS"
    echo "============================================================"

    for D in \
        "${OHPC_PUBLIC_ROOT:-/opt/ohpc/pub}/databases" \
        "${OHPC_PUBLIC_ROOT:-/opt/ohpc/pub}/database" \
        "${OHPC_PUBLIC_ROOT:-/opt/ohpc/pub}/db" \
        /scratch/global/databases
    do
        if [[ -d "$D" ]]
        then
            echo "--- $D ---"

            find "$D" \
                -maxdepth 3 \
                \( \
                    -type d \
                    -o \
                    -type f \
                \) \
                \( \
                    -iname '*kraken*' \
                    -o \
                    -iname '*bracken*' \
                    -o \
                    -iname '*pluspf*' \
                    -o \
                    -iname '*metaphlan*' \
                    -o \
                    -iname '*kaiju*' \
                    -o \
                    -iname '*centrifuge*' \
                \) \
                -print \
                2>/dev/null \
            | head -n 200
        fi
    done

    echo
    echo "============================================================"
    echo "5. HOST-FREE FASTQ AUDIT"
    echo "============================================================"

    mapfile -t R1S < <(
        find "${ROOT}/04_host_removed" \
            -type f \
            -name '*_1.hostfree.fastq.gz' \
        | sort
    )

    echo "R1_files=${#R1S[@]}"

    PAIRS=0
    MISSING_R2=0

    for R1 in "${R1S[@]}"
    do
        R2="${R1/_1.hostfree.fastq.gz/_2.hostfree.fastq.gz}"

        if [[ -s "$R2" ]]
        then
            PAIRS=$((PAIRS + 1))
        else
            echo "MISSING_R2: $R2"
            MISSING_R2=$((MISSING_R2 + 1))
        fi
    done

    echo "complete_pairs=$PAIRS"
    echo "missing_R2=$MISSING_R2"

    echo
    echo "============================================================"
    echo "6. SAMPLE LIST + FILE SIZES"
    echo "============================================================"

    printf "sample\tR1_GB\tR2_GB\ttotal_GB\n"

    for R1 in "${R1S[@]}"
    do
        R2="${R1/_1.hostfree.fastq.gz/_2.hostfree.fastq.gz}"

        SAMPLE="$(
            basename "$R1" _1.hostfree.fastq.gz
        )"

        if [[ ! -s "$R2" ]]
        then
            continue
        fi

        B1=$(stat -c '%s' "$R1")
        B2=$(stat -c '%s' "$R2")

        python3 - \
            "$SAMPLE" \
            "$B1" \
            "$B2" <<'PY'
import sys

sample = sys.argv[1]
b1 = int(sys.argv[2])
b2 = int(sys.argv[3])

gb = 1024 ** 3

print(
    f"{sample}\t"
    f"{b1/gb:.3f}\t"
    f"{b2/gb:.3f}\t"
    f"{(b1+b2)/gb:.3f}"
)
PY
    done

    echo
    echo "============================================================"
    echo "7. TOTAL HOST-FREE SIZE"
    echo "============================================================"

    python3 - "${ROOT}/04_host_removed" <<'PY'
from pathlib import Path
import sys

root = Path(sys.argv[1])

files = sorted(
    root.rglob("*.hostfree.fastq.gz")
)

total = sum(
    f.stat().st_size
    for f in files
)

print(f"files={len(files)}")
print(f"total_GB={total / 1024**3:.3f}")
print(f"total_TB={total / 1024**4:.3f}")
PY

    echo
    echo "============================================================"
    echo "8. FILESYSTEM SPACE"
    echo "============================================================"

    df -h \
        "$ROOT" \
        /scratch/global \
        2>/dev/null \
    | awk '!seen[$1]++'

    echo
    echo "============================================================"
    echo "97A AUDIT COMPLETE"
    echo "============================================================"

} | tee "$REPORT"

printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/97A_COMPLETE.ok"
