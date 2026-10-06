#!/bin/bash

#SBATCH --job-name=93e2c0_AH_audit
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=93e2c0_AH_audit_%j.out
#SBATCH --error=93e2c0_AH_audit_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

BASE="${ROOT}/95_mobilome/93E2_alienhunter_final18"
PRED="${BASE}/predictions"

OUT="${ROOT}/95_mobilome/93E2C0_alienhunter_audit"

mkdir -p "$OUT"

echo "============================================================"
echo "93E2C0 - GLOBAL ALIENHUNTER AUDIT"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

python3 - \
    "$BASE" \
    "$PRED" \
    "$OUT" <<'PY'

import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

base = Path(sys.argv[1])
pred = Path(sys.argv[2])
out = Path(sys.argv[3])

manifest = base / "93E2A_manifest.tsv"

with manifest.open() as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

if len(rows) != 271:
    raise RuntimeError(
        f"Manifest esperado=271 observado={len(rows)}"
    )

n_marker = 0
n_output = 0
n_summary = 0

feature_keys = Counter()
location_patterns = Counter()
qualifier_keys = Counter()

examples = []
missing = []

for r in rows:

    iid = r["input_id"]

    d = pred / iid

    marker = d / "93E2B_COMPLETE.ok"
    embl = d / f"{iid}.alienhunter.embl"
    summary = d / f"{iid}.task_summary.tsv"

    n_marker += marker.is_file()
    n_output += embl.is_file()
    n_summary += summary.is_file()

    if not (
        marker.is_file()
        and embl.is_file()
        and summary.is_file()
    ):
        missing.append({
            "array_index": r["array_index"],
            "input_id": iid,
            "MAG": r["MAG"],
            "contig": r["contig"],
            "marker": int(marker.is_file()),
            "output": int(embl.is_file()),
            "summary": int(summary.is_file()),
        })
        continue

    text = embl.read_text(
        errors="replace"
    )

    for line in text.splitlines():

        if not line.startswith("FT"):
            continue

        # EMBL feature line:
        # FT   feature_key        location
        m = re.match(
            r"^FT\s{3}(\S+)\s+(.+?)\s*$",
            line
        )

        if m:
            key = m.group(1)
            loc = m.group(2).strip()

            feature_keys[key] += 1

            if re.fullmatch(r"\d+\.\.\d+", loc):
                location_patterns["start..end"] += 1
            elif re.fullmatch(r"complement\(\d+\.\.\d+\)", loc):
                location_patterns["complement(start..end)"] += 1
            else:
                location_patterns[f"OTHER:{loc[:80]}"] += 1

            if len(examples) < 50:
                examples.append({
                    "input_id": iid,
                    "MAG": r["MAG"],
                    "contig": r["contig"],
                    "feature_key": key,
                    "location": loc,
                    "raw_line": line,
                })

        mq = re.match(
            r'^FT\s+/([^=\s]+)',
            line.strip()
        )

        if mq:
            qualifier_keys[mq.group(1)] += 1


with (
    out / "93E2C0_global_summary.tsv"
).open("w") as fh:

    fh.write("metric\tvalue\n")
    fh.write(f"manifest_entries\t{len(rows)}\n")
    fh.write(f"complete_markers\t{n_marker}\n")
    fh.write(f"alienhunter_outputs\t{n_output}\n")
    fh.write(f"task_summaries\t{n_summary}\n")
    fh.write(f"missing_any_component\t{len(missing)}\n")
    fh.write(f"distinct_feature_keys\t{len(feature_keys)}\n")
    fh.write(f"distinct_qualifier_keys\t{len(qualifier_keys)}\n")


with (
    out / "93E2C0_feature_keys.tsv"
).open("w") as fh:

    fh.write("feature_key\tcount\n")

    for key, n in feature_keys.most_common():
        fh.write(f"{key}\t{n}\n")


with (
    out / "93E2C0_location_patterns.tsv"
).open("w") as fh:

    fh.write("location_pattern\tcount\n")

    for key, n in location_patterns.most_common():
        fh.write(f"{key}\t{n}\n")


with (
    out / "93E2C0_qualifier_keys.tsv"
).open("w") as fh:

    fh.write("qualifier_key\tcount\n")

    for key, n in qualifier_keys.most_common():
        fh.write(f"{key}\t{n}\n")


with (
    out / "93E2C0_examples.tsv"
).open(
    "w",
    newline=""
) as fh:

    fields = [
        "input_id",
        "MAG",
        "contig",
        "feature_key",
        "location",
        "raw_line",
    ]

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(examples)


with (
    out / "93E2C0_missing.tsv"
).open(
    "w",
    newline=""
) as fh:

    fields = [
        "array_index",
        "input_id",
        "MAG",
        "contig",
        "marker",
        "output",
        "summary",
    ]

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(missing)


if len(rows) != 271:
    raise RuntimeError("Manifest no contiene 271 entradas")

if n_marker != 271:
    raise RuntimeError(
        f"Markers esperados=271 observados={n_marker}"
    )

if n_output != 271:
    raise RuntimeError(
        f"Outputs esperados=271 observados={n_output}"
    )

if n_summary != 271:
    raise RuntimeError(
        f"Summaries esperados=271 observados={n_summary}"
    )

if missing:
    raise RuntimeError(
        f"Hay {len(missing)} inputs incompletos"
    )

print("GLOBAL_271_VALIDATION=PASS")

PY


printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"${SLURM_JOB_ID}" \
> "${OUT}/93E2C0_COMPLETE.ok"


echo
echo "============================================================"
echo "GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/93E2C0_global_summary.tsv"

echo
echo "============================================================"
echo "FEATURE KEYS"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/93E2C0_feature_keys.tsv"

echo
echo "============================================================"
echo "QUALIFIER KEYS"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/93E2C0_qualifier_keys.tsv"

echo
echo "============================================================"
echo "LOCATION PATTERNS"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/93E2C0_location_patterns.tsv"

echo
echo "============================================================"
echo "93E2C0 FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"

