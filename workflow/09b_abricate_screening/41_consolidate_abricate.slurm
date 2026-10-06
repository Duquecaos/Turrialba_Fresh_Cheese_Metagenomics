#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00

#SBATCH --job-name=abricate_consolidate
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

HITS="${ROOT}/24_abricate_rep18/abricate_high_confidence_90_90.tsv"
MASTER="${ROOT}/22_final_representative_mags/representative_mags_master.tsv"

OUT="${ROOT}/25_abricate_consolidated"

mkdir -p "${OUT}"

echo "============================================================"
echo "CONSOLIDACIÓN ABRICATE ≥90/90"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

python3 - \
"${HITS}" \
"${MASTER}" \
"${OUT}" <<'PY'

import csv
import sys
from pathlib import Path
from collections import defaultdict

hits_file = Path(sys.argv[1])
master_file = Path(sys.argv[2])
outdir = Path(sys.argv[3])

# ------------------------------------------------------------
# Taxonomía
# ------------------------------------------------------------

taxonomy = {}

with open(master_file, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for r in reader:
        taxonomy[r["representative_MAG"]] = {
            "genus": r["genus"],
            "species": r["species"],
            "classification": r["classification"],
        }

# ------------------------------------------------------------
# Hits
# ------------------------------------------------------------

hits = []

with open(hits_file, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    for i, r in enumerate(reader):

        start = int(r["START"])
        end = int(r["END"])

        if start > end:
            start, end = end, start

        hits.append({
            "id": i,
            "MAG": r["MAG"],
            "contig": r["SEQUENCE"],
            "start": start,
            "end": end,
            "strand": r.get("STRAND", ""),
            "gene": r.get("GENE", ""),
            "db": r.get("source_db", ""),
            "type": r.get("evidence_type", ""),
            "identity": float(r.get("%IDENTITY", 0) or 0),
            "coverage": float(r.get("%COVERAGE", 0) or 0),
            "accession": r.get("ACCESSION", ""),
            "product": r.get("PRODUCT", ""),
            "resistance": r.get("RESISTANCE", ""),
        })

if not hits:
    raise SystemExit("ERROR: no hay hits ≥90/90")

# ------------------------------------------------------------
# Union-Find
# ------------------------------------------------------------

parent = list(range(len(hits)))

def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[rb] = ra

def overlap_fraction(a, b):

    overlap = max(
        0,
        min(a["end"], b["end"]) -
        max(a["start"], b["start"]) + 1
    )

    if overlap == 0:
        return 0.0

    la = a["end"] - a["start"] + 1
    lb = b["end"] - b["start"] + 1

    return overlap / min(la, lb)

# Solo comparar hits del mismo MAG, contig y tipo de evidencia
groups = defaultdict(list)

for h in hits:
    key = (h["MAG"], h["contig"], h["type"])
    groups[key].append(h)

for group in groups.values():

    for i in range(len(group)):
        for j in range(i + 1, len(group)):

            if overlap_fraction(group[i], group[j]) >= 0.80:
                union(group[i]["id"], group[j]["id"])

# ------------------------------------------------------------
# Construir loci
# ------------------------------------------------------------

loci = defaultdict(list)

for h in hits:
    loci[find(h["id"])].append(h)

rows = []

for locus_number, members in enumerate(
    sorted(
        loci.values(),
        key=lambda x: (
            x[0]["MAG"],
            x[0]["contig"],
            min(y["start"] for y in x)
        )
    ),
    start=1
):

    mag = members[0]["MAG"]
    evidence_type = members[0]["type"]

    start = min(x["start"] for x in members)
    end = max(x["end"] for x in members)

    genes = sorted({
        x["gene"] for x in members
        if x["gene"]
    })

    dbs = sorted({
        x["db"] for x in members
        if x["db"]
    })

    accessions = sorted({
        x["accession"] for x in members
        if x["accession"]
    })

    resistance = sorted({
        x["resistance"] for x in members
        if x["resistance"]
    })

    products = sorted({
        x["product"] for x in members
        if x["product"]
    })

    tx = taxonomy.get(mag, {})

    rows.append({
        "locus_id": f"LOCUS_{locus_number:04d}",
        "MAG": mag,
        "genus": tx.get("genus", ""),
        "species": tx.get("species", ""),
        "evidence_type": evidence_type,
        "contig": members[0]["contig"],
        "start": start,
        "end": end,
        "length_bp": end - start + 1,
        "strand": ";".join(sorted({
            x["strand"] for x in members
            if x["strand"]
        })),
        "gene_labels": ";".join(genes),
        "databases": ";".join(dbs),
        "n_databases": len(dbs),
        "n_raw_hits": len(members),
        "max_identity": max(x["identity"] for x in members),
        "max_coverage": max(x["coverage"] for x in members),
        "accessions": ";".join(accessions),
        "resistance": ";".join(resistance),
        "products": " | ".join(products),
    })

fields = [
    "locus_id",
    "MAG",
    "genus",
    "species",
    "evidence_type",
    "contig",
    "start",
    "end",
    "length_bp",
    "strand",
    "gene_labels",
    "databases",
    "n_databases",
    "n_raw_hits",
    "max_identity",
    "max_coverage",
    "accessions",
    "resistance",
    "products",
]

all_file = outdir / "abricate_nonredundant_loci_90_90.tsv"

with open(all_file, "w", newline="") as out:
    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(rows)

# ------------------------------------------------------------
# Separar por tipo
# ------------------------------------------------------------

for category, filename in [
    ("AMR", "amr_loci_90_90.tsv"),
    ("virulence", "virulence_loci_90_90.tsv"),
    ("plasmid_replicon", "plasmid_replicon_loci_90_90.tsv"),
]:

    subset = [
        r for r in rows
        if r["evidence_type"] == category
    ]

    with open(outdir / filename, "w", newline="") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(subset)

# ------------------------------------------------------------
# Resumen por MAG
# ------------------------------------------------------------

summary = defaultdict(lambda: defaultdict(int))

for r in rows:
    summary[r["MAG"]][r["evidence_type"]] += 1

summary_file = outdir / "abricate_nonredundant_summary_by_MAG.tsv"

with open(summary_file, "w") as out:
    out.write(
        "MAG\tgenus\tspecies\t"
        "AMR_loci\tvirulence_loci\tplasmid_replicon_loci\n"
    )

    for mag in sorted(taxonomy):

        tx = taxonomy[mag]

        out.write(
            f"{mag}\t"
            f"{tx['genus']}\t"
            f"{tx['species']}\t"
            f"{summary[mag]['AMR']}\t"
            f"{summary[mag]['virulence']}\t"
            f"{summary[mag]['plasmid_replicon']}\n"
        )

print(f"Hits brutos ≥90/90: {len(hits)}")
print(f"Loci no redundantes: {len(rows)}")

for category in ["AMR", "virulence", "plasmid_replicon"]:
    n = sum(
        1 for r in rows
        if r["evidence_type"] == category
    )
    print(f"{category}: {n}")

PY

echo
echo "===== RESUMEN POR MAG ====="

column -t -s $'\t' \
"${OUT}/abricate_nonredundant_summary_by_MAG.tsv"

echo
echo "===== LOCI AMR CON APOYO ≥2 BASES ====="

awk -F '\t' '
NR == 1 || ($5 == "AMR" && $13 >= 2)
' "${OUT}/abricate_nonredundant_loci_90_90.tsv" \
| column -t -s $'\t'

echo
echo "============================================================"
echo "ABRICATE CONSOLIDADO COMPLETADO"
echo "Fin: $(date)"
echo "============================================================"
