#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

#SBATCH --job-name=antismash_summary
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

RESULTS="${ROOT}/26_antismash_rep18/results"
OUT="${ROOT}/27_antismash_summary"

ANTISMASH="${ANTISMASH:-/opt/ohpc/pub/containers/BIO/antismash-7.1.0.sif}"

mkdir -p "${OUT}"

echo "============================================================"
echo "RESUMEN antiSMASH - 18 MAGs"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "============================================================"

if [[ ! -d "${RESULTS}" ]]; then
    echo "ERROR: falta ${RESULTS}" >&2
    exit 1
fi

NREGIONS=$(find "${RESULTS}" \
    -mindepth 2 \
    -maxdepth 2 \
    -type f \
    -name '*.region*.gbk' \
    | wc -l)

echo "Regiones GenBank encontradas: ${NREGIONS}"

if [[ "${NREGIONS}" -ne 66 ]]; then
    echo "ERROR: se esperaban 66 regiones y se encontraron ${NREGIONS}" >&2
    exit 1
fi

echo
echo "===== BIOPYTHON ====="

apptainer exec \
    --cleanenv \
    "${ANTISMASH}" \
    python3 -c 'import Bio; print("Biopython", Bio.__version__)'

echo
echo "============================================================"
echo "EXTRAYENDO REGIONES, PRODUCTOS Y GENES"
echo "============================================================"

srun --exact -n 1 -c 1 \
apptainer exec \
    --cleanenv \
    --bind "${ROOT}:/project:ro" \
    --bind "${OUT}:/out" \
    "${ANTISMASH}" \
    python3 - <<'PY'

import csv
import re
from pathlib import Path
from collections import Counter, defaultdict

from Bio import SeqIO

RESULTS = Path("/project/26_antismash_rep18/results")
MASTER = Path(
    "/project/22_final_representative_mags/"
    "representative_mags_master.tsv"
)
OUT = Path("/out")


# ============================================================
# Funciones
# ============================================================

def first(q, key, default=""):
    value = q.get(key, [])
    if isinstance(value, list):
        return ";".join(str(x) for x in value)
    return str(value) if value else default


def find_column(fieldnames, candidates):
    if not fieldnames:
        return None

    normalized = {
        x.lower().strip(): x
        for x in fieldnames
    }

    for candidate in candidates:
        if candidate.lower() in normalized:
            return normalized[candidate.lower()]

    return None


def yesno(value):
    value = str(value).strip().lower()
    return "YES" if value in {
        "true", "yes", "1", "y"
    } else "NO"


# RiPP amplio:
# deliberadamente sensible para no perder candidatos.
RIPP_TERMS = [
    "ripp",
    "lanthi",
    "lantibiotic",
    "lassopeptide",
    "lasso peptide",
    "sactipeptide",
    "sacti",
    "thiopeptide",
    "glycocin",
    "cyanobactin",
    "microviridin",
    "proteusin",
    "ranthipeptide",
    "thioamitide",
    "linear azol",
    "lap",
    "bacterioc",
    "lactococcin",
    "subtilosin",
    "thuricin",
]

# Señales más directamente relacionadas con bacteriocinas
# o péptidos antimicrobianos de alto interés.
BACTERIOCIN_TERMS = [
    "bacterioc",
    "lactococcin",
    "lantibiotic",
    "lanthipeptide",
    "glycocin",
    "subtilosin",
    "thuricin",
    "nisin",
    "pediocin",
    "enterocin",
    "leucocin",
    "sakacin",
    "carnobacteriocin",
    "carnocin",
    "plantaricin",
]


def contains_any(text, terms):
    t = text.lower()
    return any(term in t for term in terms)


# ============================================================
# Taxonomía
# ============================================================

taxonomy = {}

if MASTER.exists():
    with MASTER.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        mag_col = find_column(
            reader.fieldnames,
            [
                "representative_MAG",
                "representative_mag",
                "MAG",
                "mag",
            ],
        )

        genus_col = find_column(
            reader.fieldnames,
            ["genus", "gtdb_genus"]
        )

        species_col = find_column(
            reader.fieldnames,
            ["species", "gtdb_species"]
        )

        class_col = find_column(
            reader.fieldnames,
            [
                "classification",
                "taxonomy",
                "gtdb_taxonomy",
            ],
        )

        if mag_col:
            for r in reader:
                mag = r.get(mag_col, "").strip()

                if not mag:
                    continue

                taxonomy[mag] = {
                    "genus":
                        r.get(genus_col, "")
                        if genus_col else "",
                    "species":
                        r.get(species_col, "")
                        if species_col else "",
                    "taxonomy":
                        r.get(class_col, "")
                        if class_col else "",
                }


# ============================================================
# Directorios de MAGs
# ============================================================

mag_dirs = sorted([
    x for x in RESULTS.iterdir()
    if x.is_dir()
])

if len(mag_dirs) != 18:
    raise SystemExit(
        f"ERROR: se esperaban 18 directorios MAG; "
        f"se encontraron {len(mag_dirs)}"
    )


# ============================================================
# Analizar regiones
# ============================================================

regions = []
candidate_genes = []

for mag_dir in mag_dirs:

    mag = mag_dir.name
    tx = taxonomy.get(mag, {})

    region_files = sorted(
        mag_dir.glob("*.region*.gbk")
    )

    for gbk in region_files:

        record = SeqIO.read(str(gbk), "genbank")

        region_features = [
            f for f in record.features
            if f.type == "region"
        ]

        if not region_features:
            raise RuntimeError(
                f"No se encontró feature region en {gbk}"
            )

        region = region_features[0]
        q = region.qualifiers

        region_number = first(q, "region_number")

        if not region_number:
            m = re.search(
                r"\.region(\d+)\.gbk$",
                gbk.name
            )
            region_number = (
                m.group(1) if m else ""
            )

        products = first(q, "product")
        rules = first(q, "rules")
        contig_edge = yesno(
            first(q, "contig_edge")
        )

        start = int(region.location.start) + 1
        end = int(region.location.end)

        cds_features = [
            f for f in record.features
            if f.type == "CDS"
        ]

        biosynthetic_cds = 0
        direct_gene_signals = []
        ripp_gene_signals = []

        for cds in cds_features:

            cq = cds.qualifiers

            gene_kind = first(
                cq,
                "gene_kind"
            )

            if "biosynthetic" in gene_kind.lower():
                biosynthetic_cds += 1

            gene = first(cq, "gene")
            locus_tag = first(cq, "locus_tag")
            protein_id = first(cq, "protein_id")
            product = first(cq, "product")
            gene_functions = first(
                cq,
                "gene_functions"
            )
            sec_met = first(
                cq,
                "sec_met"
            )

            combined = " | ".join([
                gene,
                locus_tag,
                protein_id,
                product,
                gene_functions,
                sec_met,
            ])

            is_ripp_gene = contains_any(
                combined,
                RIPP_TERMS
            )

            is_bacteriocin_gene = contains_any(
                combined,
                BACTERIOCIN_TERMS
            )

            if is_ripp_gene or is_bacteriocin_gene:

                cds_start = (
                    int(cds.location.start) + 1
                )
                cds_end = int(cds.location.end)

                strand = cds.location.strand

                if strand == 1:
                    strand_s = "+"
                elif strand == -1:
                    strand_s = "-"
                else:
                    strand_s = ""

                reason = []

                if is_ripp_gene:
                    reason.append("RiPP_keyword")

                if is_bacteriocin_gene:
                    reason.append(
                        "bacteriocin_keyword"
                    )

                signal = (
                    gene or
                    product or
                    locus_tag or
                    protein_id or
                    "unlabelled_CDS"
                )

                if is_ripp_gene:
                    ripp_gene_signals.append(
                        signal
                    )

                if is_bacteriocin_gene:
                    direct_gene_signals.append(
                        signal
                    )

                candidate_genes.append({
                    "MAG": mag,
                    "genus": tx.get(
                        "genus", ""
                    ),
                    "species": tx.get(
                        "species", ""
                    ),
                    "contig": record.id,
                    "region_number":
                        region_number,
                    "region_products":
                        products,
                    "contig_edge":
                        contig_edge,
                    "CDS_start": cds_start,
                    "CDS_end": cds_end,
                    "strand": strand_s,
                    "gene": gene,
                    "locus_tag": locus_tag,
                    "protein_id": protein_id,
                    "product": product,
                    "gene_kind": gene_kind,
                    "gene_functions":
                        gene_functions,
                    "sec_met": sec_met,
                    "signal_reason":
                        ";".join(reason),
                })

        combined_region = (
            products + " | " +
            rules + " | " +
            " | ".join(ripp_gene_signals)
        )

        combined_bact = (
            products + " | " +
            rules + " | " +
            " | ".join(
                direct_gene_signals
            )
        )

        ripp_flag = (
            "YES"
            if contains_any(
                combined_region,
                RIPP_TERMS
            )
            else "NO"
        )

        bact_flag = (
            "YES"
            if contains_any(
                combined_bact,
                BACTERIOCIN_TERMS
            )
            else "NO"
        )

        regions.append({
            "MAG": mag,
            "genus": tx.get(
                "genus", ""
            ),
            "species": tx.get(
                "species", ""
            ),
            "taxonomy": tx.get(
                "taxonomy", ""
            ),
            "contig": record.id,
            "region_number":
                region_number,
            "start": start,
            "end": end,
            "length_bp": end - start + 1,
            "contig_edge":
                contig_edge,
            "products": products,
            "rules": rules,
            "n_CDS": len(cds_features),
            "n_biosynthetic_CDS":
                biosynthetic_cds,
            "RiPP_related":
                ripp_flag,
            "bacteriocin_priority":
                bact_flag,
            "RiPP_gene_signals":
                ";".join(sorted(
                    set(ripp_gene_signals)
                )),
            "bacteriocin_gene_signals":
                ";".join(sorted(
                    set(direct_gene_signals)
                )),
            "region_file": gbk.name,
        })


# ============================================================
# Validación
# ============================================================

if len(regions) != 66:
    raise SystemExit(
        f"ERROR: se extrajeron {len(regions)} "
        "regiones; se esperaban 66."
    )


# ============================================================
# Escribir inventario
# ============================================================

region_fields = [
    "MAG",
    "genus",
    "species",
    "taxonomy",
    "contig",
    "region_number",
    "start",
    "end",
    "length_bp",
    "contig_edge",
    "products",
    "rules",
    "n_CDS",
    "n_biosynthetic_CDS",
    "RiPP_related",
    "bacteriocin_priority",
    "RiPP_gene_signals",
    "bacteriocin_gene_signals",
    "region_file",
]

with (
    OUT / "antismash_bgc_inventory.tsv"
).open("w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=region_fields,
        delimiter="\t",
    )

    writer.writeheader()
    writer.writerows(regions)


# ============================================================
# RiPP candidates
# ============================================================

ripp = [
    r for r in regions
    if r["RiPP_related"] == "YES"
]

with (
    OUT / "antismash_ripp_candidates.tsv"
).open("w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=region_fields,
        delimiter="\t",
    )

    writer.writeheader()
    writer.writerows(ripp)


# ============================================================
# Bacteriocin priority
# ============================================================

bact = [
    r for r in regions
    if r["bacteriocin_priority"] == "YES"
]

with (
    OUT / "antismash_bacteriocin_priority.tsv"
).open("w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=region_fields,
        delimiter="\t",
    )

    writer.writeheader()
    writer.writerows(bact)


# ============================================================
# Candidate genes
# ============================================================

gene_fields = [
    "MAG",
    "genus",
    "species",
    "contig",
    "region_number",
    "region_products",
    "contig_edge",
    "CDS_start",
    "CDS_end",
    "strand",
    "gene",
    "locus_tag",
    "protein_id",
    "product",
    "gene_kind",
    "gene_functions",
    "sec_met",
    "signal_reason",
]

with (
    OUT / "antismash_candidate_genes.tsv"
).open("w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=gene_fields,
        delimiter="\t",
    )

    writer.writeheader()
    writer.writerows(candidate_genes)


# ============================================================
# Resumen por MAG
# ============================================================

by_mag = defaultdict(list)

for r in regions:
    by_mag[r["MAG"]].append(r)

summary_rows = []

for mag_dir in mag_dirs:

    mag = mag_dir.name
    tx = taxonomy.get(mag, {})
    rows = by_mag.get(mag, [])

    products = []

    for r in rows:
        if r["products"]:
            products.extend(
                x.strip()
                for x in r["products"].split(";")
                if x.strip()
            )

    summary_rows.append({
        "MAG": mag,
        "genus":
            tx.get("genus", ""),
        "species":
            tx.get("species", ""),
        "total_BGC_regions":
            len(rows),
        "non_edge_regions":
            sum(
                r["contig_edge"] == "NO"
                for r in rows
            ),
        "contig_edge_regions":
            sum(
                r["contig_edge"] == "YES"
                for r in rows
            ),
        "RiPP_related_regions":
            sum(
                r["RiPP_related"] == "YES"
                for r in rows
            ),
        "bacteriocin_priority_regions":
            sum(
                r["bacteriocin_priority"]
                == "YES"
                for r in rows
            ),
        "unique_products":
            ";".join(sorted(
                set(products)
            )),
    })

summary_fields = [
    "MAG",
    "genus",
    "species",
    "total_BGC_regions",
    "non_edge_regions",
    "contig_edge_regions",
    "RiPP_related_regions",
    "bacteriocin_priority_regions",
    "unique_products",
]

with (
    OUT / "antismash_summary_by_MAG.tsv"
).open("w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=summary_fields,
        delimiter="\t",
    )

    writer.writeheader()
    writer.writerows(summary_rows)


# ============================================================
# Conteo de productos
# ============================================================

product_counts = Counter()

for r in regions:

    if not r["products"]:
        product_counts["UNCLASSIFIED"] += 1
        continue

    products = [
        x.strip()
        for x in r["products"].split(";")
        if x.strip()
    ]

    for product in products:
        product_counts[product] += 1

with (
    OUT / "antismash_product_counts.tsv"
).open("w", newline="") as fh:

    writer = csv.writer(
        fh,
        delimiter="\t"
    )

    writer.writerow([
        "BGC_product",
        "n_regions"
    ])

    for product, n in sorted(
        product_counts.items(),
        key=lambda x: (-x[1], x[0])
    ):
        writer.writerow([
            product,
            n
        ])


# ============================================================
# Resultados
# ============================================================

print(
    f"Regiones totales: {len(regions)}"
)

print(
    "Regiones en borde de contig: "
    f"{sum(r['contig_edge'] == 'YES' for r in regions)}"
)

print(
    "Regiones no asociadas a borde: "
    f"{sum(r['contig_edge'] == 'NO' for r in regions)}"
)

print(
    "Regiones RiPP-related: "
    f"{len(ripp)}"
)

print(
    "Regiones bacteriocin-priority: "
    f"{len(bact)}"
)

print(
    "CDS con señales RiPP/bacteriocina: "
    f"{len(candidate_genes)}"
)

PY


echo
echo "============================================================"
echo "RESUMEN POR MAG"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/antismash_summary_by_MAG.tsv"

echo
echo "============================================================"
echo "TIPOS DE BGC"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/antismash_product_counts.tsv"

echo
echo "============================================================"
echo "CANDIDATOS RiPP / BACTERIOCINAS"
echo "============================================================"

cat "${OUT}/antismash_ripp_candidates.tsv" \
| column -t -s $'\t'

echo
echo "============================================================"
echo "ANTISMASH SUMMARY COMPLETADO"
echo "Resultados: ${OUT}"
echo "Fin: $(date)"
echo "============================================================"
