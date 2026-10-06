#!/bin/bash

#SBATCH --job-name=91d_resist18
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=91d_resist18_%j.out
#SBATCH --error=91d_resist18_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

OUT="${ROOT}/78_resistome_final18_integrated"

PROT="${ROOT}/23_eggnog_rep18/results/representative18.emapper.genepred.fasta"
NUC="${ROOT}/22_final_representative_mags/representative18_traceable.fasta"

ABR_AMR="${ROOT}/25_abricate_consolidated/amr_loci_90_90.tsv"
ABR_SUM="${ROOT}/25_abricate_consolidated/abricate_nonredundant_summary_by_MAG.tsv"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"
SIF="${SIF:-/opt/ohpc/pub/containers/BIO/ncbi-amrfinderplus-4.0.3.sif}"

DB_HOST="${ROOT}/databases/amrfinderplus/latest"
DB_CONT="/work/databases/amrfinderplus/latest"

mkdir -p "$OUT"

GFF="${OUT}/91D_final18_prodigal_coordinates.gff3"
AMR_COMBINED="${OUT}/91D_amrfinder_combined_plus.tsv"

echo "======================================================================"
echo "PASO 91D - RESISTOMA FINAL18 INTEGRADO"
echo "Inicio: $(date)"
echo "Host: $(hostname)"
echo "Job: ${SLURM_JOB_ID}"
echo "======================================================================"

for F in \
    "$PROT" \
    "$NUC" \
    "$ABR_AMR" \
    "$ABR_SUM" \
    "$SIF"
do
    [[ -s "$F" ]] || {
        echo "ERROR: falta archivo: $F" >&2
        exit 2
    }
done

[[ -d "$DB_HOST" ]] || {
    echo "ERROR: falta database AMRFinderPlus: $DB_HOST" >&2
    exit 3
}

echo
echo "======================================================================"
echo "1. GENERANDO GFF A PARTIR DE CABECERAS PRODIGAL"
echo "======================================================================"

python3 - <<'PY'
import os
import re
from pathlib import Path

root = Path(
    f"/scratch/global/{os.environ['USER']}/Shotgun_MAGs_Turrialba"
)

prot = (
    root
    / "23_eggnog_rep18/results"
    / "representative18.emapper.genepred.fasta"
)

nuc = (
    root
    / "22_final_representative_mags"
    / "representative18_traceable.fasta"
)

out = (
    root
    / "78_resistome_final18_integrated"
    / "91D_final18_prodigal_coordinates.gff3"
)

# ----------------------------------------------------------
# Contigs presentes en FASTA nucleotídico
# ----------------------------------------------------------

contigs = set()

with nuc.open() as fh:
    for line in fh:
        if line.startswith(">"):
            contigs.add(
                line[1:].split()[0]
            )

if not contigs:
    raise RuntimeError(
        "No se encontraron contigs en FASTA nucleotídico."
    )

# ----------------------------------------------------------
# Prodigal FASTA:
#
# >protein_id # start # end # strand # metadata
# ----------------------------------------------------------

pattern = re.compile(
    r"^>(\S+)\s+#\s+"
    r"(\d+)\s+#\s+"
    r"(\d+)\s+#\s+"
    r"(-?1)\s+#\s*(.*)$"
)

n_proteins = 0
missing_contigs = []
bad_headers = []
seen = set()

with prot.open() as inp, out.open("w") as gff:

    gff.write("##gff-version 3\n")

    for line in inp:

        if not line.startswith(">"):
            continue

        n_proteins += 1

        m = pattern.match(line.rstrip())

        if not m:
            bad_headers.append(
                line.rstrip()
            )
            continue

        protein_id = m.group(1)
        start = int(m.group(2))
        end = int(m.group(3))
        strand_num = int(m.group(4))

        if protein_id in seen:
            raise RuntimeError(
                f"Protein ID duplicado: {protein_id}"
            )

        seen.add(protein_id)

        if "_" not in protein_id:
            raise RuntimeError(
                f"No se puede derivar contig de {protein_id}"
            )

        # El último _N corresponde al número de ORF
        contig = protein_id.rsplit("_", 1)[0]

        if contig not in contigs:
            missing_contigs.append(
                (protein_id, contig)
            )
            continue

        strand = "+" if strand_num == 1 else "-"

        lo = min(start, end)
        hi = max(start, end)

        attrs = (
            f"ID={protein_id};"
            f"Name={protein_id}"
        )

        gff.write(
            "\t".join([
                contig,
                "Prodigal",
                "CDS",
                str(lo),
                str(hi),
                ".",
                strand,
                "0",
                attrs,
            ])
            + "\n"
        )

print("n_contigs_nucleotide:", len(contigs))
print("n_proteins:", n_proteins)
print("bad_headers:", len(bad_headers))
print("missing_contigs:", len(missing_contigs))

if bad_headers:
    print("\nPrimeras cabeceras problemáticas:")
    for x in bad_headers[:10]:
        print(x)

if missing_contigs:
    print("\nPrimeros contigs ausentes:")
    for x in missing_contigs[:10]:
        print(*x, sep="\t")

if bad_headers:
    raise RuntimeError(
        "Hay cabeceras Prodigal no interpretables."
    )

if missing_contigs:
    raise RuntimeError(
        "Hay proteínas cuyo contig no está en el FASTA nucleotídico."
    )

if n_proteins != len(seen):
    raise RuntimeError(
        "Número de proteínas no coincide con IDs únicos."
    )

print("GFF GENERADO CORRECTAMENTE:", out)
PY

echo
echo "--- GFF ---"
head -n 5 "$GFF"

echo
echo -n "CDS en GFF: "
awk '$3=="CDS"{n++} END{print n+0}' "$GFF"

echo
echo "======================================================================"
echo "2. AMRFINDERPLUS COMBINADO"
echo "======================================================================"

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    "$APPTAINER" exec \
    --bind "${ROOT}:/work" \
    "$SIF" \
    amrfinder \
        --protein "/work/23_eggnog_rep18/results/representative18.emapper.genepred.fasta" \
        --nucleotide "/work/22_final_representative_mags/representative18_traceable.fasta" \
        --gff "/work/78_resistome_final18_integrated/91D_final18_prodigal_coordinates.gff3" \
        --annotation_format standard \
        --plus \
        --database "$DB_CONT" \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --output "/work/78_resistome_final18_integrated/91D_amrfinder_combined_plus.tsv"

[[ -s "$AMR_COMBINED" ]] || {
    echo "ERROR: no se generó salida combinada AMRFinderPlus." >&2
    exit 4
}

echo
echo "======================================================================"
echo "3. INTEGRANDO AMRFINDERPLUS + ABRICATE 90/90"
echo "======================================================================"

python3 - <<'PY'
import csv
import os
import re
from pathlib import Path
from collections import defaultdict, Counter

root = Path(
    f"/scratch/global/{os.environ['USER']}/Shotgun_MAGs_Turrialba"
)

base = root / "78_resistome_final18_integrated"

amr_file = base / "91D_amrfinder_combined_plus.tsv"

abr_file = (
    root
    / "25_abricate_consolidated"
    / "amr_loci_90_90.tsv"
)

tax_file = (
    root
    / "25_abricate_consolidated"
    / "abricate_nonredundant_summary_by_MAG.tsv"
)

out_loci = base / "91D_integrated_AMR_loci.tsv"
out_mag = base / "91D_AMR_summary_by_MAG.tsv"
out_class = base / "91D_AMR_summary_by_class.tsv"
out_stress = base / "91D_amrfinder_STRESS.tsv"
out_vir = base / "91D_amrfinder_VIRULENCE.tsv"
out_nonamr = base / "91D_amrfinder_non_AMR.tsv"
out_summary = base / "91D_resistome_summary.tsv"

# ----------------------------------------------------------
# Helpers
# ----------------------------------------------------------

def clean(x):
    if x is None:
        return ""
    x = str(x).strip()
    if x.upper() == "NA":
        return ""
    return x

def split_semicolon(x):
    return [
        y.strip()
        for y in clean(x).split(";")
        if y.strip()
    ]

def uniq(values):
    result = []
    seen = set()

    for x in values:
        x = clean(x)
        if not x:
            continue
        if x not in seen:
            seen.add(x)
            result.append(x)

    return result

def join(values):
    return ";".join(
        sorted(uniq(values))
    )

def mag_from_contig(contig):
    contig = clean(contig)

    if "___" not in contig:
        raise RuntimeError(
            f"No se puede derivar MAG de contig: {contig}"
        )

    return contig.split("___", 1)[0]

def overlap_fraction(a1, a2, b1, b2):
    overlap = max(
        0,
        min(a2, b2) -
        max(a1, b1) + 1
    )

    if overlap == 0:
        return 0.0

    la = a2 - a1 + 1
    lb = b2 - b1 + 1

    return overlap / min(la, lb)

def norm_symbol(x):
    return re.sub(
        r"[^a-z0-9]+",
        "",
        clean(x).lower()
    )

# ----------------------------------------------------------
# Taxonomía
# ----------------------------------------------------------

taxonomy = {}

with tax_file.open(newline="") as fh:
    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in r:
        taxonomy[row["MAG"]] = {
            "genus": clean(row.get("genus")),
            "species": clean(row.get("species")),
        }

# ----------------------------------------------------------
# AMRFinder combinado
# ----------------------------------------------------------

with amr_file.open(newline="") as fh:
    r = csv.DictReader(
        fh,
        delimiter="\t"
    )
    amrf_all = list(r)
    amrf_fields = list(r.fieldnames)

required = {
    "Contig id",
    "Start",
    "Stop",
    "Element symbol",
    "Element name",
    "Scope",
    "Type",
    "Subtype",
    "Class",
    "Subclass",
    "Method",
    "% Coverage of reference",
    "% Identity to reference",
}

missing = required - set(amrf_fields)

if missing:
    raise RuntimeError(
        "Faltan campos AMRFinder: "
        + ", ".join(sorted(missing))
    )

# Guardar no-AMR por tipo
non_amr = [
    row
    for row in amrf_all
    if clean(row.get("Type")).upper() != "AMR"
]

for path, rows, target_type in [
    (
        out_stress,
        [
            x for x in amrf_all
            if clean(x.get("Type")).upper() == "STRESS"
        ],
        "STRESS",
    ),
    (
        out_vir,
        [
            x for x in amrf_all
            if clean(x.get("Type")).upper() == "VIRULENCE"
        ],
        "VIRULENCE",
    ),
    (
        out_nonamr,
        non_amr,
        "NON_AMR",
    ),
]:

    with path.open("w", newline="") as out:
        w = csv.DictWriter(
            out,
            fieldnames=amrf_fields,
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(rows)

# Convertir hits AMR a formato interno
records_amrf = []

for i, row in enumerate(amrf_all, start=1):

    if clean(row.get("Type")).upper() != "AMR":
        continue

    contig = clean(row["Contig id"])

    start = int(row["Start"])
    end = int(row["Stop"])

    lo = min(start, end)
    hi = max(start, end)

    records_amrf.append({
        "id": f"AMRF_{i:05d}",
        "source": "AMRFinderPlus",
        "MAG": mag_from_contig(contig),
        "contig": contig,
        "start": lo,
        "end": hi,
        "strand": clean(row.get("Strand")),
        "symbol": clean(row.get("Element symbol")),
        "name": clean(row.get("Element name")),
        "scope": clean(row.get("Scope")),
        "class": clean(row.get("Class")),
        "subclass": clean(row.get("Subclass")),
        "method": clean(row.get("Method")),
        "coverage": clean(
            row.get("% Coverage of reference")
        ),
        "identity": clean(
            row.get("% Identity to reference")
        ),
        "accession": clean(
            row.get("Closest reference accession")
        ),
        "protein_id": clean(
            row.get("Protein id")
        ),
        "hmm": clean(
            row.get("HMM accession")
        ),
        "raw": row,
    })

# ----------------------------------------------------------
# ABRicate 90/90
# ----------------------------------------------------------

records_abr = []

with abr_file.open(newline="") as fh:

    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in r:

        records_abr.append({
            "id": clean(row["locus_id"]),
            "source": "ABRicate_90_90",
            "MAG": clean(row["MAG"]),
            "contig": clean(row["contig"]),
            "start": int(row["start"]),
            "end": int(row["end"]),
            "strand": clean(row["strand"]),
            "symbols": split_semicolon(
                row.get("gene_labels")
            ),
            "databases": split_semicolon(
                row.get("databases")
            ),
            "n_databases": clean(
                row.get("n_databases")
            ),
            "n_raw_hits": clean(
                row.get("n_raw_hits")
            ),
            "coverage": clean(
                row.get("max_coverage")
            ),
            "identity": clean(
                row.get("max_identity")
            ),
            "accessions": split_semicolon(
                row.get("accessions")
            ),
            "resistance": split_semicolon(
                row.get("resistance")
            ),
            "products": clean(
                row.get("products")
            ),
            "raw": row,
        })

# ----------------------------------------------------------
# Unión solo entre AMRFinder y ABRicate por locus físico
# ----------------------------------------------------------

all_records = records_amrf + records_abr

parent = {
    x["id"]: x["id"]
    for x in all_records
}

def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

def union(a, b):
    ra = find(a)
    rb = find(b)

    if ra != rb:
        parent[rb] = ra

abr_by_contig = defaultdict(list)

for x in records_abr:
    abr_by_contig[
        (x["MAG"], x["contig"])
    ].append(x)

for a in records_amrf:

    for b in abr_by_contig.get(
        (a["MAG"], a["contig"]),
        []
    ):

        ov = overlap_fraction(
            a["start"],
            a["end"],
            b["start"],
            b["end"],
        )

        if ov >= 0.80:
            union(
                a["id"],
                b["id"]
            )

clusters = defaultdict(list)

for x in all_records:
    clusters[find(x["id"])].append(x)

# ----------------------------------------------------------
# Construir loci integrados
# ----------------------------------------------------------

rows = []

def symbols_from_members(members):

    xs = []

    for x in members:

        if x["source"] == "AMRFinderPlus":
            xs.append(x["symbol"])
        else:
            xs.extend(x["symbols"])

    return uniq(xs)

def symbol_concordance(amrf, abr):

    a = {
        norm_symbol(x["symbol"])
        for x in amrf
        if norm_symbol(x["symbol"])
    }

    b = {
        norm_symbol(y)
        for x in abr
        for y in x["symbols"]
        if norm_symbol(y)
    }

    if not a or not b:
        return "not_assessable"

    if a & b:
        return "concordant_gene_label"

    # permitir relación de prefijo/subcadena
    # útil para familias/variantes de alelos
    for x in a:
        for y in b:
            if len(x) >= 4 and len(y) >= 4:
                if x in y or y in x:
                    return (
                        "compatible_family_or_allele_label"
                    )

    return (
        "positional_overlap_annotation_difference"
    )

ordered_clusters = sorted(
    clusters.values(),
    key=lambda members: (
        members[0]["MAG"],
        members[0]["contig"],
        min(x["start"] for x in members),
    )
)

for number, members in enumerate(
    ordered_clusters,
    start=1
):

    amrf = [
        x for x in members
        if x["source"] == "AMRFinderPlus"
    ]

    abr = [
        x for x in members
        if x["source"] == "ABRicate_90_90"
    ]

    mag = members[0]["MAG"]
    contig = members[0]["contig"]

    if amrf and abr:
        support = (
            "AMRFinderPlus_and_ABRicate90_90"
        )
        annotation = symbol_concordance(
            amrf,
            abr
        )

    elif amrf:
        support = "AMRFinderPlus_only"
        annotation = (
            "AMRFinderPlus_only"
        )

    else:
        support = "ABRicate90_90_only"
        annotation = (
            "ABRicate90_90_only"
        )

    tx = taxonomy.get(
        mag,
        {
            "genus": "",
            "species": "",
        }
    )

    scopes = uniq(
        x["scope"]
        for x in amrf
    )

    if "core" in scopes:
        amrf_scope_summary = "core"
    elif "plus" in scopes:
        amrf_scope_summary = "plus"
    else:
        amrf_scope_summary = ""

    rows.append({
        "integrated_locus_id":
            f"AMRLOC_{number:04d}",

        "MAG": mag,
        "genus": tx["genus"],
        "species": tx["species"],

        "contig": contig,

        "start":
            min(x["start"] for x in members),

        "end":
            max(x["end"] for x in members),

        "length_bp":
            max(x["end"] for x in members)
            - min(x["start"] for x in members)
            + 1,

        "evidence_support":
            support,

        "annotation_concordance":
            annotation,

        "gene_symbols":
            join(symbols_from_members(members)),

        "amrfinder_symbols":
            join(
                x["symbol"]
                for x in amrf
            ),

        "abricate_gene_labels":
            join(
                y
                for x in abr
                for y in x["symbols"]
            ),

        "amrfinder_scope":
            amrf_scope_summary,

        "amrfinder_class":
            join(
                x["class"]
                for x in amrf
            ),

        "amrfinder_subclass":
            join(
                x["subclass"]
                for x in amrf
            ),

        "amrfinder_methods":
            join(
                x["method"]
                for x in amrf
            ),

        "amrfinder_max_identity":
            max(
                [
                    float(x["identity"])
                    for x in amrf
                    if clean(x["identity"])
                ],
                default="",
            ),

        "amrfinder_max_coverage":
            max(
                [
                    float(x["coverage"])
                    for x in amrf
                    if clean(x["coverage"])
                ],
                default="",
            ),

        "amrfinder_reference_accessions":
            join(
                x["accession"]
                for x in amrf
            ),

        "amrfinder_HMMs":
            join(
                x["hmm"]
                for x in amrf
            ),

        "abricate_databases":
            join(
                y
                for x in abr
                for y in x["databases"]
            ),

        "abricate_max_identity":
            max(
                [
                    float(x["identity"])
                    for x in abr
                    if clean(x["identity"])
                ],
                default="",
            ),

        "abricate_max_coverage":
            max(
                [
                    float(x["coverage"])
                    for x in abr
                    if clean(x["coverage"])
                ],
                default="",
            ),

        "abricate_resistance_labels":
            join(
                y
                for x in abr
                for y in x["resistance"]
            ),

        "n_amrfinder_records":
            len(amrf),

        "n_abricate_loci":
            len(abr),
    })

# ----------------------------------------------------------
# Tabla integrada
# ----------------------------------------------------------

fields = [
    "integrated_locus_id",
    "MAG",
    "genus",
    "species",
    "contig",
    "start",
    "end",
    "length_bp",
    "evidence_support",
    "annotation_concordance",
    "gene_symbols",
    "amrfinder_symbols",
    "abricate_gene_labels",
    "amrfinder_scope",
    "amrfinder_class",
    "amrfinder_subclass",
    "amrfinder_methods",
    "amrfinder_max_identity",
    "amrfinder_max_coverage",
    "amrfinder_reference_accessions",
    "amrfinder_HMMs",
    "abricate_databases",
    "abricate_max_identity",
    "abricate_max_coverage",
    "abricate_resistance_labels",
    "n_amrfinder_records",
    "n_abricate_loci",
]

with out_loci.open("w", newline="") as out:

    w = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t",
    )

    w.writeheader()
    w.writerows(rows)

# ----------------------------------------------------------
# Resumen por MAG
# ----------------------------------------------------------

mag_counts = defaultdict(
    lambda: Counter()
)

for r in rows:

    x = mag_counts[r["MAG"]]

    x["total"] += 1
    x[r["evidence_support"]] += 1

    for c in split_semicolon(
        r["amrfinder_class"]
    ):
        x[f"class::{c}"] += 1

mag_fields = [
    "MAG",
    "genus",
    "species",
    "AMR_loci_total",
    "supported_by_both",
    "AMRFinderPlus_only",
    "ABRicate90_90_only",
    "AMR_classes",
]

with out_mag.open("w", newline="") as out:

    w = csv.DictWriter(
        out,
        fieldnames=mag_fields,
        delimiter="\t",
    )

    w.writeheader()

    for mag in sorted(taxonomy):

        x = mag_counts[mag]
        tx = taxonomy[mag]

        classes = sorted(
            k.split("::", 1)[1]
            for k in x
            if k.startswith("class::")
        )

        w.writerow({
            "MAG": mag,
            "genus": tx["genus"],
            "species": tx["species"],
            "AMR_loci_total":
                x["total"],
            "supported_by_both":
                x[
                    "AMRFinderPlus_and_ABRicate90_90"
                ],
            "AMRFinderPlus_only":
                x["AMRFinderPlus_only"],
            "ABRicate90_90_only":
                x["ABRicate90_90_only"],
            "AMR_classes":
                ";".join(classes),
        })

# ----------------------------------------------------------
# Resumen por clase AMRFinder
# ----------------------------------------------------------

class_counts = Counter()

for r in rows:
    for c in split_semicolon(
        r["amrfinder_class"]
    ):
        class_counts[c] += 1

with out_class.open("w", newline="") as out:

    w = csv.writer(
        out,
        delimiter="\t",
    )

    w.writerow([
        "AMR_class",
        "integrated_loci",
    ])

    for c, n in sorted(
        class_counts.items(),
        key=lambda x: (-x[1], x[0])
    ):
        w.writerow([c, n])

# ----------------------------------------------------------
# Resumen general
# ----------------------------------------------------------

support_counts = Counter(
    r["evidence_support"]
    for r in rows
)

annotation_counts = Counter(
    r["annotation_concordance"]
    for r in rows
)

type_counts = Counter(
    clean(row.get("Type"))
    for row in amrf_all
)

scope_counts = Counter(
    clean(row.get("Scope"))
    for row in amrf_all
)

metrics = [
    (
        "amrfinder_combined_records_all_types",
        len(amrf_all),
    ),
    (
        "amrfinder_combined_AMR_records",
        len(records_amrf),
    ),
    (
        "amrfinder_combined_non_AMR_records",
        len(non_amr),
    ),
    (
        "abricate_AMR_loci_90_90",
        len(records_abr),
    ),
    (
        "integrated_unique_AMR_loci",
        len(rows),
    ),
    (
        "integrated_supported_by_both",
        support_counts[
            "AMRFinderPlus_and_ABRicate90_90"
        ],
    ),
    (
        "integrated_AMRFinderPlus_only",
        support_counts[
            "AMRFinderPlus_only"
        ],
    ),
    (
        "integrated_ABRicate90_90_only",
        support_counts[
            "ABRicate90_90_only"
        ],
    ),
]

for k, v in sorted(type_counts.items()):
    metrics.append(
        (f"amrfinder_type_{k}", v)
    )

for k, v in sorted(scope_counts.items()):
    metrics.append(
        (f"amrfinder_scope_{k}", v)
    )

for k, v in sorted(annotation_counts.items()):
    metrics.append(
        (f"annotation_{k}", v)
    )

with out_summary.open("w") as out:
    out.write("metric\tvalue\n")

    for k, v in metrics:
        out.write(
            f"{k}\t{v}\n"
        )

print("=" * 78)
print("INTEGRACION COMPLETADA")
print("=" * 78)

for k, v in metrics:
    print(k, v, sep="\t")
PY

echo
echo "======================================================================"
echo "4. PROVENANCE"
echo "======================================================================"

{
    echo -e "field\tvalue"
    echo -e "analysis\tfinal18_integrated_resistome"
    echo -e "AMRFinderPlus_version\t4.0.3"
    echo -e "AMRFinderPlus_database\t2025-07-16.1"
    echo -e "AMRFinder_mode\tcombined_protein_nucleotide_GFF"
    echo -e "AMRFinder_plus\ttrue"
    echo -e "AMRFinder_organism\tNOT_USED"
    echo -e "ABRicate_threshold\t90_identity_90_coverage"
    echo -e "ABRicate_databases\tcard;megares;ncbi;resfinder"
    echo -e "integration_overlap_threshold\t0.80_of_shorter_interval"
    echo -e "date\t$(date -Iseconds)"
} > "${OUT}/91D_run_metadata.tsv"

sha256sum \
    "$PROT" \
    "$NUC" \
    "$ABR_AMR" \
    "$AMR_COMBINED" \
    > "${OUT}/91D_input_output_sha256.tsv"

echo
echo "======================================================================"
echo "5. VALIDACION DE SALIDAS"
echo "======================================================================"

for F in \
    "${OUT}/91D_final18_prodigal_coordinates.gff3" \
    "${OUT}/91D_amrfinder_combined_plus.tsv" \
    "${OUT}/91D_integrated_AMR_loci.tsv" \
    "${OUT}/91D_AMR_summary_by_MAG.tsv" \
    "${OUT}/91D_AMR_summary_by_class.tsv" \
    "${OUT}/91D_amrfinder_STRESS.tsv" \
    "${OUT}/91D_amrfinder_VIRULENCE.tsv" \
    "${OUT}/91D_amrfinder_non_AMR.tsv" \
    "${OUT}/91D_resistome_summary.tsv" \
    "${OUT}/91D_run_metadata.tsv"
do
    [[ -s "$F" ]] || {
        echo "ERROR: salida vacía/ausente: $F" >&2
        exit 5
    }

    echo "OK: $(basename "$F")"
done

echo
echo "======================================================================"
echo "RESUMEN"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/91D_resistome_summary.tsv"

echo
echo "======================================================================"
echo "PASO 91D FINALIZÓ CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "======================================================================"
