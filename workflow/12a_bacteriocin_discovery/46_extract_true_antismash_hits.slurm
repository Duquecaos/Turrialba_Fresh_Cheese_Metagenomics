#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

#SBATCH --job-name=antismash_truehits
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

RESULTS="${ROOT}/26_antismash_rep18/results"
MASTER="${ROOT}/22_final_representative_mags/representative_mags_master.tsv"

OUT="${ROOT}/28_antismash_true_hits"

ANTISMASH="${ANTISMASH:-/opt/ohpc/pub/containers/BIO/antismash-7.1.0.sif}"

mkdir -p "${OUT}"

echo "============================================================"
echo "EXTRACCIÓN DE HITS REALES antiSMASH"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"

if [[ ! -d "${RESULTS}" ]]; then
    echo "ERROR: falta ${RESULTS}" >&2
    exit 1
fi

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

def find_col(fields, options):
    mapping = {
        x.lower().strip(): x
        for x in (fields or [])
    }

    for x in options:
        if x.lower() in mapping:
            return mapping[x.lower()]

    return None


def qvalues(feature, key):
    val = feature.qualifiers.get(key, [])

    if isinstance(val, str):
        return [val]

    return [str(x) for x in val]


def qfirst(feature, key):
    vals = qvalues(feature, key)
    return ";".join(vals)


def strand_symbol(feature):
    s = feature.location.strand

    if s == 1:
        return "+"
    if s == -1:
        return "-"
    return ""


# Perfiles directamente relacionados con bacteriocinas.
# Aquí se clasifica el PERFIL REAL observado, no la regla completa.
DIRECT_BACTERIOCIN = [
    "lactococcin",
    "lactococcin_972",
    "bacteriocin_ii",
    "bacteriociic_cy",
    "bacteriocin_iic",
    "bacteriocin_iid",
    "bacteriocin_iii",
    "bacteriocin_iii",
    "bacteriocin_iii",
    "bacteriocin_iii",
    "bacteriocin_iii",
    "bacteriocin_iii",
    "bacteriocin_ii",
    "lcng-beta",
    "strepbact",
    "linocin_m18",
    "subtilosin",
    "thuricin",
    "nisin",
    "pediocin",
    "enterocin",
    "leucocin",
    "sakacin",
    "plantaricin",
    "carnobacteriocin",
]


def is_direct_bacteriocin(profile):
    p = profile.lower().replace(" ", "_")

    return any(term in p for term in DIRECT_BACTERIOCIN)


# Tipos RiPP que antiSMASH puede reconocer de forma específica
RIPP_PRODUCTS = {
    "ripp-like",
    "ras-ripp",
    "lanthipeptide-class-i",
    "lanthipeptide-class-ii",
    "lanthipeptide-class-iii",
    "lanthipeptide-class-iv",
    "lassopeptide",
    "thiopeptide",
    "ranthipeptide",
    "microviridin",
    "glycocin",
    "cyanobactin",
    "sactipeptide",
    "proteusin",
    "thioamitide",
}


# ============================================================
# Taxonomía
# ============================================================

taxonomy = {}

if MASTER.exists():

    with MASTER.open(newline="") as fh:
        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        magcol = find_col(
            reader.fieldnames,
            [
                "representative_MAG",
                "representative_mag",
                "MAG",
            ]
        )

        gencol = find_col(
            reader.fieldnames,
            ["genus", "gtdb_genus"]
        )

        spcol = find_col(
            reader.fieldnames,
            ["species", "gtdb_species"]
        )

        for row in reader:

            if not magcol:
                continue

            mag = row.get(magcol, "").strip()

            if not mag:
                continue

            taxonomy[mag] = {
                "genus":
                    row.get(gencol, "")
                    if gencol else "",
                "species":
                    row.get(spcol, "")
                    if spcol else "",
            }


# ============================================================
# Extraer hits
# ============================================================

rows = []
direct_rows = []

nt_fasta = []
aa_fasta = []

region_count = 0

for magdir in sorted(RESULTS.iterdir()):

    if not magdir.is_dir():
        continue

    mag = magdir.name
    tx = taxonomy.get(mag, {})

    for gbk in sorted(
        magdir.glob("*.region*.gbk")
    ):

        region_count += 1

        rec = SeqIO.read(
            str(gbk),
            "genbank"
        )

        rfeatures = [
            f for f in rec.features
            if f.type == "region"
        ]

        if not rfeatures:
            continue

        rf = rfeatures[0]

        products = qfirst(
            rf,
            "product"
        )

        edge_raw = qfirst(
            rf,
            "contig_edge"
        ).lower()

        edge = (
            "YES"
            if edge_raw in {
                "true", "yes", "1"
            }
            else "NO"
        )

        m = re.search(
            r"\.region(\d+)\.gbk$",
            gbk.name
        )

        region_no = (
            m.group(1)
            if m else ""
        )

        for f in rec.features:

            if f.type != "CDS":
                continue

            gene_functions = qvalues(
                f,
                "gene_functions"
            )

            if not gene_functions:
                continue

            locus = qfirst(
                f,
                "locus_tag"
            )

            gene = qfirst(
                f,
                "gene"
            )

            product = qfirst(
                f,
                "product"
            )

            gene_kind = qfirst(
                f,
                "gene_kind"
            )

            translation = qfirst(
                f,
                "translation"
            )

            start = (
                int(f.location.start) + 1
            )

            end = int(
                f.location.end
            )

            for annotation in gene_functions:

                if (
                    "rule-based-clusters"
                    not in annotation
                ):
                    continue

                # Ejemplo:
                # biosynthetic (rule-based-clusters)
                # RiPP-like: Lactococcin

                match = re.search(
                    r"^(.+?) "
                    r"\(rule-based-clusters\) "
                    r"([^:]+):\s*(.+)$",
                    annotation
                )

                if not match:
                    continue

                role = match.group(1).strip()
                cluster_type = (
                    match.group(2).strip()
                )
                profile = (
                    match.group(3).strip()
                )

                direct = (
                    "YES"
                    if is_direct_bacteriocin(
                        profile
                    )
                    else "NO"
                )

                ripp = (
                    "YES"
                    if cluster_type.lower()
                    in RIPP_PRODUCTS
                    or "ripp"
                    in cluster_type.lower()
                    or "lanthi"
                    in cluster_type.lower()
                    or "lasso"
                    in cluster_type.lower()
                    or "thiopeptide"
                    in cluster_type.lower()
                    else "NO"
                )

                row = {
                    "MAG": mag,
                    "genus":
                        tx.get("genus", ""),
                    "species":
                        tx.get("species", ""),
                    "contig": rec.id,
                    "region":
                        region_no,
                    "region_products":
                        products,
                    "contig_edge":
                        edge,
                    "CDS_start":
                        start,
                    "CDS_end":
                        end,
                    "CDS_length_nt":
                        end - start + 1,
                    "strand":
                        strand_symbol(f),
                    "locus_tag":
                        locus,
                    "gene":
                        gene,
                    "product":
                        product,
                    "gene_kind":
                        gene_kind,
                    "rule_role":
                        role,
                    "cluster_type":
                        cluster_type,
                    "detected_profile":
                        profile,
                    "RiPP_hit":
                        ripp,
                    "direct_bacteriocin_hit":
                        direct,
                    "region_file":
                        gbk.name,
                }

                rows.append(row)

                if direct == "YES":

                    direct_rows.append(
                        row.copy()
                    )

                    seq_id = (
                        f"{mag}|{rec.id}|"
                        f"region{region_no}|"
                        f"{locus}|{profile}|"
                        f"edge={edge}"
                    )

                    nt_seq = f.extract(
                        rec.seq
                    )

                    nt_fasta.append(
                        f">{seq_id}\n"
                        f"{str(nt_seq)}\n"
                    )

                    if translation:
                        aa_fasta.append(
                            f">{seq_id}\n"
                            f"{translation}\n"
                        )


# ============================================================
# Validar regiones
# ============================================================

if region_count != 66:
    raise SystemExit(
        f"ERROR: se esperaban 66 regiones "
        f"y se procesaron {region_count}"
    )


# ============================================================
# Escritura
# ============================================================

fields = [
    "MAG",
    "genus",
    "species",
    "contig",
    "region",
    "region_products",
    "contig_edge",
    "CDS_start",
    "CDS_end",
    "CDS_length_nt",
    "strand",
    "locus_tag",
    "gene",
    "product",
    "gene_kind",
    "rule_role",
    "cluster_type",
    "detected_profile",
    "RiPP_hit",
    "direct_bacteriocin_hit",
    "region_file",
]


def write_tsv(path, data):

    with path.open(
        "w",
        newline=""
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(data)


write_tsv(
    OUT / "antismash_true_rule_hits.tsv",
    rows
)

write_tsv(
    OUT / "antismash_true_bacteriocin_hits.tsv",
    direct_rows
)

ripp_rows = [
    x for x in rows
    if x["RiPP_hit"] == "YES"
]

write_tsv(
    OUT / "antismash_true_ripp_hits.tsv",
    ripp_rows
)


with (
    OUT /
    "antismash_bacteriocin_candidates_nt.fasta"
).open("w") as fh:

    fh.writelines(nt_fasta)


with (
    OUT /
    "antismash_bacteriocin_candidates_aa.fasta"
).open("w") as fh:

    fh.writelines(aa_fasta)


# ============================================================
# Resumen por MAG
# ============================================================

all_mags = sorted(
    x.name
    for x in RESULTS.iterdir()
    if x.is_dir()
)

by_mag = defaultdict(list)

for x in rows:
    by_mag[x["MAG"]].append(x)

summary_fields = [
    "MAG",
    "genus",
    "species",
    "true_rule_hits",
    "true_RiPP_hits",
    "direct_bacteriocin_hits",
    "direct_bacteriocin_non_edge",
    "direct_profiles",
]

summary = []

for mag in all_mags:

    rr = by_mag.get(
        mag,
        []
    )

    direct = [
        x for x in rr
        if x["direct_bacteriocin_hit"]
        == "YES"
    ]

    profiles = sorted({
        x["detected_profile"]
        for x in direct
    })

    tx = taxonomy.get(
        mag,
        {}
    )

    summary.append({
        "MAG": mag,
        "genus":
            tx.get("genus", ""),
        "species":
            tx.get("species", ""),
        "true_rule_hits":
            len(rr),
        "true_RiPP_hits":
            sum(
                x["RiPP_hit"] == "YES"
                for x in rr
            ),
        "direct_bacteriocin_hits":
            len(direct),
        "direct_bacteriocin_non_edge":
            sum(
                x["contig_edge"] == "NO"
                for x in direct
            ),
        "direct_profiles":
            ";".join(profiles),
    })


with (
    OUT /
    "antismash_true_hits_summary_by_MAG.tsv"
).open("w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=summary_fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(summary)


# ============================================================
# Perfil global
# ============================================================

profile_counts = Counter(
    (
        x["cluster_type"],
        x["detected_profile"],
        x["direct_bacteriocin_hit"],
    )
    for x in rows
)

with (
    OUT /
    "antismash_detected_profile_counts.tsv"
).open("w", newline="") as fh:

    writer = csv.writer(
        fh,
        delimiter="\t"
    )

    writer.writerow([
        "cluster_type",
        "detected_profile",
        "direct_bacteriocin",
        "n_hits",
    ])

    for (
        cluster_type,
        profile,
        direct
    ), n in sorted(
        profile_counts.items(),
        key=lambda x: (
            -x[1],
            x[0][0],
            x[0][1],
        )
    ):

        writer.writerow([
            cluster_type,
            profile,
            direct,
            n,
        ])


print(
    f"Regiones procesadas: {region_count}"
)

print(
    f"Anotaciones rule-based reales: "
    f"{len(rows)}"
)

print(
    f"Hits RiPP reales: "
    f"{len(ripp_rows)}"
)

print(
    f"Hits directos de perfiles "
    f"bacteriocina: {len(direct_rows)}"
)

print(
    "Hits bacteriocina fuera de "
    "contig_edge: "
    f"{sum(x['contig_edge'] == 'NO' for x in direct_rows)}"
)

PY


echo
echo "============================================================"
echo "RESUMEN POR MAG"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/antismash_true_hits_summary_by_MAG.tsv"

echo
echo "============================================================"
echo "PERFILES DE BACTERIOCINAS REALMENTE DETECTADOS"
echo "============================================================"

if [[ $(wc -l < \
"${OUT}/antismash_true_bacteriocin_hits.tsv") -gt 1 ]]; then

    column -t -s $'\t' \
    "${OUT}/antismash_true_bacteriocin_hits.tsv"

else

    echo "No se detectaron perfiles bacteriocina directos."

fi

echo
echo "============================================================"
echo "PERFILES antiSMASH MÁS FRECUENTES"
echo "============================================================"

head -n 30 \
"${OUT}/antismash_detected_profile_counts.tsv" \
| column -t -s $'\t'

echo
echo "============================================================"
echo "ANTISMASH TRUE HITS COMPLETADO"
echo "Resultados: ${OUT}"
echo "Fin: $(date)"
echo "============================================================"
