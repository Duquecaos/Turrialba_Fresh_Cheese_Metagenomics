#!/bin/bash

#SBATCH --job-name=91e2_amrcons
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=91e2_amrcons_%j.out
#SBATCH --error=91e2_amrcons_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

INROOT="${ROOT}/79_resistome_coassemblies"
OUT="${ROOT}/80_resistome_community_consolidated"

FINAL18_AMRF="${ROOT}/78_resistome_final18_integrated/91D_amrfinder_combined_plus.tsv"
FINAL18_INT="${ROOT}/78_resistome_final18_integrated/91D_integrated_AMR_loci.tsv"

FINAL18_FAA="${ROOT}/23_eggnog_rep18/results/representative18.emapper.genepred.fasta"
FINAL18_FNA="${ROOT}/22_final_representative_mags/representative18_traceable.fasta"

mkdir -p "$OUT"

echo "======================================================================"
echo "PASO 91E2 - CONSOLIDACION RESISTOMA COMUNITARIO"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "======================================================================"

for G in L1 L2 L3 M1 M2 M3
do

    for F in \
        "${INROOT}/${G}/${G}.prodigal.faa" \
        "${INROOT}/${G}/${G}.amrfinder_plus.tsv" \
        "${INROOT}/${G}/${G}.AMR.tsv"
    do
        [[ -s "$F" ]] || {
            echo "ERROR: falta salida 91E1: $F" >&2
            exit 2
        }
    done

done

for F in \
    "$FINAL18_AMRF" \
    "$FINAL18_INT" \
    "$FINAL18_FAA" \
    "$FINAL18_FNA"
do
    [[ -s "$F" ]] || {
        echo "ERROR: falta referencia final18: $F" >&2
        exit 3
    }
done

srun --exact \
    -n 1 \
    -c 1 \
    python3 - <<'PY'
import csv
import os
import re
import hashlib
from pathlib import Path
from collections import Counter, defaultdict

root = Path(
    f"/scratch/global/{os.environ['USER']}/Shotgun_MAGs_Turrialba"
)

inroot = root / "79_resistome_coassemblies"
out = root / "80_resistome_community_consolidated"

groups = [
    "L1",
    "L2",
    "L3",
    "M1",
    "M2",
    "M3",
]

# ============================================================
# Helpers
# ============================================================

def clean(x):

    if x is None:
        return ""

    x = str(x).strip()

    if x.upper() == "NA":
        return ""

    return x


def normalize_symbol(x):

    return re.sub(
        r"[^a-z0-9]+",
        "",
        clean(x).lower()
    )


def sha256_seq(seq):

    return hashlib.sha256(
        seq.encode()
    ).hexdigest()


def revcomp(seq):

    table = str.maketrans(
        "ACGTNacgtn",
        "TGCANtgcan"
    )

    return seq.translate(table)[::-1]


def read_fasta(path):

    seqs = {}

    current = None
    buf = []

    with path.open() as fh:

        for line in fh:

            line = line.rstrip()

            if not line:
                continue

            if line.startswith(">"):

                if current is not None:
                    seqs[current] = "".join(buf)

                current = line[1:].split()[0]
                buf = []

            else:
                buf.append(line)

    if current is not None:
        seqs[current] = "".join(buf)

    return seqs


def split_symbols(x):

    return [
        y.strip()
        for y in clean(x).split(";")
        if y.strip()
    ]


# ============================================================
# FINAL18: símbolos y secuencias conocidas
# ============================================================

final18_symbols = set()

final18_integrated = (
    root
    / "78_resistome_final18_integrated"
    / "91D_integrated_AMR_loci.tsv"
)

with final18_integrated.open() as fh:

    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in r:

        for symbol in split_symbols(
            row.get("gene_symbols", "")
        ):

            n = normalize_symbol(symbol)

            if n:
                final18_symbols.add(n)


# ------------------------------------------------------------
# FASTA final18
# ------------------------------------------------------------

final18_proteins = read_fasta(
    root
    / "23_eggnog_rep18/results"
    / "representative18.emapper.genepred.fasta"
)

final18_nucl = read_fasta(
    root
    / "22_final_representative_mags"
    / "representative18_traceable.fasta"
)

final18_amrf = (
    root
    / "78_resistome_final18_integrated"
    / "91D_amrfinder_combined_plus.tsv"
)

final18_protein_hashes = set()
final18_nucleotide_hashes = set()

with final18_amrf.open() as fh:

    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in r:

        if clean(row.get("Type")) != "AMR":
            continue

        pid = clean(
            row.get("Protein id")
        )

        if pid and pid in final18_proteins:

            final18_protein_hashes.add(
                sha256_seq(
                    final18_proteins[pid]
                    .replace("*", "")
                    .upper()
                )
            )

        contig = clean(
            row.get("Contig id")
        )

        if (
            contig
            and contig in final18_nucl
            and clean(row.get("Start"))
            and clean(row.get("Stop"))
        ):

            start = int(row["Start"])
            stop = int(row["Stop"])

            lo = min(start, stop)
            hi = max(start, stop)

            seq = final18_nucl[contig][
                lo - 1:hi
            ].upper()

            if clean(row.get("Strand")) == "-":
                seq = revcomp(seq)

            final18_nucleotide_hashes.add(
                sha256_seq(seq)
            )


# ============================================================
# Coensamblajes
# ============================================================

all_rows = []
all_plus_rows = []

for group in groups:

    producer = group[0]
    biological_unit = group[1:]

    base = inroot / group

    assembly = read_fasta(
        root
        / "05_assembly"
        / group
        / "final.contigs.fa"
    )

    proteins = read_fasta(
        base
        / f"{group}.prodigal.faa"
    )

    amrf_file = (
        base
        / f"{group}.amrfinder_plus.tsv"
    )

    with amrf_file.open() as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        rows = list(reader)

    for row in rows:

        row2 = dict(row)

        row2["coassembly"] = group
        row2["producer"] = producer
        row2["biological_unit"] = biological_unit

        all_plus_rows.append(row2)

        if clean(row.get("Type")) != "AMR":
            continue

        contig = clean(
            row.get("Contig id")
        )

        pid = clean(
            row.get("Protein id")
        )

        start = int(row["Start"])
        stop = int(row["Stop"])

        lo = min(start, stop)
        hi = max(start, stop)

        strand = clean(
            row.get("Strand")
        )

        # ----------------------------------------------------
        # Protein SHA256
        # ----------------------------------------------------

        protein_hash = ""

        if pid and pid in proteins:

            pseq = (
                proteins[pid]
                .replace("*", "")
                .upper()
            )

            protein_hash = sha256_seq(
                pseq
            )

        # ----------------------------------------------------
        # Nucleotide SHA256
        # ----------------------------------------------------

        if contig not in assembly:
            raise RuntimeError(
                f"{group}: contig no encontrado: {contig}"
            )

        nseq = assembly[contig][
            lo - 1:hi
        ].upper()

        if strand == "-":
            nseq = revcomp(nseq)

        nucleotide_hash = sha256_seq(
            nseq
        )

        symbol = clean(
            row.get("Element symbol")
        )

        norm = normalize_symbol(symbol)

        exact_protein = (
            bool(protein_hash)
            and protein_hash
            in final18_protein_hashes
        )

        exact_nucleotide = (
            nucleotide_hash
            in final18_nucleotide_hashes
        )

        same_symbol = (
            bool(norm)
            and norm
            in final18_symbols
        )

        if exact_protein or exact_nucleotide:

            representation = (
                "exact_sequence_match_final18"
            )

        elif same_symbol:

            representation = (
                "same_AMR_symbol_final18_nonexact"
            )

        else:

            representation = (
                "additional_AMR_symbol_relative_to_final18"
            )

        all_rows.append({
            "community_locus_id": "",
            "coassembly": group,
            "producer": producer,
            "biological_unit": biological_unit,
            "contig": contig,
            "start": lo,
            "stop": hi,
            "strand": strand,
            "protein_id": pid,
            "element_symbol": symbol,
            "element_name":
                clean(row.get("Element name")),
            "scope":
                clean(row.get("Scope")),
            "class":
                clean(row.get("Class")),
            "subclass":
                clean(row.get("Subclass")),
            "method":
                clean(row.get("Method")),
            "coverage_reference":
                clean(
                    row.get(
                        "% Coverage of reference"
                    )
                ),
            "identity_reference":
                clean(
                    row.get(
                        "% Identity to reference"
                    )
                ),
            "closest_reference_accession":
                clean(
                    row.get(
                        "Closest reference accession"
                    )
                ),
            "HMM_accession":
                clean(
                    row.get(
                        "HMM accession"
                    )
                ),
            "protein_sha256":
                protein_hash,
            "nucleotide_sha256":
                nucleotide_hash,
            "exact_final18_protein_match":
                int(exact_protein),
            "exact_final18_nucleotide_match":
                int(exact_nucleotide),
            "symbol_seen_in_final18":
                int(same_symbol),
            "representation_vs_final18":
                representation,
        })


# ============================================================
# Orden estable + IDs
# ============================================================

all_rows.sort(
    key=lambda x: (
        x["coassembly"],
        x["contig"],
        int(x["start"]),
        int(x["stop"]),
        x["element_symbol"],
    )
)

for i, row in enumerate(
    all_rows,
    start=1
):
    row["community_locus_id"] = (
        f"COMMAMR_{i:05d}"
    )


# ============================================================
# Tabla maestra
# ============================================================

master_fields = [
    "community_locus_id",
    "coassembly",
    "producer",
    "biological_unit",
    "contig",
    "start",
    "stop",
    "strand",
    "protein_id",
    "element_symbol",
    "element_name",
    "scope",
    "class",
    "subclass",
    "method",
    "coverage_reference",
    "identity_reference",
    "closest_reference_accession",
    "HMM_accession",
    "protein_sha256",
    "nucleotide_sha256",
    "exact_final18_protein_match",
    "exact_final18_nucleotide_match",
    "symbol_seen_in_final18",
    "representation_vs_final18",
]

master = (
    out
    / "91E_community_AMR_loci.tsv"
)

with master.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=master_fields,
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(all_rows)


# ============================================================
# Todos los registros --plus
# ============================================================

plus_fields = (
    ["coassembly",
     "producer",
     "biological_unit"]
    +
    [
        x for x in
        all_plus_rows[0].keys()
        if x not in {
            "coassembly",
            "producer",
            "biological_unit"
        }
    ]
)

with (
    out
    / "91E_all_AMRFinderPlus_records.tsv"
).open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=plus_fields,
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(all_plus_rows)


# ============================================================
# Clusters exactos proteína
# ============================================================

protein_clusters = defaultdict(list)

for row in all_rows:

    if row["protein_sha256"]:

        protein_clusters[
            row["protein_sha256"]
        ].append(row)

protein_cluster_file = (
    out
    / "91E_exact_protein_clusters.tsv"
)

with protein_cluster_file.open(
    "w",
    newline=""
) as fh:

    fields = [
        "protein_cluster_id",
        "protein_sha256",
        "n_records",
        "n_coassemblies",
        "coassemblies",
        "producers",
        "symbols",
        "classes",
        "representation_vs_final18",
    ]

    w = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    w.writeheader()

    ordered = sorted(
        protein_clusters.items(),
        key=lambda x: (
            -len(x[1]),
            x[0]
        )
    )

    for i, (h, members) in enumerate(
        ordered,
        start=1
    ):

        w.writerow({
            "protein_cluster_id":
                f"COMM_PROT_{i:05d}",
            "protein_sha256":
                h,
            "n_records":
                len(members),
            "n_coassemblies":
                len({
                    x["coassembly"]
                    for x in members
                }),
            "coassemblies":
                ";".join(sorted({
                    x["coassembly"]
                    for x in members
                })),
            "producers":
                ";".join(sorted({
                    x["producer"]
                    for x in members
                })),
            "symbols":
                ";".join(sorted({
                    x["element_symbol"]
                    for x in members
                    if x["element_symbol"]
                })),
            "classes":
                ";".join(sorted({
                    x["class"]
                    for x in members
                    if x["class"]
                })),
            "representation_vs_final18":
                ";".join(sorted({
                    x["representation_vs_final18"]
                    for x in members
                })),
        })


# ============================================================
# Clusters exactos nucleótido
# ============================================================

nucleotide_clusters = defaultdict(list)

for row in all_rows:

    nucleotide_clusters[
        row["nucleotide_sha256"]
    ].append(row)

nuc_cluster_file = (
    out
    / "91E_exact_nucleotide_clusters.tsv"
)

with nuc_cluster_file.open(
    "w",
    newline=""
) as fh:

    fields = [
        "nucleotide_cluster_id",
        "nucleotide_sha256",
        "n_records",
        "n_coassemblies",
        "coassemblies",
        "producers",
        "symbols",
        "classes",
        "representation_vs_final18",
    ]

    w = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    w.writeheader()

    ordered = sorted(
        nucleotide_clusters.items(),
        key=lambda x: (
            -len(x[1]),
            x[0]
        )
    )

    for i, (h, members) in enumerate(
        ordered,
        start=1
    ):

        w.writerow({
            "nucleotide_cluster_id":
                f"COMM_NUC_{i:05d}",
            "nucleotide_sha256":
                h,
            "n_records":
                len(members),
            "n_coassemblies":
                len({
                    x["coassembly"]
                    for x in members
                }),
            "coassemblies":
                ";".join(sorted({
                    x["coassembly"]
                    for x in members
                })),
            "producers":
                ";".join(sorted({
                    x["producer"]
                    for x in members
                })),
            "symbols":
                ";".join(sorted({
                    x["element_symbol"]
                    for x in members
                    if x["element_symbol"]
                })),
            "classes":
                ";".join(sorted({
                    x["class"]
                    for x in members
                    if x["class"]
                })),
            "representation_vs_final18":
                ";".join(sorted({
                    x["representation_vs_final18"]
                    for x in members
                })),
        })


# ============================================================
# Resumen por coensamblaje
# ============================================================

by_group = defaultdict(list)

for row in all_rows:
    by_group[row["coassembly"]].append(row)

summary_group = (
    out
    / "91E_summary_by_coassembly.tsv"
)

with summary_group.open(
    "w",
    newline=""
) as fh:

    fields = [
        "coassembly",
        "producer",
        "biological_unit",
        "AMR_records",
        "unique_symbols",
        "unique_classes",
        "exact_sequence_match_final18",
        "same_symbol_final18_nonexact",
        "additional_symbol_relative_to_final18",
    ]

    w = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    w.writeheader()

    for group in groups:

        members = by_group[group]

        c = Counter(
            x["representation_vs_final18"]
            for x in members
        )

        w.writerow({
            "coassembly":
                group,
            "producer":
                group[0],
            "biological_unit":
                group[1:],
            "AMR_records":
                len(members),
            "unique_symbols":
                len({
                    x["element_symbol"]
                    for x in members
                    if x["element_symbol"]
                }),
            "unique_classes":
                len({
                    x["class"]
                    for x in members
                    if x["class"]
                }),
            "exact_sequence_match_final18":
                c[
                    "exact_sequence_match_final18"
                ],
            "same_symbol_final18_nonexact":
                c[
                    "same_AMR_symbol_final18_nonexact"
                ],
            "additional_symbol_relative_to_final18":
                c[
                    "additional_AMR_symbol_relative_to_final18"
                ],
        })


# ============================================================
# Resumen global
# ============================================================

type_counts = Counter(
    clean(x.get("Type"))
    for x in all_plus_rows
)

repr_counts = Counter(
    x["representation_vs_final18"]
    for x in all_rows
)

class_counts = Counter(
    x["class"]
    for x in all_rows
    if x["class"]
)

symbol_set = {
    x["element_symbol"]
    for x in all_rows
    if x["element_symbol"]
}

summary = [
    (
        "coassemblies_analyzed",
        len(groups)
    ),
    (
        "AMRFinderPlus_records_all_types",
        len(all_plus_rows)
    ),
    (
        "community_AMR_records",
        len(all_rows)
    ),
    (
        "community_unique_AMR_symbols",
        len(symbol_set)
    ),
    (
        "community_unique_exact_protein_sequences",
        len(protein_clusters)
    ),
    (
        "community_unique_exact_nucleotide_sequences",
        len(nucleotide_clusters)
    ),
    (
        "exact_sequence_match_final18",
        repr_counts[
            "exact_sequence_match_final18"
        ]
    ),
    (
        "same_AMR_symbol_final18_nonexact",
        repr_counts[
            "same_AMR_symbol_final18_nonexact"
        ]
    ),
    (
        "additional_AMR_symbol_relative_to_final18",
        repr_counts[
            "additional_AMR_symbol_relative_to_final18"
        ]
    ),
]

for k, v in sorted(
    type_counts.items()
):
    summary.append(
        (
            f"AMRFinder_type_{k}",
            v
        )
    )

summary_file = (
    out
    / "91E_community_resistome_summary.tsv"
)

with summary_file.open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    for k, v in summary:

        fh.write(
            f"{k}\t{v}\n"
        )


# ============================================================
# Resumen clases
# ============================================================

with (
    out
    / "91E_AMR_classes.tsv"
).open("w") as fh:

    fh.write(
        "AMR_class\trecords\n"
    )

    for cls, n in sorted(
        class_counts.items(),
        key=lambda x: (
            -x[1],
            x[0]
        )
    ):

        fh.write(
            f"{cls}\t{n}\n"
        )


# ============================================================
# Reporte consola
# ============================================================

print("=" * 78)
print("PASO 91E2 COMPLETADO")
print("=" * 78)

for k, v in summary:
    print(k, v, sep="\t")

print()
print("AMR classes:")

for cls, n in sorted(
    class_counts.items(),
    key=lambda x: (
        -x[1],
        x[0]
    )
):
    print(cls, n, sep="\t")

print()
print("AMR symbols:")
print(
    ";".join(
        sorted(symbol_set)
    )
)
PY

echo
echo "======================================================================"
echo "VALIDACION"
echo "======================================================================"

for F in \
    "${OUT}/91E_community_AMR_loci.tsv" \
    "${OUT}/91E_all_AMRFinderPlus_records.tsv" \
    "${OUT}/91E_exact_protein_clusters.tsv" \
    "${OUT}/91E_exact_nucleotide_clusters.tsv" \
    "${OUT}/91E_summary_by_coassembly.tsv" \
    "${OUT}/91E_community_resistome_summary.tsv" \
    "${OUT}/91E_AMR_classes.tsv"
do

    [[ -s "$F" ]] || {
        echo "ERROR: salida ausente/vacía: $F" >&2
        exit 4
    }

    echo "OK: $(basename "$F")"

done

{
    echo -e "field\tvalue"
    echo -e "analysis\tcommunity_resistome_consolidation"
    echo -e "source_coassemblies\tL1;L2;L3;M1;M2;M3"
    echo -e "AMRFinderPlus_version\t4.0.3"
    echo -e "AMRFinderPlus_database\t2025-07-16.1"
    echo -e "comparison_final18\texact_protein_and_nucleotide_SHA256_plus_gene_symbol"
    echo -e "date\t$(date -Iseconds)"
} > "${OUT}/91E_run_metadata.tsv"

echo
echo "======================================================================"
echo "RESUMEN GLOBAL"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/91E_community_resistome_summary.tsv"

echo
echo "======================================================================"
echo "RESUMEN POR COENSAMBLAJE"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/91E_summary_by_coassembly.tsv"

echo
echo "======================================================================"
echo "PASO 91E2 FINALIZÓ CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "======================================================================"
