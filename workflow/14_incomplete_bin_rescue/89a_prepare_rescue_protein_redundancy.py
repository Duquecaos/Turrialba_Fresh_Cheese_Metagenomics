#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import csv
import hashlib
import os
import sys


USER = os.environ["USER"]

ROOT = Path(f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba")

IN87 = ROOT / "73_incomplete_functional_rescue"
IN88 = ROOT / "74_incomplete_functional_evidence_context"
MASTER86 = (
    ROOT
    / "72_incomplete_eggnog_consolidated"
    / "gene_function_inventory_all_763306.tsv"
)
FAA_INCOMPLETE = (
    ROOT
    / "70_incomplete_gene_catalog"
    / "incomplete_lt90_new.prodigal.faa"
)

FINAL_INV = (
    ROOT
    / "53_metabolic_annotation_consolidated"
    / "gene_function_inventory_all_predicted.tsv"
)
FINAL_FAA = (
    ROOT
    / "23_eggnog_rep18"
    / "results"
    / "representative18.emapper.genepred.fasta"
)

OUT = ROOT / "75_incomplete_functional_sequence_redundancy"
OUT.mkdir(parents=True, exist_ok=True)

OUT_FAA = OUT / "selected_rescue_functional_proteins.faa"
OUT_INV = OUT / "89A_selected_rescue_functional_genes.tsv"
OUT_MAP = OUT / "final18_gene_to_MAG.tsv"
OUT_SUM = OUT / "89A_preparation_summary.tsv"


# ------------------------------------------------------------------
# Utilidades
# ------------------------------------------------------------------

def split_ids(value):
    if value is None:
        return []
    value = value.strip()
    if not value or value in {"-", "NA", "None"}:
        return []
    return [
        x.strip()
        for x in value.split(";")
        if x.strip() and x.strip() not in {"-", "NA", "None"}
    ]


def read_tsv(path):
    with path.open() as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def fasta_records(path):
    ident = None
    seq = []

    with path.open() as fh:
        for line in fh:
            line = line.rstrip()

            if not line:
                continue

            if line.startswith(">"):
                if ident is not None:
                    yield ident, "".join(seq).rstrip("*")

                ident = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())

    if ident is not None:
        yield ident, "".join(seq).rstrip("*")


def sha256_seq(seq):
    return hashlib.sha256(seq.encode()).hexdigest()


# ------------------------------------------------------------------
# Tablas funcionales del Paso 87
# ------------------------------------------------------------------

table_configs = {
    "87A_cheese_metabolic_marker_hits.tsv": {
        "id_fields": ["gene_id"],
        "label_fields": ["block", "marker"],
    },
    "87B_cheese_metabolic_systems_on_contigs.tsv": {
        "id_fields": ["gene_ids"],
        "label_fields": ["system", "status"],
    },
    "87C_biogenic_amine_candidates.tsv": {
        "id_fields": ["gene_id"],
        "label_fields": [
            "potential_product",
            "evidence_level",
            "context_support",
        ],
    },
    "87D_polyamine_arg_stress_genes.tsv": {
        "id_fields": ["gene_id"],
        "label_fields": ["block", "marker"],
    },
    "87E_Opp_Dpp_components.tsv": {
        "id_fields": ["gene_id"],
        "label_fields": ["system", "component"],
    },
    "87F_Opp_Dpp_clusters.tsv": {
        "id_fields": ["genes"],
        "label_fields": ["system", "cluster_status"],
    },
    "87G_POT_Dtp_transporters.tsv": {
        "id_fields": ["gene_id"],
        "label_fields": ["block", "marker"],
    },
    "87H_key_peptidases.tsv": {
        "id_fields": ["gene_id"],
        "label_fields": ["peptidase_type"],
    },
    "87I_surface_proteinase_CEP_candidates.tsv": {
        "id_fields": ["gene_id"],
        "label_fields": ["rescue_tier_87"],
    },
    "87J_curated_lipolysis_candidates.tsv": {
        "id_fields": ["gene_id"],
        "label_fields": ["curated_lipolysis_tier"],
    },
    "87K_curated_amino_acid_aroma_markers.tsv": {
        "id_fields": ["gene_id"],
        "label_fields": [
            "marker",
            "curated_aroma_tier",
        ],
    },
    "87L_EPS_capsule_loci.tsv": {
        "id_fields": [
            "marker_gene_ids",
            "glycosyltransferase_gene_ids",
        ],
        "label_fields": ["curated_tier"],
    },
    "87M_stress_systems_on_contigs.tsv": {
        "id_fields": ["gene_ids"],
        "label_fields": [
            "block",
            "system",
            "status",
        ],
    },
}


gene_tables = defaultdict(set)
gene_labels = defaultdict(set)

for filename, cfg in table_configs.items():

    path = IN87 / filename

    if not path.exists():
        raise SystemExit(f"ERROR: falta {path}")

    for row in read_tsv(path):

        label_parts = []

        for col in cfg["label_fields"]:
            val = (row.get(col) or "").strip()

            if val and val not in {"-", "NA", "None"}:
                label_parts.append(f"{col}={val}")

        label = "|".join(label_parts)

        ids = []

        for field in cfg["id_fields"]:
            ids.extend(split_ids(row.get(field, "")))

        for gene_id in ids:
            gene_tables[gene_id].add(filename)

            if label:
                gene_labels[gene_id].add(label)


selected = set(gene_tables)

if not selected:
    raise SystemExit(
        "ERROR: no se recuperaron genes funcionales desde las tablas 87."
    )


# ------------------------------------------------------------------
# Master incompleto 86c
# ------------------------------------------------------------------

master = {}

with MASTER86.open() as fh:
    r = csv.DictReader(fh, delimiter="\t")

    for row in r:
        gid = row["gene_id"]

        if gid in selected:
            master[gid] = row


missing_master = sorted(selected - set(master))

if missing_master:
    print(
        f"ERROR: {len(missing_master)} genes seleccionados no están "
        "en el master 86c.",
        file=sys.stderr,
    )

    for x in missing_master[:20]:
        print(x, file=sys.stderr)

    raise SystemExit(2)


# ------------------------------------------------------------------
# Proteínas incompletas seleccionadas
# ------------------------------------------------------------------

selected_seqs = {}

for ident, seq in fasta_records(FAA_INCOMPLETE):

    if ident in selected:
        selected_seqs[ident] = seq


missing_faa = sorted(selected - set(selected_seqs))

if missing_faa:
    print(
        f"ERROR: faltan {len(missing_faa)} proteínas seleccionadas "
        "en el FAA incompleto.",
        file=sys.stderr,
    )

    for x in missing_faa[:20]:
        print(x, file=sys.stderr)

    raise SystemExit(2)


# ------------------------------------------------------------------
# Inventario final18
# ------------------------------------------------------------------

final_gene_to_mag = {}

with FINAL_INV.open() as fh:
    r = csv.DictReader(fh, delimiter="\t")

    for row in r:
        final_gene_to_mag[row["gene_id"]] = row["MAG"]


# ------------------------------------------------------------------
# Proteoma final18 + hashes
# ------------------------------------------------------------------

final_seq_hash = defaultdict(list)
final_fasta_ids = set()

for ident, seq in fasta_records(FINAL_FAA):

    final_fasta_ids.add(ident)
    final_seq_hash[sha256_seq(seq)].append(ident)


overlap = len(final_fasta_ids & set(final_gene_to_mag))

if overlap < 0.95 * len(final_gene_to_mag):

    print(
        "ERROR: los IDs del FASTA final18 no corresponden de forma "
        "suficiente con el inventario funcional.",
        file=sys.stderr,
    )

    print(
        f"Inventario final18: {len(final_gene_to_mag)}",
        file=sys.stderr,
    )

    print(
        f"IDs FASTA final18: {len(final_fasta_ids)}",
        file=sys.stderr,
    )

    print(
        f"Solapamiento: {overlap}",
        file=sys.stderr,
    )

    raise SystemExit(2)


# ------------------------------------------------------------------
# FASTA funcional rescatado
# ------------------------------------------------------------------

with OUT_FAA.open("w") as out:

    for gid in sorted(selected):
        seq = selected_seqs[gid]

        out.write(f">{gid}\n")

        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")


# ------------------------------------------------------------------
# Mapa final18
# ------------------------------------------------------------------

with OUT_MAP.open("w", newline="") as out:

    w = csv.writer(out, delimiter="\t")

    w.writerow([
        "gene_id",
        "MAG",
    ])

    for gid in sorted(final_gene_to_mag):
        w.writerow([
            gid,
            final_gene_to_mag[gid],
        ])


# ------------------------------------------------------------------
# Inventario genes seleccionados
# ------------------------------------------------------------------

metadata_cols = [
    "contig_id",
    "aa_length",
    "partial",
    "gene_evidence_scope",
    "strongest_incomplete_context",
    "n_bin_memberships",
    "coassemblies",
    "producers",
    "binners",
    "quality_groups",
    "min_member_bin_completeness",
    "max_member_bin_completeness",
    "min_member_bin_contamination",
    "max_member_bin_contamination",
    "eggnog_Preferred_name",
    "eggnog_Description",
    "eggnog_EC",
    "eggnog_KEGG_ko",
    "eggnog_PFAMs",
]

fieldnames = [
    "gene_id",
    *metadata_cols,
    "evidence_tables",
    "functional_labels",
    "exact_final18_protein_match_n",
    "exact_final18_gene_ids",
    "exact_final18_MAGs",
]

n_exact = 0

with OUT_INV.open("w", newline="") as out:

    w = csv.DictWriter(
        out,
        delimiter="\t",
        fieldnames=fieldnames,
        lineterminator="\n",
    )

    w.writeheader()

    for gid in sorted(selected):

        row = master[gid]
        seq_hash = sha256_seq(selected_seqs[gid])

        exact_ids = sorted(final_seq_hash.get(seq_hash, []))

        exact_mags = sorted({
            final_gene_to_mag[x]
            for x in exact_ids
            if x in final_gene_to_mag
        })

        if exact_ids:
            n_exact += 1

        rec = {
            "gene_id": gid,
            **{
                c: row.get(c, "")
                for c in metadata_cols
            },
            "evidence_tables":
                ";".join(sorted(gene_tables[gid])),
            "functional_labels":
                ";".join(sorted(gene_labels[gid])),
            "exact_final18_protein_match_n":
                len(exact_ids),
            "exact_final18_gene_ids":
                ";".join(exact_ids),
            "exact_final18_MAGs":
                ";".join(exact_mags),
        }

        w.writerow(rec)


# ------------------------------------------------------------------
# Resumen
# ------------------------------------------------------------------

contigs = {
    master[x]["contig_id"]
    for x in selected
}

with OUT_SUM.open("w", newline="") as out:

    w = csv.writer(out, delimiter="\t")

    w.writerow(["metric", "value"])

    w.writerow([
        "selected_functional_rescue_genes",
        len(selected),
    ])

    w.writerow([
        "selected_functional_rescue_contigs",
        len(contigs),
    ])

    w.writerow([
        "selected_proteins_written",
        len(selected_seqs),
    ])

    w.writerow([
        "final18_inventory_genes",
        len(final_gene_to_mag),
    ])

    w.writerow([
        "final18_fasta_proteins",
        len(final_fasta_ids),
    ])

    w.writerow([
        "final18_inventory_fasta_ID_overlap",
        overlap,
    ])

    w.writerow([
        "rescue_genes_exactly_identical_to_final18_protein",
        n_exact,
    ])


print("=" * 72)
print("PASO 89a - PREPARACION")
print("=" * 72)
print(f"Genes funcionales rescatados seleccionados : {len(selected)}")
print(f"Contigs representados                      : {len(contigs)}")
print(f"Proteínas escritas                         : {len(selected_seqs)}")
print(f"Proteínas final18                          : {len(final_fasta_ids)}")
print(f"Exactamente idénticas a final18            : {n_exact}")
print()
print(f"Salida: {OUT}")
print()
print("PASO 89a PREPARADO CORRECTAMENTE.")
print("ES SEGURO SALIR.")
