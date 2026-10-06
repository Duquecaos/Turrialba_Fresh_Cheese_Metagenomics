#!/usr/bin/env python3

from pathlib import Path
from collections import Counter, defaultdict
import csv
import os
import re
import statistics
import sys

USER = os.environ["USER"]

ROOT = Path(
    f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"
)

CONTIG_CATALOG = (
    ROOT
    / "69_incomplete_unique_contigs"
    / "unique_contig_catalog_lt90.tsv"
)

INPUT_FASTA = (
    ROOT
    / "69_incomplete_unique_contigs"
    / "unique_contigs_lt90_not_in_final18.fna"
)

OUT = (
    ROOT
    / "70_incomplete_gene_catalog"
)

PROTEINS = (
    OUT
    / "incomplete_lt90_new.prodigal.faa"
)

GENES = (
    OUT
    / "incomplete_lt90_new.prodigal.fna"
)

GFF = (
    OUT
    / "incomplete_lt90_new.prodigal.gff"
)


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def read_fasta(path):
    name = None
    seq = []

    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)

                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())

        if name is not None:
            yield name, "".join(seq)


def parse_attrs(text):
    d = {}

    for item in text.split(";"):
        item = item.strip()

        if not item:
            continue

        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v

    return d


def write_tsv(path, rows, fields):
    with path.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(rows)


# ============================================================
# VALIDAR INPUTS
# ============================================================

for f in [
    CONTIG_CATALOG,
    INPUT_FASTA,
    PROTEINS,
    GENES,
    GFF
]:
    if not f.exists() or f.stat().st_size == 0:
        die(f"Archivo requerido faltante/vacío: {f}")


# ============================================================
# CONTIGS DE ENTRADA
# ============================================================

input_contigs = {}

for cid, seq in read_fasta(INPUT_FASTA):
    if cid in input_contigs:
        die(f"Contig duplicado en FASTA: {cid}")

    input_contigs[cid] = len(seq)


if len(input_contigs) != 140276:
    die(
        "Se esperaban 140276 contigs nuevos; "
        f"se encontraron {len(input_contigs)}"
    )


# ============================================================
# CATÁLOGO DE CONTIGS
# ============================================================

catalog = {}

with CONTIG_CATALOG.open(
    "r",
    encoding="utf-8-sig",
    newline=""
) as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for r in reader:
        cid = r["unique_contig_id"]

        if cid in input_contigs:
            catalog[cid] = r


if len(catalog) != len(input_contigs):
    missing = sorted(
        set(input_contigs)
        - set(catalog)
    )

    die(
        "Hay contigs del FASTA sin contexto en catálogo. "
        f"N={len(missing)}. Ejemplos={missing[:10]}"
    )


# ============================================================
# PROTEÍNAS
# ============================================================

protein_lengths = {}
gene_to_contig_from_faa = {}

pattern = re.compile(
    r"^(ICONTIG\d{9})_(\d+)$"
)


for gid, seq in read_fasta(PROTEINS):

    if gid in protein_lengths:
        die(f"Proteína duplicada: {gid}")

    m = pattern.match(gid)

    if not m:
        die(
            f"ID de proteína no reconocido: {gid}"
        )

    cid = m.group(1)

    if cid not in input_contigs:
        die(
            f"Proteína {gid} apunta a contig inexistente {cid}"
        )

    clean_seq = seq.rstrip("*")

    protein_lengths[gid] = len(clean_seq)

    gene_to_contig_from_faa[gid] = cid


# ============================================================
# GENES NUCLEOTÍDICOS
# ============================================================

gene_nt_lengths = {}

for gid, seq in read_fasta(GENES):

    if gid in gene_nt_lengths:
        die(f"Gen duplicado en nucleotide FASTA: {gid}")

    gene_nt_lengths[gid] = len(seq)


if set(gene_nt_lengths) != set(protein_lengths):

    only_protein = (
        set(protein_lengths)
        - set(gene_nt_lengths)
    )

    only_nt = (
        set(gene_nt_lengths)
        - set(protein_lengths)
    )

    die(
        "FASTA de proteína y genes no coinciden. "
        f"solo_proteina={len(only_protein)}; "
        f"solo_nucleotido={len(only_nt)}"
    )


# ============================================================
# PARSEAR GFF
# ============================================================

gff_rows = []
counter_by_contig = defaultdict(int)


with GFF.open(
    "r",
    encoding="utf-8"
) as fh:

    for line in fh:

        if not line.strip():
            continue

        if line.startswith("#"):
            continue

        parts = line.rstrip("\n").split("\t")

        if len(parts) != 9:
            die(
                "Línea GFF inválida: "
                + line[:200]
            )

        seqid = parts[0]
        feature_type = parts[2]

        if feature_type != "CDS":
            continue

        if seqid not in input_contigs:
            die(
                f"GFF refiere contig desconocido: {seqid}"
            )

        counter_by_contig[seqid] += 1
        local_n = counter_by_contig[seqid]

        gid = f"{seqid}_{local_n}"

        attrs = parse_attrs(parts[8])

        if gid not in protein_lengths:
            die(
                f"Gen GFF no encontrado en proteína: {gid}"
            )

        if (
            gene_to_contig_from_faa[gid]
            != seqid
        ):
            die(
                f"Inconsistencia gen/contig para {gid}"
            )

        start = int(parts[3])
        end = int(parts[4])

        partial = attrs.get(
            "partial",
            ""
        )

        contig_info = catalog[seqid]

        strongest_context = contig_info[
            "strongest_incomplete_context"
        ]

        if strongest_context == "present_in_50_89_bin":
            evidence_scope = (
                "bin_context_50_89_positive_evidence"
            )
        else:
            evidence_scope = (
                "lt50_contig_gene_centric_only"
            )

        row = {
            "gene_id":
                gid,

            "contig_id":
                seqid,

            "gene_number_on_contig":
                local_n,

            "start":
                start,

            "end":
                end,

            "strand":
                parts[6],

            "nt_length":
                gene_nt_lengths[gid],

            "aa_length":
                protein_lengths[gid],

            "partial":
                partial,

            "start_type":
                attrs.get(
                    "start_type",
                    ""
                ),

            "rbs_motif":
                attrs.get(
                    "rbs_motif",
                    ""
                ),

            "rbs_spacer":
                attrs.get(
                    "rbs_spacer",
                    ""
                ),

            "gc_cont":
                attrs.get(
                    "gc_cont",
                    ""
                ),

            "conf":
                attrs.get(
                    "conf",
                    ""
                ),

            "score":
                attrs.get(
                    "score",
                    ""
                ),

            "cscore":
                attrs.get(
                    "cscore",
                    ""
                ),

            "sscore":
                attrs.get(
                    "sscore",
                    ""
                ),

            "contig_length_bp":
                input_contigs[seqid],

            "n_bin_memberships":
                contig_info[
                    "n_bin_memberships"
                ],

            "coassemblies":
                contig_info[
                    "coassemblies"
                ],

            "producers":
                contig_info[
                    "producers"
                ],

            "binners":
                contig_info[
                    "binners"
                ],

            "quality_groups":
                contig_info[
                    "quality_groups"
                ],

            "strongest_incomplete_context":
                strongest_context,

            "min_member_bin_completeness":
                contig_info[
                    "min_member_bin_completeness"
                ],

            "max_member_bin_completeness":
                contig_info[
                    "max_member_bin_completeness"
                ],

            "min_member_bin_contamination":
                contig_info[
                    "min_member_bin_contamination"
                ],

            "max_member_bin_contamination":
                contig_info[
                    "max_member_bin_contamination"
                ],

            "gene_evidence_scope":
                evidence_scope,
        }

        gff_rows.append(row)


if len(gff_rows) != len(protein_lengths):
    die(
        "Número de CDS del GFF no coincide con proteínas: "
        f"GFF={len(gff_rows)} "
        f"FAA={len(protein_lengths)}"
    )


if {
    r["gene_id"]
    for r in gff_rows
} != set(protein_lengths):

    die(
        "IDs de GFF y proteínas no coinciden."
    )


# ============================================================
# TABLA GEN -> CONTIG
# ============================================================

fields = [
    "gene_id",
    "contig_id",
    "gene_number_on_contig",
    "start",
    "end",
    "strand",
    "nt_length",
    "aa_length",
    "partial",
    "start_type",
    "rbs_motif",
    "rbs_spacer",
    "gc_cont",
    "conf",
    "score",
    "cscore",
    "sscore",
    "contig_length_bp",
    "n_bin_memberships",
    "coassemblies",
    "producers",
    "binners",
    "quality_groups",
    "strongest_incomplete_context",
    "min_member_bin_completeness",
    "max_member_bin_completeness",
    "min_member_bin_contamination",
    "max_member_bin_contamination",
    "gene_evidence_scope",
]


write_tsv(
    OUT / "gene_to_contig_context.tsv",
    gff_rows,
    fields
)


# ============================================================
# CONTIGS SIN GENES
# ============================================================

contigs_with_genes = set(
    r["contig_id"]
    for r in gff_rows
)

contigs_without = sorted(
    set(input_contigs)
    - contigs_with_genes
)


without_rows = []

for cid in contigs_without:

    r = catalog[cid]

    without_rows.append({
        "contig_id":
            cid,

        "length_bp":
            input_contigs[cid],

        "strongest_incomplete_context":
            r[
                "strongest_incomplete_context"
            ],

        "coassemblies":
            r[
                "coassemblies"
            ],

        "binners":
            r[
                "binners"
            ],
    })


write_tsv(
    OUT / "contigs_without_predicted_genes.tsv",
    without_rows,
    [
        "contig_id",
        "length_bp",
        "strongest_incomplete_context",
        "coassemblies",
        "binners",
    ]
)


# ============================================================
# RESUMEN
# ============================================================

partial_counter = Counter(
    r["partial"]
    for r in gff_rows
)

scope_counter = Counter(
    r["gene_evidence_scope"]
    for r in gff_rows
)


aa_lengths = list(
    protein_lengths.values()
)

total_input_bp = sum(
    input_contigs.values()
)


summary = [
    (
        "input_unique_new_contigs",
        len(input_contigs)
    ),
    (
        "input_bp",
        total_input_bp
    ),
    (
        "contigs_with_predicted_genes",
        len(contigs_with_genes)
    ),
    (
        "contigs_without_predicted_genes",
        len(contigs_without)
    ),
    (
        "predicted_genes",
        len(gff_rows)
    ),
    (
        "predicted_proteins",
        len(protein_lengths)
    ),
    (
        "complete_genes_partial_00",
        partial_counter.get("00", 0)
    ),
    (
        "partial_5prime_10",
        partial_counter.get("10", 0)
    ),
    (
        "partial_3prime_01",
        partial_counter.get("01", 0)
    ),
    (
        "partial_both_11",
        partial_counter.get("11", 0)
    ),
    (
        "genes_with_50_89_context",
        scope_counter.get(
            "bin_context_50_89_positive_evidence",
            0
        )
    ),
    (
        "genes_lt50_only",
        scope_counter.get(
            "lt50_contig_gene_centric_only",
            0
        )
    ),
    (
        "median_protein_length_aa",
        round(
            statistics.median(
                aa_lengths
            ),
            2
        ) if aa_lengths else 0
    ),
    (
        "mean_protein_length_aa",
        round(
            statistics.mean(
                aa_lengths
            ),
            2
        ) if aa_lengths else 0
    ),
    (
        "genes_per_Mb",
        round(
            len(gff_rows)
            / (
                total_input_bp
                / 1_000_000
            ),
            3
        )
    ),
]


with (
    OUT
    / "gene_catalog_summary.tsv"
).open(
    "w",
    encoding="utf-8",
    newline=""
) as fh:

    writer = csv.writer(
        fh,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow(
        [
            "metric",
            "value"
        ]
    )

    writer.writerows(
        summary
    )


# ============================================================
# README
# ============================================================

readme = """
PASO 86a - CATALOGO GENICO DE CONTIGS INCOMPLETOS
==================================================

Entrada:
140276 contigs únicos no representados exactamente en los
18 MAGs finales.

Predicción:
Prodigal 2.6.3, modo metagenómico (-p meta).

Este paso NO ejecuta eggNOG.

Archivos:
- incomplete_lt90_new.prodigal.faa
- incomplete_lt90_new.prodigal.fna
- incomplete_lt90_new.prodigal.gff
- gene_to_contig_context.tsv
- contigs_without_predicted_genes.tsv
- gene_catalog_summary.tsv

Interpretación:
Los genes presentes en contigs asociados con bins de 50-89.99%
pueden conservar contexto de bin para evidencia positiva.

Los genes encontrados únicamente en bins <50% se interpretarán
principalmente a nivel gen/contig.

La ausencia de un gen en un bin incompleto NO se interpretará
como ausencia biológica.

La predicción Prodigal no será la única estrategia usada para
buscar precursores pequeños de bacteriocinas/RiPP.
"""

with (
    OUT / "README_step86a.txt"
).open(
    "w",
    encoding="utf-8"
) as fh:

    fh.write(
        readme.strip()
        + "\n"
    )


print("=" * 68)
print("PASO 86a - POSTPROCESAMIENTO COMPLETADO")
print("=" * 68)

for metric, value in summary:
    print(
        f"{metric:42s} {value}"
    )

print()

print(
    f"Salida: {OUT}"
)

print(
    "PASO 86a FINALIZÓ CORRECTAMENTE; ES SEGURO SALIR."
)
