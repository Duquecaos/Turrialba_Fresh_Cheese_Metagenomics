#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import csv
import os
import sys

USER = os.environ["USER"]

ROOT = Path(
    f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"
)

EGGBASE = (
    ROOT
    / "71_incomplete_eggnog"
)

GENE_CONTEXT = (
    ROOT
    / "70_incomplete_gene_catalog"
    / "gene_to_contig_context.tsv"
)

OUT = (
    ROOT
    / "72_incomplete_eggnog_consolidated"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)

EXPECTED_GENES = 763306
EXPECTED_CHUNKS = 32
EXPECTED_ANNOTATIONS = 698054
EXPECTED_SEEDS = 698067


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def read_emapper_annotations(path):

    header = None

    with path.open(
        "r",
        encoding="utf-8",
        errors="replace"
    ) as fh:

        for line in fh:

            line = line.rstrip("\n")

            if not line:
                continue

            if line.startswith("#query\t"):

                header = line[1:].split("\t")
                continue

            if line.startswith("#"):
                continue

            if header is None:
                die(
                    f"No se encontró header #query en {path}"
                )

            values = line.split("\t")

            if len(values) < len(header):
                values += [""] * (
                    len(header)
                    - len(values)
                )

            if len(values) > len(header):
                die(
                    f"Más campos que header en {path}"
                )

            yield dict(
                zip(
                    header,
                    values
                )
            )


def read_seed_ids(path):

    with path.open(
        "r",
        encoding="utf-8",
        errors="replace"
    ) as fh:

        for line in fh:

            if not line.strip():
                continue

            if line.startswith("#"):
                continue

            yield line.split(
                "\t",
                1
            )[0]


def write_tsv(path, rows, fields):

    with path.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as fh:

        w = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore"
        )

        w.writeheader()
        w.writerows(rows)


# ============================================================
# VALIDAR 32 CHUNKS
# ============================================================

chunk_qc = []

annotation_by_gene = {}

seed_ids = set()

annotation_header = None


for i in range(
    1,
    EXPECTED_CHUNKS + 1
):

    chunk = f"chunk_{i:02d}"

    result = (
        EGGBASE
        / "results"
        / chunk
    )

    completed = (
        result
        / "COMPLETED.tsv"
    )

    annotations = (
        result
        / f"{chunk}.emapper.annotations"
    )

    seeds = (
        result
        / f"{chunk}.emapper.seed_orthologs"
    )

    hits = (
        result
        / f"{chunk}.emapper.hits"
    )

    for f in [
        completed,
        annotations,
        seeds,
        hits
    ]:

        if (
            not f.exists()
            or f.stat().st_size == 0
        ):
            die(
                f"Falta archivo requerido: {f}"
            )

    chunk_annotations = 0
    chunk_seeds = 0

    for row in read_emapper_annotations(
        annotations
    ):

        gid = row["query"]

        if gid in annotation_by_gene:
            die(
                f"Gene duplicado entre chunks: {gid}"
            )

        annotation_by_gene[gid] = row
        chunk_annotations += 1

        current_header = list(
            row.keys()
        )

        if annotation_header is None:
            annotation_header = current_header

        elif current_header != annotation_header:
            die(
                f"Header eggNOG inconsistente en {chunk}"
            )

    for gid in read_seed_ids(
        seeds
    ):

        if gid in seed_ids:
            die(
                f"Seed gene duplicado entre chunks: {gid}"
            )

        seed_ids.add(gid)
        chunk_seeds += 1

    chunk_qc.append({
        "chunk":
            chunk,
        "annotation_records":
            chunk_annotations,
        "seed_records":
            chunk_seeds,
        "seed_minus_annotation":
            (
                chunk_seeds
                - chunk_annotations
            ),
    })


if len(annotation_by_gene) != EXPECTED_ANNOTATIONS:
    die(
        f"Anotaciones esperadas={EXPECTED_ANNOTATIONS}; "
        f"observadas={len(annotation_by_gene)}"
    )


if len(seed_ids) != EXPECTED_SEEDS:
    die(
        f"Seeds esperados={EXPECTED_SEEDS}; "
        f"observados={len(seed_ids)}"
    )


# ============================================================
# CARGAR CONTEXTO DE LOS 763306 GENES
# ============================================================

if (
    not GENE_CONTEXT.exists()
    or GENE_CONTEXT.stat().st_size == 0
):
    die(
        f"Falta contexto génico: {GENE_CONTEXT}"
    )


with GENE_CONTEXT.open(
    "r",
    encoding="utf-8-sig",
    newline=""
) as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    context_fields = reader.fieldnames

    context_rows = list(
        reader
    )


if len(context_rows) != EXPECTED_GENES:
    die(
        f"Genes contexto esperados={EXPECTED_GENES}; "
        f"observados={len(context_rows)}"
    )


context_ids = [
    r["gene_id"]
    for r in context_rows
]


if len(context_ids) != len(set(context_ids)):
    die(
        "Hay gene_id duplicados en gene_to_contig_context.tsv"
    )


context_set = set(
    context_ids
)


# ============================================================
# QC IDs EGGNOG vs CATÁLOGO
# ============================================================

annotation_ids = set(
    annotation_by_gene
)

unknown_annotations = (
    annotation_ids
    - context_set
)

unknown_seeds = (
    seed_ids
    - context_set
)


if unknown_annotations:
    die(
        "Hay anotaciones eggNOG no presentes "
        f"en catálogo Prodigal: {len(unknown_annotations)}"
    )


if unknown_seeds:
    die(
        "Hay seeds no presentes "
        f"en catálogo Prodigal: {len(unknown_seeds)}"
    )


seed_without_annotation = sorted(
    seed_ids
    - annotation_ids
)

annotation_without_seed = sorted(
    annotation_ids
    - seed_ids
)


# ============================================================
# ARCHIVO CONSOLIDADO DE ANOTACIONES
# ============================================================

raw_fields = annotation_header

raw_rows = [
    annotation_by_gene[gid]
    for gid in sorted(
        annotation_by_gene
    )
]


write_tsv(
    OUT
    / "incomplete_lt90_eggnog_annotations_698054.tsv",
    raw_rows,
    raw_fields
)


# ============================================================
# INVENTARIO COMPLETO DE 763306 GENES
# ============================================================

annotation_fields_no_query = [
    x
    for x in annotation_header
    if x != "query"
]


output_annotation_fields = [
    f"eggnog_{x}"
    for x in annotation_fields_no_query
]


inventory_fields = (
    context_fields
    + [
        "has_seed_ortholog",
        "has_eggnog_annotation",
        "eggnog_status",
    ]
    + output_annotation_fields
)


inventory_path = (
    OUT
    / "gene_function_inventory_all_763306.tsv"
)


status_counter = Counter()
scope_counter = Counter()
partial_counter = Counter()


with inventory_path.open(
    "w",
    encoding="utf-8",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=inventory_fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore"
    )

    writer.writeheader()

    for context in context_rows:

        gid = context["gene_id"]

        has_seed = (
            gid in seed_ids
        )

        has_ann = (
            gid in annotation_by_gene
        )

        if has_ann:

            status = (
                "eggnog_annotated"
            )

        elif has_seed:

            status = (
                "seed_only_no_final_annotation"
            )

        else:

            status = (
                "no_seed_no_annotation"
            )

        outrow = dict(
            context
        )

        outrow[
            "has_seed_ortholog"
        ] = int(
            has_seed
        )

        outrow[
            "has_eggnog_annotation"
        ] = int(
            has_ann
        )

        outrow[
            "eggnog_status"
        ] = status

        if has_ann:

            ann = annotation_by_gene[
                gid
            ]

            for field in (
                annotation_fields_no_query
            ):

                outrow[
                    f"eggnog_{field}"
                ] = ann.get(
                    field,
                    ""
                )

        else:

            for field in (
                annotation_fields_no_query
            ):

                outrow[
                    f"eggnog_{field}"
                ] = ""

        writer.writerow(
            outrow
        )

        status_counter[
            status
        ] += 1

        scope_counter[
            (
                context[
                    "gene_evidence_scope"
                ],
                status
            )
        ] += 1

        partial_counter[
            (
                context[
                    "partial"
                ],
                status
            )
        ] += 1


# ============================================================
# IDs ESPECIALES
# ============================================================

write_tsv(
    OUT
    / "seed_without_final_annotation.tsv",
    [
        {"gene_id": x}
        for x in seed_without_annotation
    ],
    ["gene_id"]
)


write_tsv(
    OUT
    / "annotation_without_seed.tsv",
    [
        {"gene_id": x}
        for x in annotation_without_seed
    ],
    ["gene_id"]
)


genes_without_annotation = [
    {
        "gene_id":
            r["gene_id"],

        "contig_id":
            r["contig_id"],

        "partial":
            r["partial"],

        "aa_length":
            r["aa_length"],

        "gene_evidence_scope":
            r["gene_evidence_scope"],

        "strongest_incomplete_context":
            r[
                "strongest_incomplete_context"
            ],

        "max_member_bin_completeness":
            r[
                "max_member_bin_completeness"
            ],

        "max_member_bin_contamination":
            r[
                "max_member_bin_contamination"
            ],

        "eggnog_status":
            (
                "seed_only_no_final_annotation"
                if r["gene_id"] in seed_ids
                else "no_seed_no_annotation"
            ),
    }
    for r in context_rows
    if r["gene_id"]
    not in annotation_ids
]


write_tsv(
    OUT
    / "genes_without_final_eggnog_annotation.tsv",
    genes_without_annotation,
    [
        "gene_id",
        "contig_id",
        "partial",
        "aa_length",
        "gene_evidence_scope",
        "strongest_incomplete_context",
        "max_member_bin_completeness",
        "max_member_bin_contamination",
        "eggnog_status",
    ]
)


# ============================================================
# QC POR CHUNK
# ============================================================

write_tsv(
    OUT
    / "QC_eggnog_by_chunk.tsv",
    chunk_qc,
    [
        "chunk",
        "annotation_records",
        "seed_records",
        "seed_minus_annotation",
    ]
)


# ============================================================
# QC POR NIVEL DE EVIDENCIA
# ============================================================

scope_rows = []

for (
    scope,
    status
), n in sorted(
    scope_counter.items()
):

    scope_rows.append({
        "gene_evidence_scope":
            scope,
        "eggnog_status":
            status,
        "n_genes":
            n,
    })


write_tsv(
    OUT
    / "QC_eggnog_by_gene_evidence_scope.tsv",
    scope_rows,
    [
        "gene_evidence_scope",
        "eggnog_status",
        "n_genes",
    ]
)


partial_rows = []

for (
    partial,
    status
), n in sorted(
    partial_counter.items()
):

    partial_rows.append({
        "partial":
            partial,
        "eggnog_status":
            status,
        "n_genes":
            n,
    })


write_tsv(
    OUT
    / "QC_eggnog_by_partial_status.tsv",
    partial_rows,
    [
        "partial",
        "eggnog_status",
        "n_genes",
    ]
)


# ============================================================
# RESUMEN GLOBAL
# ============================================================

n_annotated = status_counter[
    "eggnog_annotated"
]

n_seed_only = status_counter[
    "seed_only_no_final_annotation"
]

n_no_seed = status_counter[
    "no_seed_no_annotation"
]

n_without_annotation = (
    n_seed_only
    + n_no_seed
)


summary = [
    (
        "total_predicted_genes",
        EXPECTED_GENES
    ),
    (
        "genes_with_seed_ortholog",
        len(seed_ids)
    ),
    (
        "genes_with_final_eggnog_annotation",
        n_annotated
    ),
    (
        "seed_only_no_final_annotation",
        n_seed_only
    ),
    (
        "no_seed_no_annotation",
        n_no_seed
    ),
    (
        "genes_without_final_annotation",
        n_without_annotation
    ),
    (
        "annotation_rate_pct",
        round(
            100
            * n_annotated
            / EXPECTED_GENES,
            4
        )
    ),
    (
        "seed_rate_pct",
        round(
            100
            * len(seed_ids)
            / EXPECTED_GENES,
            4
        )
    ),
    (
        "annotation_without_seed",
        len(
            annotation_without_seed
        )
    ),
    (
        "chunks_integrated",
        EXPECTED_CHUNKS
    ),
]


with (
    OUT
    / "eggnog_consolidation_summary.tsv"
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
PASO 86c - CONSOLIDACION EGGNOG DE BINS INCOMPLETOS
====================================================

Entrada:
763306 proteínas predichas por Prodigal sobre los contigs
únicos nuevos de bins <90% de completitud.

eggNOG:
- eggNOG-mapper 2.1.12
- eggNOG DB 5.0.2
- DIAMOND
- 32 chunks

Este paso integra:
gene -> contig -> bin/contexto -> eggNOG

Estados:
1. eggnog_annotated
   Seed ortholog + anotación final eggNOG.

2. seed_only_no_final_annotation
   Tiene seed ortholog pero no registro final en annotations.

3. no_seed_no_annotation
   Sin seed ortholog ni anotación final.

Los genes sin anotación se mantienen explícitamente en
gene_function_inventory_all_763306.tsv.

La ausencia de anotación eggNOG NO equivale a ausencia de función.

Los genes de bins incompletos no se interpretarán mediante
ausencia de marcadores como evidencia biológica negativa.

Las bacteriocinas/RiPP se analizarán además mediante una rama
independiente del catálogo Prodigal/eggNOG.
"""

with (
    OUT
    / "README_step86c.txt"
).open(
    "w",
    encoding="utf-8"
) as fh:

    fh.write(
        readme.strip()
        + "\n"
    )


# ============================================================
# FINAL
# ============================================================

print("=" * 70)
print("PASO 86c COMPLETADO")
print("=" * 70)

for metric, value in summary:
    print(
        f"{metric:42s} {value}"
    )

print()

print(
    f"Salida: {OUT}"
)

print(
    "PASO 86c FINALIZÓ CORRECTAMENTE; ES SEGURO SALIR."
)
