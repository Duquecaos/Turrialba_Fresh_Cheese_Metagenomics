#!/usr/bin/env python3

from pathlib import Path
import os
import sys
import csv
import math

USER = os.environ["USER"]

ROOT = Path(
    f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"
)

INPUT = (
    ROOT
    / "70_incomplete_gene_catalog"
    / "incomplete_lt90_new.prodigal.faa"
)

OUT = (
    ROOT
    / "71_incomplete_eggnog"
)

CHUNKS = OUT / "chunks"

N_CHUNKS = 32

EXPECTED_PROTEINS = 763306

OUT.mkdir(
    parents=True,
    exist_ok=True
)

CHUNKS.mkdir(
    parents=True,
    exist_ok=True
)


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def read_fasta(path):
    name = None
    header = None
    seq = []

    with path.open("r", encoding="utf-8") as fh:

        for line in fh:

            line = line.rstrip("\n")

            if line.startswith(">"):

                if header is not None:
                    yield name, header, seq

                header = line
                name = line[1:].split()[0]
                seq = []

            else:

                seq.append(line)

        if header is not None:
            yield name, header, seq


if not INPUT.exists() or INPUT.stat().st_size == 0:
    die(f"Input faltante/vacío: {INPUT}")


records = list(
    read_fasta(INPUT)
)

n = len(records)

if n != EXPECTED_PROTEINS:
    die(
        f"Se esperaban {EXPECTED_PROTEINS} proteínas; "
        f"se encontraron {n}"
    )


ids = [
    r[0]
    for r in records
]

if len(ids) != len(set(ids)):
    die("Hay IDs de proteína duplicados.")


base = n // N_CHUNKS
remainder = n % N_CHUNKS

manifest = []

start = 0

for i in range(1, N_CHUNKS + 1):

    size = base + (
        1 if i <= remainder else 0
    )

    stop = start + size

    subset = records[start:stop]

    chunk_id = f"chunk_{i:02d}"

    outfile = (
        CHUNKS
        / f"{chunk_id}.faa"
    )

    with outfile.open(
        "w",
        encoding="utf-8"
    ) as fh:

        for name, header, seq_lines in subset:

            fh.write(
                header + "\n"
            )

            for seqline in seq_lines:
                fh.write(
                    seqline + "\n"
                )

    manifest.append({
        "array_index": i,
        "chunk_id": chunk_id,
        "protein_count": len(subset),
        "first_protein": subset[0][0],
        "last_protein": subset[-1][0],
        "fasta": str(outfile),
        "size_bytes": outfile.stat().st_size,
    })

    start = stop


if start != n:
    die(
        f"Error interno de partición: "
        f"{start} != {n}"
    )


manifest_path = (
    OUT / "chunk_manifest.tsv"
)

with manifest_path.open(
    "w",
    encoding="utf-8",
    newline=""
) as fh:

    fields = [
        "array_index",
        "chunk_id",
        "protein_count",
        "first_protein",
        "last_protein",
        "fasta",
        "size_bytes",
    ]

    w = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    w.writeheader()
    w.writerows(manifest)


total_check = sum(
    x["protein_count"]
    for x in manifest
)

if total_check != EXPECTED_PROTEINS:
    die(
        f"QC total falló: {total_check}"
    )


print("=" * 64)
print("PASO 86b-1 COMPLETADO")
print("=" * 64)

print(
    f"Proteínas totales:       {n}"
)

print(
    f"Número de chunks:        {N_CHUNKS}"
)

print(
    f"Proteínas/chunk mínimo:  "
    f"{min(x['protein_count'] for x in manifest)}"
)

print(
    f"Proteínas/chunk máximo:  "
    f"{max(x['protein_count'] for x in manifest)}"
)

print(
    f"Manifiesto: {manifest_path}"
)

print(
    "PASO 86b-1 FINALIZÓ CORRECTAMENTE; ES SEGURO SALIR."
)
