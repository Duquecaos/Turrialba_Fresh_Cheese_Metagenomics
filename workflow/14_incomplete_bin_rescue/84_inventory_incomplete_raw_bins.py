#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict, Counter
import csv
import gzip
import hashlib
import math
import os
import sys

USER = os.environ["USER"]

ROOT = Path(
    f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"
)

OUT = ROOT / "68_incomplete_bin_inventory"
OUT.mkdir(parents=True, exist_ok=True)

BINROOT = ROOT / "08_binning"

FINAL18_DIR = (
    ROOT
    / "22_final_representative_mags"
    / "fastas"
)

EXPECTED_RAW_BINS = 1006
COASSEMBLIES = ["L1", "L2", "L3", "M1", "M2", "M3"]


# ============================================================
# HELPERS
# ============================================================

def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def normalize_fasta_name(name):
    x = name
    changed = True

    while changed:
        changed = False
        for suffix in [".gz", ".fasta", ".fna", ".fa"]:
            if x.lower().endswith(suffix):
                x = x[:-len(suffix)]
                changed = True

    return x


def is_fasta(path):
    n = path.name.lower()

    return (
        n.endswith(".fa")
        or n.endswith(".fna")
        or n.endswith(".fasta")
        or n.endswith(".fa.gz")
        or n.endswith(".fna.gz")
        or n.endswith(".fasta.gz")
    )


def open_text(path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("r")


def fasta_stats(path):

    lengths = []
    contig_tokens = []
    gc = 0
    total = 0

    seq = []

    def process_sequence(parts):
        nonlocal gc, total

        if not parts:
            return

        s = "".join(parts).upper()

        if len(s) == 0:
            return

        L = len(s)
        lengths.append(L)
        total += L
        gc += s.count("G") + s.count("C")

        digest = hashlib.sha256(
            s.encode()
        ).hexdigest()

        contig_tokens.append(
            f"{L}:{digest}"
        )

    with open_text(path) as fh:

        for line in fh:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                process_sequence(seq)
                seq = []
            else:
                seq.append(line)

        process_sequence(seq)

    if not lengths:
        raise ValueError(
            f"FASTA sin secuencias: {path}"
        )

    lengths_sorted = sorted(
        lengths,
        reverse=True
    )

    half = total / 2.0
    cumulative = 0
    n50 = 0

    for L in lengths_sorted:
        cumulative += L

        if cumulative >= half:
            n50 = L
            break

    content_hash = hashlib.sha256(
        "\n".join(
            sorted(contig_tokens)
        ).encode()
    ).hexdigest()

    return {
        "n_contigs": len(lengths),
        "total_bp": total,
        "n50_bp": n50,
        "max_contig_bp": max(lengths),
        "gc_pct": (100.0 * gc / total),
        "sequence_content_sha256": content_hash,
    }


def read_tsv(path):

    if not path.exists():
        die(
            f"No existe tabla requerida: {path}"
        )

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        return list(reader)


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


def quality_group(comp):

    if comp >= 90:
        return "A_ge90"

    if comp >= 50:
        return "B_50_89.99"

    return "C_lt50"


def contamination_group(cont):

    if cont <= 5:
        return "contam_le5"

    if cont <= 10:
        return "contam_gt5_le10"

    if cont <= 20:
        return "contam_gt10_le20"

    return "contam_gt20"


def screening_mode(comp, cont):

    if comp >= 90:
        return "reference_ge90_not_new_target"

    if comp >= 50 and cont <= 10:
        return "incomplete_bin_positive_gene_evidence"

    if comp >= 50 and cont > 10:
        return "positive_gene_evidence_taxonomic_caution"

    if comp < 50 and cont <= 10:
        return "fragmentary_positive_gene_only"

    return "fragmentary_high_contam_contig_level_only"


# ============================================================
# INDEXAR FASTA DE LOS TRES BINNERS
# ============================================================

fasta_index = defaultdict(
    lambda: defaultdict(list)
)

source_dirs = {}

for coassembly in COASSEMBLIES:

    source_dirs[
        ("metabat2", coassembly)
    ] = (
        BINROOT
        / coassembly
        / "metabat2"
    )

    source_dirs[
        ("concoct", coassembly)
    ] = (
        BINROOT
        / coassembly
        / "concoct"
        / "fasta_bins"
    )

    source_dirs[
        ("maxbin2", coassembly)
    ] = (
        BINROOT
        / coassembly
        / "maxbin2"
    )


for key, directory in source_dirs.items():

    if not directory.exists():
        die(
            f"No existe directorio de bins: {directory}"
        )

    for p in directory.rglob("*"):

        if (
            p.is_file()
            and is_fasta(p)
        ):
            stem = normalize_fasta_name(
                p.name
            )

            fasta_index[key][stem].append(
                p
            )


# ============================================================
# CARGAR CHECKM2
# ============================================================

quality_rows = []


# ----------------------------
# MetaBAT2
# ----------------------------

metabat_table = (
    ROOT
    / "09_checkm2"
    / "checkm2_all_bins.tsv"
)

for r in read_tsv(metabat_table):

    quality_rows.append({
        "binner": "metabat2",
        "coassembly": r["coassembly"],
        "bin": r["bin"],
        "completeness": float(
            r["completeness"]
        ),
        "contamination": float(
            r["contamination"]
        ),
        "quality_source": str(
            metabat_table
        ),
    })


# ----------------------------
# CONCOCT + MaxBin2
# ----------------------------

for binner in [
    "concoct",
    "maxbin2"
]:

    for coassembly in COASSEMBLIES:

        if binner == "concoct":

            table = (
                ROOT
                / "11_checkm2_multibinner"
                / "concoct"
                / coassembly
                / (
                    f"{coassembly}_"
                    f"concoct_checkm2_summary.tsv"
                )
            )

        else:

            table = (
                ROOT
                / "11_checkm2_multibinner"
                / "maxbin2"
                / coassembly
                / (
                    f"{coassembly}_"
                    f"maxbin2_checkm2_summary.tsv"
                )
            )

        for r in read_tsv(table):

            quality_rows.append({
                "binner": binner,
                "coassembly": r["coassembly"],
                "bin": r["bin"],
                "completeness": float(
                    r["completeness"]
                ),
                "contamination": float(
                    r["contamination"]
                ),
                "quality_source": str(
                    table
                ),
            })


# ============================================================
# QC DE UNIVERSO
# ============================================================

if len(quality_rows) != EXPECTED_RAW_BINS:

    die(
        f"Se esperaban {EXPECTED_RAW_BINS} "
        f"bins brutos pero se recuperaron "
        f"{len(quality_rows)}."
    )


keys = [
    (
        r["binner"],
        r["coassembly"],
        normalize_fasta_name(
            r["bin"]
        )
    )
    for r in quality_rows
]

if len(keys) != len(set(keys)):

    die(
        "Hay claves duplicadas "
        "(binner, coassembly, bin)."
    )


# ============================================================
# HASH DE LOS 18 MAGs FINALES
# ============================================================

final18_hash_to_names = defaultdict(list)

if not FINAL18_DIR.exists():
    die(
        f"No existe: {FINAL18_DIR}"
    )


for p in sorted(
    FINAL18_DIR.iterdir()
):

    if (
        p.is_file()
        and is_fasta(p)
    ):

        stats = fasta_stats(p)

        final18_hash_to_names[
            stats[
                "sequence_content_sha256"
            ]
        ].append(
            normalize_fasta_name(
                p.name
            )
        )


if sum(
    len(v)
    for v in final18_hash_to_names.values()
) != 18:

    die(
        "No se recuperaron exactamente "
        "18 FASTA finales."
    )


# ============================================================
# INVENTARIO
# ============================================================

inventory = []
missing = []
ambiguous = []


for i, r in enumerate(
    quality_rows,
    start=1
):

    binner = r["binner"]
    coassembly = r["coassembly"]

    stem = normalize_fasta_name(
        r["bin"]
    )

    matches = fasta_index[
        (binner, coassembly)
    ].get(
        stem,
        []
    )

    if len(matches) == 0:

        missing.append({
            **r,
            "bin_stem": stem
        })

        fasta_path = ""
        stats = {}

    elif len(matches) > 1:

        ambiguous.append({
            **r,
            "bin_stem": stem,
            "matches": ";".join(
                str(x)
                for x in matches
            )
        })

        fasta_path = ""
        stats = {}

    else:

        fasta_path = str(
            matches[0]
        )

        stats = fasta_stats(
            matches[0]
        )

    comp = r["completeness"]
    cont = r["contamination"]

    qgroup = quality_group(
        comp
    )

    cgroup = contamination_group(
        cont
    )

    mode = screening_mode(
        comp,
        cont
    )

    seqhash = stats.get(
        "sequence_content_sha256",
        ""
    )

    final_matches = (
        final18_hash_to_names.get(
            seqhash,
            []
        )
        if seqhash
        else []
    )

    inventory.append({

        "raw_bin_uid":
            f"{coassembly}__{binner}__{stem}",

        "coassembly":
            coassembly,

        "producer":
            coassembly[0],

        "binner":
            binner,

        "bin":
            r["bin"],

        "bin_stem":
            stem,

        "completeness":
            comp,

        "contamination":
            cont,

        "quality_group":
            qgroup,

        "contamination_group":
            cgroup,

        "screen_target_lt90":
            int(
                comp < 90
            ),

        "screening_mode":
            mode,

        "fasta_found":
            int(
                len(matches) == 1
            ),

        "fasta_path":
            fasta_path,

        "n_contigs":
            stats.get(
                "n_contigs",
                ""
            ),

        "total_bp":
            stats.get(
                "total_bp",
                ""
            ),

        "n50_bp":
            stats.get(
                "n50_bp",
                ""
            ),

        "max_contig_bp":
            stats.get(
                "max_contig_bp",
                ""
            ),

        "gc_pct":
            (
                f"{stats['gc_pct']:.6f}"
                if "gc_pct" in stats
                else ""
            ),

        "sequence_content_sha256":
            seqhash,

        "exact_final18_match":
            ";".join(
                final_matches
            ),

        "quality_source":
            r["quality_source"],
    })


# ============================================================
# DUPLICADOS EXACTOS DE CONTENIDO
# ============================================================

hash_groups = defaultdict(list)

for r in inventory:

    h = r[
        "sequence_content_sha256"
    ]

    if h:
        hash_groups[h].append(
            r
        )


duplicate_hashes = {
    h: rows
    for h, rows in hash_groups.items()
    if len(rows) > 1
}


dup_id_by_hash = {}

for n, h in enumerate(
    sorted(
        duplicate_hashes
    ),
    start=1
):

    dup_id_by_hash[h] = (
        f"EXACTDUP{n:04d}"
    )


for r in inventory:

    h = r[
        "sequence_content_sha256"
    ]

    if h in dup_id_by_hash:

        r[
            "exact_duplicate_group"
        ] = dup_id_by_hash[h]

        r[
            "exact_duplicate_n"
        ] = len(
            duplicate_hashes[h]
        )

    else:

        r[
            "exact_duplicate_group"
        ] = ""

        r[
            "exact_duplicate_n"
        ] = 1


duplicate_rows = []

for h, rows in sorted(
    duplicate_hashes.items()
):

    group_id = dup_id_by_hash[h]

    for r in rows:

        duplicate_rows.append({
            "exact_duplicate_group":
                group_id,
            "n_members":
                len(rows),
            "sequence_content_sha256":
                h,
            "raw_bin_uid":
                r["raw_bin_uid"],
            "coassembly":
                r["coassembly"],
            "binner":
                r["binner"],
            "bin":
                r["bin"],
            "completeness":
                r["completeness"],
            "contamination":
                r["contamination"],
            "fasta_path":
                r["fasta_path"],
        })


# ============================================================
# ELEGIR REPRESENTANTE ÚNICO PARA SCREENING
# sólo elimina duplicados EXACTOS de secuencia
# ============================================================

targets = [
    r
    for r in inventory
    if r[
        "screen_target_lt90"
    ] == 1
]


targets_by_hash = defaultdict(list)

for r in targets:

    targets_by_hash[
        r[
            "sequence_content_sha256"
        ]
    ].append(
        r
    )


for h, rows in targets_by_hash.items():

    chosen = sorted(
        rows,
        key=lambda x: (
            -float(
                x["completeness"]
            ),
            float(
                x["contamination"]
            ),
            x["raw_bin_uid"]
        )
    )[0]

    for r in rows:

        r[
            "unique_sequence_screen_representative"
        ] = int(
            r["raw_bin_uid"]
            ==
            chosen[
                "raw_bin_uid"
            ]
        )


for r in inventory:

    if (
        "unique_sequence_screen_representative"
        not in r
    ):

        r[
            "unique_sequence_screen_representative"
        ] = 0


# ============================================================
# ESCRITURA
# ============================================================

inventory_fields = [
    "raw_bin_uid",
    "coassembly",
    "producer",
    "binner",
    "bin",
    "bin_stem",
    "completeness",
    "contamination",
    "quality_group",
    "contamination_group",
    "screen_target_lt90",
    "screening_mode",
    "fasta_found",
    "fasta_path",
    "n_contigs",
    "total_bp",
    "n50_bp",
    "max_contig_bp",
    "gc_pct",
    "sequence_content_sha256",
    "exact_duplicate_group",
    "exact_duplicate_n",
    "unique_sequence_screen_representative",
    "exact_final18_match",
    "quality_source",
]


write_tsv(
    OUT / "raw_bin_inventory_1006.tsv",
    inventory,
    inventory_fields
)


medium = [
    r
    for r in inventory
    if (
        50
        <= float(
            r["completeness"]
        )
        < 90
    )
]


low = [
    r
    for r in inventory
    if float(
        r["completeness"]
    ) < 50
]


write_tsv(
    OUT / "target_bins_50_89.tsv",
    medium,
    inventory_fields
)


write_tsv(
    OUT / "target_bins_lt50.tsv",
    low,
    inventory_fields
)


all_lt90 = medium + low

write_tsv(
    OUT / "target_bins_lt90_all.tsv",
    all_lt90,
    inventory_fields
)


unique_targets = [
    r
    for r in all_lt90
    if r[
        "unique_sequence_screen_representative"
    ] == 1
]


write_tsv(
    OUT / "target_unique_sequence_bins_lt90.tsv",
    unique_targets,
    inventory_fields
)


write_tsv(
    OUT / "missing_fasta.tsv",
    missing,
    [
        "binner",
        "coassembly",
        "bin",
        "bin_stem",
        "completeness",
        "contamination",
        "quality_source",
    ]
)


write_tsv(
    OUT / "ambiguous_fasta_mapping.tsv",
    ambiguous,
    [
        "binner",
        "coassembly",
        "bin",
        "bin_stem",
        "completeness",
        "contamination",
        "matches",
        "quality_source",
    ]
)


write_tsv(
    OUT / "exact_sequence_duplicate_groups.tsv",
    duplicate_rows,
    [
        "exact_duplicate_group",
        "n_members",
        "sequence_content_sha256",
        "raw_bin_uid",
        "coassembly",
        "binner",
        "bin",
        "completeness",
        "contamination",
        "fasta_path",
    ]
)


# ============================================================
# RESÚMENES
# ============================================================

summary_counter = Counter()

for r in inventory:

    summary_counter[
        (
            r["quality_group"],
            r["contamination_group"]
        )
    ] += 1


summary_rows = []

for (
    qgroup,
    cgroup
), n in sorted(
    summary_counter.items()
):

    summary_rows.append({
        "quality_group":
            qgroup,
        "contamination_group":
            cgroup,
        "n_bins":
            n,
    })


write_tsv(
    OUT / "quality_stratum_summary.tsv",
    summary_rows,
    [
        "quality_group",
        "contamination_group",
        "n_bins",
    ]
)


bc_counter = Counter()

for r in inventory:

    bc_counter[
        (
            r["binner"],
            r["coassembly"],
            r["quality_group"]
        )
    ] += 1


bc_rows = []

for (
    binner,
    coassembly,
    qgroup
), n in sorted(
    bc_counter.items()
):

    bc_rows.append({
        "binner":
            binner,
        "coassembly":
            coassembly,
        "quality_group":
            qgroup,
        "n_bins":
            n,
    })


write_tsv(
    OUT
    / "quality_summary_by_binner_coassembly.tsv",
    bc_rows,
    [
        "binner",
        "coassembly",
        "quality_group",
        "n_bins",
    ]
)


# ============================================================
# README
# ============================================================

readme = f"""
PASO 84 - INVENTARIO DE BINS INCOMPLETOS
========================================

Universo auditado:
{len(inventory)} bins brutos evaluados por CheckM2.

Fuentes:
- MetaBAT2: 09_checkm2/checkm2_all_bins.tsv
- CONCOCT: 11_checkm2_multibinner/concoct
- MaxBin2: 11_checkm2_multibinner/maxbin2

Clasificación:
A_ge90       : completeness >= 90
B_50_89.99   : 50 <= completeness < 90
C_lt50       : completeness < 50

Interpretación:
- Bins >=90 se conservan como referencia, pero no son un nuevo
  objetivo del screening.
- 50-89.99%: una detección positiva de genes puede interpretarse
  en contexto del bin, modulada por contaminación.
- <50%: las ausencias no son interpretables; las señales positivas
  se consideran principalmente evidencia gene/contig-centric.
- Contaminación >10% reduce la confianza de atribución taxonómica.

El archivo target_unique_sequence_bins_lt90.tsv elimina únicamente
duplicados EXACTOS de contenido de secuencia. No colapsa bins
parecidos ni bins del mismo organismo.

No se ha realizado todavía ninguna nueva anotación funcional.
"""

with (
    OUT
    / "README_step84.txt"
).open(
    "w",
    encoding="utf-8"
) as fh:

    fh.write(
        readme.strip()
        + "\n"
    )


# ============================================================
# QC FINAL
# ============================================================

n_hq = sum(
    float(
        r["completeness"]
    ) >= 90
    for r in inventory
)

n_medium = len(
    medium
)

n_low = len(
    low
)

n_lt90 = len(
    all_lt90
)

n_unique_lt90 = len(
    unique_targets
)

n_final_exact = sum(
    bool(
        r["exact_final18_match"]
    )
    for r in inventory
)


print(
    "=" * 64
)

print(
    "PASO 84 COMPLETADO"
)

print(
    "=" * 64
)

print(
    f"Bins brutos auditados:                 {len(inventory)}"
)

print(
    f"Completeness >=90%:                    {n_hq}"
)

print(
    f"Completeness 50-89.99%:                {n_medium}"
)

print(
    f"Completeness <50%:                     {n_low}"
)

print(
    f"Total objetivos <90%:                  {n_lt90}"
)

print(
    f"Objetivos <90% únicos por secuencia:   {n_unique_lt90}"
)

print(
    f"Grupos duplicados exactos:             {len(duplicate_hashes)}"
)

print(
    f"Bins con match exacto a final18:       {n_final_exact}"
)

print(
    f"FASTA faltantes:                       {len(missing)}"
)

print(
    f"Mapeos FASTA ambiguos:                 {len(ambiguous)}"
)

print()

print(
    f"Salida: {OUT}"
)

if missing or ambiguous:

    print(
        "PASO 84 TERMINÓ CON PROBLEMAS DE MAPEADO DE FASTA."
    )

    print(
        "REVISAR missing_fasta.tsv / ambiguous_fasta_mapping.tsv."
    )

    sys.exit(2)

print(
    "PASO 84 FINALIZÓ CORRECTAMENTE; ES SEGURO SALIR."
)
