#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import csv
import gzip
import hashlib
import os
import sys


USER = os.environ["USER"]

ROOT = Path(
    f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"
)

INV = (
    ROOT
    / "68_incomplete_bin_inventory"
    / "target_bins_lt90_all.tsv"
)

FINAL18 = (
    ROOT
    / "22_final_representative_mags"
    / "fastas"
)

OUT = (
    ROOT
    / "69_incomplete_unique_contigs"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# HELPERS
# ============================================================

def die(msg):
    print(
        f"ERROR: {msg}",
        file=sys.stderr
    )
    sys.exit(1)


def open_text(path):

    path = Path(path)

    if path.name.lower().endswith(".gz"):
        return gzip.open(
            path,
            "rt"
        )

    return path.open(
        "r",
        encoding="utf-8"
    )


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


def revcomp(seq):

    table = str.maketrans(
        "ACGTNacgtn",
        "TGCANtgcan"
    )

    return seq.translate(
        table
    )[::-1]


def canonical_sequence(seq):

    seq = (
        seq
        .replace(" ", "")
        .replace("\n", "")
        .upper()
    )

    rc = revcomp(seq)

    if rc < seq:
        return rc

    return seq


def seq_hash(seq):

    return hashlib.sha256(
        seq.encode()
    ).hexdigest()


def read_fasta(path):

    name = None
    seq = []

    with open_text(path) as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if name is not None:

                    yield (
                        name,
                        "".join(seq)
                    )

                name = line[1:].split()[0]
                seq = []

            else:

                seq.append(line)

        if name is not None:

            yield (
                name,
                "".join(seq)
            )


def write_fasta_record(
    fh,
    name,
    seq,
    width=80
):

    fh.write(
        f">{name}\n"
    )

    for i in range(
        0,
        len(seq),
        width
    ):

        fh.write(
            seq[
                i:i+width
            ]
            + "\n"
        )


# ============================================================
# CARGAR INVENTARIO DE 815 BINS
# ============================================================

if not INV.exists():
    die(
        f"No existe: {INV}"
    )


with INV.open(
    "r",
    encoding="utf-8-sig",
    newline=""
) as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    bins = list(
        reader
    )


if len(bins) != 815:

    die(
        f"Se esperaban 815 bins <90%; "
        f"se encontraron {len(bins)}"
    )


required = {
    "raw_bin_uid",
    "coassembly",
    "producer",
    "binner",
    "completeness",
    "contamination",
    "quality_group",
    "contamination_group",
    "fasta_path",
}


missing_cols = (
    required
    - set(
        bins[0].keys()
    )
)

if missing_cols:

    die(
        "Faltan columnas: "
        + ", ".join(
            sorted(
                missing_cols
            )
        )
    )


# ============================================================
# INDEXAR CONTIGS DE LOS 18 MAGs FINALES
# ============================================================

final_hash_members = defaultdict(
    set
)


final_fasta_files = sorted(
    [
        p
        for p in FINAL18.iterdir()
        if (
            p.is_file()
            and is_fasta(p)
        )
    ]
)


if len(
    final_fasta_files
) != 18:

    die(
        f"Se esperaban 18 FASTA finales; "
        f"se encontraron {len(final_fasta_files)}"
    )


for fasta in final_fasta_files:

    mag = fasta.name

    for header, seq in read_fasta(
        fasta
    ):

        canon = canonical_sequence(
            seq
        )

        h = seq_hash(
            canon
        )

        final_hash_members[h].add(
            mag
        )


# ============================================================
# CATÁLOGO ÚNICO
# ============================================================

hash_to_id = {}

meta = {}

bin_contig_rows = []

raw_memberships = 0

catalog_all = (
    OUT
    / "unique_contigs_lt90_all.fna"
)

catalog_new = (
    OUT
    / "unique_contigs_lt90_not_in_final18.fna"
)


with catalog_all.open(
    "w",
    encoding="utf-8"
) as all_fh, catalog_new.open(
    "w",
    encoding="utf-8"
) as new_fh:

    for bidx, b in enumerate(
        bins,
        start=1
    ):

        fasta = Path(
            b["fasta_path"]
        )

        if not fasta.exists():
            die(
                f"FASTA no existe: {fasta}"
            )

        n_bin_contigs = 0

        for original_header, seq in read_fasta(
            fasta
        ):

            raw_memberships += 1
            n_bin_contigs += 1

            canon = canonical_sequence(
                seq
            )

            h = seq_hash(
                canon
            )

            if h not in hash_to_id:

                uid = (
                    f"ICONTIG"
                    f"{len(hash_to_id)+1:09d}"
                )

                hash_to_id[h] = uid

                gc = (
                    canon.count("G")
                    + canon.count("C")
                )

                final_matches = sorted(
                    final_hash_members.get(
                        h,
                        set()
                    )
                )

                meta[h] = {
                    "unique_contig_id":
                        uid,
                    "sequence_sha256":
                        h,
                    "length_bp":
                        len(canon),
                    "gc_pct":
                        (
                            100.0
                            * gc
                            / len(canon)
                            if canon
                            else 0
                        ),
                    "represented_in_final18":
                        int(
                            len(
                                final_matches
                            ) > 0
                        ),
                    "final18_MAGs":
                        set(
                            final_matches
                        ),
                    "raw_bin_uids":
                        set(),
                    "coassemblies":
                        set(),
                    "producers":
                        set(),
                    "binners":
                        set(),
                    "quality_groups":
                        set(),
                    "contamination_groups":
                        set(),
                    "completeness_values":
                        [],
                    "contamination_values":
                        [],
                }

                write_fasta_record(
                    all_fh,
                    uid,
                    canon
                )

                if not final_matches:

                    write_fasta_record(
                        new_fh,
                        uid,
                        canon
                    )

            uid = hash_to_id[h]

            m = meta[h]

            m[
                "raw_bin_uids"
            ].add(
                b[
                    "raw_bin_uid"
                ]
            )

            m[
                "coassemblies"
            ].add(
                b[
                    "coassembly"
                ]
            )

            m[
                "producers"
            ].add(
                b[
                    "producer"
                ]
            )

            m[
                "binners"
            ].add(
                b[
                    "binner"
                ]
            )

            m[
                "quality_groups"
            ].add(
                b[
                    "quality_group"
                ]
            )

            m[
                "contamination_groups"
            ].add(
                b[
                    "contamination_group"
                ]
            )

            m[
                "completeness_values"
            ].append(
                float(
                    b[
                        "completeness"
                    ]
                )
            )

            m[
                "contamination_values"
            ].append(
                float(
                    b[
                        "contamination"
                    ]
                )
            )

            bin_contig_rows.append({
                "raw_bin_uid":
                    b[
                        "raw_bin_uid"
                    ],
                "coassembly":
                    b[
                        "coassembly"
                    ],
                "producer":
                    b[
                        "producer"
                    ],
                "binner":
                    b[
                        "binner"
                    ],
                "bin_completeness":
                    b[
                        "completeness"
                    ],
                "bin_contamination":
                    b[
                        "contamination"
                    ],
                "bin_quality_group":
                    b[
                        "quality_group"
                    ],
                "unique_contig_id":
                    uid,
                "original_contig_header":
                    original_header,
                "sequence_sha256":
                    h,
                "length_bp":
                    len(canon),
                "represented_in_final18":
                    int(
                        h in
                        final_hash_members
                    ),
            })


# ============================================================
# TABLA MAESTRA DE CONTIGS
# ============================================================

contig_rows = []


for h, m in sorted(
    meta.items(),
    key=lambda kv:
        kv[1][
            "unique_contig_id"
        ]
):

    comps = m[
        "completeness_values"
    ]

    conts = m[
        "contamination_values"
    ]

    qgroups = m[
        "quality_groups"
    ]

    if (
        "B_50_89.99"
        in qgroups
    ):
        strongest_context = (
            "present_in_50_89_bin"
        )

    else:
        strongest_context = (
            "lt50_only"
        )

    if (
        m[
            "represented_in_final18"
        ] == 1
    ):
        annotation_priority = (
            "already_represented_in_final18"
        )

    elif strongest_context == (
        "present_in_50_89_bin"
    ):
        annotation_priority = (
            "NEW_medium_completion_context"
        )

    else:
        annotation_priority = (
            "NEW_lt50_gene_centric"
        )

    contig_rows.append({

        "unique_contig_id":
            m[
                "unique_contig_id"
            ],

        "sequence_sha256":
            h,

        "length_bp":
            m[
                "length_bp"
            ],

        "gc_pct":
            f"{m['gc_pct']:.6f}",

        "n_bin_memberships":
            len(
                m[
                    "raw_bin_uids"
                ]
            ),

        "n_coassemblies":
            len(
                m[
                    "coassemblies"
                ]
            ),

        "coassemblies":
            ";".join(
                sorted(
                    m[
                        "coassemblies"
                    ]
                )
            ),

        "producers":
            ";".join(
                sorted(
                    m[
                        "producers"
                    ]
                )
            ),

        "binners":
            ";".join(
                sorted(
                    m[
                        "binners"
                    ]
                )
            ),

        "quality_groups":
            ";".join(
                sorted(
                    qgroups
                )
            ),

        "strongest_incomplete_context":
            strongest_context,

        "min_member_bin_completeness":
            min(
                comps
            ),

        "max_member_bin_completeness":
            max(
                comps
            ),

        "min_member_bin_contamination":
            min(
                conts
            ),

        "max_member_bin_contamination":
            max(
                conts
            ),

        "represented_in_final18":
            m[
                "represented_in_final18"
            ],

        "final18_MAGs":
            ";".join(
                sorted(
                    m[
                        "final18_MAGs"
                    ]
                )
            ),

        "annotation_priority":
            annotation_priority,

        "member_raw_bin_uids":
            ";".join(
                sorted(
                    m[
                        "raw_bin_uids"
                    ]
                )
            ),
    })


# ============================================================
# ESCRITURA TSV
# ============================================================

def write_tsv(
    path,
    rows,
    fields
):

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

        w.writerows(
            rows
        )


write_tsv(
    OUT
    / "unique_contig_catalog_lt90.tsv",
    contig_rows,
    list(
        contig_rows[0].keys()
    )
)


write_tsv(
    OUT
    / "bin_to_unique_contig_membership.tsv",
    bin_contig_rows,
    list(
        bin_contig_rows[0].keys()
    )
)


new_medium = [
    r
    for r in contig_rows
    if r[
        "annotation_priority"
    ] ==
    "NEW_medium_completion_context"
]


new_low = [
    r
    for r in contig_rows
    if r[
        "annotation_priority"
    ] ==
    "NEW_lt50_gene_centric"
]


already_final = [
    r
    for r in contig_rows
    if r[
        "represented_in_final18"
    ] == 1
]


write_tsv(
    OUT
    / "new_contigs_medium_context.tsv",
    new_medium,
    list(
        contig_rows[0].keys()
    )
)


write_tsv(
    OUT
    / "new_contigs_lt50_only.tsv",
    new_low,
    list(
        contig_rows[0].keys()
    )
)


write_tsv(
    OUT
    / "contigs_already_in_final18.tsv",
    already_final,
    list(
        contig_rows[0].keys()
    )
)


# ============================================================
# RESUMEN
# ============================================================

n_unique = len(
    contig_rows
)

n_new = (
    len(
        new_medium
    )
    +
    len(
        new_low
    )
)

total_unique_bp = sum(
    int(
        r[
            "length_bp"
        ]
    )
    for r in contig_rows
)

new_bp = sum(
    int(
        r[
            "length_bp"
        ]
    )
    for r in contig_rows
    if r[
        "represented_in_final18"
    ] == 0
)

multi_bin = sum(
    int(
        r[
            "n_bin_memberships"
        ]
    ) > 1
    for r in contig_rows
)


summary = [
    (
        "input_bins_lt90",
        len(
            bins
        )
    ),
    (
        "raw_bin_contig_memberships",
        raw_memberships
    ),
    (
        "unique_contigs_lt90",
        n_unique
    ),
    (
        "unique_contigs_shared_across_multiple_bins",
        multi_bin
    ),
    (
        "unique_contigs_already_in_final18",
        len(
            already_final
        )
    ),
    (
        "unique_new_contigs_not_in_final18",
        n_new
    ),
    (
        "new_contigs_with_50_89_context",
        len(
            new_medium
        )
    ),
    (
        "new_contigs_lt50_only",
        len(
            new_low
        )
    ),
    (
        "total_unique_bp",
        total_unique_bp
    ),
    (
        "new_bp_not_in_final18",
        new_bp
    ),
]


with (
    OUT
    / "contig_catalog_summary.tsv"
).open(
    "w",
    encoding="utf-8",
    newline=""
) as fh:

    w = csv.writer(
        fh,
        delimiter="\t",
        lineterminator="\n"
    )

    w.writerow(
        [
            "metric",
            "value"
        ]
    )

    w.writerows(
        summary
    )


# ============================================================
# README
# ============================================================

readme = """
PASO 85 - CATALOGO UNICO DE CONTIGS DE BINS <90%
=================================================

Este paso desduplica por SECUENCIA EXACTA de contig, considerando
tambien la reversa complementaria como la misma secuencia.

Mantiene:
- bin de procedencia
- coensamblaje
- productor
- binner
- completitud/contaminacion del bin
- numero de bins que contienen cada contig
- coincidencia exacta con contigs de los 18 MAGs finales

Archivos principales:

unique_contigs_lt90_all.fna
    Todos los contigs unicos presentes en los 815 bins <90%.

unique_contigs_lt90_not_in_final18.fna
    Contigs unicos que NO estan ya representados exactamente
    en los 18 MAGs finales.

unique_contig_catalog_lt90.tsv
    Catalogo y procedencia de cada contig.

new_contigs_medium_context.tsv
    Contigs nuevos presentes al menos en un bin 50-89.99%.

new_contigs_lt50_only.tsv
    Contigs nuevos encontrados solo en bins <50%.

contigs_already_in_final18.tsv
    Contigs de bins incompletos que ya estan representados
    en los MAGs finales.

El siguiente paso debe predecir/anotar genes unicamente sobre
unique_contigs_lt90_not_in_final18.fna.
"""

with (
    OUT
    / "README_step85.txt"
).open(
    "w",
    encoding="utf-8"
) as fh:

    fh.write(
        readme.strip()
        + "\n"
    )


print(
    "=" * 64
)

print(
    "PASO 85 COMPLETADO"
)

print(
    "=" * 64
)

for metric, value in summary:
    print(
        f"{metric:46s} {value}"
    )

print()

print(
    f"Salida: {OUT}"
)

print(
    "PASO 85 FINALIZÓ CORRECTAMENTE; ES SEGURO SALIR."
)
