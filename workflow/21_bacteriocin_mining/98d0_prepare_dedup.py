#!/usr/bin/env python3

import csv
import hashlib
import sys
from pathlib import Path
from collections import defaultdict


if len(sys.argv) != 10:
    raise SystemExit(
        "Usage: 98d0_prepare_dedup.py "
        "<protein_synthesis.tsv> "
        "<contig_synthesis.tsv> "
        "<98C2_context.tsv> "
        "<positive_proteins.faa> "
        "<master_candidates.faa> "
        "<unique_contigs.fna> "
        "<ATTRLOC_loci.tsv> "
        "<antismash_rep18_results_dir> "
        "<OUT>"
    )


PROTEIN_TSV = Path(sys.argv[1])
CONTIG_TSV = Path(sys.argv[2])
C2_TSV = Path(sys.argv[3])
POSITIVE_FAA = Path(sys.argv[4])
MASTER_FAA = Path(sys.argv[5])
UNIQUE_FNA = Path(sys.argv[6])
ATTRLOC_TSV = Path(sys.argv[7])
ASROOT = Path(sys.argv[8])
OUT = Path(sys.argv[9])

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# Helpers
# ============================================================

def read_tsv(path):

    with path.open(
        newline=""
    ) as fh:

        return list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


def write_tsv(
    path,
    rows,
    fields=None
):

    with path.open(
        "w",
        newline=""
    ) as fh:

        if fields is None:

            if rows:
                fields = list(
                    rows[0].keys()
                )

            else:
                fields = [
                    "status"
                ]


        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore"
        )

        writer.writeheader()

        for row in rows:
            writer.writerow(
                row
            )


def read_fasta(path):

    records = {}

    current = None
    seq = []

    def save():

        if current is not None:

            records[
                current
            ] = "".join(
                seq
            ).upper()

    with path.open(
        errors="replace"
    ) as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                save()

                current = (
                    line[1:]
                    .split()[0]
                )

                seq = []

            else:

                seq.append(
                    line
                )

    save()

    return records


def revcomp(seq):

    table = str.maketrans(
        "ACGTRYMKBDHVNacgtrymkbdhvn",
        "TGCAYRKMVHDBNtgcayrkmvhdbn"
    )

    return seq.translate(
        table
    )[::-1]


def parse_genbank_records(path):

    records = {}

    aliases = set()
    seq = []
    in_origin = False


    def save():

        nonlocal aliases, seq

        if not aliases:
            return

        sequence = "".join(
            seq
        ).upper()

        if sequence:

            for alias in aliases:
                records[
                    alias
                ] = sequence


    with path.open(
        errors="replace"
    ) as fh:

        for line in fh:

            if line.startswith("LOCUS"):

                save()

                aliases = set()
                seq = []
                in_origin = False

                fields = line.split()

                if len(fields) >= 2:
                    aliases.add(
                        fields[1]
                    )


            elif line.startswith("ACCESSION"):

                fields = line.split()[1:]

                aliases.update(
                    fields
                )


            elif line.startswith("VERSION"):

                fields = line.split()

                if len(fields) >= 2:

                    aliases.add(
                        fields[1]
                    )

                    aliases.add(
                        fields[1]
                        .split(
                            ".",
                            1
                        )[0]
                    )


            elif line.startswith("ORIGIN"):

                in_origin = True


            elif line.startswith("//"):

                save()

                aliases = set()
                seq = []
                in_origin = False


            elif in_origin:

                bases = "".join(
                    ch
                    for ch in line
                    if ch.upper()
                    in {
                        "A",
                        "C",
                        "G",
                        "T",
                        "N",
                    }
                )

                seq.append(
                    bases
                )


    save()

    return records


# ============================================================
# Inputs
# ============================================================

protein_rows = read_tsv(
    PROTEIN_TSV
)

contig_rows = read_tsv(
    CONTIG_TSV
)

c2_rows = read_tsv(
    C2_TSV
)

attrloc_rows = read_tsv(
    ATTRLOC_TSV
)

positive_proteins = read_fasta(
    POSITIVE_FAA
)

master_proteins = read_fasta(
    MASTER_FAA
)

unique_contigs = read_fasta(
    UNIQUE_FNA
)


if len(contig_rows) != 40:

    raise RuntimeError(
        f"Expected 40 primary contigs; "
        f"observed {len(contig_rows)}"
    )


if len(attrloc_rows) != 7:

    raise RuntimeError(
        f"Expected 7 ATTRLOC loci; "
        f"observed {len(attrloc_rows)}"
    )


# ============================================================
# ATTRLOC map for known protein hypotheses
# ============================================================

ATTRLOC_MAP = {
    "L2_petauri_k141_179388_lactococcin_like":
        "ATTRLOC001",

    "M2_lactis_k141_64662_orf00015":
        "ATTRLOC002;ATTRLOC006",

    "L3_laudensis_k141_46847_orf00027":
        "ATTRLOC003",

    "L3_laudensis_k141_46847_orf00030":
        "ATTRLOC003",

    "L3_laudensis_k141_84984_sORF2":
        "ATTRLOC004",

    "M2_lactis_k141_25580_orf00005":
        "ATTRLOC005",

    "M2_lactis_k141_64662_orf00009":
        "ATTRLOC006",

    "M2_lactis_k141_78158_orf00001":
        "ATTRLOC007",
}


# ============================================================
# 1. Classify the 40 contig contexts
# ============================================================

context_class = {}

context_class_rows = []


for row in contig_rows:

    cid = row[
        "contig_id"
    ]

    linked = row.get(
        "linked_existing_ATTRLOC",
        ""
    )

    specialized = int(
        row.get(
            "specialized_context_support",
            "0"
        )
        or 0
    )


    if linked:

        klass = (
            "recurrence_linked_context"
        )

    elif specialized == 1:

        klass = (
            "additional_specialized_context"
        )

    else:

        klass = (
            "additional_profile_only_context"
        )


    context_class[
        cid
    ] = klass


    context_class_rows.append({
        **row,

        "98D_context_class":
            klass,
    })


write_tsv(
    OUT
    / "98D0_context_classification.tsv",
    context_class_rows
)


# ============================================================
# 2. Protein -> context metadata
# ============================================================

additional_protein_rows = []

colocated_nonhomology_rows = []


for row in protein_rows:

    pid = row[
        "protein_id"
    ]

    cid = row[
        "contig_id"
    ]

    klass = context_class[
        cid
    ]


    if klass.startswith(
        "additional_"
    ):

        additional_protein_rows.append({
            **row,

            "98D_context_class":
                klass,
        })


    elif (
        klass
        ==
        "recurrence_linked_context"
        and
        not row.get(
            "best_query_candidate",
            ""
        )
    ):

        colocated_nonhomology_rows.append({
            **row,

            "98D_context_class":
                klass,
        })


write_tsv(
    OUT
    / "98D0_additional_candidate_proteins.tsv",
    additional_protein_rows
)


write_tsv(
    OUT
    / "98D0_nonhomology_proteins_colocated_with_recurrence.tsv",
    colocated_nonhomology_rows
)


additional_ids = {
    row[
        "protein_id"
    ]
    for row in additional_protein_rows
}


# ============================================================
# 3. Extract additional proteins
# ============================================================

missing_protein = (
    additional_ids
    -
    set(
        positive_proteins
    )
)


if missing_protein:

    raise RuntimeError(
        "Missing additional proteins: "
        + ",".join(
            sorted(
                missing_protein
            )
        )
    )


with (
    OUT
    / "98D0_additional_candidate_proteins.faa"
).open("w") as fh:

    for pid in sorted(
        additional_ids
    ):

        seq = positive_proteins[
            pid
        ].replace(
            "*",
            ""
        )

        fh.write(
            f">{pid}\n{seq}\n"
        )


# ============================================================
# 4. Coordinates for CDS extraction
# ============================================================

c2_by_protein = {
    row[
        "protein_id"
    ]:
        row
    for row in c2_rows
}


cds_metadata = []


with (
    OUT
    / "98D0_additional_candidate_CDS.fna"
).open("w") as out_fh:

    for i, pid in enumerate(
        sorted(
            additional_ids
        ),
        start=1
    ):

        row = c2_by_protein.get(
            pid
        )

        if row is None:

            raise RuntimeError(
                f"No 98C2 coordinates for {pid}"
            )


        cid = row[
            "contig_id"
        ]


        if cid not in unique_contigs:

            raise RuntimeError(
                f"Missing contig sequence: {cid}"
            )


        start = int(
            row[
                "gene_start"
            ]
        )

        end = int(
            row[
                "gene_end"
            ]
        )


        strand = str(
            row[
                "gene_strand"
            ]
        )


        left = min(
            start,
            end
        )

        right = max(
            start,
            end
        )


        seq = unique_contigs[
            cid
        ][
            left - 1:right
        ]


        if strand in {
            "-1",
            "-",
        }:

            seq = revcomp(
                seq
            )


        safe_id = (
            f"ADDCDS{i:03d}"
        )


        aa_len = len(
            positive_proteins[
                pid
            ].replace(
                "*",
                ""
            )
        )


        expected_no_stop = (
            aa_len * 3
        )

        expected_with_stop = (
            aa_len * 3
            + 3
        )


        if len(seq) == expected_no_stop:

            length_relation = (
                "aa_x3"
            )

        elif len(seq) == expected_with_stop:

            length_relation = (
                "aa_x3_plus_stop"
            )

        else:

            length_relation = (
                "other_or_partial"
            )


        out_fh.write(
            f">{safe_id}\n{seq}\n"
        )


        cds_metadata.append({
            "panel_id":
                safe_id,

            "protein_id":
                pid,

            "contig_id":
                cid,

            "context_class":
                context_class[
                    cid
                ],

            "gene_start":
                start,

            "gene_end":
                end,

            "strand":
                strand,

            "aa_length":
                aa_len,

            "cds_length_nt":
                len(
                    seq
                ),

            "length_relation":
                length_relation,
        })


write_tsv(
    OUT
    / "98D0_additional_CDS_metadata.tsv",
    cds_metadata
)


# ============================================================
# 5. Protein redundancy panel:
#       8 established master sequences
#       + additional proteins
# ============================================================

protein_panel_meta = []


with (
    OUT
    / "98D0_protein_redundancy_panel.faa"
).open("w") as fh:

    for i, (
        original_id,
        seq
    ) in enumerate(
        sorted(
            master_proteins.items()
        ),
        start=1
    ):

        safe = (
            f"KNOWNPROT{i:03d}"
        )

        root = (
            original_id.split(
                "|",
                1
            )[0]
        )


        fh.write(
            f">{safe}\n"
            f"{seq.replace('*', '')}\n"
        )


        protein_panel_meta.append({
            "panel_id":
                safe,

            "source_type":
                "existing_master_candidate",

            "source_id":
                original_id,

            "source_root":
                root,

            "linked_ATTRLOC":
                ATTRLOC_MAP.get(
                    root,
                    ""
                ),

            "context_id":
                "",
        })


    offset = len(
        protein_panel_meta
    )


    for j, row in enumerate(
        sorted(
            additional_protein_rows,
            key=lambda x:
                x[
                    "protein_id"
                ]
        ),
        start=1
    ):

        pid = row[
            "protein_id"
        ]

        safe = (
            f"ADDPROT{j:03d}"
        )

        seq = (
            positive_proteins[
                pid
            ]
            .replace(
                "*",
                ""
            )
        )


        fh.write(
            f">{safe}\n{seq}\n"
        )


        protein_panel_meta.append({
            "panel_id":
                safe,

            "source_type":
                "additional_candidate_protein",

            "source_id":
                pid,

            "source_root":
                pid,

            "linked_ATTRLOC":
                "",

            "context_id":
                row[
                    "contig_id"
                ],
        })


write_tsv(
    OUT
    / "98D0_protein_redundancy_panel_metadata.tsv",
    protein_panel_meta
)


# ============================================================
# 6. Extract the seven established ATTRLOC contigs
#    from antiSMASH whole-MAG GenBank outputs
# ============================================================

known_context_sequences = {}

known_context_meta = []


for row in attrloc_rows:

    locus = row[
        "locus_id"
    ]

    mag = row[
        "MAG"
    ]

    contig = row[
        "contig"
    ]


    gbk = (
        ASROOT
        / mag
        / f"{mag}.gbk"
    )


    if not gbk.is_file():

        raise RuntimeError(
            f"antiSMASH combined GBK absent: {gbk}"
        )


    records = parse_genbank_records(
        gbk
    )


    if contig not in records:

        raise RuntimeError(
            f"Contig {contig} not found in {gbk}"
        )


    seq = records[
        contig
    ]


    panel_id = (
        f"KNOWNCTX__{locus}"
    )


    known_context_sequences[
        panel_id
    ] = seq


    known_context_meta.append({
        "panel_id":
            panel_id,

        "source_type":
            "existing_ATTRLOC_context",

        "source_id":
            locus,

        "MAG":
            mag,

        "original_contig":
            contig,

        "context_class":
            "existing_ATTRLOC",

        "length_bp":
            len(
                seq
            ),

        "specialized_context_support":
            "NA",

        "producer":
            "",
    })


# ============================================================
# 7. Additional context panel
# ============================================================

additional_context_rows = [
    row
    for row in context_class_rows
    if row[
        "98D_context_class"
    ].startswith(
        "additional_"
    )
]


if len(
    additional_context_rows
) != 28:

    raise RuntimeError(
        f"Expected 28 additional contexts; "
        f"observed {len(additional_context_rows)}"
    )


context_panel_meta = list(
    known_context_meta
)


with (
    OUT
    / "98D0_context_redundancy_panel.fna"
).open("w") as fh:

    for panel_id, seq in sorted(
        known_context_sequences.items()
    ):

        fh.write(
            f">{panel_id}\n"
        )

        for i in range(
            0,
            len(seq),
            80
        ):

            fh.write(
                seq[
                    i:i+80
                ]
                + "\n"
            )


    for row in sorted(
        additional_context_rows,
        key=lambda x:
            x[
                "contig_id"
            ]
    ):

        cid = row[
            "contig_id"
        ]


        if cid not in unique_contigs:

            raise RuntimeError(
                f"Additional contig absent from FASTA: {cid}"
            )


        seq = unique_contigs[
            cid
        ]


        panel_id = (
            f"ADDCTX__{cid}"
        )


        fh.write(
            f">{panel_id}\n"
        )

        for i in range(
            0,
            len(seq),
            80
        ):

            fh.write(
                seq[
                    i:i+80
                ]
                + "\n"
            )


        context_panel_meta.append({
            "panel_id":
                panel_id,

            "source_type":
                "additional_context",

            "source_id":
                cid,

            "MAG":
                "",

            "original_contig":
                cid,

            "context_class":
                row[
                    "98D_context_class"
                ],

            "length_bp":
                len(
                    seq
                ),

            "specialized_context_support":
                row[
                    "specialized_context_support"
                ],

            "producer":
                row[
                    "producer"
                ],
        })


write_tsv(
    OUT
    / "98D0_context_redundancy_panel_metadata.tsv",
    context_panel_meta
)


# ============================================================
# 8. Preparation summary
# ============================================================

class_counts = defaultdict(
    int
)


for row in context_class_rows:

    class_counts[
        row[
            "98D_context_class"
        ]
    ] += 1


with (
    OUT
    / "98D0_preparation_summary.tsv"
).open(
    "w",
    newline=""
) as fh:

    w = csv.writer(
        fh,
        delimiter="\t",
        lineterminator="\n"
    )

    w.writerow([
        "metric",
        "value"
    ])

    metrics = [
        (
            "primary_contexts",
            len(
                contig_rows
            )
        ),
        (
            "recurrence_linked_contexts",
            class_counts[
                "recurrence_linked_context"
            ]
        ),
        (
            "additional_specialized_contexts",
            class_counts[
                "additional_specialized_context"
            ]
        ),
        (
            "additional_profile_only_contexts",
            class_counts[
                "additional_profile_only_context"
            ]
        ),
        (
            "additional_contexts_total",
            len(
                additional_context_rows
            )
        ),
        (
            "additional_candidate_proteins",
            len(
                additional_protein_rows
            )
        ),
        (
            "nonhomology_proteins_colocated_with_recurrence",
            len(
                colocated_nonhomology_rows
            )
        ),
        (
            "known_master_protein_sequences",
            len(
                master_proteins
            )
        ),
        (
            "known_ATTRLOC_contexts_extracted",
            len(
                known_context_sequences
            )
        ),
        (
            "protein_panel_sequences",
            len(
                protein_panel_meta
            )
        ),
        (
            "context_panel_sequences",
            len(
                context_panel_meta
            )
        ),
        (
            "reads_remapped",
            "NO"
        ),
        (
            "new_bacteriocin_prediction",
            "NO"
        ),
    ]


    w.writerows(
        metrics
    )


print(
    f"ADDITIONAL_CONTEXTS={len(additional_context_rows)}"
)

print(
    f"ADDITIONAL_PROTEINS={len(additional_protein_rows)}"
)

print(
    f"KNOWN_ATTRLOC_CONTEXTS={len(known_context_sequences)}"
)

print(
    "98D0_PREPARE=PASS"
)
