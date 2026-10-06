#!/usr/bin/env python3

import csv
import hashlib
import re
import sys
from pathlib import Path
from collections import defaultdict


if len(sys.argv) != 6:
    raise SystemExit(
        "Usage: 98d1ar_build_locus_windows.py "
        "<D0_mapping_manifest.tsv> "
        "<98C2_positive_context.tsv> "
        "<unique_contigs.fna> "
        "<established_ATTRLOC_reference.fna> "
        "<OUT>"
    )


MANIFEST = Path(sys.argv[1])
C2 = Path(sys.argv[2])
CONTIG_FASTA = Path(sys.argv[3])
KNOWN_FASTA = Path(sys.argv[4])
OUT = Path(sys.argv[5])

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# Helpers
# ============================================================

def read_tsv(path):
    with path.open(newline="") as fh:
        return list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


def read_fasta_full_headers(path):

    records = []

    header = None
    seq = []

    def save():
        if header is not None:
            records.append(
                (
                    header,
                    "".join(seq).upper()
                )
            )

    with path.open(errors="replace") as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                save()

                header = line[1:]
                seq = []

            else:

                seq.append(line)

    save()

    return records


def read_fasta_ids(path):

    return {
        header.split()[0]: seq
        for header, seq
        in read_fasta_full_headers(path)
    }


def revcomp(seq):

    return seq.translate(
        str.maketrans(
            "ACGTRYMKBDHVN",
            "TGCAYRKMVHDBN"
        )
    )[::-1]


def canonical_hash(seq):

    canonical = min(
        seq,
        revcomp(seq)
    )

    return hashlib.sha256(
        canonical.encode()
    ).hexdigest()


# ============================================================
# Inputs
# ============================================================

manifest = read_tsv(
    MANIFEST
)

c2 = read_tsv(
    C2
)

contigs = read_fasta_ids(
    CONTIG_FASTA
)


# ============================================================
# 1. Recover established ATTRLOC reference
# ============================================================

known = {}


for header, seq in read_fasta_full_headers(
    KNOWN_FASTA
):

    match = re.search(
        r"ATTRLOC\d{3}",
        header
    )

    if not match:
        continue

    locus = match.group(0)

    if locus in known:

        raise RuntimeError(
            f"Duplicate established reference for {locus}"
        )

    known[locus] = seq


expected = {
    f"ATTRLOC{i:03d}"
    for i in range(1, 8)
}


missing = (
    expected
    -
    set(known)
)


extra = (
    set(known)
    -
    expected
)


if missing or extra:

    raise RuntimeError(
        f"Established ATTRLOC reference mismatch. "
        f"missing={sorted(missing)}, extra={sorted(extra)}"
    )


# ============================================================
# 2. Select 23 deduplicated additional representatives
# ============================================================

allowed_roles = {
    "primary_additional_mapping_representative",
    "exploratory_additional_mapping_representative",
}


selected = [
    row
    for row in manifest
    if row[
        "mapping_role"
    ] in allowed_roles
]


if len(selected) != 23:

    raise RuntimeError(
        f"Expected 23 additional representatives; "
        f"found {len(selected)}"
    )


# ============================================================
# 3. Primary candidate coordinates from 98C2
# ============================================================

primary_by_contig = defaultdict(
    list
)


for row in c2:

    if str(
        row.get(
            "primary_for_specialized_followup",
            ""
        )
    ) != "1":
        continue

    primary_by_contig[
        row[
            "contig_id"
        ]
    ].append(
        row
    )


# ============================================================
# 4. Build records
# ============================================================

records = []


# Established contexts: use ORIGINAL mapping reference unchanged.
for locus in sorted(known):

    seq = known[
        locus
    ]

    records.append({
        "reference_id":
            locus,

        "source_id":
            locus,

        "reference_class":
            "existing_ATTRLOC",

        "mapping_role":
            "existing_ATTRLOC_context",

        "context_group":
            locus,

        "producer":
            "",

        "full_contig_length":
            len(seq),

        "window_start":
            1,

        "window_end":
            len(seq),

        "reference_length":
            len(seq),

        "candidate_protein_ids":
            "",

        "candidate_gene_span_start":
            "",

        "candidate_gene_span_end":
            "",

        "candidate_evidence_classes":
            "",

        "sequence":
            seq,
    })


FLANK = 5000


for row in sorted(
    selected,
    key=lambda x:
        x[
            "contig_id"
        ]
):

    cid = row[
        "contig_id"
    ]


    if cid not in contigs:

        raise RuntimeError(
            f"Missing contig sequence: {cid}"
        )


    candidates = primary_by_contig.get(
        cid,
        []
    )


    if not candidates:

        raise RuntimeError(
            f"No primary candidate coordinates for {cid}"
        )


    starts = []
    ends = []


    for candidate in candidates:

        a = int(
            candidate[
                "gene_start"
            ]
        )

        b = int(
            candidate[
                "gene_end"
            ]
        )

        starts.append(
            min(a, b)
        )

        ends.append(
            max(a, b)
        )


    gene_left = min(
        starts
    )

    gene_right = max(
        ends
    )


    full_seq = contigs[
        cid
    ]

    full_len = len(
        full_seq
    )


    window_start = max(
        1,
        gene_left - FLANK
    )

    window_end = min(
        full_len,
        gene_right + FLANK
    )


    seq = full_seq[
        window_start - 1:
        window_end
    ]


    if (
        row[
            "mapping_role"
        ]
        ==
        "primary_additional_mapping_representative"
    ):

        ref_class = (
            "additional_specialized"
        )

    else:

        ref_class = (
            "additional_exploratory"
        )


    records.append({
        "reference_id":
            f"ADDLocus_{cid}",

        "source_id":
            cid,

        "reference_class":
            ref_class,

        "mapping_role":
            row[
                "mapping_role"
            ],

        "context_group":
            row[
                "context_group"
            ],

        "producer":
            row[
                "producer"
            ],

        "full_contig_length":
            full_len,

        "window_start":
            window_start,

        "window_end":
            window_end,

        "reference_length":
            len(seq),

        "candidate_protein_ids":
            ";".join(
                c[
                    "protein_id"
                ]
                for c in candidates
            ),

        "candidate_gene_span_start":
            gene_left,

        "candidate_gene_span_end":
            gene_right,

        "candidate_evidence_classes":
            ";".join(
                sorted({
                    c[
                        "direct_evidence_class"
                    ]
                    for c in candidates
                })
            ),

        "sequence":
            seq,
    })


# ============================================================
# 5. Validate 30 references
# ============================================================

if len(records) != 30:

    raise RuntimeError(
        f"Expected 30 references; found {len(records)}"
    )


ids = [
    r[
        "reference_id"
    ]
    for r in records
]


if len(ids) != len(set(ids)):

    raise RuntimeError(
        "Duplicate reference IDs"
    )


# Exact orientation-insensitive duplicates should not remain.
hash_groups = defaultdict(
    list
)


for row in records:

    hash_groups[
        canonical_hash(
            row[
                "sequence"
            ]
        )
    ].append(
        row[
            "reference_id"
        ]
    )


duplicates = [
    members
    for members in hash_groups.values()
    if len(members) > 1
]


if duplicates:

    raise RuntimeError(
        "Exact reference duplicates remain: "
        + repr(duplicates)
    )


# ============================================================
# 6. FASTA
# ============================================================

REF = (
    OUT
    / "98D1AR_competitive_locus_windows_30.fna"
)


with REF.open("w") as fh:

    for row in records:

        fh.write(
            f">{row['reference_id']}\n"
        )

        seq = row[
            "sequence"
        ]

        for i in range(
            0,
            len(seq),
            80
        ):

            fh.write(
                seq[i:i+80]
                + "\n"
            )


# ============================================================
# 7. Metadata
# ============================================================

fields = [
    "reference_id",
    "source_id",
    "reference_class",
    "mapping_role",
    "context_group",
    "producer",
    "full_contig_length",
    "window_start",
    "window_end",
    "reference_length",
    "candidate_protein_ids",
    "candidate_gene_span_start",
    "candidate_gene_span_end",
    "candidate_evidence_classes",
]


with (
    OUT
    / "98D1AR_competitive_locus_windows_30_metadata.tsv"
).open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()

    for row in records:

        writer.writerow({
            key:
                row[
                    key
                ]
            for key in fields
        })


# ============================================================
# 8. Sample manifest
# ============================================================

samples = [
    "L1_1", "L1_2", "L1_3",
    "L2_1", "L2_2", "L2_3",
    "L3_1", "L3_2", "L3_3",
    "M1_1", "M1_2", "M1_3",
    "M2_1", "M2_2", "M2_3",
    "M3_1", "M3_2", "M3_3",
]


with (
    OUT
    / "98D1AR_sample_manifest.tsv"
).open(
    "w",
    newline=""
) as fh:

    writer = csv.writer(
        fh,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow([
        "array_index",
        "sample",
        "producer_code",
        "biological_unit",
        "time_code",
    ])


    for idx, sample in enumerate(
        samples
    ):

        left, time_code = (
            sample.split("_")
        )

        writer.writerow([
            idx,
            sample,
            sample[0],
            left[1:],
            time_code,
        ])


# ============================================================
# 9. Summary
# ============================================================

additional_lengths = [
    int(
        r[
            "reference_length"
        ]
    )
    for r in records
    if r[
        "reference_class"
    ].startswith(
        "additional_"
    )
]


known_lengths = [
    int(
        r[
            "reference_length"
        ]
    )
    for r in records
    if r[
        "reference_class"
    ]
    ==
    "existing_ATTRLOC"
]


with (
    OUT
    / "98D1AR_reference_summary.tsv"
).open(
    "w",
    newline=""
) as fh:

    writer = csv.writer(
        fh,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow([
        "metric",
        "value"
    ])

    summary = [
        (
            "competitive_reference_sequences",
            30
        ),
        (
            "existing_ATTRLOC_contexts",
            7
        ),
        (
            "additional_specialized_windows",
            20
        ),
        (
            "additional_exploratory_windows",
            3
        ),
        (
            "window_flank_bp",
            5000
        ),
        (
            "known_reference_min_bp",
            min(
                known_lengths
            )
        ),
        (
            "known_reference_max_bp",
            max(
                known_lengths
            )
        ),
        (
            "additional_window_min_bp",
            min(
                additional_lengths
            )
        ),
        (
            "additional_window_max_bp",
            max(
                additional_lengths
            )
        ),
        (
            "exact_duplicate_reference_sequences",
            0
        ),
        (
            "samples",
            18
        ),
        (
            "reference_design",
            "established_ATTRLOC_local_contexts_plus_candidate_centered_additional_windows"
        ),
        (
            "reads_mapped",
            "NO"
        ),
        (
            "amplicon_data_used",
            "NO"
        ),
        (
            "supersedes_reference_from_job_133413_for_mapping",
            "YES"
        ),
        (
            "next_step",
            "98D1B_competitive_mapping_18_samples"
        ),
    ]

    writer.writerows(
        summary
    )


print(
    "ESTABLISHED_ATTRLOC=7"
)

print(
    "ADDITIONAL_SPECIALIZED=20"
)

print(
    "ADDITIONAL_EXPLORATORY=3"
)

print(
    "TOTAL_REFERENCE_SEQUENCES=30"
)

print(
    f"ADDITIONAL_WINDOW_MIN={min(additional_lengths)}"
)

print(
    f"ADDITIONAL_WINDOW_MAX={max(additional_lengths)}"
)

print(
    "98D1AR_PREPARE=PASS"
)
