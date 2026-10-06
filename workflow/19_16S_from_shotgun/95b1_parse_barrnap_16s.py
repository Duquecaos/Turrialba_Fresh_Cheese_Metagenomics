#!/usr/bin/env python3

import csv
import re
import sys
from pathlib import Path


source_type = sys.argv[1]
source_id = sys.argv[2]
coassembly = sys.argv[3]
fasta = Path(sys.argv[4])
gff = Path(sys.argv[5])
outdir = Path(sys.argv[6])

outdir.mkdir(
    parents=True,
    exist_ok=True
)


def revcomp(seq):
    table = str.maketrans(
        "ACGTRYMKBDHVNacgtrymkbdhvn",
        "TGCAYRKMVHDBNtgcayrkmvhdbn"
    )
    return seq.translate(table)[::-1]


# ==============================================================
# Parse Barrnap GFF
# ==============================================================

features_16s = []
total_features = 0

with gff.open(errors="replace") as fh:

    for line in fh:

        if not line.strip() or line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 9:
            continue

        (
            seqid,
            source,
            feature_type,
            start,
            end,
            score,
            strand,
            phase,
            attributes,
        ) = fields

        total_features += 1

        text = attributes.lower()

        is_16s = (
            "16s_rrna" in text
            or "16s ribosomal rna" in text
        )

        if not is_16s:
            continue

        s = int(start)
        e = int(end)

        if s > e:
            s, e = e, s

        partial_attr = (
            "YES"
            if "partial" in text
            else "NO"
        )

        features_16s.append({
            "seqid": seqid,
            "source": source,
            "feature_type": feature_type,
            "start": s,
            "end": e,
            "score": score,
            "strand": strand,
            "attributes": attributes,
            "barrnap_partial_attribute_present":
                partial_attr,
        })


# ==============================================================
# Load only contigs containing a 16S prediction
# ==============================================================

wanted = {
    x["seqid"]
    for x in features_16s
}

sequences = {}

current_id = None
chunks = []

with fasta.open() as fh:

    for line in fh:

        if line.startswith(">"):

            if (
                current_id is not None
                and current_id in wanted
            ):
                sequences[current_id] = "".join(chunks)

            current_id = (
                line[1:]
                .strip()
                .split()[0]
            )

            chunks = []

        else:
            chunks.append(
                line.strip()
            )

    if (
        current_id is not None
        and current_id in wanted
    ):
        sequences[current_id] = "".join(chunks)


missing_contigs = sorted(
    wanted - set(sequences)
)

if missing_contigs:
    raise RuntimeError(
        "Contigs 16S ausentes del FASTA: "
        + ",".join(missing_contigs)
    )


# ==============================================================
# Extract sequences and annotate descriptive completeness
# ==============================================================

records = []

for i, hit in enumerate(
    features_16s,
    start=1
):

    seq = sequences[
        hit["seqid"]
    ]

    contig_len = len(seq)

    s = hit["start"]
    e = hit["end"]

    if not (
        1 <= s <= e <= contig_len
    ):
        raise RuntimeError(
            f"Coordenadas inválidas: "
            f"{hit['seqid']}:{s}-{e} "
            f"len={contig_len}"
        )

    subseq = seq[
        s - 1:e
    ]

    if hit["strand"] == "-":
        subseq = revcomp(
            subseq
        )

    length_bp = len(subseq)

    touches_left = (
        s == 1
    )

    touches_right = (
        e == contig_len
    )

    touches_edge = (
        touches_left
        or touches_right
    )

    if length_bp >= 1300:
        length_class = (
            "near_full_length_candidate_ge1300"
        )

    elif length_bp >= 800:
        length_class = (
            "intermediate_fragment_800_1299"
        )

    else:
        length_class = (
            "short_fragment_lt800"
        )

    locus_id = (
        f"SSU16S_{source_type}_"
        f"{source_id}_{i:03d}"
    )

    records.append({
        "locus_id": locus_id,
        "source_type": source_type,
        "source_id": source_id,
        "coassembly": coassembly,
        "contig": hit["seqid"],
        "contig_length_bp": contig_len,
        "start": s,
        "end": e,
        "strand": hit["strand"],
        "length_bp": length_bp,
        "barrnap_score": hit["score"],
        "touches_left_edge":
            "YES" if touches_left else "NO",
        "touches_right_edge":
            "YES" if touches_right else "NO",
        "edge_truncation_possible":
            "YES" if touches_edge else "NO",
        "barrnap_partial_attribute_present":
            hit[
                "barrnap_partial_attribute_present"
            ],
        "length_class": length_class,
        "attributes": hit["attributes"],
        "sequence": subseq,
    })


# ==============================================================
# Write table
# ==============================================================

table = (
    outdir
    / f"{source_type}__{source_id}.16S.tsv"
)

fields = [
    "locus_id",
    "source_type",
    "source_id",
    "coassembly",
    "contig",
    "contig_length_bp",
    "start",
    "end",
    "strand",
    "length_bp",
    "barrnap_score",
    "touches_left_edge",
    "touches_right_edge",
    "edge_truncation_possible",
    "barrnap_partial_attribute_present",
    "length_class",
    "attributes",
]

with table.open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        extrasaction="ignore",
    )

    writer.writeheader()
    writer.writerows(records)


# ==============================================================
# FASTA
# ==============================================================

fasta_out = (
    outdir
    / f"{source_type}__{source_id}.16S.fasta"
)

with fasta_out.open("w") as fh:

    for r in records:

        fh.write(
            ">"
            + r["locus_id"]
            + " "
            + f"contig={r['contig']} "
            + f"coords={r['start']}-{r['end']} "
            + f"strand={r['strand']} "
            + f"length={r['length_bp']} "
            + f"source={source_type}:{source_id}"
            + "\n"
        )

        seq = r["sequence"]

        for j in range(
            0,
            len(seq),
            80
        ):
            fh.write(
                seq[j:j+80]
                + "\n"
            )


# ==============================================================
# Summary
# ==============================================================

summary = (
    outdir
    / f"{source_type}__{source_id}.summary.tsv"
)

with summary.open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    fh.write(
        f"source_type\t{source_type}\n"
    )

    fh.write(
        f"source_id\t{source_id}\n"
    )

    fh.write(
        f"coassembly\t{coassembly}\n"
    )

    fh.write(
        f"barrnap_total_features\t"
        f"{total_features}\n"
    )

    fh.write(
        f"barrnap_16S_features\t"
        f"{len(records)}\n"
    )

    fh.write(
        "16S_near_full_ge1300\t"
        f"{sum(r['length_bp'] >= 1300 for r in records)}\n"
    )

    fh.write(
        "16S_800_1299\t"
        f"{sum(800 <= r['length_bp'] < 1300 for r in records)}\n"
    )

    fh.write(
        "16S_lt800\t"
        f"{sum(r['length_bp'] < 800 for r in records)}\n"
    )

    fh.write(
        "16S_edge_touching\t"
        f"{sum(r['edge_truncation_possible']=='YES' for r in records)}\n"
    )

    fh.write(
        "16S_internal\t"
        f"{sum(r['edge_truncation_possible']=='NO' for r in records)}\n"
    )

    fh.write(
        "16S_with_barrnap_partial_attribute\t"
        f"{sum(r['barrnap_partial_attribute_present']=='YES' for r in records)}\n"
    )


print(
    f"SOURCE={source_type}:{source_id}"
)

print(
    f"BARRNAP_FEATURES={total_features}"
)

print(
    f"SSU16S_FEATURES={len(records)}"
)

print(
    "95B1_PARSE=PASS"
)
