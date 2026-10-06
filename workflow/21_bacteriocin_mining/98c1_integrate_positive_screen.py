#!/usr/bin/env python3

import csv
import re
import sys
from pathlib import Path
from collections import defaultdict


if len(sys.argv) != 8:
    raise SystemExit(
        "Usage: 98c1_integrate_positive_screen.py "
        "<catalog.faa> <eggnog.tsv> "
        "<ga.tblout> <noga.tblout> "
        "<blast.tsv> <master.faa> <OUT>"
    )


CATALOG = Path(sys.argv[1])
EGGNOG = Path(sys.argv[2])
GA_TBL = Path(sys.argv[3])
NOGA_TBL = Path(sys.argv[4])
BLAST = Path(sys.argv[5])
MASTER = Path(sys.argv[6])
OUT = Path(sys.argv[7])

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# Helpers
# ============================================================

def parse_hmmer_tbl(path, threshold_mode):

    hits = []

    if not path.is_file():
        return hits

    with path.open(
        errors="replace"
    ) as fh:

        for line in fh:

            if (
                not line.strip()
                or line.startswith("#")
            ):
                continue

            f = line.split()

            if len(f) < 18:
                continue

            # hmmscan --tblout:
            # target(profile), accession,
            # query(sequence), accession,
            # full_E, score, bias, ...
            model = f[0]
            protein = f[2]

            try:
                evalue = float(f[4])
                score = float(f[5])
            except ValueError:
                continue

            hits.append({
                "protein_id":
                    protein,

                "HMM":
                    model,

                "hmm_evalue":
                    evalue,

                "hmm_score":
                    score,

                "hmm_threshold_mode":
                    threshold_mode,
            })

    return hits


def parse_fasta_headers_for_targets(path, wanted):

    metadata = {}

    current_id = None
    seq = []

    def save():

        if (
            current_id is None
            or current_id not in wanted
        ):
            return

        header = metadata[
            current_id
        ][
            "raw_header"
        ]

        metadata[
            current_id
        ][
            "aa_length"
        ] = len(
            "".join(seq)
            .replace("*", "")
        )

    with path.open(
        errors="replace"
    ) as fh:

        for line in fh:

            line = line.rstrip("\n")

            if line.startswith(">"):

                save()

                header = line[1:]

                current_id = (
                    header.split()[0]
                )

                seq = []

                if current_id in wanted:

                    contig = re.sub(
                        r"_\d+$",
                        "",
                        current_id
                    )

                    partial = ""
                    start_type = ""
                    strand = ""
                    gene_start = ""
                    gene_end = ""

                    m = re.search(
                        r"#\s*(\d+)\s*#\s*(\d+)\s*#\s*(-?1)\s*#",
                        header
                    )

                    if m:
                        gene_start = m.group(1)
                        gene_end = m.group(2)
                        strand = m.group(3)

                    pm = re.search(
                        r"partial=([01]{2})",
                        header
                    )

                    if pm:
                        partial = pm.group(1)

                    sm = re.search(
                        r"start_type=([^;]+)",
                        header
                    )

                    if sm:
                        start_type = sm.group(1)

                    metadata[
                        current_id
                    ] = {
                        "protein_id":
                            current_id,

                        "contig_id":
                            contig,

                        "gene_start":
                            gene_start,

                        "gene_end":
                            gene_end,

                        "strand":
                            strand,

                        "partial":
                            partial,

                        "start_type":
                            start_type,

                        "raw_header":
                            header,

                        "aa_length":
                            0,
                    }

            else:

                if (
                    current_id is not None
                    and current_id in wanted
                ):
                    seq.append(
                        line.strip()
                    )

    save()

    return metadata


def extract_positive_fasta(
    source,
    wanted,
    destination
):

    write = False

    with source.open(
        errors="replace"
    ) as inp, destination.open(
        "w"
    ) as out:

        for line in inp:

            if line.startswith(">"):

                pid = (
                    line[1:]
                    .split()[0]
                )

                write = (
                    pid in wanted
                )

            if write:
                out.write(line)


# ============================================================
# 1. HMM hits
# ============================================================

ga_hits = parse_hmmer_tbl(
    GA_TBL,
    "GA_curated"
)

noga_hits = parse_hmmer_tbl(
    NOGA_TBL,
    "E1e-5_noGA"
)


hmm_by_protein = defaultdict(list)

for hit in ga_hits + noga_hits:

    hmm_by_protein[
        hit[
            "protein_id"
        ]
    ].append(hit)


# ============================================================
# 2. BLASTP-short
# ============================================================

blast_hits = []

blast_by_protein = defaultdict(list)


with BLAST.open(
    errors="replace"
) as fh:

    for line in fh:

        if not line.strip():
            continue

        f = line.rstrip("\n").split("\t")

        if len(f) < 12:
            continue

        (
            qseqid,
            sseqid,
            pident,
            length,
            qlen,
            slen,
            qstart,
            qend,
            sstart,
            send,
            evalue,
            bitscore,
        ) = f[:12]

        pident = float(pident)
        aln_len = int(length)
        qlen = int(qlen)
        slen = int(slen)
        evalue = float(evalue)
        bitscore = float(bitscore)

        query_span = abs(
            int(qend) - int(qstart)
        ) + 1

        qcov = (
            100.0 * query_span / qlen
            if qlen
            else 0.0
        )

        qcov = min(
            100.0,
            qcov
        )


        if (
            pident == 100.0
            and qcov == 100.0
            and slen == qlen
        ):

            klass = (
                "exact_full_length"
            )

        elif (
            pident == 100.0
            and qcov == 100.0
        ):

            klass = (
                "exact_query_contained"
            )

        elif (
            pident >= 90.0
            and qcov >= 90.0
            and bitscore >= 25.0
        ):

            klass = (
                "near_exact"
            )

        elif (
            pident >= 70.0
            and qcov >= 80.0
            and bitscore >= 25.0
            and evalue <= 1e-2
        ):

            klass = (
                "strong_short_protein_homology"
            )

        else:

            # Keep permissive raw BLAST output on disk,
            # but do not promote weaker matches.
            continue


        row = {
            "query_candidate":
                qseqid,

            "protein_id":
                sseqid,

            "pident":
                pident,

            "alignment_length":
                aln_len,

            "query_length":
                qlen,

            "subject_length":
                slen,

            "query_coverage_pct":
                qcov,

            "evalue":
                evalue,

            "bitscore":
                bitscore,

            "blast_class":
                klass,
        }

        blast_hits.append(row)

        blast_by_protein[
            sseqid
        ].append(row)


# ============================================================
# 3. Union of positive proteins
# ============================================================

positive_ids = set(
    hmm_by_protein
) | set(
    blast_by_protein
)


metadata = parse_fasta_headers_for_targets(
    CATALOG,
    positive_ids
)


# ============================================================
# 4. eggNOG only for positive proteins
# ============================================================

eggnog = {}


with EGGNOG.open(
    errors="replace"
) as fh:

    header = None

    for line in fh:

        line = line.rstrip("\n")

        if not line:
            continue

        if line.startswith("#"):
            continue

        fields = line.split("\t")

        if header is None:

            header = fields
            continue

        if len(fields) != len(header):
            continue

        query = fields[0]

        if query not in positive_ids:
            continue

        row = dict(
            zip(
                header,
                fields
            )
        )

        eggnog[
            query
        ] = row


# ============================================================
# 5. Integrated positive table
# ============================================================

positive_rows = []


for protein in sorted(
    positive_ids
):

    hh = hmm_by_protein.get(
        protein,
        []
    )

    bb = blast_by_protein.get(
        protein,
        []
    )

    ga_models = sorted({
        h["HMM"]
        for h in hh
        if h[
            "hmm_threshold_mode"
        ] == "GA_curated"
    })

    noga_models = sorted({
        h["HMM"]
        for h in hh
        if h[
            "hmm_threshold_mode"
        ] == "E1e-5_noGA"
    })


    blast_classes = sorted({
        b[
            "blast_class"
        ]
        for b in bb
    })


    if bb and ga_models:

        screen_class = (
            "sequence_homology_plus_GA_HMM"
        )

    elif bb and noga_models:

        screen_class = (
            "sequence_homology_plus_supplementary_HMM"
        )

    elif bb:

        screen_class = (
            "sequence_homology_only"
        )

    elif ga_models:

        screen_class = (
            "GA_HMM_profile_only"
        )

    else:

        screen_class = (
            "supplementary_noGA_HMM_only"
        )


    m = metadata.get(
        protein,
        {}
    )

    e = eggnog.get(
        protein,
        {}
    )


    best_blast = None

    if bb:

        best_blast = sorted(
            bb,
            key=lambda x: (
                -x[
                    "query_coverage_pct"
                ],
                -x[
                    "pident"
                ],
                x[
                    "evalue"
                ],
                -x[
                    "bitscore"
                ],
            )
        )[0]


    positive_rows.append({
        "protein_id":
            protein,

        "contig_id":
            m.get(
                "contig_id",
                re.sub(
                    r"_\d+$",
                    "",
                    protein
                )
            ),

        "aa_length":
            m.get(
                "aa_length",
                ""
            ),

        "gene_start":
            m.get(
                "gene_start",
                ""
            ),

        "gene_end":
            m.get(
                "gene_end",
                ""
            ),

        "strand":
            m.get(
                "strand",
                ""
            ),

        "partial":
            m.get(
                "partial",
                ""
            ),

        "start_type":
            m.get(
                "start_type",
                ""
            ),

        "GA_HMM_models":
            ";".join(
                ga_models
            ),

        "supplementary_HMM_models":
            ";".join(
                noga_models
            ),

        "HMM_hit_count":
            len(hh),

        "BLAST_hit_count":
            len(bb),

        "BLAST_classes":
            ";".join(
                blast_classes
            ),

        "best_query_candidate":
            (
                best_blast[
                    "query_candidate"
                ]
                if best_blast
                else ""
            ),

        "best_BLAST_pident":
            (
                f"{best_blast['pident']:.3f}"
                if best_blast
                else ""
            ),

        "best_BLAST_qcov_pct":
            (
                f"{best_blast['query_coverage_pct']:.3f}"
                if best_blast
                else ""
            ),

        "best_BLAST_evalue":
            (
                best_blast[
                    "evalue"
                ]
                if best_blast
                else ""
            ),

        "best_BLAST_bitscore":
            (
                best_blast[
                    "bitscore"
                ]
                if best_blast
                else ""
            ),

        "screen_class":
            screen_class,

        "eggnog_annotated":
            int(
                protein in eggnog
            ),

        "eggnog_COG_category":
            e.get(
                "COG_category",
                ""
            ),

        "eggnog_Description":
            e.get(
                "Description",
                ""
            ),

        "eggnog_Preferred_name":
            e.get(
                "Preferred_name",
                ""
            ),

        "eggnog_PFAMs":
            e.get(
                "PFAMs",
                ""
            ),
    })


# ============================================================
# 6. Output tables
# ============================================================

def write_dict_table(
    path,
    rows
):

    if not rows:

        path.write_text(
            ""
        )
        return

    fields = list(
        rows[0].keys()
    )

    with path.open(
        "w",
        newline=""
    ) as fh:

        w = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n"
        )

        w.writeheader()
        w.writerows(rows)


write_dict_table(
    OUT
    / "98C1_positive_proteins.tsv",
    positive_rows
)

write_dict_table(
    OUT
    / "98C1_promoted_BLAST_hits.tsv",
    blast_hits
)

write_dict_table(
    OUT
    / "98C1_HMM_hits.tsv",
    ga_hits + noga_hits
)


# ============================================================
# 7. Positive protein FASTA
# ============================================================

extract_positive_fasta(
    CATALOG,
    positive_ids,
    OUT
    / "98C1_positive_proteins.faa"
)


# ============================================================
# 8. Summary
# ============================================================

ga_ids = {
    h[
        "protein_id"
    ]
    for h in ga_hits
}

noga_ids = {
    h[
        "protein_id"
    ]
    for h in noga_hits
}

blast_ids = set(
    blast_by_protein
)

both_hmm_blast = (
    set(
        hmm_by_protein
    )
    &
    blast_ids
)


class_counts = defaultdict(
    int
)

for row in positive_rows:

    class_counts[
        row[
            "screen_class"
        ]
    ] += 1


with (
    OUT
    / "98C1_global_summary.tsv"
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
            "GA_HMM_hits",
            len(
                ga_hits
            )
        ),
        (
            "GA_HMM_unique_proteins",
            len(
                ga_ids
            )
        ),
        (
            "noGA_HMM_hits",
            len(
                noga_hits
            )
        ),
        (
            "noGA_HMM_unique_proteins",
            len(
                noga_ids
            )
        ),
        (
            "promoted_BLAST_hits",
            len(
                blast_hits
            )
        ),
        (
            "BLAST_unique_proteins",
            len(
                blast_ids
            )
        ),
        (
            "HMM_and_BLAST_unique_proteins",
            len(
                both_hmm_blast
            )
        ),
        (
            "positive_union_unique_proteins",
            len(
                positive_ids
            )
        ),
        (
            "positive_union_unique_contigs",
            len({
                r[
                    "contig_id"
                ]
                for r in positive_rows
            })
        ),
        (
            "positive_proteins_with_eggNOG",
            sum(
                int(
                    r[
                        "eggnog_annotated"
                    ]
                )
                for r in positive_rows
            )
        ),
    ]


    for key in sorted(
        class_counts
    ):

        metrics.append(
            (
                f"screen_class__{key}",
                class_counts[
                    key
                ]
            )
        )


    metrics += [
        (
            "bacteriocin_function_confirmed",
            "NO"
        ),
        (
            "new_candidate_screen_performed",
            "YES"
        ),
        (
            "reads_remapped",
            "NO"
        ),
        (
            "amplicon_data_used",
            "NO"
        ),
        (
            "next_step",
            "98C2_positive_contig_and_bin_context"
        ),
    ]


    w.writerows(
        metrics
    )


print(
    f"GA_HMM_UNIQUE={len(ga_ids)}"
)

print(
    f"NOGA_HMM_UNIQUE={len(noga_ids)}"
)

print(
    f"BLAST_UNIQUE={len(blast_ids)}"
)

print(
    f"UNION={len(positive_ids)}"
)

print(
    "98C1_INTEGRATION=PASS"
)
