#!/usr/bin/env python3

import csv
import hashlib
import sys
from pathlib import Path
from collections import defaultdict, Counter


if len(sys.argv) != 2:
    raise SystemExit(
        "Usage: 98d0_parse_dedup.py <OUT>"
    )


OUT = Path(sys.argv[1])


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

            fields = (
                list(
                    rows[0].keys()
                )
                if rows
                else ["status"]
            )


        w = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore"
        )

        w.writeheader()

        for row in rows:
            w.writerow(
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

    return seq.translate(
        str.maketrans(
            "ACGTRYMKBDHVN",
            "TGCAYRKMVHDBN"
        )
    )[::-1]


class DSU:

    def __init__(self, members):

        self.parent = {
            x: x
            for x in members
        }


    def find(self, x):

        while self.parent[x] != x:

            self.parent[x] = (
                self.parent[
                    self.parent[x]
                ]
            )

            x = self.parent[x]

        return x


    def union(self, a, b):

        ra = self.find(a)
        rb = self.find(b)

        if ra == rb:
            return

        if ra < rb:
            self.parent[rb] = ra
        else:
            self.parent[ra] = rb


def make_groups(
    members,
    exact_keys,
    near_pairs
):

    dsu = DSU(
        members
    )

    by_exact = defaultdict(
        list
    )


    for member, key in exact_keys.items():

        by_exact[
            key
        ].append(
            member
        )


    for values in by_exact.values():

        first = values[0]

        for x in values[1:]:

            dsu.union(
                first,
                x
            )


    for a, b in near_pairs:

        dsu.union(
            a,
            b
        )


    groups = defaultdict(
        list
    )


    for m in members:

        groups[
            dsu.find(
                m
            )
        ].append(
            m
        )


    ordered = sorted(
        (
            sorted(
                values
            )
            for values in groups.values()
        ),
        key=lambda x:
            x[0]
    )


    return ordered


# ============================================================
# Inputs
# ============================================================

protein_meta = read_tsv(
    OUT
    / "98D0_protein_redundancy_panel_metadata.tsv"
)

context_meta = read_tsv(
    OUT
    / "98D0_context_redundancy_panel_metadata.tsv"
)

cds_meta = read_tsv(
    OUT
    / "98D0_additional_CDS_metadata.tsv"
)

protein_fasta = read_fasta(
    OUT
    / "98D0_protein_redundancy_panel.faa"
)

context_fasta = read_fasta(
    OUT
    / "98D0_context_redundancy_panel.fna"
)

cds_fasta = read_fasta(
    OUT
    / "98D0_additional_candidate_CDS.fna"
)


protein_meta_by_id = {
    r[
        "panel_id"
    ]:
        r
    for r in protein_meta
}

context_meta_by_id = {
    r[
        "panel_id"
    ]:
        r
    for r in context_meta
}

cds_meta_by_id = {
    r[
        "panel_id"
    ]:
        r
    for r in cds_meta
}


# ============================================================
# 1. Protein exact hashes
# ============================================================

protein_exact_keys = {
    pid:
        hashlib.sha256(
            seq.encode()
        ).hexdigest()
    for pid, seq in protein_fasta.items()
}


# ============================================================
# 2. Parse protein BLAST
# ============================================================

protein_pairs_best = {}


with (
    OUT
    / "98D0_protein_all_vs_all.tsv"
).open(
    errors="replace"
) as fh:

    for line in fh:

        if not line.strip():
            continue

        f = line.rstrip(
            "\n"
        ).split(
            "\t"
        )


        if len(f) < 12:
            continue


        (
            q,
            s,
            pident,
            aln_len,
            qlen,
            slen,
            qstart,
            qend,
            sstart,
            send,
            evalue,
            bitscore,
        ) = f[:12]


        if q == s:
            continue


        if (
            q not in protein_fasta
            or
            s not in protein_fasta
        ):
            continue


        pident = float(
            pident
        )

        qlen = int(
            qlen
        )

        slen = int(
            slen
        )

        qspan = (
            abs(
                int(qend)
                -
                int(qstart)
            )
            + 1
        )

        sspan = (
            abs(
                int(send)
                -
                int(sstart)
            )
            + 1
        )

        qcov = min(
            100.0,
            100.0
            * qspan
            / qlen
        )

        scov = min(
            100.0,
            100.0
            * sspan
            / slen
        )


        pair = tuple(
            sorted(
                [
                    q,
                    s,
                ]
            )
        )


        record = {
            "protein_a":
                pair[0],

            "protein_b":
                pair[1],

            "pident":
                pident,

            "qcov_pct":
                qcov,

            "scov_pct":
                scov,

            "evalue":
                float(
                    evalue
                ),

            "bitscore":
                float(
                    bitscore
                ),
        }


        old = protein_pairs_best.get(
            pair
        )


        if (
            old is None
            or
            record[
                "bitscore"
            ]
            >
            old[
                "bitscore"
            ]
        ):

            protein_pairs_best[
                pair
            ] = record


protein_pair_rows = []

protein_near_edges = []


for pair, row in sorted(
    protein_pairs_best.items()
):

    a = pair[0]
    b = pair[1]


    exact = (
        protein_exact_keys[
            a
        ]
        ==
        protein_exact_keys[
            b
        ]
    )


    if exact:

        klass = (
            "exact_sequence_duplicate"
        )

        protein_near_edges.append(
            pair
        )


    elif (
        row[
            "pident"
        ] >= 95.0
        and
        row[
            "qcov_pct"
        ] >= 90.0
        and
        row[
            "scov_pct"
        ] >= 90.0
    ):

        klass = (
            "near_redundant_sequence"
        )

        protein_near_edges.append(
            pair
        )


    elif (
        row[
            "pident"
        ] >= 80.0
        and
        row[
            "qcov_pct"
        ] >= 80.0
        and
        row[
            "scov_pct"
        ] >= 80.0
    ):

        klass = (
            "related_sequence_not_collapsed"
        )


    else:

        continue


    protein_pair_rows.append({
        **row,

        "source_type_a":
            protein_meta_by_id[
                a
            ][
                "source_type"
            ],

        "source_id_a":
            protein_meta_by_id[
                a
            ][
                "source_id"
            ],

        "source_type_b":
            protein_meta_by_id[
                b
            ][
                "source_type"
            ],

        "source_id_b":
            protein_meta_by_id[
                b
            ][
                "source_id"
            ],

        "redundancy_class":
            klass,
    })


write_tsv(
    OUT
    / "98D0_protein_similarity_pairs.tsv",
    protein_pair_rows
)


# ============================================================
# 3. Protein redundancy groups
# ============================================================

protein_groups = make_groups(
    list(
        protein_fasta
    ),
    protein_exact_keys,
    protein_near_edges
)


protein_group_rows = []

protein_group_for_member = {}


for i, members in enumerate(
    protein_groups,
    start=1
):

    gid = (
        f"PROTGRP{i:03d}"
    )


    known = [
        x
        for x in members
        if protein_meta_by_id[
            x
        ][
            "source_type"
        ]
        ==
        "existing_master_candidate"
    ]


    additional = [
        x
        for x in members
        if protein_meta_by_id[
            x
        ][
            "source_type"
        ]
        ==
        "additional_candidate_protein"
    ]


    if (
        known
        and
        additional
    ):

        interpretation = (
            "additional_sequence_near_existing_candidate"
        )

    elif len(
        additional
    ) > 1:

        interpretation = (
            "additional_sequence_redundancy_group"
        )

    elif additional:

        interpretation = (
            "additional_sequence_unique_at_operational_threshold"
        )

    else:

        interpretation = (
            "existing_candidate_only"
        )


    for member in members:

        protein_group_for_member[
            member
        ] = gid


    protein_group_rows.append({
        "protein_group":
            gid,

        "member_count":
            len(
                members
            ),

        "members":
            ";".join(
                members
            ),

        "existing_master_members":
            ";".join(
                protein_meta_by_id[
                    x
                ][
                    "source_id"
                ]
                for x in known
            ),

        "additional_members":
            ";".join(
                protein_meta_by_id[
                    x
                ][
                    "source_id"
                ]
                for x in additional
            ),

        "interpretation":
            interpretation,
    })


write_tsv(
    OUT
    / "98D0_protein_redundancy_groups.tsv",
    protein_group_rows
)


# ============================================================
# 4. CDS exact clusters
# ============================================================

cds_exact_keys = {
    cid:
        hashlib.sha256(
            seq.encode()
        ).hexdigest()
    for cid, seq in cds_fasta.items()
}


by_cds_hash = defaultdict(
    list
)


for cid, key in cds_exact_keys.items():

    by_cds_hash[
        key
    ].append(
        cid
    )


cds_cluster_rows = []


for i, members in enumerate(
    sorted(
        (
            sorted(
                x
            )
            for x in by_cds_hash.values()
        ),
        key=lambda x:
            x[0]
    ),
    start=1
):

    gid = (
        f"CDSCLUST{i:03d}"
    )


    cds_cluster_rows.append({
        "CDS_cluster":
            gid,

        "member_count":
            len(
                members
            ),

        "panel_ids":
            ";".join(
                members
            ),

        "protein_ids":
            ";".join(
                cds_meta_by_id[
                    x
                ][
                    "protein_id"
                ]
                for x in members
            ),

        "contigs":
            ";".join(
                cds_meta_by_id[
                    x
                ][
                    "contig_id"
                ]
                for x in members
            ),

        "exact_nucleotide_sequence":
            "YES",
    })


write_tsv(
    OUT
    / "98D0_CDS_exact_clusters.tsv",
    cds_cluster_rows
)


# ============================================================
# 5. Context exact hashes, orientation-insensitive
# ============================================================

context_exact_keys = {}


for cid, seq in context_fasta.items():

    canonical = min(
        seq,
        revcomp(
            seq
        )
    )

    context_exact_keys[
        cid
    ] = hashlib.sha256(
        canonical.encode()
    ).hexdigest()


# ============================================================
# 6. Context BLAST
# ============================================================

context_pairs_best = {}


with (
    OUT
    / "98D0_context_all_vs_all.tsv"
).open(
    errors="replace"
) as fh:

    for line in fh:

        if not line.strip():
            continue


        f = line.rstrip(
            "\n"
        ).split(
            "\t"
        )


        if len(f) < 12:
            continue


        (
            q,
            s,
            pident,
            aln_len,
            qlen,
            slen,
            qstart,
            qend,
            sstart,
            send,
            evalue,
            bitscore,
        ) = f[:12]


        if q == s:
            continue


        if (
            q not in context_fasta
            or
            s not in context_fasta
        ):
            continue


        pident = float(
            pident
        )

        aln_len = int(
            aln_len
        )

        qlen = int(
            qlen
        )

        slen = int(
            slen
        )


        shorter_cov = min(
            100.0,
            100.0
            * aln_len
            / min(
                qlen,
                slen
            )
        )


        pair = tuple(
            sorted(
                [
                    q,
                    s,
                ]
            )
        )


        record = {
            "context_a":
                pair[0],

            "context_b":
                pair[1],

            "pident":
                pident,

            "aligned_length":
                aln_len,

            "length_a":
                len(
                    context_fasta[
                        pair[0]
                    ]
                ),

            "length_b":
                len(
                    context_fasta[
                        pair[1]
                    ]
                ),

            "aligned_fraction_shorter_pct":
                shorter_cov,

            "evalue":
                float(
                    evalue
                ),

            "bitscore":
                float(
                    bitscore
                ),
        }


        old = context_pairs_best.get(
            pair
        )


        if (
            old is None
            or
            record[
                "bitscore"
            ]
            >
            old[
                "bitscore"
            ]
        ):

            context_pairs_best[
                pair
            ] = record


context_pair_rows = []

context_near_edges = []


for pair, row in sorted(
    context_pairs_best.items()
):

    a = pair[0]
    b = pair[1]


    exact = (
        context_exact_keys[
            a
        ]
        ==
        context_exact_keys[
            b
        ]
    )


    if exact:

        klass = (
            "exact_context_duplicate"
        )

        context_near_edges.append(
            pair
        )


    elif (
        row[
            "pident"
        ] >= 95.0
        and
        row[
            "aligned_fraction_shorter_pct"
        ] >= 85.0
    ):

        klass = (
            "near_redundant_context"
        )

        context_near_edges.append(
            pair
        )


    elif (
        row[
            "pident"
        ] >= 90.0
        and
        row[
            "aligned_fraction_shorter_pct"
        ] >= 50.0
    ):

        klass = (
            "related_context_not_collapsed"
        )


    else:

        continue


    context_pair_rows.append({
        **row,

        "source_type_a":
            context_meta_by_id[
                a
            ][
                "source_type"
            ],

        "source_id_a":
            context_meta_by_id[
                a
            ][
                "source_id"
            ],

        "source_type_b":
            context_meta_by_id[
                b
            ][
                "source_type"
            ],

        "source_id_b":
            context_meta_by_id[
                b
            ][
                "source_id"
            ],

        "redundancy_class":
            klass,
    })


write_tsv(
    OUT
    / "98D0_context_similarity_pairs.tsv",
    context_pair_rows
)


# ============================================================
# 7. Context groups
# ============================================================

context_groups = make_groups(
    list(
        context_fasta
    ),
    context_exact_keys,
    context_near_edges
)


context_group_rows = []

context_group_for_member = {}

context_group_rep = {}


for i, members in enumerate(
    context_groups,
    start=1
):

    gid = (
        f"CTXGRP{i:03d}"
    )


    known = [
        x
        for x in members
        if context_meta_by_id[
            x
        ][
            "source_type"
        ]
        ==
        "existing_ATTRLOC_context"
    ]


    additional = [
        x
        for x in members
        if context_meta_by_id[
            x
        ][
            "source_type"
        ]
        ==
        "additional_context"
    ]


    if known:

        representative = sorted(
            known
        )[0]


        if additional:

            interpretation = (
                "additional_context_near_existing_ATTRLOC"
            )

        else:

            interpretation = (
                "existing_ATTRLOC_only"
            )


    elif additional:

        # Prefer specialized contexts over profile-only,
        # then longest context, then stable lexical ID.
        representative = sorted(
            additional,
            key=lambda x: (
                0
                if context_meta_by_id[
                    x
                ][
                    "context_class"
                ]
                ==
                "additional_specialized_context"
                else 1,

                -int(
                    context_meta_by_id[
                        x
                    ][
                        "length_bp"
                    ]
                ),

                x,
            )
        )[0]


        if len(
            additional
        ) > 1:

            interpretation = (
                "additional_context_redundancy_group"
            )

        else:

            interpretation = (
                "additional_context_unique_at_operational_threshold"
            )


    else:

        representative = members[0]

        interpretation = (
            "unresolved"
        )


    for member in members:

        context_group_for_member[
            member
        ] = gid

        context_group_rep[
            member
        ] = representative


    context_group_rows.append({
        "context_group":
            gid,

        "member_count":
            len(
                members
            ),

        "members":
            ";".join(
                members
            ),

        "existing_ATTRLOC_members":
            ";".join(
                context_meta_by_id[
                    x
                ][
                    "source_id"
                ]
                for x in known
            ),

        "additional_context_members":
            ";".join(
                context_meta_by_id[
                    x
                ][
                    "source_id"
                ]
                for x in additional
            ),

        "representative_panel_id":
            representative,

        "representative_source_id":
            context_meta_by_id[
                representative
            ][
                "source_id"
            ],

        "interpretation":
            interpretation,
    })


write_tsv(
    OUT
    / "98D0_context_redundancy_groups.tsv",
    context_group_rows
)


# ============================================================
# 8. Mapping manifest for additional contexts
# ============================================================

mapping_rows = []


for panel_id, meta in sorted(
    context_meta_by_id.items()
):

    if meta[
        "source_type"
    ] != "additional_context":
        continue


    gid = context_group_for_member[
        panel_id
    ]

    representative = context_group_rep[
        panel_id
    ]

    representative_meta = (
        context_meta_by_id[
            representative
        ]
    )


    group_row = next(
        x
        for x in context_group_rows
        if x[
            "context_group"
        ] == gid
    )


    group_has_existing = bool(
        group_row[
            "existing_ATTRLOC_members"
        ]
    )


    if group_has_existing:

        mapping_role = (
            "covered_by_existing_ATTRLOC_context_group"
        )


    elif panel_id == representative:

        if (
            meta[
                "context_class"
            ]
            ==
            "additional_specialized_context"
        ):

            mapping_role = (
                "primary_additional_mapping_representative"
            )

        else:

            mapping_role = (
                "exploratory_additional_mapping_representative"
            )


    else:

        mapping_role = (
            "redundancy_group_member_not_primary_representative"
        )


    mapping_rows.append({
        "contig_id":
            meta[
                "source_id"
            ],

        "panel_id":
            panel_id,

        "context_class":
            meta[
                "context_class"
            ],

        "producer":
            meta[
                "producer"
            ],

        "length_bp":
            meta[
                "length_bp"
            ],

        "context_group":
            gid,

        "group_interpretation":
            group_row[
                "interpretation"
            ],

        "group_existing_ATTRLOC_members":
            group_row[
                "existing_ATTRLOC_members"
            ],

        "group_representative_source_id":
            representative_meta[
                "source_id"
            ],

        "mapping_role":
            mapping_role,
    })


write_tsv(
    OUT
    / "98D0_additional_context_mapping_manifest.tsv",
    mapping_rows
)


# ============================================================
# 9. Additional protein manifest
# ============================================================

protein_manifest_rows = []


for panel_id, meta in sorted(
    protein_meta_by_id.items()
):

    if (
        meta[
            "source_type"
        ]
        !=
        "additional_candidate_protein"
    ):
        continue


    gid = protein_group_for_member[
        panel_id
    ]


    group = next(
        x
        for x in protein_group_rows
        if x[
            "protein_group"
        ] == gid
    )


    protein_manifest_rows.append({
        "protein_id":
            meta[
                "source_id"
            ],

        "contig_id":
            meta[
                "context_id"
            ],

        "protein_group":
            gid,

        "group_interpretation":
            group[
                "interpretation"
            ],

        "near_existing_master_candidates":
            group[
                "existing_master_members"
            ],
    })


write_tsv(
    OUT
    / "98D0_additional_protein_manifest.tsv",
    protein_manifest_rows
)


# ============================================================
# 10. Global summary
# ============================================================

mapping_roles = Counter(
    row[
        "mapping_role"
    ]
    for row in mapping_rows
)


additional_context_near_existing = sum(
    row[
        "group_interpretation"
    ]
    ==
    "additional_context_near_existing_ATTRLOC"
    for row in mapping_rows
)


additional_protein_near_existing = sum(
    row[
        "group_interpretation"
    ]
    ==
    "additional_sequence_near_existing_candidate"
    for row in protein_manifest_rows
)


multi_additional_context_groups = sum(
    (
        row[
            "interpretation"
        ]
        ==
        "additional_context_redundancy_group"
    )
    for row in context_group_rows
)


multi_cds_clusters = sum(
    int(
        row[
            "member_count"
        ]
    ) > 1
    for row in cds_cluster_rows
)


with (
    OUT
    / "98D0_global_summary.tsv"
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
            "known_master_proteins",
            sum(
                r[
                    "source_type"
                ]
                ==
                "existing_master_candidate"
                for r in protein_meta
            )
        ),
        (
            "additional_candidate_proteins",
            len(
                protein_manifest_rows
            )
        ),
        (
            "additional_protein_groups_operational",
            len({
                r[
                    "protein_group"
                ]
                for r in protein_manifest_rows
            })
        ),
        (
            "additional_proteins_near_existing_master_candidate",
            additional_protein_near_existing
        ),
        (
            "additional_CDS_sequences",
            len(
                cds_meta
            )
        ),
        (
            "additional_CDS_exact_clusters",
            len(
                cds_cluster_rows
            )
        ),
        (
            "additional_CDS_exact_multimember_clusters",
            multi_cds_clusters
        ),
        (
            "known_ATTRLOC_contexts",
            sum(
                r[
                    "source_type"
                ]
                ==
                "existing_ATTRLOC_context"
                for r in context_meta
            )
        ),
        (
            "additional_contexts",
            len(
                mapping_rows
            )
        ),
        (
            "additional_contexts_near_existing_ATTRLOC_context",
            additional_context_near_existing
        ),
        (
            "additional_only_multimember_context_groups",
            multi_additional_context_groups
        ),
        (
            "primary_additional_mapping_representatives",
            mapping_roles[
                "primary_additional_mapping_representative"
            ]
        ),
        (
            "exploratory_additional_mapping_representatives",
            mapping_roles[
                "exploratory_additional_mapping_representative"
            ]
        ),
        (
            "additional_contexts_covered_by_existing_ATTRLOC_group",
            mapping_roles[
                "covered_by_existing_ATTRLOC_context_group"
            ]
        ),
        (
            "additional_context_redundant_group_members",
            mapping_roles[
                "redundancy_group_member_not_primary_representative"
            ]
        ),
        (
            "protein_near_redundancy_threshold",
            "PID>=95%;coverage_both>=90%"
        ),
        (
            "context_near_redundancy_threshold",
            "PID>=95%;aligned_fraction_shorter>=85%"
        ),
        (
            "context_homology_group_equivalent_to_same_biological_locus",
            "NO"
        ),
        (
            "additional_context_equivalent_to_novel_bacteriocin",
            "NO"
        ),
        (
            "reads_remapped",
            "NO"
        ),
        (
            "98D0_status",
            "COMPLETE"
        ),
        (
            "next_step",
            "98D1_optional_competitive_read_mapping_of_deduplicated_context_groups"
        ),
    ]


    w.writerows(
        metrics
    )


print(
    f"ADDITIONAL_PROTEINS={len(protein_manifest_rows)}"
)

print(
    f"ADDITIONAL_CONTEXTS={len(mapping_rows)}"
)

print(
    f"PRIMARY_MAPPING_REPS="
    f"{mapping_roles['primary_additional_mapping_representative']}"
)

print(
    "98D0_PARSE=PASS"
)
