#!/usr/bin/env python3

import csv
import hashlib
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path(sys.argv[1])

BASE = (
    ROOT
    / "97_16S_shotgun"
    / "95B1_barrnap"
)

MANIFEST = (
    BASE
    / "95B1_input_manifest.tsv"
)

AUDIT = (
    ROOT
    / "97_16S_shotgun"
    / "95A_audit"
)

FINAL18_MANIFEST = (
    AUDIT
    / "95A_final18.tsv"
)

OUT = (
    ROOT
    / "97_16S_shotgun"
    / "95B2_consolidated"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# Helpers
# ============================================================

def read_tsv(path):

    if not path.is_file():
        raise RuntimeError(
            f"Archivo faltante: {path}"
        )

    with path.open(
        newline="",
        errors="replace"
    ) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        rows = list(reader)
        fields = reader.fieldnames or []

    return rows, fields


def write_tsv(path, rows, fields):

    with path.open(
        "w",
        newline=""
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path):

    sequences = {}

    if not path.is_file():
        raise RuntimeError(
            f"FASTA faltante: {path}"
        )

    current = None
    chunks = []

    with path.open() as fh:

        for line in fh:

            line = line.rstrip("\n")

            if line.startswith(">"):

                if current is not None:
                    sequences[current] = (
                        "".join(chunks).upper()
                    )

                current = (
                    line[1:]
                    .split()[0]
                )

                chunks = []

            else:
                chunks.append(
                    line.strip()
                )

        if current is not None:
            sequences[current] = (
                "".join(chunks).upper()
            )

    return sequences


def revcomp(seq):

    table = str.maketrans(
        "ACGTRYMKBDHVN",
        "TGCAYRKMVHDBN"
    )

    return seq.translate(
        table
    )[::-1]


def canonical_sequence(seq):

    seq = seq.upper()
    rc = revcomp(seq)

    return min(
        seq,
        rc
    )


def integer(value):

    return int(float(value))


def overlap(a1, a2, b1, b2):

    return (
        a1 <= b2
        and b1 <= a2
    )


def overlap_length(
    a1,
    a2,
    b1,
    b2
):

    if not overlap(
        a1,
        a2,
        b1,
        b2
    ):
        return 0

    return (
        min(a2, b2)
        - max(a1, b1)
        + 1
    )


def unique_join(values):

    return ";".join(
        sorted(
            {
                str(x).strip()
                for x in values
                if str(x).strip()
            }
        )
    )


# ============================================================
# 1. Read manifest and all Barrnap outputs
# ============================================================

manifest, _ = read_tsv(
    MANIFEST
)

if len(manifest) != 24:
    raise RuntimeError(
        f"Manifest esperado=24 observado={len(manifest)}"
    )


coassembly_rows = []
mag_rows = []

coassembly_seq = {}
mag_seq = {}


for entry in manifest:

    source_type = entry["source_type"]
    source_id = entry["source_id"]

    task_id = (
        f"{source_type}__{source_id}"
    )

    taskdir = (
        BASE
        / "results"
        / task_id
    )

    table = (
        taskdir
        / f"{task_id}.16S.tsv"
    )

    fasta = (
        taskdir
        / f"{task_id}.16S.fasta"
    )

    marker = (
        taskdir
        / "95B1_COMPLETE.ok"
    )

    if not marker.is_file():
        raise RuntimeError(
            f"Marker faltante: {marker}"
        )

    rows, _ = read_tsv(
        table
    )

    seqs = read_fasta(
        fasta
    )

    for row in rows:

        locus = row["locus_id"]

        if locus not in seqs:
            raise RuntimeError(
                f"Secuencia faltante para {locus}"
            )

        seq = seqs[locus]

        expected_length = integer(
            row["length_bp"]
        )

        if len(seq) != expected_length:
            raise RuntimeError(
                f"{locus}: TSV length={expected_length}, "
                f"FASTA length={len(seq)}"
            )

        if source_type == "coassembly":

            coassembly_rows.append(
                dict(row)
            )

            coassembly_seq[locus] = seq

        elif source_type == "MAG":

            mag_rows.append(
                dict(row)
            )

            mag_seq[locus] = seq

        else:
            raise RuntimeError(
                f"source_type inesperado: {source_type}"
            )


if len(coassembly_rows) != 158:
    raise RuntimeError(
        "Loci 16S coassembly esperado=158 "
        f"observado={len(coassembly_rows)}"
    )

if len(mag_rows) != 6:
    raise RuntimeError(
        "Loci 16S MAG esperado=6 "
        f"observado={len(mag_rows)}"
    )


# ============================================================
# 2. Exact final18 contig membership
# ============================================================

final18_entries, _ = read_tsv(
    FINAL18_MANIFEST
)

if len(final18_entries) != 18:
    raise RuntimeError(
        f"Final18 esperado=18 observado={len(final18_entries)}"
    )


mag_contig_membership = defaultdict(list)


for entry in final18_entries:

    mag = entry["MAG"]
    coassembly = entry["coassembly"]
    fasta = Path(
        entry["fasta"]
    )

    if not fasta.is_file():
        raise RuntimeError(
            f"FASTA MAG faltante: {fasta}"
        )

    with fasta.open() as fh:

        for line in fh:

            if not line.startswith(">"):
                continue

            contig = (
                line[1:]
                .strip()
                .split()[0]
            )

            mag_contig_membership[
                (
                    coassembly,
                    contig,
                )
            ].append(
                mag
            )


# ============================================================
# 3. Annotate canonical coassembly loci
# ============================================================

canonical_rows = []


for row in coassembly_rows:

    locus = row["locus_id"]

    key = (
        row["coassembly"],
        row["contig"],
    )

    mags = sorted(
        set(
            mag_contig_membership.get(
                key,
                []
            )
        )
    )

    seq = coassembly_seq[
        locus
    ]

    seq_canonical = canonical_sequence(
        seq
    )

    sha = hashlib.sha256(
        seq_canonical.encode()
    ).hexdigest()

    r = dict(row)

    r.update({
        "n_final18_MAGs_exact":
            len(mags),

        "final18_MAGs_exact":
            ";".join(mags),

        "exact_final18_contig_membership":
            "YES"
            if mags
            else "NO",

        "orientation_insensitive_sequence_sha256":
            sha,

        "sequence_length_verified":
            "YES",

        "primary_catalog_role":
            "coassembly_16S_primary_locus",
    })

    canonical_rows.append(
        r
    )


canonical_rows.sort(
    key=lambda x: (
        x["coassembly"],
        x["contig"],
        integer(x["start"]),
        integer(x["end"]),
    )
)


# ============================================================
# 4. Cross-check MAG Barrnap calls against coassembly calls
# ============================================================

coassembly_by_contig = defaultdict(list)

for row in canonical_rows:

    coassembly_by_contig[
        (
            row["coassembly"],
            row["contig"],
        )
    ].append(row)


mag_crosscheck = []


for mag_hit in mag_rows:

    key = (
        mag_hit["coassembly"],
        mag_hit["contig"],
    )

    ms = integer(
        mag_hit["start"]
    )

    me = integer(
        mag_hit["end"]
    )

    candidates = []

    for co_hit in coassembly_by_contig.get(
        key,
        []
    ):

        cs = integer(
            co_hit["start"]
        )

        ce = integer(
            co_hit["end"]
        )

        ov = overlap_length(
            ms,
            me,
            cs,
            ce
        )

        if ov <= 0:
            continue

        candidates.append(
            (
                ov,
                co_hit
            )
        )

    candidates.sort(
        key=lambda x: (
            -x[0],
            x[1]["locus_id"],
        )
    )

    if not candidates:

        raise RuntimeError(
            "MAG 16S sin correspondiente coassembly: "
            f"{mag_hit['locus_id']} "
            f"{mag_hit['coassembly']}:{mag_hit['contig']}"
        )

    best_overlap, best = (
        candidates[0]
    )

    mag_sequence = mag_seq[
        mag_hit["locus_id"]
    ]

    co_sequence = coassembly_seq[
        best["locus_id"]
    ]

    same_interval = (
        integer(mag_hit["start"])
        == integer(best["start"])
        and integer(mag_hit["end"])
        == integer(best["end"])
        and mag_hit["strand"]
        == best["strand"]
    )

    exact_sequence = (
        canonical_sequence(mag_sequence)
        == canonical_sequence(co_sequence)
    )

    mag_crosscheck.append({
        "MAG_locus_id":
            mag_hit["locus_id"],

        "MAG":
            mag_hit["source_id"],

        "coassembly":
            mag_hit["coassembly"],

        "contig":
            mag_hit["contig"],

        "MAG_start":
            mag_hit["start"],

        "MAG_end":
            mag_hit["end"],

        "MAG_strand":
            mag_hit["strand"],

        "MAG_length_bp":
            mag_hit["length_bp"],

        "coassembly_locus_id":
            best["locus_id"],

        "coassembly_start":
            best["start"],

        "coassembly_end":
            best["end"],

        "coassembly_strand":
            best["strand"],

        "coassembly_length_bp":
            best["length_bp"],

        "overlap_bp":
            best_overlap,

        "same_interval_and_strand":
            "YES"
            if same_interval
            else "NO",

        "orientation_insensitive_exact_sequence":
            "YES"
            if exact_sequence
            else "NO",

        "interpretation":
            "MAG_call_is_secondary_view_of_coassembly_locus",
    })


if len(mag_crosscheck) != 6:
    raise RuntimeError(
        f"MAG crosscheck esperado=6 "
        f"observado={len(mag_crosscheck)}"
    )


# ============================================================
# 5. Exact sequence clustering of PRIMARY coassembly loci
# ============================================================

clusters = defaultdict(list)


for row in canonical_rows:

    clusters[
        row[
            "orientation_insensitive_sequence_sha256"
        ]
    ].append(row)


cluster_records = []
cluster_member_records = []
representative_sequences = {}


sorted_clusters = sorted(
    clusters.items(),
    key=lambda x: (
        min(
            r["locus_id"]
            for r in x[1]
        )
    )
)


for idx, (
    sha,
    members,
) in enumerate(
    sorted_clusters,
    start=1
):

    cluster_id = (
        f"SSU16S_EXACT{idx:04d}"
    )

    members = sorted(
        members,
        key=lambda x: x["locus_id"]
    )

    representative = members[0]

    rep_locus = (
        representative[
            "locus_id"
        ]
    )

    representative_sequences[
        cluster_id
    ] = coassembly_seq[
        rep_locus
    ]

    coassemblies = sorted(
        {
            r["coassembly"]
            for r in members
        }
    )

    producers = sorted(
        {
            r["coassembly"][0]
            for r in members
        }
    )

    mags = sorted(
        {
            mag
            for r in members
            for mag in r[
                "final18_MAGs_exact"
            ].split(";")
            if mag
        }
    )

    lengths = {
        integer(
            r["length_bp"]
        )
        for r in members
    }

    if len(lengths) != 1:
        raise RuntimeError(
            f"{cluster_id}: exact cluster "
            "con longitudes distintas"
        )

    cluster_records.append({
        "exact_cluster_id":
            cluster_id,

        "sequence_sha256":
            sha,

        "representative_locus_id":
            rep_locus,

        "sequence_length_bp":
            next(iter(lengths)),

        "n_loci":
            len(members),

        "n_coassemblies":
            len(coassemblies),

        "coassemblies":
            ";".join(coassemblies),

        "n_producers":
            len(producers),

        "producers":
            ";".join(producers),

        "both_producers":
            "YES"
            if len(producers) == 2
            else "NO",

        "n_final18_MAGs":
            len(mags),

        "final18_MAGs":
            ";".join(mags),

        "recurrence_interpretation":
            "exact_sequence_recurrence_not_prevalence",
    })

    for r in members:

        cluster_member_records.append({
            "exact_cluster_id":
                cluster_id,

            "locus_id":
                r["locus_id"],

            "coassembly":
                r["coassembly"],

            "contig":
                r["contig"],

            "start":
                r["start"],

            "end":
                r["end"],

            "strand":
                r["strand"],

            "length_bp":
                r["length_bp"],

            "length_class":
                r["length_class"],

            "edge_truncation_possible":
                r[
                    "edge_truncation_possible"
                ],

            "final18_MAGs_exact":
                r[
                    "final18_MAGs_exact"
                ],
        })


cluster_lookup = {
    r["locus_id"]:
    r["exact_cluster_id"]
    for r in cluster_member_records
}


for row in canonical_rows:

    row[
        "exact_sequence_cluster_id"
    ] = cluster_lookup[
        row["locus_id"]
    ]


# ============================================================
# 6. Write primary catalog
# ============================================================

canonical_fields = list(
    canonical_rows[0].keys()
)

write_tsv(
    OUT
    / "95B2_coassembly_16S_canonical.tsv",
    canonical_rows,
    canonical_fields
)


# ============================================================
# 7. Write all primary sequences
# ============================================================

with (
    OUT
    / "95B2_coassembly_16S_all.fasta"
).open("w") as fh:

    for row in canonical_rows:

        locus = row["locus_id"]
        seq = coassembly_seq[locus]

        fh.write(
            ">"
            + locus
            + " "
            + f"coassembly={row['coassembly']} "
            + f"contig={row['contig']} "
            + f"coords={row['start']}-{row['end']} "
            + f"strand={row['strand']} "
            + f"length={row['length_bp']} "
            + f"class={row['length_class']}"
            + "\n"
        )

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
# 8. Cluster outputs
# ============================================================

cluster_fields = [
    "exact_cluster_id",
    "sequence_sha256",
    "representative_locus_id",
    "sequence_length_bp",
    "n_loci",
    "n_coassemblies",
    "coassemblies",
    "n_producers",
    "producers",
    "both_producers",
    "n_final18_MAGs",
    "final18_MAGs",
    "recurrence_interpretation",
]

write_tsv(
    OUT
    / "95B2_exact_sequence_clusters.tsv",
    cluster_records,
    cluster_fields
)


member_fields = [
    "exact_cluster_id",
    "locus_id",
    "coassembly",
    "contig",
    "start",
    "end",
    "strand",
    "length_bp",
    "length_class",
    "edge_truncation_possible",
    "final18_MAGs_exact",
]

write_tsv(
    OUT
    / "95B2_exact_sequence_cluster_members.tsv",
    cluster_member_records,
    member_fields
)


with (
    OUT
    / "95B2_unique_exact_representatives.fasta"
).open("w") as fh:

    for cluster in cluster_records:

        cid = cluster[
            "exact_cluster_id"
        ]

        seq = representative_sequences[
            cid
        ]

        fh.write(
            ">"
            + cid
            + " "
            + f"representative={cluster['representative_locus_id']} "
            + f"length={cluster['sequence_length_bp']} "
            + f"n_loci={cluster['n_loci']}"
            + "\n"
        )

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
# 9. MAG cross-check
# ============================================================

mag_cross_fields = [
    "MAG_locus_id",
    "MAG",
    "coassembly",
    "contig",
    "MAG_start",
    "MAG_end",
    "MAG_strand",
    "MAG_length_bp",
    "coassembly_locus_id",
    "coassembly_start",
    "coassembly_end",
    "coassembly_strand",
    "coassembly_length_bp",
    "overlap_bp",
    "same_interval_and_strand",
    "orientation_insensitive_exact_sequence",
    "interpretation",
]

write_tsv(
    OUT
    / "95B2_MAG_16S_crosscheck.tsv",
    mag_crosscheck,
    mag_cross_fields
)


# ============================================================
# 10. Summary by coassembly
# ============================================================

summary_coassembly = []


for group in [
    "L1",
    "L2",
    "L3",
    "M1",
    "M2",
    "M3",
]:

    rr = [
        r
        for r in canonical_rows
        if r["coassembly"] == group
    ]

    summary_coassembly.append({
        "coassembly":
            group,

        "16S_loci":
            len(rr),

        "near_full_ge1300":
            sum(
                integer(r["length_bp"])
                >= 1300
                for r in rr
            ),

        "fragment_800_1299":
            sum(
                800
                <= integer(r["length_bp"])
                < 1300
                for r in rr
            ),

        "fragment_lt800":
            sum(
                integer(r["length_bp"])
                < 800
                for r in rr
            ),

        "edge_touching":
            sum(
                r[
                    "edge_truncation_possible"
                ]
                == "YES"
                for r in rr
            ),

        "exact_final18_contig_membership":
            sum(
                r[
                    "exact_final18_contig_membership"
                ]
                == "YES"
                for r in rr
            ),

        "unique_exact_sequence_clusters":
            len(
                {
                    r[
                        "exact_sequence_cluster_id"
                    ]
                    for r in rr
                }
            ),
    })


write_tsv(
    OUT
    / "95B2_summary_by_coassembly.tsv",
    summary_coassembly,
    list(
        summary_coassembly[0].keys()
    )
)


# ============================================================
# 11. Summary by MAG
# ============================================================

mag_summary = []


for entry in final18_entries:

    mag = entry["MAG"]

    accepted_calls = [
        r
        for r in mag_rows
        if r["source_id"] == mag
    ]

    primary_loci = [
        r
        for r in canonical_rows
        if mag in r[
            "final18_MAGs_exact"
        ].split(";")
    ]

    mag_summary.append({
        "MAG":
            mag,

        "coassembly":
            entry["coassembly"],

        "Barrnap_16S_calls_in_MAG":
            len(accepted_calls),

        "primary_coassembly_16S_loci_on_MAG_contigs":
            len(primary_loci),

        "near_full_ge1300_in_MAG":
            sum(
                integer(r["length_bp"])
                >= 1300
                for r in accepted_calls
            ),

        "fragment_lt1300_in_MAG":
            sum(
                integer(r["length_bp"])
                < 1300
                for r in accepted_calls
            ),

        "absence_interpretation":
            "negative_not_taxon_absence",
    })


write_tsv(
    OUT
    / "95B2_summary_by_MAG.tsv",
    mag_summary,
    list(
        mag_summary[0].keys()
    )
)


# ============================================================
# 12. Global summary
# ============================================================

assigned_primary = [
    r
    for r in canonical_rows
    if r[
        "exact_final18_contig_membership"
    ] == "YES"
]

multimember_clusters = [
    r
    for r in cluster_records
    if integer(r["n_loci"]) > 1
]

cross_coassembly_clusters = [
    r
    for r in cluster_records
    if integer(r["n_coassemblies"]) > 1
]

both_producer_clusters = [
    r
    for r in cluster_records
    if r["both_producers"] == "YES"
]


metrics = [
    (
        "primary_coassembly_16S_loci",
        len(canonical_rows)
    ),
    (
        "primary_16S_near_full_ge1300",
        sum(
            integer(r["length_bp"]) >= 1300
            for r in canonical_rows
        )
    ),
    (
        "primary_16S_800_1299",
        sum(
            800 <= integer(r["length_bp"]) < 1300
            for r in canonical_rows
        )
    ),
    (
        "primary_16S_lt800",
        sum(
            integer(r["length_bp"]) < 800
            for r in canonical_rows
        )
    ),
    (
        "primary_16S_edge_touching",
        sum(
            r["edge_truncation_possible"] == "YES"
            for r in canonical_rows
        )
    ),
    (
        "primary_16S_on_final18_MAG_contigs",
        len(assigned_primary)
    ),
    (
        "MAG_Barrnap_16S_calls",
        len(mag_rows)
    ),
    (
        "MAGs_with_accepted_16S",
        len(
            {
                r["source_id"]
                for r in mag_rows
            }
        )
    ),
    (
        "MAG_calls_crosschecked_to_coassembly",
        len(mag_crosscheck)
    ),
    (
        "MAG_calls_same_interval_and_strand",
        sum(
            r["same_interval_and_strand"] == "YES"
            for r in mag_crosscheck
        )
    ),
    (
        "MAG_calls_exact_sequence_to_coassembly",
        sum(
            r[
                "orientation_insensitive_exact_sequence"
            ] == "YES"
            for r in mag_crosscheck
        )
    ),
    (
        "orientation_insensitive_exact_sequence_clusters",
        len(cluster_records)
    ),
    (
        "exact_multimember_clusters",
        len(multimember_clusters)
    ),
    (
        "exact_cross_coassembly_clusters",
        len(cross_coassembly_clusters)
    ),
    (
        "exact_both_producer_clusters",
        len(both_producer_clusters)
    ),
]


with (
    OUT
    / "95B2_global_summary.tsv"
).open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    for key, value in metrics:
        fh.write(
            f"{key}\t{value}\n"
        )


# ============================================================
# 13. Methodological scope
# ============================================================

with (
    OUT
    / "95B2_methodological_scope.tsv"
).open("w") as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        "primary_catalog\tcoassembly_Barrnap_16S\n"
    )

    fh.write(
        "MAG_calls_added_as_independent_loci\tNO\n"
    )

    fh.write(
        "MAG_association_method\t"
        "exact_coassembly_contig_membership\n"
    )

    fh.write(
        "MAG_Barrnap_calls_role\tsecondary_crosscheck\n"
    )

    fh.write(
        "near_full_definition\tlength_ge1300bp_descriptive_only\n"
    )

    fh.write(
        "near_full_equals_complete_16S\tNO\n"
    )

    fh.write(
        "edge_touching_equals_complete_boundary\tNO\n"
    )

    fh.write(
        "Barrnap_rejected_ultrashort_predictions_in_primary_catalog\tNO\n"
    )

    fh.write(
        "exact_sequence_cluster_orientation\tstrand_insensitive_SHA256\n"
    )

    fh.write(
        "exact_sequence_recurrence_equals_prevalence\tNO\n"
    )

    fh.write(
        "assembled_16S_locus_count_equals_taxon_abundance\tNO\n"
    )

    fh.write(
        "taxonomy_target\tunique_exact_sequence_representatives\n"
    )


# ============================================================
# 14. Final validation
# ============================================================

required = [
    "95B2_coassembly_16S_canonical.tsv",
    "95B2_coassembly_16S_all.fasta",
    "95B2_exact_sequence_clusters.tsv",
    "95B2_exact_sequence_cluster_members.tsv",
    "95B2_unique_exact_representatives.fasta",
    "95B2_MAG_16S_crosscheck.tsv",
    "95B2_summary_by_coassembly.tsv",
    "95B2_summary_by_MAG.tsv",
    "95B2_global_summary.tsv",
    "95B2_methodological_scope.tsv",
]


for name in required:

    if not (
        OUT / name
    ).is_file():

        raise RuntimeError(
            f"Salida faltante: {name}"
        )


print(
    f"PRIMARY_16S={len(canonical_rows)}"
)

print(
    f"MAG_16S_CALLS={len(mag_rows)}"
)

print(
    f"EXACT_SEQUENCE_CLUSTERS={len(cluster_records)}"
)

print(
    f"PRIMARY_ON_FINAL18_CONTIGS={len(assigned_primary)}"
)

print(
    "95B2_CONSOLIDATION=PASS"
)

