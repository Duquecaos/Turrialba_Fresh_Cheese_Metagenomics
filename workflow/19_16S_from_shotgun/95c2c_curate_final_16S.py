#!/usr/bin/env python3

import csv
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path(sys.argv[1])

C2B = (
    ROOT
    / "97_16S_shotgun"
    / "95C2B_final_bacterial_catalog"
)

B2 = (
    ROOT
    / "97_16S_shotgun"
    / "95B2_consolidated"
)

OUT = (
    ROOT
    / "97_16S_shotgun"
    / "95C2C_final_curated_bacterial_catalog"
)

REP = (
    C2B
    / "95C2B_exact_representatives_final_QC.tsv"
)

PRIMARY = (
    C2B
    / "95C2B_primary_158_loci_final_QC.tsv"
)

FASTA = (
    B2
    / "95B2_unique_exact_representatives.fasta"
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

        return (
            list(reader),
            reader.fieldnames or []
        )


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

    seqs = {}

    current = None
    chunks = []

    with path.open() as fh:

        for line in fh:

            line = line.rstrip("\n")

            if line.startswith(">"):

                if current is not None:
                    seqs[current] = "".join(
                        chunks
                    ).upper()

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
            seqs[current] = "".join(
                chunks
            ).upper()

    return seqs


def write_fasta(path, ids, seqs):

    with path.open("w") as fh:

        for seqid in ids:

            if seqid not in seqs:
                raise RuntimeError(
                    f"Secuencia faltante: {seqid}"
                )

            seq = seqs[seqid]

            fh.write(
                f">{seqid}\n"
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


def yes(value):

    return (
        str(value)
        .strip()
        .lower()
        in {
            "yes",
            "true",
            "1",
            "y",
        }
    )


def record_text(row):

    fields = [
        "deepest_informative_LCA_taxon",
        "best_subject_tax_label",
        "best_subject_SILVA_lineage",
        "SILVA_near_best_LCA_lineage",
        "SILVA_near_best_LCA_terminal",
    ]

    return " ; ".join(
        str(
            row.get(
                field,
                ""
            )
        )
        for field in fields
    ).lower()


def known_nonbacterial_flag(row):

    text = record_text(
        row
    )

    terms = [
        "gossypium",
        "eukaryota",
        "mitochondria",
        "mitochondrion",
        "chloroplast",
        "plastid",
    ]

    hits = [
        x
        for x in terms
        if x in text
    ]

    return hits


# ============================================================
# Load
# ============================================================

rep, rep_fields = read_tsv(
    REP
)

primary, primary_fields = read_tsv(
    PRIMARY
)

seqs = read_fasta(
    FASTA
)


if len(rep) != 148:
    raise RuntimeError(
        f"Representatives={len(rep)} expected=148"
    )

if len(primary) != 158:
    raise RuntimeError(
        f"Primary={len(primary)} expected=158"
    )

if len(seqs) != 148:
    raise RuntimeError(
        f"FASTA={len(seqs)} expected=148"
    )


# ============================================================
# Explicit curated exclusion
# ============================================================

final_rows = []
newly_excluded = []


for row in rep:

    old_target = yes(
        row[
            "revised_bacterial_quantification_target"
        ]
    )

    text = record_text(
        row
    )

    # --------------------------------------------------------
    # Curated biological correction:
    # Gossypium hirsutum is a eukaryotic plant and therefore
    # cannot be retained in the bacterial 16S catalog.
    # --------------------------------------------------------

    curated_gossypium = (
        "gossypium" in text
    )

    reasons = []

    if curated_gossypium:
        reasons.append(
            "curated_eukaryotic_taxon_Gossypium"
        )


    final_target = (
        old_target
        and not curated_gossypium
    )


    rec = dict(row)

    rec.update({
        "manual_curated_eukaryote_exclusion":
            "YES"
            if curated_gossypium
            else "NO",

        "manual_curated_exclusion_reason":
            ";".join(reasons),

        "final_bacterial_quantification_target_v2":
            "YES"
            if final_target
            else "NO",

        "final_taxonomy_reporting_eligible_v2":
            "YES"
            if final_target
            else "NO",

        "changed_from_95C2B_target":
            "YES"
            if old_target != final_target
            else "NO",

        "final_catalog_status":
            (
                "accepted_bacterial_16S"
                if final_target
                else "excluded_from_primary_bacterial_catalog"
            ),
    })

    final_rows.append(
        rec
    )

    if (
        old_target
        and not final_target
    ):
        newly_excluded.append(
            rec
        )


final_fields = list(
    final_rows[0].keys()
)


write_tsv(
    OUT
    / "95C2C_all_representatives_final_QC.tsv",
    final_rows,
    final_fields
)


write_tsv(
    OUT
    / "95C2C_newly_excluded_curated.tsv",
    newly_excluded,
    final_fields
)


# There must be at least one correction because 95C2B
# explicitly reported Gossypium in the accepted table.

if len(newly_excluded) < 1:
    raise RuntimeError(
        "No se encontró el Gossypium que debía ser excluido"
    )


# ============================================================
# Final accepted catalog
# ============================================================

accepted = [
    x
    for x in final_rows
    if x[
        "final_bacterial_quantification_target_v2"
    ] == "YES"
]


write_tsv(
    OUT
    / "95C2C_final_bacterial_targets.tsv",
    accepted,
    final_fields
)


accepted_ids = [
    x["exact_cluster_id"]
    for x in accepted
]


write_fasta(
    OUT
    / "95C2C_final_bacterial_16S.fasta",
    accepted_ids,
    seqs
)


# ============================================================
# Propagate to 158 original loci
# ============================================================

final_by_cluster = {
    x["exact_cluster_id"]: x
    for x in final_rows
}


primary_final = []


for row in primary:

    cluster = (
        row[
            "taxonomy_exact_cluster_id"
        ]
    )

    if cluster not in final_by_cluster:
        raise RuntimeError(
            f"Cluster faltante: {cluster}"
        )

    q = final_by_cluster[
        cluster
    ]

    rec = dict(row)

    rec.update({
        "manual_curated_eukaryote_exclusion":
            q[
                "manual_curated_eukaryote_exclusion"
            ],

        "manual_curated_exclusion_reason":
            q[
                "manual_curated_exclusion_reason"
            ],

        "final_bacterial_quantification_target_v2":
            q[
                "final_bacterial_quantification_target_v2"
            ],

        "final_catalog_status":
            q[
                "final_catalog_status"
            ],
    })

    primary_final.append(
        rec
    )


primary_final_fields = list(
    primary_final[0].keys()
)


write_tsv(
    OUT
    / "95C2C_primary_158_loci_final_QC.tsv",
    primary_final,
    primary_final_fields
)


# ============================================================
# Final taxon summary
# ============================================================

groups = defaultdict(list)


for row in accepted:

    taxon = (
        row.get(
            "deepest_informative_LCA_taxon",
            ""
        )
        or "UNRESOLVED"
    )

    groups[taxon].append(
        row
    )


summary = []


for taxon, rows in groups.items():

    coassemblies = sorted({
        x
        for r in rows
        for x in r[
            "coassemblies"
        ].split(";")
        if x
    })

    producers = sorted({
        x
        for r in rows
        for x in r[
            "producers"
        ].split(";")
        if x
    })

    summary.append({
        "deepest_informative_LCA_taxon":
            taxon,

        "exact_sequence_clusters":
            len(rows),

        "coassemblies":
            ";".join(
                coassemblies
            ),

        "n_coassemblies":
            len(coassemblies),

        "producers":
            ";".join(
                producers
            ),

        "n_producers":
            len(producers),

        "A_high_similarity_ge1300":
            sum(
                r[
                    "revised_resolution_tier"
                ]
                ==
                "A_high_similarity_ge1300"
                for r in rows
            ),

        "B_moderate_similarity":
            sum(
                r[
                    "revised_resolution_tier"
                ]
                ==
                "B_moderate_similarity"
                for r in rows
            ),

        "C_fragmentary_or_broad":
            sum(
                r[
                    "revised_resolution_tier"
                ]
                ==
                "C_fragmentary_or_broad"
                for r in rows
            ),
    })


summary.sort(
    key=lambda x: (
        -x[
            "exact_sequence_clusters"
        ],
        x[
            "deepest_informative_LCA_taxon"
        ],
    )
)


write_tsv(
    OUT
    / "95C2C_final_bacterial_taxon_summary.tsv",
    summary,
    [
        "deepest_informative_LCA_taxon",
        "exact_sequence_clusters",
        "coassemblies",
        "n_coassemblies",
        "producers",
        "n_producers",
        "A_high_similarity_ge1300",
        "B_moderate_similarity",
        "C_fragmentary_or_broad",
    ]
)


# ============================================================
# Known non-bacterial flag audit
# ============================================================

remaining_flags = []


for row in accepted:

    flags = known_nonbacterial_flag(
        row
    )

    if not flags:
        continue

    rec = dict(row)

    rec[
        "remaining_known_nonbacterial_flags"
    ] = ";".join(flags)

    remaining_flags.append(
        rec
    )


flag_fields = (
    final_fields
    + [
        "remaining_known_nonbacterial_flags"
    ]
)


write_tsv(
    OUT
    / "95C2C_remaining_known_nonbacterial_flags.tsv",
    remaining_flags,
    flag_fields
)


if remaining_flags:
    raise RuntimeError(
        "Persisten señales eucariotas/organelarias conocidas "
        f"en targets finales: {len(remaining_flags)}"
    )


# ============================================================
# Global summary
# ============================================================

old_targets = sum(
    yes(
        x[
            "revised_bacterial_quantification_target"
        ]
    )
    for x in final_rows
)

primary_targets = sum(
    x[
        "final_bacterial_quantification_target_v2"
    ]
    == "YES"
    for x in primary_final
)


metrics = [
    (
        "exact_sequence_clusters_total",
        len(final_rows)
    ),
    (
        "95C2B_targets_before_manual_curation",
        old_targets
    ),
    (
        "curated_clusters_removed",
        len(newly_excluded)
    ),
    (
        "final_bacterial_16S_clusters",
        len(accepted)
    ),
    (
        "primary_16S_loci_total",
        len(primary_final)
    ),
    (
        "primary_final_bacterial_16S_loci",
        primary_targets
    ),
    (
        "remaining_known_nonbacterial_flags",
        len(remaining_flags)
    ),
]


with (
    OUT
    / "95C2C_global_summary.tsv"
).open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    for key, value in metrics:

        fh.write(
            f"{key}\t{value}\n"
        )


# ============================================================
# Scope
# ============================================================

with (
    OUT
    / "95C2C_methodological_scope.tsv"
).open("w") as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        "purpose\tfinal_curated_bacterial_16S_catalog\n"
    )

    fh.write(
        "manual_curated_exclusion\t"
        "Gossypium_hirsutum_eukaryotic_plant\n"
    )

    fh.write(
        "Barrnap_rerun\tNO\n"
    )

    fh.write(
        "SILVA_BLAST_rerun\tNO\n"
    )

    fh.write(
        "nonbacterial_sequences_in_primary_quantification\tNO\n"
    )

    fh.write(
        "low_resolution_bacterial_sequences_in_primary_quantification\tNO\n"
    )

    fh.write(
        "final_catalog_unit\torientation_insensitive_exact_16S_sequence_cluster\n"
    )

    fh.write(
        "cluster_count_equals_abundance\tNO\n"
    )


# ============================================================
# Validate
# ============================================================

required = [
    "95C2C_all_representatives_final_QC.tsv",
    "95C2C_newly_excluded_curated.tsv",
    "95C2C_final_bacterial_targets.tsv",
    "95C2C_final_bacterial_16S.fasta",
    "95C2C_primary_158_loci_final_QC.tsv",
    "95C2C_final_bacterial_taxon_summary.tsv",
    "95C2C_remaining_known_nonbacterial_flags.tsv",
    "95C2C_global_summary.tsv",
    "95C2C_methodological_scope.tsv",
]


for name in required:

    if not (
        OUT / name
    ).is_file():

        raise RuntimeError(
            f"Salida faltante: {name}"
        )


print(
    f"TARGETS_BEFORE={old_targets}"
)

print(
    f"CURATED_REMOVED={len(newly_excluded)}"
)

print(
    f"FINAL_BACTERIAL_CLUSTERS={len(accepted)}"
)

print(
    f"PRIMARY_FINAL_BACTERIAL_LOCI={primary_targets}"
)

print(
    "95C2C_FINAL_CATALOG=PASS"
)

