#!/usr/bin/env python3

import csv
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(sys.argv[1])

C2 = (
    ROOT
    / "97_16S_shotgun"
    / "95C2_bacterial_QC"
)

B2 = (
    ROOT
    / "97_16S_shotgun"
    / "95B2_consolidated"
)

OUT = (
    ROOT
    / "97_16S_shotgun"
    / "95C2B_final_bacterial_catalog"
)

REP = (
    C2
    / "95C2_exact_representatives_QC.tsv"
)

PRIMARY = (
    C2
    / "95C2_primary_158_loci_QC.tsv"
)

FASTA = (
    B2
    / "95B2_unique_exact_representatives.fasta"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


def read_tsv(path):

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


def write_tsv(
    path,
    rows,
    fields
):

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


def write_fasta(
    path,
    ids,
    seqs
):

    with path.open("w") as fh:

        for seqid in ids:

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


def euk_tax_label(value):

    v = (
        value or ""
    ).strip().lower()

    return (
        v == "d__eukaryota"
        or v.startswith(
            "d__eukaryota;"
        )
    )


def organelle_keyword(value):

    v = (
        value or ""
    ).lower()

    return any(
        x in v
        for x in (
            "mitochondria",
            "mitochondrion",
            "chloroplast",
            "plastid",
        )
    )


def starts_eukaryota(value):

    v = (
        value or ""
    ).strip().lower()

    return (
        v == "eukaryota"
        or v.startswith(
            "eukaryota;"
        )
    )


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
# Revised representative classification
# ============================================================

revised = []

newly_excluded = []


for row in rep:

    cluster = (
        row[
            "exact_cluster_id"
        ]
    )

    best_label = (
        row.get(
            "best_subject_tax_label",
            ""
        )
    )

    best_lineage = (
        row.get(
            "best_subject_SILVA_lineage",
            ""
        )
    )

    lca_lineage = (
        row.get(
            "SILVA_near_best_LCA_lineage",
            ""
        )
    )

    old_context = (
        row[
            "final_lineage_context"
        ]
    )

    old_target = yes(
        row[
            "bacterial_quantification_target"
        ]
    )


    flag_taxmap_euk = (
        euk_tax_label(
            best_label
        )
    )

    flag_best_euk = (
        starts_eukaryota(
            best_lineage
        )
    )

    flag_lca_euk = (
        starts_eukaryota(
            lca_lineage
        )
    )

    flag_organelle = (
        organelle_keyword(
            best_lineage
        )
        or
        organelle_keyword(
            lca_lineage
        )
    )


    exclusion_reasons = []

    if flag_taxmap_euk:
        exclusion_reasons.append(
            "best_subject_taxmap_d__Eukaryota"
        )

    if flag_best_euk:
        exclusion_reasons.append(
            "best_subject_lineage_Eukaryota"
        )

    if flag_lca_euk:
        exclusion_reasons.append(
            "near_best_LCA_Eukaryota"
        )

    if flag_organelle:
        exclusion_reasons.append(
            "organelle_keyword"
        )


    revised_excluded = bool(
        exclusion_reasons
    )


    if revised_excluded:

        revised_context = (
            "eukaryotic_or_organelle"
        )

        revised_tier = (
            "EXCLUDED_nonbacterial_or_organelle"
        )

        revised_target = False

        reportable = False

    else:

        revised_context = old_context

        revised_tier = (
            row[
                "final_resolution_tier"
            ]
        )

        revised_target = (
            revised_context
            == "bacterial_nonorganelle"
            and revised_tier
            in {
                "A_high_similarity_ge1300",
                "B_moderate_similarity",
                "C_fragmentary_or_broad",
            }
        )

        reportable = (
            revised_target
        )


    rec = dict(row)

    rec.update({
        "filter_flag_taxmap_Eukaryota":
            "YES"
            if flag_taxmap_euk
            else "NO",

        "filter_flag_best_lineage_Eukaryota":
            "YES"
            if flag_best_euk
            else "NO",

        "filter_flag_LCA_Eukaryota":
            "YES"
            if flag_lca_euk
            else "NO",

        "filter_flag_organelle_keyword":
            "YES"
            if flag_organelle
            else "NO",

        "final_exclusion_reasons":
            ";".join(
                exclusion_reasons
            ),

        "revised_lineage_context":
            revised_context,

        "revised_resolution_tier":
            revised_tier,

        "revised_taxonomy_reporting_eligible":
            "YES"
            if reportable
            else "NO",

        "revised_bacterial_quantification_target":
            "YES"
            if revised_target
            else "NO",

        "changed_from_95C2_target":
            "YES"
            if old_target != revised_target
            else "NO",

        "final_catalog_guardrail":
            "eukaryotic_or_organelle_sequences_excluded_from_bacterial_quantification",
    })

    revised.append(
        rec
    )

    if (
        old_target
        and not revised_target
    ):
        newly_excluded.append(
            rec
        )


revised_fields = list(
    revised[0].keys()
)


write_tsv(
    OUT
    / "95C2B_exact_representatives_final_QC.tsv",
    revised,
    revised_fields
)

write_tsv(
    OUT
    / "95C2B_newly_excluded_from_95C2.tsv",
    newly_excluded,
    revised_fields
)


# ============================================================
# Final accepted/excluded sets
# ============================================================

accepted = [
    x
    for x in revised
    if x[
        "revised_bacterial_quantification_target"
    ] == "YES"
]

excluded = [
    x
    for x in revised
    if x[
        "revised_lineage_context"
    ] == "eukaryotic_or_organelle"
]

low_bacterial = [
    x
    for x in revised
    if (
        x[
            "revised_lineage_context"
        ]
        == "bacterial_nonorganelle"
        and
        x[
            "revised_resolution_tier"
        ]
        == "D_low_resolution"
    )
]


write_tsv(
    OUT
    / "95C2B_final_bacterial_targets.tsv",
    accepted,
    revised_fields
)

write_tsv(
    OUT
    / "95C2B_final_excluded_eukaryotic_or_organelle.tsv",
    excluded,
    revised_fields
)

write_tsv(
    OUT
    / "95C2B_final_low_resolution_bacterial.tsv",
    low_bacterial,
    revised_fields
)


accepted_ids = [
    x[
        "exact_cluster_id"
    ]
    for x in accepted
]

write_fasta(
    OUT
    / "95C2B_final_bacterial_quantification_targets.fasta",
    accepted_ids,
    seqs
)


# ============================================================
# Propagate revised status to 158 primary loci
# ============================================================

rev_by_cluster = {
    x[
        "exact_cluster_id"
    ]: x
    for x in revised
}


primary_final = []

for row in primary:

    cluster = (
        row[
            "taxonomy_exact_cluster_id"
        ]
    )

    q = rev_by_cluster[
        cluster
    ]

    rec = dict(row)

    rec.update({
        "revised_lineage_context":
            q[
                "revised_lineage_context"
            ],

        "revised_resolution_tier":
            q[
                "revised_resolution_tier"
            ],

        "revised_bacterial_quantification_target":
            q[
                "revised_bacterial_quantification_target"
            ],

        "final_exclusion_reasons":
            q[
                "final_exclusion_reasons"
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
    / "95C2B_primary_158_loci_final_QC.tsv",
    primary_final,
    primary_final_fields
)


# ============================================================
# Final taxon summary
# ============================================================

tax_groups = {}

for row in accepted:

    taxon = (
        row[
            "deepest_informative_LCA_taxon"
        ]
        or "UNRESOLVED"
    )

    tax_groups.setdefault(
        taxon,
        []
    ).append(row)


summary = []

for taxon, rows in tax_groups.items():

    coassemblies = sorted({
        x
        for row in rows
        for x in row[
            "coassemblies"
        ].split(";")
        if x
    })

    producers = sorted({
        x
        for row in rows
        for x in row[
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
                x[
                    "revised_resolution_tier"
                ]
                ==
                "A_high_similarity_ge1300"
                for x in rows
            ),

        "B_moderate_similarity":
            sum(
                x[
                    "revised_resolution_tier"
                ]
                ==
                "B_moderate_similarity"
                for x in rows
            ),

        "C_fragmentary_or_broad":
            sum(
                x[
                    "revised_resolution_tier"
                ]
                ==
                "C_fragmentary_or_broad"
                for x in rows
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
    / "95C2B_final_bacterial_taxon_summary.tsv",
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
# Audit contradictions
# ============================================================

contradictions = [
    x
    for x in accepted
    if (
        euk_tax_label(
            x.get(
                "best_subject_tax_label",
                ""
            )
        )
        or
        starts_eukaryota(
            x.get(
                "best_subject_SILVA_lineage",
                ""
            )
        )
        or
        starts_eukaryota(
            x.get(
                "SILVA_near_best_LCA_lineage",
                ""
            )
        )
        or
        organelle_keyword(
            x.get(
                "best_subject_SILVA_lineage",
                ""
            )
        )
        or
        organelle_keyword(
            x.get(
                "SILVA_near_best_LCA_lineage",
                ""
            )
        )
    )
]


write_tsv(
    OUT
    / "95C2B_remaining_filter_contradictions.tsv",
    contradictions,
    revised_fields
)


if contradictions:

    raise RuntimeError(
        "Quedaron contradicciones eucariotas/organelarias "
        f"en targets: {len(contradictions)}"
    )


# ============================================================
# Summary
# ============================================================

context_counts = Counter(
    x[
        "revised_lineage_context"
    ]
    for x in revised
)


metrics = [
    (
        "exact_sequence_clusters_total",
        len(revised)
    ),
    (
        "95C2_original_quantification_targets",
        sum(
            yes(
                x[
                    "bacterial_quantification_target"
                ]
            )
            for x in revised
        )
    ),
    (
        "newly_excluded_from_95C2_targets",
        len(
            newly_excluded
        )
    ),
    (
        "final_bacterial_quantification_targets",
        len(
            accepted
        )
    ),
    (
        "final_eukaryotic_or_organelle_clusters",
        len(
            excluded
        )
    ),
    (
        "final_low_resolution_bacterial_clusters",
        len(
            low_bacterial
        )
    ),
    (
        "primary_16S_loci_total",
        len(
            primary_final
        )
    ),
    (
        "primary_final_bacterial_target_loci",
        sum(
            x[
                "revised_bacterial_quantification_target"
            ]
            == "YES"
            for x in primary_final
        )
    ),
    (
        "primary_final_excluded_eukaryotic_or_organelle_loci",
        sum(
            x[
                "revised_lineage_context"
            ]
            == "eukaryotic_or_organelle"
            for x in primary_final
        )
    ),
    (
        "remaining_filter_contradictions",
        len(
            contradictions
        )
    ),
]


with (
    OUT
    / "95C2B_global_summary.tsv"
).open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    for key, value in metrics:

        fh.write(
            f"{key}\t{value}\n"
        )


with (
    OUT
    / "95C2B_methodological_scope.tsv"
).open("w") as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        "purpose\tfinal_bacterial_16S_catalog_before_read_quantification\n"
    )

    fh.write(
        "additional_exclusion_guard\t"
        "best_subject_taxmap_d__Eukaryota\n"
    )

    fh.write(
        "lineage_guards\t"
        "Eukaryota_Mitochondria_Chloroplast_plastid\n"
    )

    fh.write(
        "rerun_BLAST_required\tNO\n"
    )

    fh.write(
        "rerun_Barrnap_required\tNO\n"
    )

    fh.write(
        "excluded_sequences_used_for_bacterial_quantification\tNO\n"
    )

    fh.write(
        "low_resolution_bacterial_used_for_primary_quantification\tNO\n"
    )

    fh.write(
        "exact_sequence_cluster_count_equals_abundance\tNO\n"
    )


print(
    f"ORIGINAL_TARGETS="
    f"{sum(yes(x['bacterial_quantification_target']) for x in revised)}"
)

print(
    f"NEWLY_EXCLUDED="
    f"{len(newly_excluded)}"
)

print(
    f"FINAL_TARGETS="
    f"{len(accepted)}"
)

print(
    f"FINAL_EUK_ORGANELLE="
    f"{len(excluded)}"
)

print(
    f"LOW_RES_BACTERIAL="
    f"{len(low_bacterial)}"
)

print(
    "95C2B_FINAL_FILTER=PASS"
)

