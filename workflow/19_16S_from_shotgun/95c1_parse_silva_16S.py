#!/usr/bin/env python3

import csv
import sys
from collections import defaultdict, Counter
from pathlib import Path


ROOT = Path(sys.argv[1])

BASE95B2 = (
    ROOT
    / "97_16S_shotgun"
    / "95B2_consolidated"
)

OUT = (
    ROOT
    / "97_16S_shotgun"
    / "95C1_SILVA_taxonomy"
)

RAW = (
    OUT
    / "95C1_SILVA1382_blast_raw.tsv"
)

TAXMAP = (
    ROOT
    / "databases/mdmcleaner"
    / "lookup_repair_20260911"
    / "silva1382_acc2taxid.tsv"
)

CLUSTERS = (
    BASE95B2
    / "95B2_exact_sequence_clusters.tsv"
)

MEMBERS = (
    BASE95B2
    / "95B2_exact_sequence_cluster_members.tsv"
)

CANONICAL = (
    BASE95B2
    / "95B2_coassembly_16S_canonical.tsv"
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
        errors="replace",
        newline=""
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


def number(x):

    try:
        return float(x)
    except Exception:
        return None


def integer(x):

    return int(float(x))


def unique_join(values):

    return ";".join(
        sorted(
            {
                str(v).strip()
                for v in values
                if str(v).strip()
            }
        )
    )


def normalize_lineage_from_title(
    sseqid,
    stitle
):
    """
    SILVA tax_silva FASTA normalmente guarda la taxonomía
    en el título de la secuencia.

    Conservamos el lineage como texto SILVA; no inventamos
    rangos a partir de su posición.
    """

    title = (
        stitle or ""
    ).strip()

    if not title:
        return ""

    # Si el accession aparece al inicio del título, retirarlo.
    if title.startswith(sseqid):

        title = title[
            len(sseqid):
        ].strip()

    # Requiere estructura taxonómica delimitada por ';'
    if ";" not in title:
        return ""

    parts = [
        x.strip()
        for x in title.split(";")
        if x.strip()
    ]

    return ";".join(parts)


def lineage_lca(lineages):
    """
    Longest common exact semicolon-delimited prefix.
    """

    parsed = [
        [
            x.strip()
            for x in lin.split(";")
            if x.strip()
        ]
        for lin in lineages
        if lin
    ]

    if not parsed:
        return ""

    n = min(
        len(x)
        for x in parsed
    )

    common = []

    for i in range(n):

        values = {
            x[i]
            for x in parsed
        }

        if len(values) != 1:
            break

        common.append(
            parsed[0][i]
        )

    return ";".join(common)


def terminal_taxon(lineage):

    if not lineage:
        return ""

    parts = [
        x.strip()
        for x in lineage.split(";")
        if x.strip()
    ]

    return (
        parts[-1]
        if parts
        else ""
    )


# ============================================================
# 1. Tax map
# ============================================================

taxmap = {}

with TAXMAP.open(
    errors="replace"
) as fh:

    for line in fh:

        line = line.rstrip("\n")

        if not line:
            continue

        parts = line.split(
            "\t"
        )

        if len(parts) < 2:
            continue

        taxmap[
            parts[0]
        ] = parts[1]


if len(taxmap) < 500000:
    raise RuntimeError(
        f"Tax map inesperadamente pequeño: {len(taxmap)}"
    )


# ============================================================
# 2. BLAST raw
# ============================================================

blast_fields = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
    "qcovs",
    "stitle",
]

blast_rows = []

with RAW.open(
    errors="replace"
) as fh:

    reader = csv.DictReader(
        fh,
        fieldnames=blast_fields,
        delimiter="\t"
    )

    for row in reader:
        blast_rows.append(
            row
        )


hits_by_query = defaultdict(list)

for row in blast_rows:

    row["pident_num"] = number(
        row["pident"]
    )

    row["bitscore_num"] = number(
        row["bitscore"]
    )

    row["qcov_num"] = number(
        row["qcovs"]
    )

    row["qlen_num"] = integer(
        row["qlen"]
    )

    row["subject_tax_label"] = (
        taxmap.get(
            row["sseqid"],
            ""
        )
    )

    row["subject_silva_lineage"] = (
        normalize_lineage_from_title(
            row["sseqid"],
            row["stitle"],
        )
    )

    hits_by_query[
        row["qseqid"]
    ].append(row)


# ============================================================
# 3. Cluster metadata
# ============================================================

clusters, _ = read_tsv(
    CLUSTERS
)

if len(clusters) != 148:
    raise RuntimeError(
        f"Clusters esperado=148 observado={len(clusters)}"
    )

cluster_meta = {
    x["exact_cluster_id"]: x
    for x in clusters
}


# ============================================================
# 4. Conservative assignments
# ============================================================

assignments = []

for cluster in clusters:

    qid = (
        cluster[
            "exact_cluster_id"
        ]
    )

    qlen = integer(
        cluster[
            "sequence_length_bp"
        ]
    )

    hits = hits_by_query.get(
        qid,
        []
    )

    if not hits:

        assignments.append({
            "exact_cluster_id": qid,
            "representative_locus_id":
                cluster[
                    "representative_locus_id"
                ],
            "query_length_bp": qlen,
            "n_BLAST_hits": 0,
            "best_sseqid": "",
            "best_percent_identity": "",
            "best_query_coverage_percent": "",
            "best_bitscore": "",
            "best_evalue": "",
            "best_subject_tax_label": "",
            "best_subject_SILVA_lineage": "",
            "near_best_hit_count": 0,
            "near_best_unique_subjects": 0,
            "near_best_unique_tax_labels": 0,
            "near_best_unique_lineages": 0,
            "SILVA_near_best_LCA_lineage": "",
            "SILVA_near_best_LCA_terminal":
                "",
            "title_lineage_parseable":
                "NO",
            "sequence_similarity_resolution_tier":
                "unresolved_no_hit",
            "n_coassemblies":
                cluster["n_coassemblies"],
            "coassemblies":
                cluster["coassemblies"],
            "n_producers":
                cluster["n_producers"],
            "producers":
                cluster["producers"],
            "both_producers":
                cluster["both_producers"],
            "final18_MAGs":
                cluster["final18_MAGs"],
            "interpretation_guardrail":
                "16S_similarity_assignment_not_organism_proof",
        })

        continue


    # Orden explícito
    hits = sorted(
        hits,
        key=lambda x: (
            -x["bitscore_num"],
            -x["pident_num"],
            -x["qcov_num"],
            x["sseqid"],
        )
    )

    best = hits[0]

    best_bitscore = (
        best["bitscore_num"]
    )

    best_pident = (
        best["pident_num"]
    )

    best_qcov = (
        best["qcov_num"]
    )


    # Near-best set:
    #   >=99% del bitscore máximo
    #   y no más de 0.5 puntos porcentuales debajo del mejor %ID.
    near_best = [
        h
        for h in hits
        if (
            h["bitscore_num"]
            >= 0.99 * best_bitscore
            and
            h["pident_num"]
            >= best_pident - 0.5
        )
    ]


    lineages = [
        h[
            "subject_silva_lineage"
        ]
        for h in near_best
        if h[
            "subject_silva_lineage"
        ]
    ]

    lca = lineage_lca(
        lineages
    )

    lca_terminal = terminal_taxon(
        lca
    )


    # --------------------------------------------------------
    # Resolution tier: descriptive/conservative.
    #
    # NO equivale automáticamente a asignación formal de rango.
    # Fragmentos cortos no reciben "species candidate".
    # --------------------------------------------------------

    if (
        qlen >= 1200
        and best_qcov >= 90
        and best_pident >= 98.7
    ):

        tier = (
            "high_similarity_near_full_candidate"
        )

    elif (
        qlen >= 800
        and best_qcov >= 90
        and best_pident >= 94.5
    ):

        tier = (
            "moderate_high_similarity_candidate"
        )

    elif (
        best_qcov >= 80
        and best_pident >= 90
    ):

        tier = (
            "broad_similarity_or_fragmentary"
        )

    else:

        tier = (
            "low_resolution"
        )


    assignments.append({
        "exact_cluster_id": qid,

        "representative_locus_id":
            cluster[
                "representative_locus_id"
            ],

        "query_length_bp":
            qlen,

        "n_BLAST_hits":
            len(hits),

        "best_sseqid":
            best["sseqid"],

        "best_percent_identity":
            best["pident"],

        "best_query_coverage_percent":
            best["qcovs"],

        "best_bitscore":
            best["bitscore"],

        "best_evalue":
            best["evalue"],

        "best_subject_tax_label":
            best[
                "subject_tax_label"
            ],

        "best_subject_SILVA_lineage":
            best[
                "subject_silva_lineage"
            ],

        "near_best_hit_count":
            len(near_best),

        "near_best_unique_subjects":
            len(
                {
                    h["sseqid"]
                    for h in near_best
                }
            ),

        "near_best_unique_tax_labels":
            len(
                {
                    h[
                        "subject_tax_label"
                    ]
                    for h in near_best
                    if h[
                        "subject_tax_label"
                    ]
                }
            ),

        "near_best_unique_lineages":
            len(
                {
                    h[
                        "subject_silva_lineage"
                    ]
                    for h in near_best
                    if h[
                        "subject_silva_lineage"
                    ]
                }
            ),

        "SILVA_near_best_LCA_lineage":
            lca,

        "SILVA_near_best_LCA_terminal":
            lca_terminal,

        "title_lineage_parseable":
            "YES"
            if lineages
            else "NO",

        "sequence_similarity_resolution_tier":
            tier,

        "n_coassemblies":
            cluster[
                "n_coassemblies"
            ],

        "coassemblies":
            cluster[
                "coassemblies"
            ],

        "n_producers":
            cluster[
                "n_producers"
            ],

        "producers":
            cluster[
                "producers"
            ],

        "both_producers":
            cluster[
                "both_producers"
            ],

        "final18_MAGs":
            cluster[
                "final18_MAGs"
            ],

        "interpretation_guardrail":
            "16S_similarity_assignment_not_organism_proof",
    })


assignment_fields = list(
    assignments[0].keys()
)

write_tsv(
    OUT
    / "95C1_exact_representative_taxonomy.tsv",
    assignments,
    assignment_fields
)


# ============================================================
# 5. Propagate assignments back to 158 primary loci
# ============================================================

members, _ = read_tsv(
    MEMBERS
)

canonical, canonical_fields = read_tsv(
    CANONICAL
)

if len(members) != 158:
    raise RuntimeError(
        f"Cluster members esperado=158 observado={len(members)}"
    )

if len(canonical) != 158:
    raise RuntimeError(
        f"Canonical esperado=158 observado={len(canonical)}"
    )


assignment_by_cluster = {
    x["exact_cluster_id"]: x
    for x in assignments
}

cluster_by_locus = {
    x["locus_id"]:
        x["exact_cluster_id"]
    for x in members
}


primary_taxonomy = []

for row in canonical:

    locus = row["locus_id"]

    cluster_id = (
        cluster_by_locus[
            locus
        ]
    )

    tax = assignment_by_cluster[
        cluster_id
    ]

    rec = dict(row)

    rec.update({
        "taxonomy_exact_cluster_id":
            cluster_id,

        "SILVA_best_sseqid":
            tax["best_sseqid"],

        "SILVA_best_percent_identity":
            tax[
                "best_percent_identity"
            ],

        "SILVA_best_query_coverage_percent":
            tax[
                "best_query_coverage_percent"
            ],

        "SILVA_best_subject_tax_label":
            tax[
                "best_subject_tax_label"
            ],

        "SILVA_best_subject_lineage":
            tax[
                "best_subject_SILVA_lineage"
            ],

        "SILVA_LCA_lineage":
            tax[
                "SILVA_near_best_LCA_lineage"
            ],

        "SILVA_LCA_terminal":
            tax[
                "SILVA_near_best_LCA_terminal"
            ],

        "sequence_similarity_resolution_tier":
            tax[
                "sequence_similarity_resolution_tier"
            ],

        "taxonomy_interpretation_guardrail":
            "16S_similarity_assignment_not_organism_proof",
    })

    primary_taxonomy.append(
        rec
    )


primary_fields = list(
    primary_taxonomy[0].keys()
)

write_tsv(
    OUT
    / "95C1_primary_158_loci_taxonomy.tsv",
    primary_taxonomy,
    primary_fields
)


# ============================================================
# 6. Resolution summary
# ============================================================

tier_counts = Counter(
    x[
        "sequence_similarity_resolution_tier"
    ]
    for x in assignments
)

tier_rows = [
    {
        "resolution_tier": k,
        "exact_sequence_clusters": v,
    }
    for k, v in sorted(
        tier_counts.items()
    )
]

write_tsv(
    OUT
    / "95C1_resolution_summary.tsv",
    tier_rows,
    [
        "resolution_tier",
        "exact_sequence_clusters",
    ]
)


# ============================================================
# 7. LCA terminal summary
# ============================================================

taxon_groups = defaultdict(list)

for row in assignments:

    taxon = (
        row[
            "SILVA_near_best_LCA_terminal"
        ]
        or
        row[
            "best_subject_tax_label"
        ]
        or
        "UNRESOLVED"
    )

    taxon_groups[
        taxon
    ].append(row)


taxon_summary = []

for taxon, rows in taxon_groups.items():

    taxon_summary.append({
        "taxon_or_label":
            taxon,

        "exact_sequence_clusters":
            len(rows),

        "coassemblies":
            unique_join(
                x["coassemblies"]
                for x in rows
            ),

        "producers":
            unique_join(
                x["producers"]
                for x in rows
            ),

        "both_producers_clusters":
            sum(
                x["both_producers"]
                == "YES"
                for x in rows
            ),

        "high_similarity_near_full":
            sum(
                x[
                    "sequence_similarity_resolution_tier"
                ]
                ==
                "high_similarity_near_full_candidate"
                for x in rows
            ),
    })


taxon_summary.sort(
    key=lambda x: (
        -x[
            "exact_sequence_clusters"
        ],
        x["taxon_or_label"],
    )
)


write_tsv(
    OUT
    / "95C1_taxon_summary.tsv",
    taxon_summary,
    list(
        taxon_summary[0].keys()
    )
)


# ============================================================
# 8. By coassembly
# ============================================================

summary_by_coassembly = []

for group in [
    "L1",
    "L2",
    "L3",
    "M1",
    "M2",
    "M3",
]:

    rr = [
        x
        for x in primary_taxonomy
        if x["coassembly"] == group
    ]

    summary_by_coassembly.append({
        "coassembly":
            group,

        "primary_16S_loci":
            len(rr),

        "with_BLAST_assignment":
            sum(
                bool(
                    x[
                        "SILVA_best_sseqid"
                    ]
                )
                for x in rr
            ),

        "high_similarity_near_full":
            sum(
                x[
                    "sequence_similarity_resolution_tier"
                ]
                ==
                "high_similarity_near_full_candidate"
                for x in rr
            ),

        "moderate_high_similarity":
            sum(
                x[
                    "sequence_similarity_resolution_tier"
                ]
                ==
                "moderate_high_similarity_candidate"
                for x in rr
            ),

        "fragmentary_or_broad":
            sum(
                x[
                    "sequence_similarity_resolution_tier"
                ]
                ==
                "broad_similarity_or_fragmentary"
                for x in rr
            ),

        "low_resolution":
            sum(
                x[
                    "sequence_similarity_resolution_tier"
                ]
                ==
                "low_resolution"
                for x in rr
            ),
    })


write_tsv(
    OUT
    / "95C1_summary_by_coassembly.tsv",
    summary_by_coassembly,
    list(
        summary_by_coassembly[0].keys()
    )
)


# ============================================================
# 9. Unresolved/ambiguous
# ============================================================

unresolved = [
    x
    for x in assignments
    if (
        x[
            "sequence_similarity_resolution_tier"
        ]
        in {
            "low_resolution",
            "unresolved_no_hit",
        }
        or
        not x[
            "SILVA_near_best_LCA_lineage"
        ]
    )
]

write_tsv(
    OUT
    / "95C1_unresolved_or_low_resolution.tsv",
    unresolved,
    assignment_fields
)


# ============================================================
# 10. Global summary
# ============================================================

metrics = [
    (
        "exact_sequence_queries",
        len(assignments)
    ),
    (
        "queries_with_any_BLAST_hit",
        sum(
            x["n_BLAST_hits"] > 0
            for x in assignments
        )
    ),
    (
        "queries_with_parseable_SILVA_lineage",
        sum(
            x[
                "title_lineage_parseable"
            ] == "YES"
            for x in assignments
        )
    ),
    (
        "high_similarity_near_full_candidate",
        tier_counts.get(
            "high_similarity_near_full_candidate",
            0
        )
    ),
    (
        "moderate_high_similarity_candidate",
        tier_counts.get(
            "moderate_high_similarity_candidate",
            0
        )
    ),
    (
        "broad_similarity_or_fragmentary",
        tier_counts.get(
            "broad_similarity_or_fragmentary",
            0
        )
    ),
    (
        "low_resolution",
        tier_counts.get(
            "low_resolution",
            0
        )
    ),
    (
        "unresolved_no_hit",
        tier_counts.get(
            "unresolved_no_hit",
            0
        )
    ),
    (
        "primary_loci_with_taxonomy_propagated",
        len(primary_taxonomy)
    ),
    (
        "exact_clusters_both_producers",
        sum(
            x["both_producers"]
            == "YES"
            for x in assignments
        )
    ),
]


with (
    OUT
    / "95C1_global_summary.tsv"
).open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    for k, v in metrics:
        fh.write(
            f"{k}\t{v}\n"
        )


# ============================================================
# 11. Methodological scope
# ============================================================

with (
    OUT
    / "95C1_methodological_scope.tsv"
).open("w") as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        "reference_database\t"
        "SILVA_138.2_SSURef_NR99_tax_silva\n"
    )

    fh.write(
        "queries\t148_orientation_insensitive_exact_representatives\n"
    )

    fh.write(
        "blast_task\tblastn\n"
    )

    fh.write(
        "max_target_sequences\t100\n"
    )

    fh.write(
        "near_best_definition\t"
        "bitscore_ge_99pct_best_and_identity_within_0.5pct_best\n"
    )

    fh.write(
        "taxonomy_consensus\t"
        "longest_common_prefix_of_parseable_SILVA_lineages_in_near_best_hits\n"
    )

    fh.write(
        "best_reference_species_equals_query_species\tNO\n"
    )

    fh.write(
        "high_similarity_near_full_equals_formal_species_assignment\tNO\n"
    )

    fh.write(
        "short_fragment_taxonomy_requires_conservative_interpretation\tYES\n"
    )

    fh.write(
        "exact_sequence_cluster_recurrence_equals_prevalence\tNO\n"
    )

    fh.write(
        "assembled_16S_count_equals_abundance\tNO\n"
    )


# ============================================================
# 12. Validation
# ============================================================

required = [
    "95C1_exact_representative_taxonomy.tsv",
    "95C1_primary_158_loci_taxonomy.tsv",
    "95C1_resolution_summary.tsv",
    "95C1_taxon_summary.tsv",
    "95C1_summary_by_coassembly.tsv",
    "95C1_unresolved_or_low_resolution.tsv",
    "95C1_global_summary.tsv",
    "95C1_methodological_scope.tsv",
]

for name in required:

    if not (
        OUT / name
    ).is_file():

        raise RuntimeError(
            f"Salida faltante: {name}"
        )


if len(assignments) != 148:
    raise RuntimeError(
        f"Assignments={len(assignments)}, expected=148"
    )

if len(primary_taxonomy) != 158:
    raise RuntimeError(
        f"Primary loci={len(primary_taxonomy)}, expected=158"
    )


print(
    f"QUERIES={len(assignments)}"
)

print(
    f"PRIMARY_LOCI={len(primary_taxonomy)}"
)

print(
    f"RAW_BLAST_HITS={len(blast_rows)}"
)

print(
    "95C1_PARSE=PASS"
)

