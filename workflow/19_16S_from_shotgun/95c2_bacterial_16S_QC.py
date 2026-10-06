#!/usr/bin/env python3

import csv
import sys
from collections import defaultdict, Counter
from pathlib import Path


ROOT = Path(sys.argv[1])

C1 = (
    ROOT
    / "97_16S_shotgun"
    / "95C1_SILVA_taxonomy"
)

B2 = (
    ROOT
    / "97_16S_shotgun"
    / "95B2_consolidated"
)

OUT = (
    ROOT
    / "97_16S_shotgun"
    / "95C2_bacterial_QC"
)

REP_TAX = (
    C1
    / "95C1_exact_representative_taxonomy.tsv"
)

PRIMARY = (
    C1
    / "95C1_primary_158_loci_taxonomy.tsv"
)

REP_FASTA = (
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


def number(value):

    try:
        return float(value)
    except Exception:
        return None


def integer(value):

    return int(float(value))


def lineage_tokens(lineage):

    return [
        x.strip()
        for x in (
            lineage or ""
        ).split(";")
        if x.strip()
    ]


def lineage_context(lineage):

    tokens = lineage_tokens(
        lineage
    )

    lower = [
        x.lower()
        for x in tokens
    ]

    if any(
        x == "eukaryota"
        for x in lower
    ):
        return (
            "eukaryotic_or_organelle"
        )

    if any(
        (
            "mitochondria" in x
            or "chloroplast" in x
            or "plastid" in x
        )
        for x in lower
    ):
        return (
            "eukaryotic_or_organelle"
        )

    if (
        tokens
        and tokens[0].lower()
        == "bacteria"
    ):
        return (
            "bacterial_nonorganelle"
        )

    if (
        tokens
        and tokens[0].lower()
        == "archaea"
    ):
        return (
            "archaeal_nonorganelle"
        )

    return "unresolved_lineage"


def is_placeholder(token):

    x = token.lower().strip()

    bad_substrings = [
        "uncultured",
        "unidentified",
        "incertae sedis",
        "environmental sample",
        "metagenome",
        "unknown",
    ]

    return any(
        bad in x
        for bad in bad_substrings
    )


def deepest_informative_taxon(
    lineage
):

    tokens = lineage_tokens(
        lineage
    )

    for token in reversed(tokens):

        if is_placeholder(
            token
        ):
            continue

        if token.lower() in {
            "bacteria",
            "archaea",
            "eukaryota",
        }:
            continue

        return token

    if tokens:
        return tokens[0]

    return ""


def lca_lineages(lineages):

    parsed = [
        lineage_tokens(x)
        for x in lineages
        if x
    ]

    if not parsed:
        return ""

    n = min(
        len(x)
        for x in parsed
    )

    common = []

    for i in range(n):

        vals = {
            x[i]
            for x in parsed
        }

        if len(vals) != 1:
            break

        common.append(
            parsed[0][i]
        )

    return ";".join(common)


def split_semicolon_values(
    values
):

    out = set()

    for value in values:

        for x in (
            value or ""
        ).split(";"):

            x = x.strip()

            if x:
                out.add(x)

    return sorted(out)


# ============================================================
# Load
# ============================================================

rep_rows, rep_fields = read_tsv(
    REP_TAX
)

primary_rows, primary_fields = read_tsv(
    PRIMARY
)

sequences = read_fasta(
    REP_FASTA
)


if len(rep_rows) != 148:
    raise RuntimeError(
        f"Representatives expected=148 "
        f"observed={len(rep_rows)}"
    )

if len(primary_rows) != 158:
    raise RuntimeError(
        f"Primary loci expected=158 "
        f"observed={len(primary_rows)}"
    )

if len(sequences) != 148:
    raise RuntimeError(
        f"FASTA sequences expected=148 "
        f"observed={len(sequences)}"
    )


# ============================================================
# QC exact representatives
# ============================================================

qc = []


for row in rep_rows:

    cluster = (
        row["exact_cluster_id"]
    )

    if cluster not in sequences:
        raise RuntimeError(
            f"FASTA missing {cluster}"
        )

    qlen = integer(
        row["query_length_bp"]
    )

    ident = number(
        row[
            "best_percent_identity"
        ]
    )

    qcov = number(
        row[
            "best_query_coverage_percent"
        ]
    )

    lineage = (
        row[
            "SILVA_near_best_LCA_lineage"
        ]
        or
        row[
            "best_subject_SILVA_lineage"
        ]
    )

    context = lineage_context(
        lineage
    )

    informative = (
        deepest_informative_taxon(
            lineage
        )
    )


    # --------------------------------------------------------
    # FINAL consistency with 95B2:
    # near-full descriptive threshold = >=1300 bp
    # --------------------------------------------------------

    if context != "bacterial_nonorganelle":

        tier = (
            "EXCLUDED_nonbacterial_or_organelle"
        )

    elif (
        ident is not None
        and qcov is not None
        and qlen >= 1300
        and ident >= 98.7
        and qcov >= 90
    ):

        tier = (
            "A_high_similarity_ge1300"
        )

    elif (
        ident is not None
        and qcov is not None
        and qlen >= 800
        and ident >= 94.5
        and qcov >= 90
    ):

        tier = (
            "B_moderate_similarity"
        )

    elif (
        ident is not None
        and qcov is not None
        and ident >= 90
        and qcov >= 80
    ):

        tier = (
            "C_fragmentary_or_broad"
        )

    else:

        tier = (
            "D_low_resolution"
        )


    taxonomy_reporting = (
        context
        == "bacterial_nonorganelle"
        and tier
        in {
            "A_high_similarity_ge1300",
            "B_moderate_similarity",
            "C_fragmentary_or_broad",
        }
    )

    quantification_target = (
        taxonomy_reporting
    )


    rec = dict(row)

    rec.update({
        "final_lineage_context":
            context,

        "deepest_informative_LCA_taxon":
            informative,

        "final_resolution_tier":
            tier,

        "taxonomy_reporting_eligible":
            "YES"
            if taxonomy_reporting
            else "NO",

        "bacterial_quantification_target":
            "YES"
            if quantification_target
            else "NO",

        "final_near_full_definition":
            "length_ge1300bp",

        "interpretation_guardrail_final":
            "sequence_similarity_not_formal_species_identification",
    })

    qc.append(rec)


qc_fields = list(
    qc[0].keys()
)

write_tsv(
    OUT
    / "95C2_exact_representatives_QC.tsv",
    qc,
    qc_fields
)


# ============================================================
# FASTA subsets
# ============================================================

quant_ids = [
    x["exact_cluster_id"]
    for x in qc
    if x[
        "bacterial_quantification_target"
    ] == "YES"
]

high_ids = [
    x["exact_cluster_id"]
    for x in qc
    if x[
        "final_resolution_tier"
    ]
    == "A_high_similarity_ge1300"
]

excluded_ids = [
    x["exact_cluster_id"]
    for x in qc
    if x[
        "final_lineage_context"
    ]
    == "eukaryotic_or_organelle"
]

low_bacterial_ids = [
    x["exact_cluster_id"]
    for x in qc
    if (
        x[
            "final_lineage_context"
        ]
        == "bacterial_nonorganelle"
        and
        x[
            "final_resolution_tier"
        ]
        == "D_low_resolution"
    )
]


write_fasta(
    OUT
    / "95C2_bacterial_quantification_targets.fasta",
    quant_ids,
    sequences
)

write_fasta(
    OUT
    / "95C2_high_confidence_ge1300.fasta",
    high_ids,
    sequences
)

write_fasta(
    OUT
    / "95C2_excluded_eukaryotic_or_organelle.fasta",
    excluded_ids,
    sequences
)

write_fasta(
    OUT
    / "95C2_low_resolution_bacterial.fasta",
    low_bacterial_ids,
    sequences
)


# ============================================================
# Propagate QC to 158 primary loci
# ============================================================

qc_by_cluster = {
    x["exact_cluster_id"]: x
    for x in qc
}

primary_qc = []


for row in primary_rows:

    cluster = (
        row[
            "taxonomy_exact_cluster_id"
        ]
    )

    if cluster not in qc_by_cluster:
        raise RuntimeError(
            f"QC cluster missing: {cluster}"
        )

    q = qc_by_cluster[
        cluster
    ]

    rec = dict(row)

    rec.update({
        "final_lineage_context":
            q[
                "final_lineage_context"
            ],

        "deepest_informative_LCA_taxon":
            q[
                "deepest_informative_LCA_taxon"
            ],

        "final_resolution_tier":
            q[
                "final_resolution_tier"
            ],

        "taxonomy_reporting_eligible":
            q[
                "taxonomy_reporting_eligible"
            ],

        "bacterial_quantification_target":
            q[
                "bacterial_quantification_target"
            ],
    })

    primary_qc.append(
        rec
    )


primary_qc_fields = list(
    primary_qc[0].keys()
)

write_tsv(
    OUT
    / "95C2_primary_158_loci_QC.tsv",
    primary_qc,
    primary_qc_fields
)


# ============================================================
# Exclusion / low-resolution tables
# ============================================================

excluded = [
    x
    for x in qc
    if x[
        "final_lineage_context"
    ]
    == "eukaryotic_or_organelle"
]

low_bacterial = [
    x
    for x in qc
    if (
        x[
            "final_lineage_context"
        ]
        == "bacterial_nonorganelle"
        and
        x[
            "final_resolution_tier"
        ]
        == "D_low_resolution"
    )
]

reportable = [
    x
    for x in qc
    if x[
        "taxonomy_reporting_eligible"
    ] == "YES"
]


write_tsv(
    OUT
    / "95C2_excluded_eukaryotic_or_organelle.tsv",
    excluded,
    qc_fields
)

write_tsv(
    OUT
    / "95C2_low_resolution_bacterial.tsv",
    low_bacterial,
    qc_fields
)

write_tsv(
    OUT
    / "95C2_taxonomically_reportable_bacterial.tsv",
    reportable,
    qc_fields
)


# ============================================================
# Taxon summary, cleaned
# ============================================================

tax_groups = defaultdict(list)


for row in reportable:

    taxon = (
        row[
            "deepest_informative_LCA_taxon"
        ]
        or "UNRESOLVED"
    )

    tax_groups[
        taxon
    ].append(row)


tax_summary = []


for taxon, rows in tax_groups.items():

    coassemblies = (
        split_semicolon_values(
            x["coassemblies"]
            for x in rows
        )
    )

    producers = (
        split_semicolon_values(
            x["producers"]
            for x in rows
        )
    )

    tax_summary.append({
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
                    "final_resolution_tier"
                ]
                ==
                "A_high_similarity_ge1300"
                for x in rows
            ),

        "B_moderate_similarity":
            sum(
                x[
                    "final_resolution_tier"
                ]
                ==
                "B_moderate_similarity"
                for x in rows
            ),

        "C_fragmentary_or_broad":
            sum(
                x[
                    "final_resolution_tier"
                ]
                ==
                "C_fragmentary_or_broad"
                for x in rows
            ),

        "interpretation":
            "assembled_exact_sequence_clusters_not_abundance",
    })


tax_summary.sort(
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
    / "95C2_clean_bacterial_taxon_summary.tsv",
    tax_summary,
    list(
        tax_summary[0].keys()
    )
)


# ============================================================
# Summary by coassembly
# ============================================================

coassembly_summary = []


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
        for x in primary_qc
        if x["coassembly"] == group
    ]

    coassembly_summary.append({
        "coassembly":
            group,

        "primary_16S_loci":
            len(rr),

        "bacterial_supported_loci":
            sum(
                x[
                    "bacterial_quantification_target"
                ]
                == "YES"
                for x in rr
            ),

        "eukaryotic_or_organelle_loci":
            sum(
                x[
                    "final_lineage_context"
                ]
                ==
                "eukaryotic_or_organelle"
                for x in rr
            ),

        "bacterial_low_resolution_loci":
            sum(
                (
                    x[
                        "final_lineage_context"
                    ]
                    ==
                    "bacterial_nonorganelle"
                    and
                    x[
                        "final_resolution_tier"
                    ]
                    ==
                    "D_low_resolution"
                )
                for x in rr
            ),

        "A_high_similarity_ge1300":
            sum(
                x[
                    "final_resolution_tier"
                ]
                ==
                "A_high_similarity_ge1300"
                for x in rr
            ),

        "B_moderate_similarity":
            sum(
                x[
                    "final_resolution_tier"
                ]
                ==
                "B_moderate_similarity"
                for x in rr
            ),

        "C_fragmentary_or_broad":
            sum(
                x[
                    "final_resolution_tier"
                ]
                ==
                "C_fragmentary_or_broad"
                for x in rr
            ),
    })


write_tsv(
    OUT
    / "95C2_summary_by_coassembly.tsv",
    coassembly_summary,
    list(
        coassembly_summary[0].keys()
    )
)


# ============================================================
# MAG-associated loci + cross-locus consistency
# ============================================================

mag_loci = defaultdict(list)


for row in primary_qc:

    mags = [
        x
        for x in (
            row.get(
                "final18_MAGs_exact",
                ""
            )
        ).split(";")
        if x
    ]

    for mag in mags:
        mag_loci[mag].append(
            row
        )


mag_summary = []


for mag, rows in sorted(
    mag_loci.items()
):

    lineages = [
        x["SILVA_LCA_lineage"]
        for x in rows
        if x["SILVA_LCA_lineage"]
    ]

    cross_lca = (
        lca_lineages(
            lineages
        )
    )

    depth = len(
        lineage_tokens(
            cross_lca
        )
    )

    taxa = sorted(
        {
            x[
                "deepest_informative_LCA_taxon"
            ]
            for x in rows
            if x[
                "deepest_informative_LCA_taxon"
            ]
        }
    )

    if len(rows) == 1:

        flag = (
            "single_16S_locus"
        )

    elif depth >= 6:

        flag = (
            "cross_locus_consistent_genus_or_deeper"
        )

    elif depth >= 5:

        flag = (
            "cross_locus_consistent_family_level"
        )

    else:

        flag = (
            "discordant_high_level_16S_lineages"
        )

    mag_summary.append({
        "MAG":
            mag,

        "n_16S_loci":
            len(rows),

        "16S_locus_ids":
            ";".join(
                x["locus_id"]
                for x in rows
            ),

        "informative_taxa":
            ";".join(taxa),

        "cross_locus_LCA":
            cross_lca,

        "cross_locus_LCA_depth":
            depth,

        "16S_consistency_flag":
            flag,

        "A_high_similarity_ge1300":
            sum(
                x[
                    "final_resolution_tier"
                ]
                ==
                "A_high_similarity_ge1300"
                for x in rows
            ),

        "bacterial_supported_for_quantification":
            sum(
                x[
                    "bacterial_quantification_target"
                ]
                == "YES"
                for x in rows
            ),

        "interpretation":
            "discordance_is_QC_signal_not_proof_of_contamination",
    })


write_tsv(
    OUT
    / "95C2_MAG_16S_consistency.tsv",
    mag_summary,
    list(
        mag_summary[0].keys()
    )
)


mag_detail = [
    x
    for x in primary_qc
    if x.get(
        "exact_final18_contig_membership",
        ""
    ) == "YES"
]

write_tsv(
    OUT
    / "95C2_MAG_associated_16S_detail.tsv",
    mag_detail,
    primary_qc_fields
)


# ============================================================
# Global summary
# ============================================================

context_counts = Counter(
    x[
        "final_lineage_context"
    ]
    for x in qc
)

tier_counts = Counter(
    x[
        "final_resolution_tier"
    ]
    for x in qc
)


metrics = [
    (
        "exact_sequence_clusters_total",
        len(qc)
    ),
    (
        "bacterial_nonorganelle_clusters",
        context_counts.get(
            "bacterial_nonorganelle",
            0
        )
    ),
    (
        "eukaryotic_or_organelle_clusters",
        context_counts.get(
            "eukaryotic_or_organelle",
            0
        )
    ),
    (
        "archaeal_nonorganelle_clusters",
        context_counts.get(
            "archaeal_nonorganelle",
            0
        )
    ),
    (
        "unresolved_lineage_clusters",
        context_counts.get(
            "unresolved_lineage",
            0
        )
    ),
    (
        "A_high_similarity_ge1300",
        tier_counts.get(
            "A_high_similarity_ge1300",
            0
        )
    ),
    (
        "B_moderate_similarity",
        tier_counts.get(
            "B_moderate_similarity",
            0
        )
    ),
    (
        "C_fragmentary_or_broad",
        tier_counts.get(
            "C_fragmentary_or_broad",
            0
        )
    ),
    (
        "D_low_resolution",
        tier_counts.get(
            "D_low_resolution",
            0
        )
    ),
    (
        "excluded_nonbacterial_or_organelle",
        tier_counts.get(
            "EXCLUDED_nonbacterial_or_organelle",
            0
        )
    ),
    (
        "bacterial_quantification_target_clusters",
        len(quant_ids)
    ),
    (
        "primary_16S_loci_total",
        len(primary_qc)
    ),
    (
        "primary_bacterial_quantification_target_loci",
        sum(
            x[
                "bacterial_quantification_target"
            ]
            == "YES"
            for x in primary_qc
        )
    ),
    (
        "primary_eukaryotic_or_organelle_loci",
        sum(
            x[
                "final_lineage_context"
            ]
            ==
            "eukaryotic_or_organelle"
            for x in primary_qc
        )
    ),
    (
        "MAG_associated_primary_loci",
        len(mag_detail)
    ),
    (
        "MAGs_with_16S",
        len(mag_summary)
    ),
    (
        "MAGs_with_discordant_high_level_16S",
        sum(
            x[
                "16S_consistency_flag"
            ]
            ==
            "discordant_high_level_16S_lineages"
            for x in mag_summary
        )
    ),
]


with (
    OUT
    / "95C2_global_summary.tsv"
).open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    for key, value in metrics:

        fh.write(
            f"{key}\t{value}\n"
        )


# ============================================================
# Methodological scope
# ============================================================

with (
    OUT
    / "95C2_methodological_scope.tsv"
).open("w") as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        "analysis_scope\tshotgun_bacterial_16S_only\n"
    )

    fh.write(
        "reference\tSILVA_138.2_SSURef_NR99\n"
    )

    fh.write(
        "nonbacterial_filter\t"
        "Eukaryota_or_lineage_containing_Mitochondria_Chloroplast_or_plastid\n"
    )

    fh.write(
        "final_near_full_threshold_bp\t1300\n"
    )

    fh.write(
        "A_threshold\t"
        "length_ge1300_identity_ge98.7_querycov_ge90\n"
    )

    fh.write(
        "B_threshold\t"
        "length_ge800_identity_ge94.5_querycov_ge90\n"
    )

    fh.write(
        "C_threshold\t"
        "identity_ge90_querycov_ge80\n"
    )

    fh.write(
        "D_low_resolution_used_for_primary_quantification\tNO\n"
    )

    fh.write(
        "organellar_sequences_used_for_bacterial_quantification\tNO\n"
    )

    fh.write(
        "deepest_taxon_source\tnear_best_SILVA_LCA\n"
    )

    fh.write(
        "species_level_label_equals_formal_species_identification\tNO\n"
    )

    fh.write(
        "MAG_16S_discordance_equals_contamination_proof\tNO\n"
    )

    fh.write(
        "exact_sequence_cluster_count_equals_abundance\tNO\n"
    )


# ============================================================
# Validation
# ============================================================

required = [
    "95C2_exact_representatives_QC.tsv",
    "95C2_primary_158_loci_QC.tsv",
    "95C2_bacterial_quantification_targets.fasta",
    "95C2_high_confidence_ge1300.fasta",
    "95C2_excluded_eukaryotic_or_organelle.fasta",
    "95C2_low_resolution_bacterial.fasta",
    "95C2_excluded_eukaryotic_or_organelle.tsv",
    "95C2_low_resolution_bacterial.tsv",
    "95C2_taxonomically_reportable_bacterial.tsv",
    "95C2_clean_bacterial_taxon_summary.tsv",
    "95C2_summary_by_coassembly.tsv",
    "95C2_MAG_16S_consistency.tsv",
    "95C2_MAG_associated_16S_detail.tsv",
    "95C2_global_summary.tsv",
    "95C2_methodological_scope.tsv",
]


for name in required:

    if not (
        OUT / name
    ).is_file():

        raise RuntimeError(
            f"Salida faltante: {name}"
        )


if len(qc) != 148:
    raise RuntimeError(
        "QC representatives != 148"
    )

if len(primary_qc) != 158:
    raise RuntimeError(
        "Primary QC != 158"
    )


print(
    f"REP_CLUSTERS={len(qc)}"
)

print(
    f"BACTERIAL_QUANT_TARGETS={len(quant_ids)}"
)

print(
    f"EXCLUDED_ORGANELLAR={len(excluded_ids)}"
)

print(
    f"LOW_RESOLUTION_BACTERIAL={len(low_bacterial_ids)}"
)

print(
    f"MAGS_WITH_16S={len(mag_summary)}"
)

print(
    "95C2_QC=PASS"
)

