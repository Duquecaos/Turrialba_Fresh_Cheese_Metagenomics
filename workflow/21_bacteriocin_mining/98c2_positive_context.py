#!/usr/bin/env python3

import csv
import re
import sys
from pathlib import Path
from collections import defaultdict, Counter


if len(sys.argv) != 7:
    raise SystemExit(
        "Usage: 98c2_positive_context.py "
        "<positive.tsv> <gene_context.tsv> "
        "<bin_membership.tsv> <contig_catalog.tsv> "
        "<eggnog.tsv> <OUT>"
    )


POSITIVE = Path(sys.argv[1])
GENES = Path(sys.argv[2])
MEMBERSHIP = Path(sys.argv[3])
CONTIGCAT = Path(sys.argv[4])
EGGNOG = Path(sys.argv[5])
OUT = Path(sys.argv[6])

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# Helpers
# ============================================================

def int_or_zero(x):
    try:
        return int(float(x))
    except Exception:
        return 0


def float_or_none(x):
    try:
        return float(x)
    except Exception:
        return None


def uniq(values):
    return list(
        dict.fromkeys(
            x
            for x in values
            if x not in (
                "",
                None
            )
        )
    )


def direct_evidence_class(row):

    screen = row.get(
        "screen_class",
        ""
    )

    blast = row.get(
        "BLAST_classes",
        ""
    )

    exact = (
        "exact_full_length" in blast
        or
        "exact_query_contained" in blast
    )

    near = (
        "near_exact" in blast
    )


    if (
        screen
        == "sequence_homology_plus_GA_HMM"
    ):

        if exact:
            return (
                "A_exact_sequence_plus_GA_HMM"
            )

        if near:
            return (
                "A_near_exact_sequence_plus_GA_HMM"
            )

        return (
            "B_strong_sequence_homology_plus_GA_HMM"
        )


    if (
        screen
        == "sequence_homology_only"
    ):

        if exact:
            return (
                "A_exact_sequence_only"
            )

        if near:
            return (
                "A_near_exact_sequence_only"
            )

        return (
            "B_sequence_homology_only"
        )


    if (
        screen
        == "GA_HMM_profile_only"
    ):

        return (
            "C_GA_HMM_profile_only"
        )


    if (
        screen
        == "supplementary_noGA_HMM_only"
    ):

        return (
            "D_supplementary_noGA_HMM_only"
        )


    return (
        "UNRESOLVED"
    )


def is_primary(row):

    return int(
        row.get(
            "screen_class",
            ""
        )
        != "supplementary_noGA_HMM_only"
    )


# ============================================================
# Auxiliary context annotation categories
#
# These are descriptive flags only.
# ============================================================

RX_SPECIFIC = re.compile(
    r"bacteriocin|"
    r"lactococcin|"
    r"pediocin|"
    r"\bnisin\b|"
    r"lantibiotic|"
    r"lanthipept|"
    r"lassopept|"
    r"thiopept|"
    r"microviridin|"
    r"ranthipept|"
    r"sactipept|"
    r"glycocin|"
    r"\bripp\b|"
    r"ribosomally synthesized",
    re.I
)

RX_IMMUNITY = re.compile(
    r"\bimmunity\b|"
    r"self[- ]protection",
    re.I
)

RX_TRANSPORT = re.compile(
    r"\bABC\b|"
    r"ATP[- ]binding cassette|"
    r"\bexporter\b|"
    r"\bsecretion\b|"
    r"\btransporter\b",
    re.I
)

RX_PROCESS = re.compile(
    r"peptidase|"
    r"protease|"
    r"leader peptidase|"
    r"dehydratase|"
    r"cyclase|"
    r"radical SAM|"
    r"\bLanB\b|"
    r"\bLanC\b|"
    r"\bLanM\b|"
    r"\bRRE\b|"
    r"methyltransferase|"
    r"glycosyltransferase",
    re.I
)

RX_REGULATION = re.compile(
    r"histidine kinase|"
    r"response regulator|"
    r"two[- ]component|"
    r"transcriptional regulator",
    re.I
)


def annotation_flags(text):

    return {
        "specific_bacteriocin_RiPP":
            int(
                bool(
                    RX_SPECIFIC.search(
                        text
                    )
                )
            ),

        "immunity":
            int(
                bool(
                    RX_IMMUNITY.search(
                        text
                    )
                )
            ),

        "transport_export":
            int(
                bool(
                    RX_TRANSPORT.search(
                        text
                    )
                )
            ),

        "processing_modification":
            int(
                bool(
                    RX_PROCESS.search(
                        text
                    )
                )
            ),

        "regulation":
            int(
                bool(
                    RX_REGULATION.search(
                        text
                    )
                )
            ),
    }


# ============================================================
# 1. Positive proteins
# ============================================================

with POSITIVE.open(
    newline=""
) as fh:

    positive_rows = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )


if not positive_rows:
    raise RuntimeError(
        "Positive table is empty"
    )


positive_by_id = {
    r["protein_id"]:
        r
    for r in positive_rows
}


positive_ids = set(
    positive_by_id
)

positive_contigs = {
    r["contig_id"]
    for r in positive_rows
}


# Add operational evidence class.
for r in positive_rows:

    r[
        "direct_evidence_class"
    ] = direct_evidence_class(
        r
    )

    r[
        "primary_for_specialized_followup"
    ] = is_primary(
        r
    )


# ============================================================
# 2. Gene context
#
# Only retain genes from positive contigs.
# ============================================================

genes_by_contig = defaultdict(
    list
)

target_gene_meta = {}


with GENES.open(
    newline=""
) as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in reader:

        contig = row.get(
            "contig_id",
            ""
        )

        if contig not in positive_contigs:
            continue


        row[
            "_gene_number"
        ] = int_or_zero(
            row.get(
                "gene_number_on_contig",
                "0"
            )
        )


        row[
            "_start"
        ] = int_or_zero(
            row.get(
                "start",
                "0"
            )
        )


        row[
            "_end"
        ] = int_or_zero(
            row.get(
                "end",
                "0"
            )
        )


        genes_by_contig[
            contig
        ].append(
            row
        )


        gene_id = row.get(
            "gene_id",
            ""
        )

        if gene_id in positive_ids:

            target_gene_meta[
                gene_id
            ] = row


for contig in genes_by_contig:

    genes_by_contig[
        contig
    ].sort(
        key=lambda x:
            x[
                "_gene_number"
            ]
    )


missing_genes = sorted(
    positive_ids
    -
    set(
        target_gene_meta
    )
)


if missing_genes:

    (
        OUT
        / "98C2_missing_positive_gene_ids.txt"
    ).write_text(
        "\n".join(
            missing_genes
        )
        + "\n"
    )

    raise RuntimeError(
        f"{len(missing_genes)} positive gene IDs "
        "not found in gene_to_contig_context.tsv"
    )


# ============================================================
# 3. Unique contig metadata
# ============================================================

contig_meta = {}


with CONTIGCAT.open(
    newline=""
) as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in reader:

        cid = row.get(
            "unique_contig_id",
            ""
        )

        if cid in positive_contigs:

            contig_meta[
                cid
            ] = row


missing_contigs = sorted(
    positive_contigs
    -
    set(
        contig_meta
    )
)


if missing_contigs:

    (
        OUT
        / "98C2_missing_positive_contigs.txt"
    ).write_text(
        "\n".join(
            missing_contigs
        )
        + "\n"
    )

    raise RuntimeError(
        f"{len(missing_contigs)} positive contigs "
        "not found in unique contig catalog"
    )


# ============================================================
# 4. Raw-bin memberships
# ============================================================

memberships_by_contig = defaultdict(
    list
)


with MEMBERSHIP.open(
    newline=""
) as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in reader:

        cid = row.get(
            "unique_contig_id",
            ""
        )

        if cid in positive_contigs:

            memberships_by_contig[
                cid
            ].append(
                row
            )


membership_rows = []


for protein in positive_rows:

    cid = protein[
        "contig_id"
    ]

    for m in memberships_by_contig.get(
        cid,
        []
    ):

        membership_rows.append({
            "protein_id":
                protein[
                    "protein_id"
                ],

            "direct_evidence_class":
                protein[
                    "direct_evidence_class"
                ],

            "primary_for_specialized_followup":
                protein[
                    "primary_for_specialized_followup"
                ],

            **m,
        })


# ============================================================
# 5. Identify neighborhood genes
#
# +/- 5 genes from the positive gene.
# ============================================================

neighbor_gene_ids = set()

neighbor_raw = []


for protein in positive_rows:

    pid = protein[
        "protein_id"
    ]

    cid = protein[
        "contig_id"
    ]

    target = target_gene_meta[
        pid
    ]

    target_num = target[
        "_gene_number"
    ]

    genes = genes_by_contig[
        cid
    ]


    for g in genes:

        relative = (
            g[
                "_gene_number"
            ]
            -
            target_num
        )


        if abs(relative) > 5:
            continue


        neighbor_gene_ids.add(
            g[
                "gene_id"
            ]
        )


        target_mid = (
            target[
                "_start"
            ]
            +
            target[
                "_end"
            ]
        ) / 2.0


        gene_mid = (
            g[
                "_start"
            ]
            +
            g[
                "_end"
            ]
        ) / 2.0


        neighbor_raw.append({
            "target_protein_id":
                pid,

            "target_direct_evidence_class":
                protein[
                    "direct_evidence_class"
                ],

            "target_primary_followup":
                protein[
                    "primary_for_specialized_followup"
                ],

            "contig_id":
                cid,

            "neighbor_gene_id":
                g[
                    "gene_id"
                ],

            "relative_gene_order":
                relative,

            "center_distance_bp":
                round(
                    gene_mid
                    -
                    target_mid,
                    1
                ),

            "neighbor_gene_number":
                g.get(
                    "gene_number_on_contig",
                    ""
                ),

            "neighbor_start":
                g.get(
                    "start",
                    ""
                ),

            "neighbor_end":
                g.get(
                    "end",
                    ""
                ),

            "neighbor_strand":
                g.get(
                    "strand",
                    ""
                ),

            "neighbor_aa_length":
                g.get(
                    "aa_length",
                    ""
                ),

            "neighbor_partial":
                g.get(
                    "partial",
                    ""
                ),

            "is_target_gene":
                int(
                    g[
                        "gene_id"
                    ]
                    == pid
                ),
        })


# ============================================================
# 6. eggNOG annotation for neighborhood genes only
# ============================================================

eggnog = {}


with EGGNOG.open(
    newline=""
) as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    fieldnames = (
        reader.fieldnames
        or []
    )


    if not fieldnames:

        raise RuntimeError(
            "eggNOG table has no header"
        )


    query_field = fieldnames[0]


    for row in reader:

        q = row.get(
            query_field,
            ""
        )

        if q in neighbor_gene_ids:

            eggnog[
                q
            ] = row


# ============================================================
# 7. Final neighborhood table
# ============================================================

neighbor_rows = []


for row in neighbor_raw:

    gid = row[
        "neighbor_gene_id"
    ]

    e = eggnog.get(
        gid,
        {}
    )


    description = e.get(
        "Description",
        ""
    )

    preferred = e.get(
        "Preferred_name",
        ""
    )

    pfams = e.get(
        "PFAMs",
        ""
    )


    annotation_text = " | ".join(
        x
        for x in [
            description,
            preferred,
            pfams,
        ]
        if x
    )


    flags = annotation_flags(
        annotation_text
    )


    neighbor_rows.append({
        **row,

        "eggnog_annotated":
            int(
                gid in eggnog
            ),

        "eggnog_COG_category":
            e.get(
                "COG_category",
                ""
            ),

        "eggnog_Preferred_name":
            preferred,

        "eggnog_Description":
            description,

        "eggnog_PFAMs":
            pfams,

        **flags,
    })


# ============================================================
# 8. Positive-protein context summary
# ============================================================

neighbors_by_target = defaultdict(
    list
)


for row in neighbor_rows:

    neighbors_by_target[
        row[
            "target_protein_id"
        ]
    ].append(
        row
    )


context_rows = []


for p in positive_rows:

    pid = p[
        "protein_id"
    ]

    cid = p[
        "contig_id"
    ]

    target = target_gene_meta[
        pid
    ]

    cm = contig_meta[
        cid
    ]

    nn = neighbors_by_target[
        pid
    ]


    nonself = [
        x
        for x in nn
        if x[
            "is_target_gene"
        ] == 0
    ]


    all_genes = genes_by_contig[
        cid
    ]

    target_index = next(
        i
        for i, g
        in enumerate(
            all_genes
        )
        if g[
            "gene_id"
        ] == pid
    )


    left_context_truncated = int(
        target_index < 5
    )

    right_context_truncated = int(
        (
            len(
                all_genes
            )
            - 1
            - target_index
        ) < 5
    )


    memberships = memberships_by_contig.get(
        cid,
        []
    )


    raw_bins = uniq([
        x.get(
            "raw_bin_uid",
            ""
        )
        for x in memberships
    ])


    original_headers = uniq([
        x.get(
            "original_contig_header",
            ""
        )
        for x in memberships
    ])


    complete_values = [
        float_or_none(
            x.get(
                "bin_completeness",
                ""
            )
        )
        for x in memberships
    ]

    complete_values = [
        x
        for x in complete_values
        if x is not None
    ]


    contam_values = [
        float_or_none(
            x.get(
                "bin_contamination",
                ""
            )
        )
        for x in memberships
    ]

    contam_values = [
        x
        for x in contam_values
        if x is not None
    ]


    represented = int_or_zero(
        cm.get(
            "represented_in_final18",
            "0"
        )
    )


    additional_scope = (
        "additional_incomplete_context_not_in_final18"
        if represented == 0
        else
        "sequence_already_represented_in_final18"
    )


    specific_count = sum(
        x[
            "specific_bacteriocin_RiPP"
        ]
        for x in nonself
    )

    immunity_count = sum(
        x[
            "immunity"
        ]
        for x in nonself
    )

    transport_count = sum(
        x[
            "transport_export"
        ]
        for x in nonself
    )

    process_count = sum(
        x[
            "processing_modification"
        ]
        for x in nonself
    )

    regulation_count = sum(
        x[
            "regulation"
        ]
        for x in nonself
    )


    if (
        specific_count > 0
        or immunity_count > 0
    ):

        context_flag = (
            "specific_RiPP_or_immunity_neighbor"
        )

    elif (
        transport_count > 0
        and process_count > 0
    ):

        context_flag = (
            "generic_transport_plus_processing_context"
        )

    elif (
        transport_count > 0
        or process_count > 0
        or regulation_count > 0
    ):

        context_flag = (
            "generic_supportive_context_feature"
        )

    else:

        context_flag = (
            "no_auxiliary_context_feature_detected"
        )


    context_rows.append({
        **p,

        "gene_number_on_contig":
            target.get(
                "gene_number_on_contig",
                ""
            ),

        "gene_start":
            target.get(
                "start",
                ""
            ),

        "gene_end":
            target.get(
                "end",
                ""
            ),

        "gene_strand":
            target.get(
                "strand",
                ""
            ),

        "contig_length_bp":
            cm.get(
                "length_bp",
                target.get(
                    "contig_length_bp",
                    ""
                )
            ),

        "contig_gc_pct":
            cm.get(
                "gc_pct",
                ""
            ),

        "contig_coassemblies":
            cm.get(
                "coassemblies",
                ""
            ),

        "contig_producers":
            cm.get(
                "producers",
                ""
            ),

        "contig_binners":
            cm.get(
                "binners",
                ""
            ),

        "contig_quality_groups":
            cm.get(
                "quality_groups",
                ""
            ),

        "strongest_incomplete_context":
            cm.get(
                "strongest_incomplete_context",
                ""
            ),

        "represented_in_final18":
            represented,

        "final18_MAGs":
            cm.get(
                "final18_MAGs",
                ""
            ),

        "additional_scope":
            additional_scope,

        "n_raw_bin_memberships":
            len(
                raw_bins
            ),

        "raw_bin_uids":
            ";".join(
                raw_bins
            ),

        "original_contig_headers":
            ";".join(
                original_headers
            ),

        "min_member_bin_completeness":
            min(
                complete_values
            )
            if complete_values
            else "",

        "max_member_bin_completeness":
            max(
                complete_values
            )
            if complete_values
            else "",

        "min_member_bin_contamination":
            min(
                contam_values
            )
            if contam_values
            else "",

        "max_member_bin_contamination":
            max(
                contam_values
            )
            if contam_values
            else "",

        "neighbor_genes_examined_excluding_target":
            len(
                nonself
            ),

        "left_context_truncated_by_contig":
            left_context_truncated,

        "right_context_truncated_by_contig":
            right_context_truncated,

        "specific_RiPP_neighbor_count":
            specific_count,

        "immunity_neighbor_count":
            immunity_count,

        "transport_export_neighbor_count":
            transport_count,

        "processing_modification_neighbor_count":
            process_count,

        "regulation_neighbor_count":
            regulation_count,

        "auxiliary_context_flag":
            context_flag,
    })


# ============================================================
# 9. Contig-level primary-target summary
# ============================================================

positives_by_contig = defaultdict(
    list
)


for row in context_rows:

    positives_by_contig[
        row[
            "contig_id"
        ]
    ].append(
        row
    )


primary_contig_rows = []
exploratory_contig_rows = []
multihit_rows = []


for cid, rows in sorted(
    positives_by_contig.items()
):

    primary = [
        r
        for r in rows
        if str(
            r[
                "primary_for_specialized_followup"
            ]
        ) == "1"
    ]


    cm = contig_meta[
        cid
    ]


    base = {
        "contig_id":
            cid,

        "positive_protein_count":
            len(
                rows
            ),

        "primary_positive_protein_count":
            len(
                primary
            ),

        "positive_protein_ids":
            ";".join(
                r[
                    "protein_id"
                ]
                for r in rows
            ),

        "primary_protein_ids":
            ";".join(
                r[
                    "protein_id"
                ]
                for r in primary
            ),

        "direct_evidence_classes":
            ";".join(
                uniq([
                    r[
                        "direct_evidence_class"
                    ]
                    for r in rows
                ])
            ),

        "best_query_candidates":
            ";".join(
                uniq([
                    r.get(
                        "best_query_candidate",
                        ""
                    )
                    for r in rows
                ])
            ),

        "GA_HMM_models":
            ";".join(
                uniq([
                    r.get(
                        "GA_HMM_models",
                        ""
                    )
                    for r in rows
                ])
            ),

        "length_bp":
            cm.get(
                "length_bp",
                ""
            ),

        "coassemblies":
            cm.get(
                "coassemblies",
                ""
            ),

        "producers":
            cm.get(
                "producers",
                ""
            ),

        "n_bin_memberships":
            cm.get(
                "n_bin_memberships",
                ""
            ),

        "strongest_incomplete_context":
            cm.get(
                "strongest_incomplete_context",
                ""
            ),

        "represented_in_final18":
            cm.get(
                "represented_in_final18",
                ""
            ),

        "final18_MAGs":
            cm.get(
                "final18_MAGs",
                ""
            ),

        "member_raw_bin_uids":
            cm.get(
                "member_raw_bin_uids",
                ""
            ),
    }


    if primary:

        primary_contig_rows.append(
            base
        )

    else:

        exploratory_contig_rows.append(
            base
        )


    if len(primary) >= 2:

        multihit_rows.append(
            base
        )


# ============================================================
# 10. Writers
# ============================================================

def write_table(
    path,
    rows,
    fieldnames=None
):

    with path.open(
        "w",
        newline=""
    ) as fh:

        if rows:

            fields = (
                fieldnames
                if fieldnames
                else list(
                    rows[0].keys()
                )
            )

        else:

            fields = (
                fieldnames
                if fieldnames
                else ["status"]
            )


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


write_table(
    OUT
    / "98C2_positive_context_summary.tsv",
    context_rows
)


write_table(
    OUT
    / "98C2_neighbor_genes.tsv",
    neighbor_rows
)


membership_fields = (
    list(
        membership_rows[0].keys()
    )
    if membership_rows
    else [
        "protein_id",
        "direct_evidence_class",
        "primary_for_specialized_followup",
        "raw_bin_uid",
        "unique_contig_id",
    ]
)


write_table(
    OUT
    / "98C2_positive_bin_memberships.tsv",
    membership_rows,
    membership_fields
)


contig_fields = [
    "contig_id",
    "positive_protein_count",
    "primary_positive_protein_count",
    "positive_protein_ids",
    "primary_protein_ids",
    "direct_evidence_classes",
    "best_query_candidates",
    "GA_HMM_models",
    "length_bp",
    "coassemblies",
    "producers",
    "n_bin_memberships",
    "strongest_incomplete_context",
    "represented_in_final18",
    "final18_MAGs",
    "member_raw_bin_uids",
]


write_table(
    OUT
    / "98C2_primary_target_contigs.tsv",
    primary_contig_rows,
    contig_fields
)


write_table(
    OUT
    / "98C2_exploratory_noGA_contigs.tsv",
    exploratory_contig_rows,
    contig_fields
)


write_table(
    OUT
    / "98C2_multihit_primary_contigs.tsv",
    multihit_rows,
    contig_fields
)


# ============================================================
# 11. Global summary
# ============================================================

primary_proteins = [
    r
    for r in context_rows
    if str(
        r[
            "primary_for_specialized_followup"
        ]
    ) == "1"
]


additional_primary_contigs = {
    r[
        "contig_id"
    ]
    for r in primary_proteins
    if str(
        r[
            "represented_in_final18"
        ]
    ) == "0"
}


represented_primary_contigs = {
    r[
        "contig_id"
    ]
    for r in primary_proteins
    if str(
        r[
            "represented_in_final18"
        ]
    ) == "1"
}


evidence_counter = Counter(
    r[
        "direct_evidence_class"
    ]
    for r in context_rows
)


context_counter = Counter(
    r[
        "auxiliary_context_flag"
    ]
    for r in context_rows
)


with (
    OUT
    / "98C2_global_summary.tsv"
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
            "positive_proteins_total",
            len(
                context_rows
            )
        ),
        (
            "positive_contigs_total",
            len(
                positives_by_contig
            )
        ),
        (
            "primary_positive_proteins",
            len(
                primary_proteins
            )
        ),
        (
            "primary_target_contigs",
            len(
                primary_contig_rows
            )
        ),
        (
            "exploratory_noGA_only_contigs",
            len(
                exploratory_contig_rows
            )
        ),
        (
            "multihit_primary_contigs",
            len(
                multihit_rows
            )
        ),
        (
            "primary_contigs_additional_to_final18",
            len(
                additional_primary_contigs
            )
        ),
        (
            "primary_contigs_already_represented_in_final18",
            len(
                represented_primary_contigs
            )
        ),
        (
            "positive_proteins_with_specific_RiPP_or_immunity_neighbor",
            sum(
                r[
                    "auxiliary_context_flag"
                ]
                ==
                "specific_RiPP_or_immunity_neighbor"
                for r in context_rows
            )
        ),
    ]


    for key, value in sorted(
        evidence_counter.items()
    ):

        metrics.append(
            (
                f"evidence_class__{key}",
                value
            )
        )


    for key, value in sorted(
        context_counter.items()
    ):

        metrics.append(
            (
                f"context_flag__{key}",
                value
            )
        )


    metrics += [
        (
            "specialized_prediction_run",
            "NO"
        ),
        (
            "reads_remapped",
            "NO"
        ),
        (
            "bacteriocin_function_confirmed",
            "NO"
        ),
        (
            "amplicon_data_used",
            "NO"
        ),
        (
            "next_step",
            "98C3_targeted_specialized_prediction_on_primary_positive_contexts"
        ),
    ]


    w.writerows(
        metrics
    )


print(
    f"POSITIVE_PROTEINS={len(context_rows)}"
)

print(
    f"PRIMARY_PROTEINS={len(primary_proteins)}"
)

print(
    f"PRIMARY_CONTIGS={len(primary_contig_rows)}"
)

print(
    f"MULTIHIT_CONTIGS={len(multihit_rows)}"
)

print(
    "98C2=PASS"
)
