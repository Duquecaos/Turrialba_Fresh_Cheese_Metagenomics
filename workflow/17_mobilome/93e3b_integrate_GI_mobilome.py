#!/usr/bin/env python3

import csv
import sys
from collections import defaultdict, Counter
from pathlib import Path


ROOT = Path(sys.argv[1])

GI_FILE = (
    ROOT
    / "95_mobilome/93E2C_alienhunter_regions"
    / "93E2C_candidate_GI_regions.tsv"
)

AMR_FILE = (
    ROOT
    / "95_mobilome/93D4_IS_ARG_integron_plasmid"
    / "93D4_AMR_nearest_IS.tsv"
)

IS_FILE = (
    ROOT
    / "95_mobilome/93D4_IS_ARG_integron_plasmid"
    / "93D4_IS_integrated.tsv"
)

INTEGRON_FILE = (
    ROOT
    / "95_mobilome/93D2_integron_ARG_plasmid"
    / "93D2_integron_structures_integrated.tsv"
)

PLASMID_FILE = (
    ROOT
    / "95_mobilome/93C1_plasmid_cargo"
    / "93C1_plasmid_candidate_integrated.tsv"
)

CTX_FILE = (
    ROOT
    / "95_mobilome/93D5_mobile_contexts"
    / "93D5_local_contexts.tsv"
)

OUT = (
    ROOT
    / "95_mobilome/93E3B_GI_mobile_integration"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ======================================================================
# Helpers
# ======================================================================

def read_tsv(path):

    if not path.is_file():
        raise RuntimeError(
            f"Archivo faltante: {path}"
        )

    with path.open() as fh:
        return list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


def require_columns(rows, required, label):

    if not rows:
        raise RuntimeError(
            f"{label}: tabla vacía"
        )

    cols = set(rows[0].keys())

    missing = sorted(
        set(required) - cols
    )

    if missing:
        raise RuntimeError(
            f"{label}: faltan columnas: "
            + ",".join(missing)
        )


def integer(value, label="value"):

    try:
        return int(float(value))
    except Exception:
        raise RuntimeError(
            f"No pude convertir {label}='{value}' a entero"
        )


def numeric(value):

    try:
        return float(value)
    except Exception:
        return None


def yes(value):

    return str(value).strip().upper() in {
        "YES",
        "TRUE",
        "1",
        "Y",
    }


def overlap(a1, a2, b1, b2):

    return (
        a1 <= b2
        and b1 <= a2
    )


def overlap_length(a1, a2, b1, b2):

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


def relation_to_GI(
    gi_start,
    gi_end,
    feature_start,
    feature_end,
):

    if not overlap(
        gi_start,
        gi_end,
        feature_start,
        feature_end
    ):
        return ""

    if (
        feature_start >= gi_start
        and feature_end <= gi_end
    ):
        return "feature_fully_inside_GI"

    if (
        gi_start >= feature_start
        and gi_end <= feature_end
    ):
        return "GI_fully_inside_feature"

    return "partial_overlap"


def unique_join(values):

    clean = sorted(
        {
            str(x).strip()
            for x in values
            if str(x).strip()
        }
    )

    return ";".join(clean)


# ======================================================================
# Load inputs
# ======================================================================

GI = read_tsv(GI_FILE)
AMR = read_tsv(AMR_FILE)
IS = read_tsv(IS_FILE)
INTEGRON = read_tsv(INTEGRON_FILE)
PLASMID = read_tsv(PLASMID_FILE)
CTX = read_tsv(CTX_FILE)


require_columns(
    GI,
    [
        "GI_candidate_id",
        "coassembly",
        "MAG",
        "contig",
        "start",
        "end",
        "length_bp",
        "edge_truncation_possible",
    ],
    "GI"
)

require_columns(
    AMR,
    [
        "community_locus_id",
        "coassembly",
        "contig",
        "AMR_start",
        "AMR_end",
        "AMR_symbol",
        "AMR_class",
    ],
    "AMR"
)

require_columns(
    IS,
    [
        "IS_id",
        "coassembly",
        "contig",
        "family",
        "IS_completeness",
        "is_start",
        "is_end",
    ],
    "IS"
)

require_columns(
    INTEGRON,
    [
        "structure_id",
        "coassembly",
        "contig",
        "integron_type",
        "structure_start",
        "structure_end",
    ],
    "INTEGRON"
)

require_columns(
    PLASMID,
    [
        "seq_name",
        "coassembly",
        "original_contig",
        "plasmid_score",
        "has_conjugation_annotation",
        "conjugation_genes_genomad",
        "replicon_supported",
        "plasmidfinder_replicons",
    ],
    "PLASMID"
)

require_columns(
    CTX,
    [
        "context_id",
        "coassembly",
        "contig",
        "context_start",
        "context_end",
        "component_signature",
        "exact_cluster_id",
    ],
    "CTX"
)


# Expected upstream counts
expected = {
    "GI": 1251,
    "AMR": 565,
    "IS": 6391,
    "INTEGRON": 375,
    "PLASMID": 13599,
    "CTX": 22,
}

observed = {
    "GI": len(GI),
    "AMR": len(AMR),
    "IS": len(IS),
    "INTEGRON": len(INTEGRON),
    "PLASMID": len(PLASMID),
    "CTX": len(CTX),
}

for key in expected:

    if observed[key] != expected[key]:

        raise RuntimeError(
            f"{key}: esperado={expected[key]} "
            f"observado={observed[key]}"
        )


# ======================================================================
# Build indices by exact coassembly + original contig
# ======================================================================

amr_by_contig = defaultdict(list)
is_by_contig = defaultdict(list)
integron_by_contig = defaultdict(list)
plasmid_by_contig = defaultdict(list)
ctx_by_contig = defaultdict(list)


for r in AMR:
    amr_by_contig[
        (
            r["coassembly"],
            r["contig"],
        )
    ].append(r)


for r in IS:
    is_by_contig[
        (
            r["coassembly"],
            r["contig"],
        )
    ].append(r)


for r in INTEGRON:
    integron_by_contig[
        (
            r["coassembly"],
            r["contig"],
        )
    ].append(r)


for r in PLASMID:
    plasmid_by_contig[
        (
            r["coassembly"],
            r["original_contig"],
        )
    ].append(r)


for r in CTX:
    ctx_by_contig[
        (
            r["coassembly"],
            r["contig"],
        )
    ].append(r)


# ======================================================================
# Link tables
# ======================================================================

amr_links = []
is_links = []
integron_links = []
context_links = []

integrated = []


for gi in GI:

    gid = gi["GI_candidate_id"]

    coassembly = gi["coassembly"]
    contig = gi["contig"]

    key = (
        coassembly,
        contig
    )

    gs = integer(
        gi["start"],
        f"{gid}.start"
    )

    ge = integer(
        gi["end"],
        f"{gid}.end"
    )

    # ------------------------------------------------------------------
    # ARG overlap
    # ------------------------------------------------------------------

    linked_amr = []

    for r in amr_by_contig.get(
        key,
        []
    ):

        fs = integer(
            r["AMR_start"],
            "AMR_start"
        )

        fe = integer(
            r["AMR_end"],
            "AMR_end"
        )

        if not overlap(
            gs,
            ge,
            fs,
            fe
        ):
            continue

        relation = relation_to_GI(
            gs,
            ge,
            fs,
            fe
        )

        ov = overlap_length(
            gs,
            ge,
            fs,
            fe
        )

        rec = {
            "GI_candidate_id": gid,
            "coassembly": coassembly,
            "MAG": gi["MAG"],
            "contig": contig,
            "GI_start": gs,
            "GI_end": ge,
            "community_locus_id": r["community_locus_id"],
            "AMR_start": fs,
            "AMR_end": fe,
            "AMR_symbol": r["AMR_symbol"],
            "AMR_class": r["AMR_class"],
            "AMR_method": r.get(
                "AMR_method",
                ""
            ),
            "overlap_bp": ov,
            "relation_to_GI": relation,
        }

        amr_links.append(rec)
        linked_amr.append(r)


    # ------------------------------------------------------------------
    # IS overlap
    # ------------------------------------------------------------------

    linked_is = []

    for r in is_by_contig.get(
        key,
        []
    ):

        fs = integer(
            r["is_start"],
            "IS_start"
        )

        fe = integer(
            r["is_end"],
            "IS_end"
        )

        if not overlap(
            gs,
            ge,
            fs,
            fe
        ):
            continue

        relation = relation_to_GI(
            gs,
            ge,
            fs,
            fe
        )

        ov = overlap_length(
            gs,
            ge,
            fs,
            fe
        )

        rec = {
            "GI_candidate_id": gid,
            "coassembly": coassembly,
            "MAG": gi["MAG"],
            "contig": contig,
            "GI_start": gs,
            "GI_end": ge,
            "IS_id": r["IS_id"],
            "IS_family": r["family"],
            "IS_type": r.get(
                "type",
                ""
            ),
            "IS_completeness": r["IS_completeness"],
            "IS_start": fs,
            "IS_end": fe,
            "overlap_bp": ov,
            "relation_to_GI": relation,
        }

        is_links.append(rec)
        linked_is.append(r)


    # ------------------------------------------------------------------
    # Integron overlap
    # ------------------------------------------------------------------

    linked_integrons = []

    for r in integron_by_contig.get(
        key,
        []
    ):

        fs = integer(
            r["structure_start"],
            "structure_start"
        )

        fe = integer(
            r["structure_end"],
            "structure_end"
        )

        if not overlap(
            gs,
            ge,
            fs,
            fe
        ):
            continue

        relation = relation_to_GI(
            gs,
            ge,
            fs,
            fe
        )

        ov = overlap_length(
            gs,
            ge,
            fs,
            fe
        )

        rec = {
            "GI_candidate_id": gid,
            "coassembly": coassembly,
            "MAG": gi["MAG"],
            "contig": contig,
            "GI_start": gs,
            "GI_end": ge,
            "structure_id": r["structure_id"],
            "integron_type": r["integron_type"],
            "structure_start": fs,
            "structure_end": fe,
            "n_attC": r.get(
                "n_attC",
                ""
            ),
            "n_AMR_direct_structure": r.get(
                "n_AMR_direct_structure",
                ""
            ),
            "AMR_symbols_direct": r.get(
                "AMR_symbols_direct",
                ""
            ),
            "overlap_bp": ov,
            "relation_to_GI": relation,
        }

        integron_links.append(rec)
        linked_integrons.append(r)


    # ------------------------------------------------------------------
    # Reconstructed ARG-IS context overlap
    # DERIVED support only; not counted as independent feature class.
    # ------------------------------------------------------------------

    linked_ctx = []

    for r in ctx_by_contig.get(
        key,
        []
    ):

        fs = integer(
            r["context_start"],
            "context_start"
        )

        fe = integer(
            r["context_end"],
            "context_end"
        )

        if not overlap(
            gs,
            ge,
            fs,
            fe
        ):
            continue

        relation = relation_to_GI(
            gs,
            ge,
            fs,
            fe
        )

        ov = overlap_length(
            gs,
            ge,
            fs,
            fe
        )

        rec = {
            "GI_candidate_id": gid,
            "coassembly": coassembly,
            "MAG": gi["MAG"],
            "contig": contig,
            "GI_start": gs,
            "GI_end": ge,
            "context_id": r["context_id"],
            "context_start": fs,
            "context_end": fe,
            "component_signature": r["component_signature"],
            "exact_cluster_id": r["exact_cluster_id"],
            "plasmid_candidate_context": r.get(
                "plasmid_candidate",
                ""
            ),
            "overlap_bp": ov,
            "relation_to_GI": relation,
        }

        context_links.append(rec)
        linked_ctx.append(r)


    # ------------------------------------------------------------------
    # Plasmid-candidate context at WHOLE-CONTIG level
    # ------------------------------------------------------------------

    plasmid_rows = plasmid_by_contig.get(
        key,
        []
    )

    plasmid_scores = [
        numeric(
            x["plasmid_score"]
        )
        for x in plasmid_rows
    ]

    plasmid_scores = [
        x
        for x in plasmid_scores
        if x is not None
    ]


    # ------------------------------------------------------------------
    # Evidence fields
    # ------------------------------------------------------------------

    has_arg = (
        len(linked_amr) > 0
    )

    has_is = (
        len(linked_is) > 0
    )

    has_integron = (
        len(linked_integrons) > 0
    )

    has_complete_integron = any(
        str(x["integron_type"]).strip().lower()
        == "complete"
        for x in linked_integrons
    )

    has_in0 = any(
        str(x["integron_type"]).strip().lower()
        == "in0"
        for x in linked_integrons
    )

    has_calin = any(
        str(x["integron_type"]).strip().lower()
        == "calin"
        for x in linked_integrons
    )

    has_mobile_marker = (
        has_is
        or has_integron
    )

    has_arg_plus_mobile = (
        has_arg
        and has_mobile_marker
    )

    on_plasmid_candidate = (
        len(plasmid_rows) > 0
    )

    has_mobile_context = (
        len(linked_ctx) > 0
    )


    if has_arg_plus_mobile:

        evidence_class = (
            "composition_plus_ARG_plus_mobile_marker"
        )

    elif has_mobile_marker:

        evidence_class = (
            "composition_plus_mobile_marker"
        )

    elif has_arg:

        evidence_class = (
            "composition_plus_ARG"
        )

    else:

        evidence_class = (
            "composition_only"
        )


    mobile_marker_types = []

    if has_is:
        mobile_marker_types.append(
            "IS"
        )

    if has_integron:
        mobile_marker_types.append(
            "integron_related_structure"
        )


    # IMPORTANT:
    # plasmid status remains a separate CONTIG-LEVEL context.
    if on_plasmid_candidate:

        plasmid_context = (
            "on_genomad_plasmid_candidate_contig"
        )

        chromosomal_GI_caution = "YES"

    else:

        plasmid_context = (
            "not_on_genomad_plasmid_candidate_contig"
        )

        chromosomal_GI_caution = "NO"


    outrow = dict(gi)

    outrow.update({
        "has_overlapping_ARG": (
            "YES" if has_arg else "NO"
        ),
        "n_overlapping_ARG_loci": len(linked_amr),
        "AMR_locus_ids": unique_join(
            x["community_locus_id"]
            for x in linked_amr
        ),
        "AMR_symbols": unique_join(
            x["AMR_symbol"]
            for x in linked_amr
        ),
        "AMR_classes": unique_join(
            x["AMR_class"]
            for x in linked_amr
        ),

        "has_overlapping_IS": (
            "YES" if has_is else "NO"
        ),
        "n_overlapping_IS": len(linked_is),
        "IS_ids": unique_join(
            x["IS_id"]
            for x in linked_is
        ),
        "IS_families": unique_join(
            x["family"]
            for x in linked_is
        ),
        "n_complete_IS": sum(
            str(x["IS_completeness"])
            .strip()
            .lower()
            == "complete"
            for x in linked_is
        ),
        "n_partial_IS": sum(
            str(x["IS_completeness"])
            .strip()
            .lower()
            == "partial"
            for x in linked_is
        ),

        "has_overlapping_integron_structure": (
            "YES"
            if has_integron
            else "NO"
        ),
        "n_overlapping_integron_structures": (
            len(linked_integrons)
        ),
        "integron_structure_ids": unique_join(
            x["structure_id"]
            for x in linked_integrons
        ),
        "integron_types": unique_join(
            x["integron_type"]
            for x in linked_integrons
        ),
        "has_complete_integron": (
            "YES"
            if has_complete_integron
            else "NO"
        ),
        "has_In0": (
            "YES"
            if has_in0
            else "NO"
        ),
        "has_CALIN": (
            "YES"
            if has_calin
            else "NO"
        ),

        "has_interval_resolved_mobile_marker": (
            "YES"
            if has_mobile_marker
            else "NO"
        ),
        "mobile_marker_types": (
            ";".join(
                mobile_marker_types
            )
        ),
        "has_ARG_plus_mobile_marker": (
            "YES"
            if has_arg_plus_mobile
            else "NO"
        ),

        "evidence_class": evidence_class,

        "plasmid_candidate_contig": (
            "YES"
            if on_plasmid_candidate
            else "NO"
        ),
        "n_genomad_plasmid_candidates_on_contig": (
            len(plasmid_rows)
        ),
        "plasmid_seq_names": unique_join(
            x["seq_name"]
            for x in plasmid_rows
        ),
        "plasmid_max_score": (
            max(plasmid_scores)
            if plasmid_scores
            else ""
        ),
        "plasmid_replicon_supported": (
            "YES"
            if any(
                yes(
                    x.get(
                        "replicon_supported",
                        ""
                    )
                )
                for x in plasmid_rows
            )
            else "NO"
        ),
        "plasmidfinder_replicons": unique_join(
            x.get(
                "plasmidfinder_replicons",
                ""
            )
            for x in plasmid_rows
        ),
        "plasmid_conjugation_annotation": (
            "YES"
            if any(
                yes(
                    x.get(
                        "has_conjugation_annotation",
                        ""
                    )
                )
                for x in plasmid_rows
            )
            else "NO"
        ),
        "plasmid_conjugation_genes": unique_join(
            x.get(
                "conjugation_genes_genomad",
                ""
            )
            for x in plasmid_rows
        ),
        "plasmid_context": plasmid_context,

        "chromosomal_GI_interpretation_caution": (
            chromosomal_GI_caution
        ),

        "has_overlapping_93D5_mobile_context": (
            "YES"
            if has_mobile_context
            else "NO"
        ),
        "n_overlapping_93D5_mobile_contexts": (
            len(linked_ctx)
        ),
        "mobile_context_ids": unique_join(
            x["context_id"]
            for x in linked_ctx
        ),
        "mobile_context_exact_clusters": unique_join(
            x["exact_cluster_id"]
            for x in linked_ctx
        ),
        "mobile_context_is_derived_evidence": (
            "YES"
        ),

        "interpretation_guardrail": (
            "candidate_GI_not_validated_HGT"
        ),
    })

    integrated.append(
        outrow
    )


# ======================================================================
# Write link tables
# ======================================================================

def write_rows(
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
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows)


write_rows(
    OUT / "93E3B_GI_ARG_links.tsv",
    amr_links,
    [
        "GI_candidate_id",
        "coassembly",
        "MAG",
        "contig",
        "GI_start",
        "GI_end",
        "community_locus_id",
        "AMR_start",
        "AMR_end",
        "AMR_symbol",
        "AMR_class",
        "AMR_method",
        "overlap_bp",
        "relation_to_GI",
    ]
)


write_rows(
    OUT / "93E3B_GI_IS_links.tsv",
    is_links,
    [
        "GI_candidate_id",
        "coassembly",
        "MAG",
        "contig",
        "GI_start",
        "GI_end",
        "IS_id",
        "IS_family",
        "IS_type",
        "IS_completeness",
        "IS_start",
        "IS_end",
        "overlap_bp",
        "relation_to_GI",
    ]
)


write_rows(
    OUT / "93E3B_GI_integron_links.tsv",
    integron_links,
    [
        "GI_candidate_id",
        "coassembly",
        "MAG",
        "contig",
        "GI_start",
        "GI_end",
        "structure_id",
        "integron_type",
        "structure_start",
        "structure_end",
        "n_attC",
        "n_AMR_direct_structure",
        "AMR_symbols_direct",
        "overlap_bp",
        "relation_to_GI",
    ]
)


write_rows(
    OUT / "93E3B_GI_mobile_context_links.tsv",
    context_links,
    [
        "GI_candidate_id",
        "coassembly",
        "MAG",
        "contig",
        "GI_start",
        "GI_end",
        "context_id",
        "context_start",
        "context_end",
        "component_signature",
        "exact_cluster_id",
        "plasmid_candidate_context",
        "overlap_bp",
        "relation_to_GI",
    ]
)


# ======================================================================
# Integrated master table
# ======================================================================

if len(integrated) != 1251:
    raise RuntimeError(
        f"Integrated GI rows={len(integrated)}, expected=1251"
    )

if len(
    {
        x["GI_candidate_id"]
        for x in integrated
    }
) != 1251:

    raise RuntimeError(
        "GI_candidate_id no es único"
    )


master_fields = list(
    integrated[0].keys()
)

write_rows(
    OUT / "93E3B_GI_mobile_integrated.tsv",
    integrated,
    master_fields
)


# ======================================================================
# Focus tables
# ======================================================================

arg_mobile = [
    x
    for x in integrated
    if x["has_ARG_plus_mobile_marker"]
    == "YES"
]

mobile_supported = [
    x
    for x in integrated
    if x[
        "has_interval_resolved_mobile_marker"
    ] == "YES"
]

complete_integron = [
    x
    for x in integrated
    if x["has_complete_integron"]
    == "YES"
]

plasmid_context_candidates = [
    x
    for x in integrated
    if x["plasmid_candidate_contig"]
    == "YES"
]

composition_only = [
    x
    for x in integrated
    if x["evidence_class"]
    == "composition_only"
]


write_rows(
    OUT / "93E3B_ARG_plus_mobile_marker_GI.tsv",
    arg_mobile,
    master_fields
)

write_rows(
    OUT / "93E3B_mobile_marker_supported_GI.tsv",
    mobile_supported,
    master_fields
)

write_rows(
    OUT / "93E3B_complete_integron_GI.tsv",
    complete_integron,
    master_fields
)

write_rows(
    OUT / "93E3B_plasmid_context_GI.tsv",
    plasmid_context_candidates,
    master_fields
)

write_rows(
    OUT / "93E3B_composition_only_GI.tsv",
    composition_only,
    master_fields
)


# ======================================================================
# Evidence-class summary
# ======================================================================

class_counts = Counter(
    x["evidence_class"]
    for x in integrated
)

with (
    OUT
    / "93E3B_evidence_class_summary.tsv"
).open("w") as fh:

    fh.write(
        "evidence_class\tcandidate_regions\n"
    )

    for cls, n in sorted(
        class_counts.items()
    ):
        fh.write(
            f"{cls}\t{n}\n"
        )


# ======================================================================
# Per MAG / coassembly summaries
# ======================================================================

def grouped_summary(
    group_field,
    filename
):

    groups = defaultdict(list)

    for x in integrated:
        groups[
            x[group_field]
        ].append(x)

    rows = []

    for group in sorted(groups):

        rr = groups[group]

        rows.append({
            group_field: group,
            "candidate_regions": len(rr),
            "with_ARG": sum(
                x["has_overlapping_ARG"]
                == "YES"
                for x in rr
            ),
            "with_IS": sum(
                x["has_overlapping_IS"]
                == "YES"
                for x in rr
            ),
            "with_integron_structure": sum(
                x[
                    "has_overlapping_integron_structure"
                ]
                == "YES"
                for x in rr
            ),
            "with_complete_integron": sum(
                x["has_complete_integron"]
                == "YES"
                for x in rr
            ),
            "with_mobile_marker": sum(
                x[
                    "has_interval_resolved_mobile_marker"
                ]
                == "YES"
                for x in rr
            ),
            "ARG_plus_mobile_marker": sum(
                x["has_ARG_plus_mobile_marker"]
                == "YES"
                for x in rr
            ),
            "on_plasmid_candidate_contig": sum(
                x["plasmid_candidate_contig"]
                == "YES"
                for x in rr
            ),
            "with_93D5_mobile_context": sum(
                x[
                    "has_overlapping_93D5_mobile_context"
                ]
                == "YES"
                for x in rr
            ),
            "composition_only": sum(
                x["evidence_class"]
                == "composition_only"
                for x in rr
            ),
            "edge_touching": sum(
                x[
                    "edge_truncation_possible"
                ]
                == "YES"
                for x in rr
            ),
        })

    fields = list(
        rows[0].keys()
    )

    write_rows(
        OUT / filename,
        rows,
        fields
    )


grouped_summary(
    "MAG",
    "93E3B_summary_by_MAG.tsv"
)

grouped_summary(
    "coassembly",
    "93E3B_summary_by_coassembly.tsv"
)


# ======================================================================
# Global summary
# ======================================================================

unique_amr = {
    x["community_locus_id"]
    for x in amr_links
}

unique_is = {
    x["IS_id"]
    for x in is_links
}

unique_integrons = {
    x["structure_id"]
    for x in integron_links
}

unique_contexts = {
    x["context_id"]
    for x in context_links
}

plasmid_contigs_screened = {
    (
        x["coassembly"],
        x["contig"],
    )
    for x in integrated
    if x["plasmid_candidate_contig"]
    == "YES"
}


metrics = [
    (
        "candidate_GI_regions",
        len(integrated)
    ),
    (
        "GI_with_overlapping_ARG",
        sum(
            x["has_overlapping_ARG"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "unique_ARG_loci_overlapping_GI",
        len(unique_amr)
    ),
    (
        "GI_with_overlapping_IS",
        sum(
            x["has_overlapping_IS"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "unique_IS_predictions_overlapping_GI",
        len(unique_is)
    ),
    (
        "GI_with_overlapping_integron_structure",
        sum(
            x[
                "has_overlapping_integron_structure"
            ]
            == "YES"
            for x in integrated
        )
    ),
    (
        "unique_integron_structures_overlapping_GI",
        len(unique_integrons)
    ),
    (
        "GI_with_complete_integron",
        sum(
            x["has_complete_integron"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "GI_with_In0",
        sum(
            x["has_In0"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "GI_with_CALIN",
        sum(
            x["has_CALIN"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "GI_with_interval_resolved_mobile_marker",
        sum(
            x[
                "has_interval_resolved_mobile_marker"
            ]
            == "YES"
            for x in integrated
        )
    ),
    (
        "GI_with_ARG_plus_mobile_marker",
        len(arg_mobile)
    ),
    (
        "GI_on_genomad_plasmid_candidate_contig",
        len(plasmid_context_candidates)
    ),
    (
        "unique_screened_contigs_also_plasmid_candidate",
        len(plasmid_contigs_screened)
    ),
    (
        "GI_on_replicon_supported_plasmid_candidate_contig",
        sum(
            x["plasmid_replicon_supported"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "GI_on_conjugation_annotated_plasmid_candidate_contig",
        sum(
            x["plasmid_conjugation_annotation"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "GI_with_93D5_mobile_context_overlap",
        sum(
            x[
                "has_overlapping_93D5_mobile_context"
            ]
            == "YES"
            for x in integrated
        )
    ),
    (
        "unique_93D5_mobile_contexts_overlapping_GI",
        len(unique_contexts)
    ),
    (
        "composition_only_GI",
        len(composition_only)
    ),
    (
        "edge_touching_GI",
        sum(
            x["edge_truncation_possible"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "ARG_GI_link_rows",
        len(amr_links)
    ),
    (
        "IS_GI_link_rows",
        len(is_links)
    ),
    (
        "integron_GI_link_rows",
        len(integron_links)
    ),
    (
        "mobile_context_GI_link_rows",
        len(context_links)
    ),
]


with (
    OUT
    / "93E3B_global_summary.tsv"
).open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    for key, value in metrics:
        fh.write(
            f"{key}\t{value}\n"
        )


# ======================================================================
# Scientific guardrails
# ======================================================================

with (
    OUT
    / "93E3B_methodological_scope.tsv"
).open("w") as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        "GI_definition\t"
        "AlienHunter_composition_based_candidate_region\n"
    )

    fh.write(
        "GI_prediction_equals_validated_genomic_island\tNO\n"
    )

    fh.write(
        "GI_prediction_equals_demonstrated_HGT\tNO\n"
    )

    fh.write(
        "ARG_overlap_equals_mobility_evidence\tNO\n"
    )

    fh.write(
        "interval_resolved_mobile_markers\tIS_and_integron_related_structures\n"
    )

    fh.write(
        "integron_types_preserved\tcomplete_In0_CALIN\n"
    )

    fh.write(
        "CALIN_or_In0_called_complete_integron\tNO\n"
    )

    fh.write(
        "plasmid_candidate_status_resolution\twhole_contig_context\n"
    )

    fh.write(
        "plasmid_candidate_equals_validated_plasmid\tNO\n"
    )

    fh.write(
        "plasmid_context_used_as_interval_overlap\tNO\n"
    )

    fh.write(
        "93D5_mobile_context_independent_evidence\tNO\n"
    )

    fh.write(
        "edge_touching_candidate_boundary_complete\tNO\n"
    )


# ======================================================================
# Validation
# ======================================================================

required_outputs = [
    "93E3B_GI_mobile_integrated.tsv",
    "93E3B_GI_ARG_links.tsv",
    "93E3B_GI_IS_links.tsv",
    "93E3B_GI_integron_links.tsv",
    "93E3B_GI_mobile_context_links.tsv",
    "93E3B_ARG_plus_mobile_marker_GI.tsv",
    "93E3B_mobile_marker_supported_GI.tsv",
    "93E3B_complete_integron_GI.tsv",
    "93E3B_plasmid_context_GI.tsv",
    "93E3B_composition_only_GI.tsv",
    "93E3B_evidence_class_summary.tsv",
    "93E3B_summary_by_MAG.tsv",
    "93E3B_summary_by_coassembly.tsv",
    "93E3B_global_summary.tsv",
    "93E3B_methodological_scope.tsv",
]


for name in required_outputs:

    path = OUT / name

    if not path.is_file():
        raise RuntimeError(
            f"Salida faltante: {path}"
        )


print(
    f"GI_CANDIDATES={len(integrated)}"
)

print(
    f"ARG_GI_LINKS={len(amr_links)}"
)

print(
    f"IS_GI_LINKS={len(is_links)}"
)

print(
    f"INTEGRON_GI_LINKS={len(integron_links)}"
)

print(
    f"ARG_PLUS_MOBILE={len(arg_mobile)}"
)

print(
    "93E3B_INTEGRATION=PASS"
)

