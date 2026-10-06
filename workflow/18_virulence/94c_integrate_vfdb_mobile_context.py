#!/usr/bin/env python3

import csv
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path(sys.argv[1])

VFDB_FILE = (
    ROOT
    / "96_virulence/94B2_VFDB_consolidated"
    / "94B2_community_high90_90.tsv"
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

GI_FILE = (
    ROOT
    / "95_mobilome/93E3B_GI_mobile_integration"
    / "93E3B_GI_mobile_integrated.tsv"
)

AMR_FILE = (
    ROOT
    / "95_mobilome/93D4_IS_ARG_integron_plasmid"
    / "93D4_AMR_nearest_IS.tsv"
)

OUT = (
    ROOT
    / "96_virulence/94C_VFDB_mobile_context"
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

    with path.open(errors="replace") as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        rows = list(reader)
        fields = reader.fieldnames or []

    return rows, fields


def require_columns(rows, required, label):

    if not rows:
        raise RuntimeError(
            f"{label}: tabla vacía"
        )

    missing = (
        set(required)
        - set(rows[0].keys())
    )

    if missing:
        raise RuntimeError(
            f"{label}: faltan columnas: "
            + ",".join(sorted(missing))
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
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows)


def integer(value):

    return int(float(value))


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


def interval_distance(a1, a2, b1, b2):
    """
    Número de bases interpuestas entre intervalos.
    overlap -> 0
    adyacentes -> 0
    """

    if overlap(
        a1,
        a2,
        b1,
        b2
    ):
        return 0

    if a2 < b1:
        return max(
            0,
            b1 - a2 - 1
        )

    return max(
        0,
        a1 - b2 - 1
    )


def overlap_bp(a1, a2, b1, b2):

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


# ======================================================================
# Load
# ======================================================================

VFDB, vf_fields = read_tsv(
    VFDB_FILE
)

IS, _ = read_tsv(
    IS_FILE
)

INTEGRON, _ = read_tsv(
    INTEGRON_FILE
)

PLASMID, _ = read_tsv(
    PLASMID_FILE
)

GI, _ = read_tsv(
    GI_FILE
)

AMR, _ = read_tsv(
    AMR_FILE
)


# ======================================================================
# Validate schemas
# ======================================================================

require_columns(
    VFDB,
    [
        "community_vfdb_locus_id",
        "coassembly",
        "SEQUENCE",
        "START",
        "END",
        "GENE",
        "ACCESSION",
        "%COVERAGE",
        "%IDENTITY",
        "assignment_level",
        "final18_MAGs_exact",
        "vf_factor_annotation",
    ],
    "VFDB"
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
        "replicon_supported",
        "plasmidfinder_replicons",
        "has_conjugation_annotation",
        "conjugation_genes_genomad",
    ],
    "PLASMID"
)

require_columns(
    GI,
    [
        "GI_candidate_id",
        "coassembly",
        "contig",
        "start",
        "end",
        "has_interval_resolved_mobile_marker",
        "has_ARG_plus_mobile_marker",
        "evidence_class",
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


expected_counts = {
    "VFDB": 51,
    "IS": 6391,
    "INTEGRON": 375,
    "PLASMID": 13599,
    "GI": 1251,
    "AMR": 565,
}

observed_counts = {
    "VFDB": len(VFDB),
    "IS": len(IS),
    "INTEGRON": len(INTEGRON),
    "PLASMID": len(PLASMID),
    "GI": len(GI),
    "AMR": len(AMR),
}

for key in expected_counts:

    if (
        observed_counts[key]
        != expected_counts[key]
    ):

        raise RuntimeError(
            f"{key}: "
            f"esperado={expected_counts[key]} "
            f"observado={observed_counts[key]}"
        )


# ======================================================================
# Indices
# ======================================================================

is_by_contig = defaultdict(list)
int_by_contig = defaultdict(list)
plasmid_by_contig = defaultdict(list)
gi_by_contig = defaultdict(list)
amr_by_contig = defaultdict(list)


for row in IS:

    is_by_contig[
        (
            row["coassembly"],
            row["contig"],
        )
    ].append(row)


for row in INTEGRON:

    int_by_contig[
        (
            row["coassembly"],
            row["contig"],
        )
    ].append(row)


for row in PLASMID:

    plasmid_by_contig[
        (
            row["coassembly"],
            row["original_contig"],
        )
    ].append(row)


for row in GI:

    gi_by_contig[
        (
            row["coassembly"],
            row["contig"],
        )
    ].append(row)


for row in AMR:

    amr_by_contig[
        (
            row["coassembly"],
            row["contig"],
        )
    ].append(row)


# ======================================================================
# Integration
# ======================================================================

integrated = []

is_links = []
integron_links = []
gi_links = []
amr_links = []


for vf in VFDB:

    locus_id = (
        vf["community_vfdb_locus_id"]
    )

    group = vf["coassembly"]
    contig = vf["SEQUENCE"]

    key = (
        group,
        contig,
    )

    vs = integer(
        vf["START"]
    )

    ve = integer(
        vf["END"]
    )

    if vs > ve:
        vs, ve = ve, vs


    # ==================================================================
    # IS
    # ==================================================================

    is_candidates = []

    for row in is_by_contig.get(
        key,
        []
    ):

        fs = integer(
            row["is_start"]
        )

        fe = integer(
            row["is_end"]
        )

        if fs > fe:
            fs, fe = fe, fs

        dist = interval_distance(
            vs,
            ve,
            fs,
            fe
        )

        ov = overlap_bp(
            vs,
            ve,
            fs,
            fe
        )

        if dist <= 5000:

            rec = {
                "community_vfdb_locus_id":
                    locus_id,
                "coassembly": group,
                "contig": contig,
                "VFDB_gene": vf["GENE"],
                "VFDB_start": vs,
                "VFDB_end": ve,
                "IS_id": row["IS_id"],
                "IS_family": row["family"],
                "IS_completeness":
                    row["IS_completeness"],
                "IS_start": fs,
                "IS_end": fe,
                "distance_bp": dist,
                "coordinate_overlap":
                    "YES"
                    if ov > 0
                    else "NO",
                "overlap_bp": ov,
            }

            is_links.append(rec)
            is_candidates.append(rec)


    nearest_is_distance = (
        min(
            x["distance_bp"]
            for x in is_candidates
        )
        if is_candidates
        else ""
    )

    overlapping_is = [
        x
        for x in is_candidates
        if x["coordinate_overlap"]
        == "YES"
    ]

    is_within1 = [
        x
        for x in is_candidates
        if x["distance_bp"] <= 1000
    ]

    is_within5 = is_candidates


    # ==================================================================
    # Integrons
    # ==================================================================

    int_candidates = []

    for row in int_by_contig.get(
        key,
        []
    ):

        fs = integer(
            row["structure_start"]
        )

        fe = integer(
            row["structure_end"]
        )

        if fs > fe:
            fs, fe = fe, fs

        dist = interval_distance(
            vs,
            ve,
            fs,
            fe
        )

        ov = overlap_bp(
            vs,
            ve,
            fs,
            fe
        )

        if dist <= 5000:

            rec = {
                "community_vfdb_locus_id":
                    locus_id,
                "coassembly": group,
                "contig": contig,
                "VFDB_gene": vf["GENE"],
                "VFDB_start": vs,
                "VFDB_end": ve,
                "structure_id":
                    row["structure_id"],
                "integron_type":
                    row["integron_type"],
                "structure_start": fs,
                "structure_end": fe,
                "distance_bp": dist,
                "coordinate_overlap":
                    "YES"
                    if ov > 0
                    else "NO",
                "overlap_bp": ov,
            }

            integron_links.append(rec)
            int_candidates.append(rec)


    nearest_int_distance = (
        min(
            x["distance_bp"]
            for x in int_candidates
        )
        if int_candidates
        else ""
    )

    overlapping_int = [
        x
        for x in int_candidates
        if x["coordinate_overlap"]
        == "YES"
    ]

    int_within1 = [
        x
        for x in int_candidates
        if x["distance_bp"] <= 1000
    ]

    int_within5 = int_candidates


    # ==================================================================
    # AlienHunter GI candidate
    # ==================================================================

    gi_candidates = []

    for row in gi_by_contig.get(
        key,
        []
    ):

        gs = integer(
            row["start"]
        )

        ge = integer(
            row["end"]
        )

        if gs > ge:
            gs, ge = ge, gs

        if not overlap(
            vs,
            ve,
            gs,
            ge
        ):
            continue

        ov = overlap_bp(
            vs,
            ve,
            gs,
            ge
        )

        vf_fully_inside = (
            vs >= gs
            and ve <= ge
        )

        rec = {
            "community_vfdb_locus_id":
                locus_id,
            "coassembly": group,
            "contig": contig,
            "VFDB_gene": vf["GENE"],
            "VFDB_start": vs,
            "VFDB_end": ve,
            "GI_candidate_id":
                row["GI_candidate_id"],
            "GI_start": gs,
            "GI_end": ge,
            "overlap_bp": ov,
            "VFDB_fully_inside_GI":
                "YES"
                if vf_fully_inside
                else "NO",
            "GI_has_mobile_marker":
                row[
                    "has_interval_resolved_mobile_marker"
                ],
            "GI_has_ARG_plus_mobile_marker":
                row[
                    "has_ARG_plus_mobile_marker"
                ],
            "GI_evidence_class":
                row["evidence_class"],
        }

        gi_links.append(rec)
        gi_candidates.append(rec)


    gi_mobile_supported = [
        x
        for x in gi_candidates
        if yes(
            x["GI_has_mobile_marker"]
        )
    ]


    # ==================================================================
    # AMR
    # ==================================================================

    amr_candidates = []

    for row in amr_by_contig.get(
        key,
        []
    ):

        fs = integer(
            row["AMR_start"]
        )

        fe = integer(
            row["AMR_end"]
        )

        if fs > fe:
            fs, fe = fe, fs

        dist = interval_distance(
            vs,
            ve,
            fs,
            fe
        )

        ov = overlap_bp(
            vs,
            ve,
            fs,
            fe
        )

        if dist <= 5000:

            rec = {
                "community_vfdb_locus_id":
                    locus_id,
                "coassembly": group,
                "contig": contig,
                "VFDB_gene": vf["GENE"],
                "VFDB_start": vs,
                "VFDB_end": ve,
                "community_AMR_locus_id":
                    row["community_locus_id"],
                "AMR_symbol":
                    row["AMR_symbol"],
                "AMR_class":
                    row["AMR_class"],
                "AMR_start": fs,
                "AMR_end": fe,
                "distance_bp": dist,
                "coordinate_overlap":
                    "YES"
                    if ov > 0
                    else "NO",
                "overlap_bp": ov,
            }

            amr_links.append(rec)
            amr_candidates.append(rec)


    nearest_amr_distance = (
        min(
            x["distance_bp"]
            for x in amr_candidates
        )
        if amr_candidates
        else ""
    )

    amr_overlap = [
        x
        for x in amr_candidates
        if x["coordinate_overlap"]
        == "YES"
    ]

    amr_within1 = [
        x
        for x in amr_candidates
        if x["distance_bp"] <= 1000
    ]

    amr_within5 = amr_candidates


    # ==================================================================
    # Plasmid contig context
    # ==================================================================

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

    on_plasmid_candidate = (
        len(plasmid_rows) > 0
    )

    plasmid_replicon_supported = any(
        yes(
            x.get(
                "replicon_supported",
                ""
            )
        )
        for x in plasmid_rows
    )

    plasmid_conjugation = any(
        yes(
            x.get(
                "has_conjugation_annotation",
                ""
            )
        )
        for x in plasmid_rows
    )


    # ==================================================================
    # Evidence dimensions
    # ==================================================================

    direct_mobile_overlap = (
        bool(overlapping_is)
        or bool(overlapping_int)
    )

    mobile_within1kb = (
        bool(is_within1)
        or bool(int_within1)
    )

    mobile_within5kb = (
        bool(is_within5)
        or bool(int_within5)
    )

    overlaps_GI = bool(
        gi_candidates
    )

    overlaps_mobile_supported_GI = bool(
        gi_mobile_supported
    )

    has_AMR_5kb = bool(
        amr_within5
    )

    selected_context_features = []

    if direct_mobile_overlap:
        selected_context_features.append(
            "interval_mobile_marker_overlap"
        )

    if mobile_within1kb:
        selected_context_features.append(
            "mobile_marker_within_1kb"
        )

    if mobile_within5kb:
        selected_context_features.append(
            "mobile_marker_within_5kb"
        )

    if overlaps_GI:
        selected_context_features.append(
            "AlienHunter_GI_overlap"
        )

    if overlaps_mobile_supported_GI:
        selected_context_features.append(
            "mobile_supported_GI_overlap"
        )

    if on_plasmid_candidate:
        selected_context_features.append(
            "genomad_plasmid_candidate_contig"
        )

    if plasmid_replicon_supported:
        selected_context_features.append(
            "PlasmidFinder_replicon_support"
        )

    if plasmid_conjugation:
        selected_context_features.append(
            "genomad_conjugation_annotation"
        )

    if has_AMR_5kb:
        selected_context_features.append(
            "AMR_within_5kb"
        )


    # ==================================================================
    # Master row
    # ==================================================================

    rec = dict(vf)

    rec.update({
        "has_IS_overlap":
            "YES"
            if overlapping_is
            else "NO",

        "n_IS_overlap":
            len(overlapping_is),

        "has_IS_within_1kb":
            "YES"
            if is_within1
            else "NO",

        "n_IS_within_1kb":
            len(is_within1),

        "has_IS_within_5kb":
            "YES"
            if is_within5
            else "NO",

        "n_IS_within_5kb":
            len(is_within5),

        "nearest_IS_distance_bp":
            nearest_is_distance,

        "nearby_IS_ids":
            unique_join(
                x["IS_id"]
                for x in is_within5
            ),

        "nearby_IS_families":
            unique_join(
                x["IS_family"]
                for x in is_within5
            ),

        "has_integron_overlap":
            "YES"
            if overlapping_int
            else "NO",

        "n_integron_overlap":
            len(overlapping_int),

        "has_integron_within_1kb":
            "YES"
            if int_within1
            else "NO",

        "has_integron_within_5kb":
            "YES"
            if int_within5
            else "NO",

        "nearest_integron_distance_bp":
            nearest_int_distance,

        "nearby_integron_ids":
            unique_join(
                x["structure_id"]
                for x in int_within5
            ),

        "nearby_integron_types":
            unique_join(
                x["integron_type"]
                for x in int_within5
            ),

        "has_direct_mobile_marker_overlap":
            "YES"
            if direct_mobile_overlap
            else "NO",

        "has_mobile_marker_within_1kb":
            "YES"
            if mobile_within1kb
            else "NO",

        "has_mobile_marker_within_5kb":
            "YES"
            if mobile_within5kb
            else "NO",

        "overlaps_AlienHunter_GI":
            "YES"
            if overlaps_GI
            else "NO",

        "n_overlapping_GI":
            len(gi_candidates),

        "overlapping_GI_ids":
            unique_join(
                x["GI_candidate_id"]
                for x in gi_candidates
            ),

        "overlaps_mobile_supported_GI":
            "YES"
            if overlaps_mobile_supported_GI
            else "NO",

        "plasmid_candidate_contig":
            "YES"
            if on_plasmid_candidate
            else "NO",

        "n_genomad_plasmid_candidates":
            len(plasmid_rows),

        "plasmid_seq_names":
            unique_join(
                x["seq_name"]
                for x in plasmid_rows
            ),

        "plasmid_max_score":
            (
                max(plasmid_scores)
                if plasmid_scores
                else ""
            ),

        "plasmid_replicon_supported":
            "YES"
            if plasmid_replicon_supported
            else "NO",

        "plasmidfinder_replicons":
            unique_join(
                x.get(
                    "plasmidfinder_replicons",
                    ""
                )
                for x in plasmid_rows
            ),

        "plasmid_conjugation_annotation":
            "YES"
            if plasmid_conjugation
            else "NO",

        "plasmid_conjugation_genes":
            unique_join(
                x.get(
                    "conjugation_genes_genomad",
                    ""
                )
                for x in plasmid_rows
            ),

        "has_AMR_overlap":
            "YES"
            if amr_overlap
            else "NO",

        "has_AMR_within_1kb":
            "YES"
            if amr_within1
            else "NO",

        "has_AMR_within_5kb":
            "YES"
            if amr_within5
            else "NO",

        "nearest_AMR_distance_bp":
            nearest_amr_distance,

        "nearby_AMR_locus_ids":
            unique_join(
                x["community_AMR_locus_id"]
                for x in amr_within5
            ),

        "nearby_AMR_symbols":
            unique_join(
                x["AMR_symbol"]
                for x in amr_within5
            ),

        "nearby_AMR_classes":
            unique_join(
                x["AMR_class"]
                for x in amr_within5
            ),

        "selected_context_features":
            ";".join(
                selected_context_features
            ),

        "has_selected_mobile_context":
            "YES"
            if (
                direct_mobile_overlap
                or mobile_within1kb
                or overlaps_mobile_supported_GI
                or on_plasmid_candidate
            )
            else "NO",

        "interpretation_guardrail":
            "VFDB_homology_not_pathogenic_phenotype",
    })

    integrated.append(rec)


# ======================================================================
# Write links
# ======================================================================

is_fields = [
    "community_vfdb_locus_id",
    "coassembly",
    "contig",
    "VFDB_gene",
    "VFDB_start",
    "VFDB_end",
    "IS_id",
    "IS_family",
    "IS_completeness",
    "IS_start",
    "IS_end",
    "distance_bp",
    "coordinate_overlap",
    "overlap_bp",
]

write_tsv(
    OUT / "94C_VFDB_IS_links_within5kb.tsv",
    is_links,
    is_fields
)


int_fields = [
    "community_vfdb_locus_id",
    "coassembly",
    "contig",
    "VFDB_gene",
    "VFDB_start",
    "VFDB_end",
    "structure_id",
    "integron_type",
    "structure_start",
    "structure_end",
    "distance_bp",
    "coordinate_overlap",
    "overlap_bp",
]

write_tsv(
    OUT / "94C_VFDB_integron_links_within5kb.tsv",
    integron_links,
    int_fields
)


gi_fields = [
    "community_vfdb_locus_id",
    "coassembly",
    "contig",
    "VFDB_gene",
    "VFDB_start",
    "VFDB_end",
    "GI_candidate_id",
    "GI_start",
    "GI_end",
    "overlap_bp",
    "VFDB_fully_inside_GI",
    "GI_has_mobile_marker",
    "GI_has_ARG_plus_mobile_marker",
    "GI_evidence_class",
]

write_tsv(
    OUT / "94C_VFDB_GI_links.tsv",
    gi_links,
    gi_fields
)


amr_fields = [
    "community_vfdb_locus_id",
    "coassembly",
    "contig",
    "VFDB_gene",
    "VFDB_start",
    "VFDB_end",
    "community_AMR_locus_id",
    "AMR_symbol",
    "AMR_class",
    "AMR_start",
    "AMR_end",
    "distance_bp",
    "coordinate_overlap",
    "overlap_bp",
]

write_tsv(
    OUT / "94C_VFDB_AMR_links_within5kb.tsv",
    amr_links,
    amr_fields
)


# ======================================================================
# Master table
# ======================================================================

if len(integrated) != 51:
    raise RuntimeError(
        f"Master VFDB esperado=51 observado={len(integrated)}"
    )

if len(
    {
        x["community_vfdb_locus_id"]
        for x in integrated
    }
) != 51:
    raise RuntimeError(
        "IDs VFDB no son únicos"
    )


master_fields = list(
    integrated[0].keys()
)

write_tsv(
    OUT / "94C_VFDB_high90_90_mobile_integrated.tsv",
    integrated,
    master_fields
)


# ======================================================================
# Focus tables
# ======================================================================

selected_mobile = [
    x
    for x in integrated
    if x["has_selected_mobile_context"]
    == "YES"
]

direct_mobile = [
    x
    for x in integrated
    if x["has_direct_mobile_marker_overlap"]
    == "YES"
]

plasmid_context = [
    x
    for x in integrated
    if x["plasmid_candidate_contig"]
    == "YES"
]

mobile_GI = [
    x
    for x in integrated
    if x["overlaps_mobile_supported_GI"]
    == "YES"
]

amr5 = [
    x
    for x in integrated
    if x["has_AMR_within_5kb"]
    == "YES"
]


for filename, rows in [
    (
        "94C_selected_mobile_context_VFDB.tsv",
        selected_mobile,
    ),
    (
        "94C_direct_mobile_marker_overlap_VFDB.tsv",
        direct_mobile,
    ),
    (
        "94C_plasmid_context_VFDB.tsv",
        plasmid_context,
    ),
    (
        "94C_mobile_supported_GI_VFDB.tsv",
        mobile_GI,
    ),
    (
        "94C_AMR_within5kb_VFDB.tsv",
        amr5,
    ),
]:

    write_tsv(
        OUT / filename,
        rows,
        master_fields
    )


# ======================================================================
# Summary by coassembly
# ======================================================================

summary_group = []

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
        for x in integrated
        if x["coassembly"] == group
    ]

    summary_group.append({
        "coassembly": group,
        "VFDB_high90_90": len(rr),
        "exact_final18_contig": sum(
            x["assignment_level"]
            == "exact_final18_contig"
            for x in rr
        ),
        "community_only": sum(
            x["assignment_level"]
            == "community_contig_only"
            for x in rr
        ),
        "IS_overlap": sum(
            x["has_IS_overlap"]
            == "YES"
            for x in rr
        ),
        "IS_within_1kb": sum(
            x["has_IS_within_1kb"]
            == "YES"
            for x in rr
        ),
        "integron_overlap": sum(
            x["has_integron_overlap"]
            == "YES"
            for x in rr
        ),
        "integron_within_1kb": sum(
            x["has_integron_within_1kb"]
            == "YES"
            for x in rr
        ),
        "AlienHunter_GI_overlap": sum(
            x["overlaps_AlienHunter_GI"]
            == "YES"
            for x in rr
        ),
        "mobile_supported_GI_overlap": sum(
            x["overlaps_mobile_supported_GI"]
            == "YES"
            for x in rr
        ),
        "plasmid_candidate_contig": sum(
            x["plasmid_candidate_contig"]
            == "YES"
            for x in rr
        ),
        "AMR_within_5kb": sum(
            x["has_AMR_within_5kb"]
            == "YES"
            for x in rr
        ),
        "selected_mobile_context": sum(
            x["has_selected_mobile_context"]
            == "YES"
            for x in rr
        ),
    })


write_tsv(
    OUT / "94C_summary_by_coassembly.tsv",
    summary_group,
    list(
        summary_group[0].keys()
    )
)


# ======================================================================
# Summary by gene/accession
# ======================================================================

gene_groups = defaultdict(list)

for row in integrated:

    gene_groups[
        (
            row["GENE"],
            row["ACCESSION"],
            row["vf_factor_annotation"],
        )
    ].append(row)


gene_summary = []

for key, rr in gene_groups.items():

    gene_summary.append({
        "GENE": key[0],
        "ACCESSION": key[1],
        "vf_factor_annotation": key[2],
        "n_hits": len(rr),
        "coassemblies": unique_join(
            x["coassembly"]
            for x in rr
        ),
        "producers": unique_join(
            x["producer"]
            for x in rr
        ),
        "exact_final18_MAGs": unique_join(
            x["final18_MAGs_exact"]
            for x in rr
        ),
        "with_IS_overlap": sum(
            x["has_IS_overlap"] == "YES"
            for x in rr
        ),
        "with_mobile_marker_within_1kb": sum(
            x["has_mobile_marker_within_1kb"]
            == "YES"
            for x in rr
        ),
        "on_plasmid_candidate_contig": sum(
            x["plasmid_candidate_contig"]
            == "YES"
            for x in rr
        ),
        "overlap_AlienHunter_GI": sum(
            x["overlaps_AlienHunter_GI"]
            == "YES"
            for x in rr
        ),
        "overlap_mobile_supported_GI": sum(
            x["overlaps_mobile_supported_GI"]
            == "YES"
            for x in rr
        ),
        "AMR_within_5kb": sum(
            x["has_AMR_within_5kb"]
            == "YES"
            for x in rr
        ),
        "selected_mobile_context": sum(
            x["has_selected_mobile_context"]
            == "YES"
            for x in rr
        ),
    })


gene_summary.sort(
    key=lambda x: (
        -x["n_hits"],
        x["GENE"],
        x["ACCESSION"],
    )
)


write_tsv(
    OUT / "94C_summary_by_gene.tsv",
    gene_summary,
    list(
        gene_summary[0].keys()
    )
)


# ======================================================================
# Global summary
# ======================================================================

metrics = [
    (
        "VFDB_high90_90_total",
        len(integrated)
    ),
    (
        "VFDB_exact_final18_contig",
        sum(
            x["assignment_level"]
            == "exact_final18_contig"
            for x in integrated
        )
    ),
    (
        "VFDB_community_only",
        sum(
            x["assignment_level"]
            == "community_contig_only"
            for x in integrated
        )
    ),
    (
        "VFDB_with_IS_overlap",
        sum(
            x["has_IS_overlap"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "VFDB_with_IS_within_1kb",
        sum(
            x["has_IS_within_1kb"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "VFDB_with_IS_within_5kb",
        sum(
            x["has_IS_within_5kb"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "VFDB_with_integron_overlap",
        sum(
            x["has_integron_overlap"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "VFDB_with_integron_within_1kb",
        sum(
            x["has_integron_within_1kb"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "VFDB_with_integron_within_5kb",
        sum(
            x["has_integron_within_5kb"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "VFDB_with_direct_mobile_marker_overlap",
        len(direct_mobile)
    ),
    (
        "VFDB_with_mobile_marker_within_1kb",
        sum(
            x["has_mobile_marker_within_1kb"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "VFDB_with_mobile_marker_within_5kb",
        sum(
            x["has_mobile_marker_within_5kb"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "VFDB_overlapping_AlienHunter_GI",
        sum(
            x["overlaps_AlienHunter_GI"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "VFDB_overlapping_mobile_supported_GI",
        len(mobile_GI)
    ),
    (
        "VFDB_on_genomad_plasmid_candidate_contig",
        len(plasmid_context)
    ),
    (
        "VFDB_on_replicon_supported_plasmid_candidate_contig",
        sum(
            x["plasmid_replicon_supported"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "VFDB_on_conjugation_annotated_plasmid_candidate_contig",
        sum(
            x["plasmid_conjugation_annotation"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "VFDB_with_AMR_overlap",
        sum(
            x["has_AMR_overlap"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "VFDB_with_AMR_within_1kb",
        sum(
            x["has_AMR_within_1kb"]
            == "YES"
            for x in integrated
        )
    ),
    (
        "VFDB_with_AMR_within_5kb",
        len(amr5)
    ),
    (
        "VFDB_with_selected_mobile_context",
        len(selected_mobile)
    ),
    (
        "IS_link_rows_within5kb",
        len(is_links)
    ),
    (
        "integron_link_rows_within5kb",
        len(integron_links)
    ),
    (
        "GI_link_rows",
        len(gi_links)
    ),
    (
        "AMR_link_rows_within5kb",
        len(amr_links)
    ),
]


with (
    OUT / "94C_global_summary.tsv"
).open("w") as fh:

    fh.write("metric\tvalue\n")

    for key, value in metrics:
        fh.write(
            f"{key}\t{value}\n"
        )


# ======================================================================
# Methodological scope
# ======================================================================

with (
    OUT / "94C_methodological_scope.tsv"
).open("w") as fh:

    fh.write("field\tvalue\n")

    fh.write(
        "VFDB_scope\tcommunity_high90_90_only\n"
    )

    fh.write(
        "VFDB_hit_equals_pathogenic_phenotype\tNO\n"
    )

    fh.write(
        "VFDB_reference_organism_equals_sample_taxonomy\tNO\n"
    )

    fh.write(
        "direct_mobile_marker_definition\t"
        "coordinate_overlap_with_IS_or_integron_related_structure\n"
    )

    fh.write(
        "proximity_context_thresholds\t"
        "1kb_and_5kb_intervening_bases\n"
    )

    fh.write(
        "proximity_equals_same_mobile_element\tNO\n"
    )

    fh.write(
        "plasmid_candidate_context_resolution\twhole_contig\n"
    )

    fh.write(
        "genomad_plasmid_candidate_equals_validated_plasmid\tNO\n"
    )

    fh.write(
        "genomad_conjugation_annotation_equals_demonstrated_conjugation\tNO\n"
    )

    fh.write(
        "AlienHunter_GI_equals_validated_genomic_island\tNO\n"
    )

    fh.write(
        "CALIN_or_In0_equals_complete_integron\tNO\n"
    )

    fh.write(
        "AMR_proximity_equals_cotransfer\tNO\n"
    )

    fh.write(
        "cross_coassembly_recurrence_equals_prevalence\tNO\n"
    )


# ======================================================================
# Final validation
# ======================================================================

required_outputs = [
    "94C_VFDB_high90_90_mobile_integrated.tsv",
    "94C_VFDB_IS_links_within5kb.tsv",
    "94C_VFDB_integron_links_within5kb.tsv",
    "94C_VFDB_GI_links.tsv",
    "94C_VFDB_AMR_links_within5kb.tsv",
    "94C_selected_mobile_context_VFDB.tsv",
    "94C_direct_mobile_marker_overlap_VFDB.tsv",
    "94C_plasmid_context_VFDB.tsv",
    "94C_mobile_supported_GI_VFDB.tsv",
    "94C_AMR_within5kb_VFDB.tsv",
    "94C_summary_by_coassembly.tsv",
    "94C_summary_by_gene.tsv",
    "94C_global_summary.tsv",
    "94C_methodological_scope.tsv",
]

for name in required_outputs:

    if not (
        OUT / name
    ).is_file():

        raise RuntimeError(
            f"Salida faltante: {name}"
        )


print(
    f"VFDB_HIGH90_90={len(integrated)}"
)

print(
    f"SELECTED_MOBILE_CONTEXT={len(selected_mobile)}"
)

print(
    f"DIRECT_MOBILE_OVERLAP={len(direct_mobile)}"
)

print(
    f"PLASMID_CONTEXT={len(plasmid_context)}"
)

print(
    f"MOBILE_GI={len(mobile_GI)}"
)

print(
    f"AMR_WITHIN5KB={len(amr5)}"
)

print(
    "94C_INTEGRATION=PASS"
)

