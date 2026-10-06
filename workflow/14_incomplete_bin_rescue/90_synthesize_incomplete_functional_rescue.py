#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict, Counter
import csv
import os


# ============================================================
# PATHS
# ============================================================

USER = os.environ["USER"]

ROOT = Path(
    f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"
)

FINAL_MET = (
    ROOT
    / "55_cheese_metabolic_capabilities"
    / "metabolic_capability_matrix_18MAGs.tsv"
)

FINAL_BA = (
    ROOT
    / "55_cheese_metabolic_capabilities"
    / "curated_biogenic_amine_matrix_18MAGs.tsv"
)

FINAL_PROT = (
    ROOT
    / "58_surface_proteinase_curated"
    / "proteolytic_capability_matrix_18MAGs_78c.tsv"
)

FINAL_REMAIN = (
    ROOT
    / "60_remaining_cheese_function_curated"
    / "remaining_cheese_functions_curated_matrix_18MAGs.tsv"
)

FINAL_MASTER = (
    ROOT
    / "61_MAG_functional_master"
    / "MAG_functional_master_18MAGs.tsv"
)

IN88 = (
    ROOT
    / "74_incomplete_functional_evidence_context"
)

IN89 = (
    ROOT
    / "75_incomplete_functional_sequence_redundancy"
)

OUT = (
    ROOT
    / "76_incomplete_functional_synthesis"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# CAPABILITIES
# ============================================================

CAPABILITIES = [
    "galactose_Leloir",
    "lactose_transport_betaGal",
    "lactose_PTS_LacEFG",
    "lactate_LDH",
    "acetate_PtaAckA",
    "formate_PFL",
    "citrate_fermentation",
    "acetoin_branch",
    "butanediol_branch",

    "histamine",
    "tyramine",
    "cadaverine",
    "putrescine",

    "Opp_transport",
    "Dpp_transport",
    "POT_Dtp",
    "key_peptidases",
    "CEP_like",

    "lipolysis",
    "amino_acid_aroma",
    "EPS_capsule",

    "acid_stress",
    "osmotic_stress",
    "oxidative_stress",
]

CAP_ORDER = {
    c: i
    for i, c in enumerate(CAPABILITIES)
}


# ============================================================
# HELPERS
# ============================================================

def clean(x):
    if x is None:
        return ""

    x = str(x).strip()

    if x in {
        "",
        "-",
        "NA",
        "N/A",
        "None",
        "nan",
    }:
        return ""

    return x


def intval(x):
    try:
        return int(float(clean(x)))
    except Exception:
        return 0


def split_ids(x):
    return [
        y.strip()
        for y in clean(x).split(";")
        if y.strip()
    ]


def read_tsv(path):
    with path.open(
        newline="",
        encoding="utf-8"
    ) as fh:
        yield from csv.DictReader(
            fh,
            delimiter="\t"
        )


def write_tsv(path, rows, fields):
    with path.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as fh:

        w = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=fields,
            lineterminator="\n",
            extrasaction="ignore",
        )

        w.writeheader()

        for row in rows:
            w.writerow(row)


def context_level(cls):

    if cls in {
        "B50_89_contam_le5",
        "B50_89_contam_gt5_le10",
    }:
        return "good_context"

    if cls == "B50_89_contam_gt10_le20":
        return "cautious_context"

    if cls == "B50_89_contam_gt20":
        return "high_contamination_context"

    if cls == "lt50_contig_gene_centric_only":
        return "fragmentary_context"

    return "unknown_context"


def integrity_complete(value):

    return clean(value) in {
        "complete_predicted_CDS",
        "all_component_CDS_complete",
    }


# ============================================================
# VALIDATE INPUTS
# ============================================================

required = [
    FINAL_MET,
    FINAL_BA,
    FINAL_PROT,
    FINAL_REMAIN,
    FINAL_MASTER,

    IN88 / "contextualized_87B_cheese_metabolic_systems_on_contigs.tsv",
    IN88 / "contextualized_87C_biogenic_amine_candidates.tsv",
    IN88 / "contextualized_87F_Opp_Dpp_clusters.tsv",
    IN88 / "contextualized_87G_POT_Dtp_transporters.tsv",
    IN88 / "contextualized_87H_key_peptidases.tsv",
    IN88 / "contextualized_87J_curated_lipolysis_candidates.tsv",
    IN88 / "contextualized_87K_curated_amino_acid_aroma_markers.tsv",
    IN88 / "contextualized_87L_EPS_capsule_loci.tsv",
    IN88 / "contextualized_87M_stress_systems_on_contigs.tsv",

    IN89 / "89E_rescue_near_identity_cluster_members.tsv",
    IN89 / "89K_best_final18_protein_match_corrected.tsv",
    IN89 / "89L_priority_CEP_sequence_redundancy_corrected.tsv",
]

for p in required:
    if not p.exists() or p.stat().st_size == 0:
        raise RuntimeError(
            f"Falta input requerido o está vacío: {p}"
        )


# ============================================================
# FINAL18: MAG -> PRODUCER
# ============================================================

mag_producer = {}

for row in read_tsv(FINAL_MASTER):

    mag = clean(row.get("MAG"))

    if not mag:
        continue

    # El productor se infiere del identificador del MAG,
    # no de representative_source.
    #
    # Ejemplos válidos:
    #   L1__L1_bin.51  -> L
    #   L2__...        -> L
    #   M1__...        -> M
    #   M3__...        -> M

    if mag.startswith(("L1__", "L2__", "L3__")):
        producer = "L"

    elif mag.startswith(("M1__", "M2__", "M3__")):
        producer = "M"

    else:
        raise RuntimeError(
            f"No se puede resolver productor desde MAG={mag!r}"
        )

    mag_producer[mag] = producer


if len(mag_producer) != 18:
    raise RuntimeError(
        "Se esperaban 18 MAGs en el master final18; "
        f"se encontraron {len(mag_producer)}."
    )


# ============================================================
# FINAL18: STATES
# strong / partial / absent
# ============================================================

final_state = defaultdict(dict)


def register(cap, mag, state):

    if cap not in CAP_ORDER:
        raise RuntimeError(
            f"Capability no reconocida: {cap}"
        )

    if state not in {
        "strong",
        "partial",
        "absent",
    }:
        raise RuntimeError(
            f"Estado inválido: {state}"
        )

    if mag in final_state[cap]:
        raise RuntimeError(
            f"Estado duplicado: {cap} / {mag}"
        )

    final_state[cap][mag] = state


# ------------------------------------------------------------
# Metabolic capabilities
# ------------------------------------------------------------

met_strong = {
    "galactose_Leloir": {
        "complete_core",
    },

    "lactose_transport_betaGal": {
        "supported",
    },

    "lactose_PTS_LacEFG": {
        "complete_marker_set",
    },

    "lactate_LDH": {
        "supported",
    },

    "acetate_PtaAckA": {
        "complete_marker_set",
    },

    "formate_PFL": {
        "complete_marker_set",
    },

    "citrate_fermentation": {
        "strong_system",
    },

    "acetoin_branch": {
        "supported_marker_pair",
    },

    "butanediol_branch": {
        "supported",
    },
}

met_partial = {
    "galactose_Leloir": {
        "partial",
    },

    "lactose_transport_betaGal": {
        "partial",
    },

    "lactose_PTS_LacEFG": {
        "partial",
    },

    "lactate_LDH": set(),

    "acetate_PtaAckA": {
        "partial",
    },

    "formate_PFL": set(),

    "citrate_fermentation": {
        "partial",
    },

    "acetoin_branch": {
        "partial",
    },

    "butanediol_branch": set(),
}

for row in read_tsv(FINAL_MET):

    mag = row["MAG"]

    for cap in met_strong:

        val = clean(row.get(cap))

        if val in met_strong[cap]:
            state = "strong"

        elif val in met_partial[cap]:
            state = "partial"

        elif val in {
            "absent",
            "not_detected",
        }:
            state = "absent"

        else:
            raise RuntimeError(
                f"Estado final18 inesperado: "
                f"{cap} / {mag} / {val!r}"
            )

        register(
            cap,
            mag,
            state
        )


# ------------------------------------------------------------
# Biogenic amines
# ------------------------------------------------------------

for row in read_tsv(FINAL_BA):

    mag = row["MAG"]

    for cap in [
        "histamine",
        "tyramine",
        "cadaverine",
        "putrescine",
    ]:

        val = clean(row.get(cap))

        if val in {
            "system_supported",
            "enzyme_plus_transport_context",
        }:
            state = "strong"

        elif val == "enzyme_only":
            state = "partial"

        elif val == "not_detected":
            state = "absent"

        else:
            raise RuntimeError(
                f"Estado BA inesperado: "
                f"{cap} / {mag} / {val!r}"
            )

        register(
            cap,
            mag,
            state
        )


# ------------------------------------------------------------
# Proteolysis
# ------------------------------------------------------------

for row in read_tsv(FINAL_PROT):

    mag = row["MAG"]

    if intval(
        row.get("Opp_complete_clusters")
    ) > 0:
        state = "strong"

    elif intval(
        row.get("Opp_partial_clusters")
    ) > 0:
        state = "partial"

    else:
        state = "absent"

    register(
        "Opp_transport",
        mag,
        state
    )

    if intval(
        row.get("Dpp_complete_clusters")
    ) > 0:
        state = "strong"

    elif intval(
        row.get("Dpp_partial_clusters")
    ) > 0:
        state = "partial"

    else:
        state = "absent"

    register(
        "Dpp_transport",
        mag,
        state
    )

    register(
        "POT_Dtp",
        mag,
        (
            "strong"
            if intval(
                row.get("POT_Dtp_genes")
            ) > 0
            else "absent"
        )
    )

    register(
        "key_peptidases",
        mag,
        (
            "strong"
            if intval(
                row.get("n_key_peptidase_types")
            ) > 0
            else "absent"
        )
    )

    if intval(
        row.get(
            "CEP_A_high_confidence_n"
        )
    ) > 0:
        state = "strong"

    elif intval(
        row.get(
            "CEP_B_probable_n"
        )
    ) > 0:
        state = "partial"

    else:
        state = "absent"

    register(
        "CEP_like",
        mag,
        state
    )


# ------------------------------------------------------------
# Remaining functions
# ------------------------------------------------------------

for row in read_tsv(FINAL_REMAIN):

    mag = row["MAG"]

    # Lipolysis
    val = clean(
        row.get(
            "lipolysis_evidence_curated"
        )
    )

    if val in {
        "specific_lipase_candidate_recovered",
        "esterase_lipase_candidate_recovered",
    }:
        state = "strong"

    elif val == "broad_GDSL_hydrolase_only":
        state = "partial"

    elif val == "no_curated_lipase_candidate":
        state = "absent"

    else:
        raise RuntimeError(
            f"Estado lipolysis inesperado: "
            f"{mag} / {val!r}"
        )

    register(
        "lipolysis",
        mag,
        state
    )

    # Amino acid aroma
    val = clean(
        row.get(
            "amino_acid_aroma_evidence_curated"
        )
    )

    if val == "not_detected":
        state = "absent"

    elif val:
        state = "strong"

    else:
        raise RuntimeError(
            f"Estado aroma vacío para {mag}"
        )

    register(
        "amino_acid_aroma",
        mag,
        state
    )

    # EPS/capsule
    val = clean(
        row.get(
            "EPS_capsule_evidence_curated"
        )
    )

    if val in {
        "mixed_EPS_capsule_like_locus_recovered",
        "capsule_biosynthesis_like_locus_recovered",
        "separate_EPS_like_and_capsule_like_loci",
        "EPS_biosynthesis_like_locus_recovered",
    }:
        state = "strong"

    elif val == "single_marker_context_only":
        state = "partial"

    elif val == "no_resolved_EPS_capsule_locus":
        state = "absent"

    else:
        raise RuntimeError(
            f"Estado EPS inesperado: "
            f"{mag} / {val!r}"
        )

    register(
        "EPS_capsule",
        mag,
        state
    )

    # Acid stress
    val = clean(
        row.get(
            "acid_stress_evidence"
        )
    )

    if val == "system_supported":
        state = "strong"

    elif val == "marker_only":
        state = "partial"

    elif val == "not_detected":
        state = "absent"

    else:
        raise RuntimeError(
            f"Estado acid stress inesperado: "
            f"{mag} / {val!r}"
        )

    register(
        "acid_stress",
        mag,
        state
    )

    # Osmotic stress
    val = clean(
        row.get(
            "osmotic_stress_evidence_curated"
        )
    )

    if val == "strong_multigene_osmoadaptation_system":
        state = "strong"

    elif val == "partial_or_single_component_only":
        state = "partial"

    else:
        raise RuntimeError(
            f"Estado osmótico inesperado: "
            f"{mag} / {val!r}"
        )

    register(
        "osmotic_stress",
        mag,
        state
    )

    # Oxidative stress
    val = clean(
        row.get(
            "oxidative_stress_evidence"
        )
    )

    if val in {
        "broad_antioxidant_repertoire",
        "supported",
    }:
        state = "strong"

    else:
        raise RuntimeError(
            f"Estado oxidative inesperado: "
            f"{mag} / {val!r}"
        )

    register(
        "oxidative_stress",
        mag,
        state
    )


# ============================================================
# FINAL18 COMPLETENESS QC
# ============================================================

for cap in CAPABILITIES:

    if len(
        final_state[cap]
    ) != 18:
        raise RuntimeError(
            f"{cap}: se esperaban 18 MAGs, "
            f"se obtuvieron {len(final_state[cap])}"
        )


# ============================================================
# SEQUENCE CLUSTERS + CORRECTED FINAL18 MATCH
# ============================================================

gene_cluster = {}

for row in read_tsv(
    IN89
    / "89E_rescue_near_identity_cluster_members.tsv"
):
    gene_cluster[
        row["gene_id"]
    ] = row["cluster_id"]


gene_match = {}

for row in read_tsv(
    IN89
    / "89K_best_final18_protein_match_corrected.tsv"
):
    gene_match[
        row["gene_id"]
    ] = row["match_class"]


MATCH_RANK = {
    "A_exact_full_length": 7,
    "B_near_identical_full_length": 6,
    "C_near_identical_query_fragment": 5,
    "D_close_homolog": 4,
    "E_homologous_or_domain_level": 3,
    "F_weak_similarity": 2,
    "G_no_detected_similarity": 1,
}


# ============================================================
# RESCUE EVIDENCE RECORDS
# ============================================================

rescue_records = []


def add_record(
    capability,
    source_table,
    row,
    functional_strength,
    genes,
    descriptor,
):

    if capability not in CAP_ORDER:
        raise RuntimeError(
            f"Capability rescue desconocida: {capability}"
        )

    genes = sorted(
        set(
            g
            for g in genes
            if g
        )
    )

    ctx = clean(
        row.get(
            "context_support_class"
        )
    )

    integ = clean(
        row.get(
            "prediction_integrity"
        )
    )

    producer = clean(
        row.get(
            "best_support_bin_producer"
        )
    )

    if producer not in {
        "L",
        "M",
    }:
        raise RuntimeError(
            f"Productor rescue no resuelto: "
            f"{capability} / {producer!r}"
        )

    rec = {
        "capability":
            capability,

        "source_table":
            source_table,

        "producer":
            producer,

        "contig_id":
            clean(
                row.get(
                    "contig_id"
                )
            ),

        "functional_strength":
            functional_strength,

        "context_support_class":
            ctx,

        "context_level":
            context_level(ctx),

        "prediction_integrity":
            integ,

        "descriptor":
            descriptor,

        "genes":
            ";".join(genes),

        "n_genes":
            len(genes),

        "high_confidence_rescue":
            int(
                functional_strength == "strong"
                and context_level(ctx) == "good_context"
                and integrity_complete(integ)
            ),
    }

    rescue_records.append(
        rec
    )


# ------------------------------------------------------------
# Metabolic systems
# ------------------------------------------------------------

MET_FILE = (
    IN88
    / "contextualized_87B_cheese_metabolic_systems_on_contigs.tsv"
)

for row in read_tsv(MET_FILE):

    cap = clean(
        row.get("system")
    )

    if cap not in met_strong:
        raise RuntimeError(
            f"Sistema metabólico rescue no esperado: {cap}"
        )

    status = clean(
        row.get("status")
    )

    functional_strength = (
        "partial"
        if "partial" in status.lower()
        else "strong"
    )

    add_record(
        capability=cap,
        source_table=MET_FILE.name,
        row=row,
        functional_strength=functional_strength,
        genes=split_ids(
            row.get("gene_ids")
        ),
        descriptor=status,
    )


# ------------------------------------------------------------
# Biogenic amines
# ------------------------------------------------------------

BA_FILE = (
    IN88
    / "contextualized_87C_biogenic_amine_candidates.tsv"
)

for row in read_tsv(BA_FILE):

    cap = clean(
        row.get(
            "potential_product"
        )
    )

    if cap not in {
        "histamine",
        "tyramine",
        "cadaverine",
        "putrescine",
    }:
        raise RuntimeError(
            f"Producto BA no esperado: {cap}"
        )

    level = clean(
        row.get(
            "evidence_level"
        )
    )

    if level in {
        "system_supported",
        "enzyme_plus_transport_context",
    }:
        strength = "strong"

    elif level == "enzyme_only":
        strength = "partial"

    else:
        strength = "weak"

    add_record(
        capability=cap,
        source_table=BA_FILE.name,
        row=row,
        functional_strength=strength,
        genes=[
            clean(
                row.get(
                    "gene_id"
                )
            )
        ],
        descriptor=(
            level
            + (
                ":"
                + clean(
                    row.get(
                        "context_support"
                    )
                )
                if clean(
                    row.get(
                        "context_support"
                    )
                )
                else ""
            )
        ),
    )


# ------------------------------------------------------------
# Opp / Dpp
# ------------------------------------------------------------

OPP_FILE = (
    IN88
    / "contextualized_87F_Opp_Dpp_clusters.tsv"
)

for row in read_tsv(OPP_FILE):

    sysname = clean(
        row.get(
            "system"
        )
    ).upper()

    if sysname == "OPP":
        cap = "Opp_transport"

    elif sysname == "DPP":
        cap = "Dpp_transport"

    else:
        raise RuntimeError(
            f"Sistema Opp/Dpp inesperado: {sysname}"
        )

    status = clean(
        row.get(
            "cluster_status"
        )
    )

    if status == "complete_named_cluster":
        strength = "strong"

    elif status == "partial_named_cluster":
        strength = "partial"

    else:
        strength = "weak"

    add_record(
        capability=cap,
        source_table=OPP_FILE.name,
        row=row,
        functional_strength=strength,
        genes=split_ids(
            row.get(
                "genes"
            )
        ),
        descriptor=status,
    )


# ------------------------------------------------------------
# POT / Dtp
# ------------------------------------------------------------

POT_FILE = (
    IN88
    / "contextualized_87G_POT_Dtp_transporters.tsv"
)

for row in read_tsv(POT_FILE):

    add_record(
        capability="POT_Dtp",
        source_table=POT_FILE.name,
        row=row,
        functional_strength="strong",
        genes=[
            clean(
                row.get(
                    "gene_id"
                )
            )
        ],
        descriptor=clean(
            row.get(
                "marker"
            )
        ),
    )


# ------------------------------------------------------------
# Key peptidases
# ------------------------------------------------------------

PEP_FILE = (
    IN88
    / "contextualized_87H_key_peptidases.tsv"
)

for row in read_tsv(PEP_FILE):

    add_record(
        capability="key_peptidases",
        source_table=PEP_FILE.name,
        row=row,
        functional_strength="strong",
        genes=[
            clean(
                row.get(
                    "gene_id"
                )
            )
        ],
        descriptor=clean(
            row.get(
                "peptidase_type"
            )
        ),
    )


# ------------------------------------------------------------
# CEP: use only corrected prioritized CEPs
# ------------------------------------------------------------

CEP_FILE = (
    IN89
    / "89L_priority_CEP_sequence_redundancy_corrected.tsv"
)

for row in read_tsv(CEP_FILE):

    tier = clean(
        row.get(
            "rescue_tier_87"
        )
    )

    if tier == (
        "A_high_confidence_CEP_like_complete_prediction"
    ):
        strength = "strong"

    elif tier == "T_truncated_CEP_like_candidate":
        strength = "partial"

    else:
        raise RuntimeError(
            f"Tier CEP inesperado en tabla prioritaria: {tier}"
        )

    add_record(
        capability="CEP_like",
        source_table=CEP_FILE.name,
        row=row,
        functional_strength=strength,
        genes=[
            clean(
                row.get(
                    "gene_id"
                )
            )
        ],
        descriptor=tier,
    )


# ------------------------------------------------------------
# Lipolysis
# ------------------------------------------------------------

LIP_FILE = (
    IN88
    / "contextualized_87J_curated_lipolysis_candidates.tsv"
)

for row in read_tsv(LIP_FILE):

    tier = clean(
        row.get(
            "curated_lipolysis_tier"
        )
    )

    if tier.startswith(
        ("A_", "B_")
    ):
        strength = "strong"

    elif tier.startswith("C_"):
        strength = "partial"

    else:
        strength = "weak"

    add_record(
        capability="lipolysis",
        source_table=LIP_FILE.name,
        row=row,
        functional_strength=strength,
        genes=[
            clean(
                row.get(
                    "gene_id"
                )
            )
        ],
        descriptor=tier,
    )


# ------------------------------------------------------------
# Aroma
# ------------------------------------------------------------

AROMA_FILE = (
    IN88
    / "contextualized_87K_curated_amino_acid_aroma_markers.tsv"
)

for row in read_tsv(AROMA_FILE):

    tier = clean(
        row.get(
            "curated_aroma_tier"
        )
    )

    discordant = (
        clean(
            row.get(
                "annotation_discordance"
            )
        ) == "1"
    )

    if tier.startswith(
        ("A_", "B_")
    ):
        strength = "strong"

    elif tier.startswith("C_"):
        strength = "partial"

    else:
        strength = "weak"

    if discordant:

        if strength == "strong":
            strength = "partial"

        elif strength == "partial":
            strength = "weak"

    add_record(
        capability="amino_acid_aroma",
        source_table=AROMA_FILE.name,
        row=row,
        functional_strength=strength,
        genes=[
            clean(
                row.get(
                    "gene_id"
                )
            )
        ],
        descriptor=tier,
    )


# ------------------------------------------------------------
# EPS/capsule
# ------------------------------------------------------------

EPS_FILE = (
    IN88
    / "contextualized_87L_EPS_capsule_loci.tsv"
)

for row in read_tsv(EPS_FILE):

    tier = clean(
        row.get(
            "curated_tier"
        )
    )

    if tier.startswith("A_"):
        strength = "strong"

    elif tier.startswith("B_"):
        strength = "partial"

    else:
        strength = "weak"

    genes = (
        split_ids(
            row.get(
                "marker_gene_ids"
            )
        )
        +
        split_ids(
            row.get(
                "glycosyltransferase_gene_ids"
            )
        )
    )

    add_record(
        capability="EPS_capsule",
        source_table=EPS_FILE.name,
        row=row,
        functional_strength=strength,
        genes=genes,
        descriptor=tier,
    )


# ------------------------------------------------------------
# Stress
# ------------------------------------------------------------

STRESS_FILE = (
    IN88
    / "contextualized_87M_stress_systems_on_contigs.tsv"
)

for row in read_tsv(STRESS_FILE):

    block = clean(
        row.get(
            "block"
        )
    )

    if block == "acid_stress":
        cap = "acid_stress"

    elif block == "osmotic_stress":
        cap = "osmotic_stress"

    elif block == "oxidative_stress":
        cap = "oxidative_stress"

    else:
        raise RuntimeError(
            f"Bloque stress inesperado: {block}"
        )

    status = clean(
        row.get(
            "status"
        )
    )

    low = status.lower()

    if cap == "acid_stress":

        strength = (
            "strong"
            if "system_supported" in low
            else "partial"
        )

    elif cap == "osmotic_stress":

        strength = (
            "strong"
            if (
                "cluster" in low
                or "complete" in low
                or "system_supported" in low
            )
            else "partial"
        )

    else:
        # Oxidative repertoire:
        # >=2 marker types = stronger multicomponent support.
        strength = (
            "strong"
            if intval(
                row.get(
                    "n_marker_types_detected"
                )
            ) >= 2
            else "partial"
        )

    add_record(
        capability=cap,
        source_table=STRESS_FILE.name,
        row=row,
        functional_strength=strength,
        genes=split_ids(
            row.get(
                "gene_ids"
            )
        ),
        descriptor=(
            clean(
                row.get(
                    "system"
                )
            )
            + ":"
            + status
        ),
    )


# ============================================================
# QC RESCUE COVERAGE
# ============================================================

caps_with_rescue = {
    x["capability"]
    for x in rescue_records
}

missing_rescue_caps = [
    c
    for c in CAPABILITIES
    if c not in caps_with_rescue
]

if missing_rescue_caps:
    print(
        "WARNING: capacidades sin evidencia rescue: "
        + ";".join(missing_rescue_caps)
    )


# ============================================================
# MAP RESCUE GENES -> CLUSTERS / FINAL18 MATCH CLASSES
# ============================================================

cap_genes = defaultdict(set)

for rec in rescue_records:

    for gene in split_ids(
        rec["genes"]
    ):
        cap_genes[
            rec["capability"]
        ].add(
            gene
        )


cap_cluster_genes = defaultdict(
    lambda: defaultdict(set)
)

for cap, genes in cap_genes.items():

    for gene in genes:

        if gene not in gene_cluster:
            raise RuntimeError(
                f"Gene rescue sin cluster Paso 89: {gene}"
            )

        if gene not in gene_match:
            raise RuntimeError(
                f"Gene rescue sin match corregido Paso 89c: {gene}"
            )

        cid = gene_cluster[gene]

        cap_cluster_genes[
            cap
        ][
            cid
        ].add(
            gene
        )


cap_cluster_class = defaultdict(dict)

for cap, clusters in cap_cluster_genes.items():

    for cid, genes in clusters.items():

        classes = [
            gene_match[g]
            for g in genes
        ]

        best = max(
            classes,
            key=lambda x: MATCH_RANK[x]
        )

        cap_cluster_class[
            cap
        ][
            cid
        ] = best


# ============================================================
# FINAL18 GLOBAL / PRODUCER SUMMARIES
# ============================================================

def final_counts(cap, producer=None):

    states = []

    for mag, state in final_state[
        cap
    ].items():

        if (
            producer is not None
            and mag_producer[mag] != producer
        ):
            continue

        states.append(
            state
        )

    c = Counter(states)

    return {
        "strong":
            c.get("strong", 0),

        "partial":
            c.get("partial", 0),

        "absent":
            c.get("absent", 0),

        "n_MAGs":
            len(states),
    }


def representation_from_counts(counts):

    if counts["strong"] > 0:
        return "represented_strong"

    if counts["partial"] > 0:
        return "represented_partial_only"

    return "not_recovered"


# ============================================================
# GLOBAL CAPABILITY TABLE
# ============================================================

global_rows = []


for cap in CAPABILITIES:

    fc = final_counts(
        cap
    )

    final_rep = representation_from_counts(
        fc
    )

    rr = [
        x
        for x in rescue_records
        if x["capability"] == cap
    ]

    n_records = len(rr)

    n_strong_func = sum(
        x[
            "functional_strength"
        ] == "strong"
        for x in rr
    )

    n_high = sum(
        int(
            x[
                "high_confidence_rescue"
            ]
        )
        for x in rr
    )

    n_good_context = sum(
        x["context_level"]
        == "good_context"
        for x in rr
    )

    n_cautious_context = sum(
        x["context_level"]
        == "cautious_context"
        for x in rr
    )

    n_weak_context = sum(
        x["context_level"]
        in {
            "high_contamination_context",
            "fragmentary_context",
        }
        for x in rr
    )

    contigs = {
        x["contig_id"]
        for x in rr
        if x["contig_id"]
    }

    clusters = cap_cluster_class[
        cap
    ]

    cluster_counts = Counter(
        clusters.values()
    )

    n_clusters = len(
        clusters
    )

    n_AB = (
        cluster_counts.get(
            "A_exact_full_length",
            0
        )
        +
        cluster_counts.get(
            "B_near_identical_full_length",
            0
        )
    )

    n_not_AB = (
        n_clusters
        - n_AB
    )

    n_without_close = sum(
        cluster_counts.get(
            cls,
            0
        )
        for cls in {
            "E_homologous_or_domain_level",
            "F_weak_similarity",
            "G_no_detected_similarity",
        }
    )

    # Diversidad de secuencia respecto de final18:
    # C es un fragmento casi idéntico y NO debe contarse como
    # diversidad proteica adicional.
    n_divergent_D_to_G = (
        cluster_counts.get(
            "D_close_homolog",
            0
        )
        + n_without_close
    )

    if n_high > 0:
        rescue_level = (
            "high_confidence_positive_rescue"
        )

    elif n_strong_func > 0:
        rescue_level = (
            "positive_rescue_context_or_integrity_limited"
        )

    elif n_records > 0:
        rescue_level = (
            "screening_or_partial_rescue_evidence"
        )

    else:
        rescue_level = (
            "no_rescue_evidence"
        )

    if final_rep == "not_recovered":

        if n_high > 0:
            impact = (
                "candidate_additional_capability_"
                "not_recovered_in_final18"
            )

            changes_global = (
                "candidate_change_requires_validation"
            )

        elif n_records > 0:
            impact = (
                "limited_signal_for_capability_"
                "not_recovered_in_final18"
            )

            changes_global = (
                "no_change_insufficient_rescue_strength"
            )

        else:
            impact = "absent_in_both"
            changes_global = "no"

    elif final_rep == "represented_partial_only":

        if n_high > 0:
            impact = (
                "rescue_strengthens_globally_"
                "partial_final18_capability"
            )

            changes_global = (
                "strengthens_existing_partial_conclusion"
            )

        else:
            impact = (
                "supplementary_evidence_for_"
                "partial_final18_capability"
            )

            changes_global = "no_new_category"

    else:

        if (
            n_divergent_D_to_G > 0
            and n_records > 0
        ):
            impact = (
                "known_capability_with_"
                "additional_sequence_diversity"
            )

        elif (
            cluster_counts.get(
                "C_near_identical_query_fragment",
                0
            ) > 0
            and n_records > 0
        ):
            impact = (
                "known_capability_with_"
                "fragmentary_near_identical_rescue"
            )

        elif n_records > 0:
            impact = (
                "known_capability_recovered_"
                "redundantly"
            )

        else:
            impact = (
                "known_final18_capability_"
                "without_rescue_signal"
            )

        changes_global = (
            "no_new_global_capability_category"
        )

    global_rows.append({
        "capability":
            cap,

        "final18_n_MAGs":
            fc["n_MAGs"],

        "final18_strong_MAGs":
            fc["strong"],

        "final18_partial_MAGs":
            fc["partial"],

        "final18_absent_MAGs":
            fc["absent"],

        "final18_global_representation":
            final_rep,

        "rescue_evidence_records":
            n_records,

        "rescue_strong_function_records":
            n_strong_func,

        "rescue_good_context_records":
            n_good_context,

        "rescue_cautious_context_records":
            n_cautious_context,

        "rescue_high_or_fragmentary_context_records":
            n_weak_context,

        "rescue_high_confidence_records":
            n_high,

        "rescue_unique_contigs":
            len(contigs),

        "rescue_supporting_genes":
            len(
                cap_genes[
                    cap
                ]
            ),

        "rescue_sequence_clusters":
            n_clusters,

        "clusters_exact_final18_A":
            cluster_counts.get(
                "A_exact_full_length",
                0
            ),

        "clusters_near_identical_final18_B":
            cluster_counts.get(
                "B_near_identical_full_length",
                0
            ),

        "clusters_fragment_near_identical_C":
            cluster_counts.get(
                "C_near_identical_query_fragment",
                0
            ),

        "clusters_close_homolog_D":
            cluster_counts.get(
                "D_close_homolog",
                0
            ),

        "clusters_homologous_domain_E":
            cluster_counts.get(
                "E_homologous_or_domain_level",
                0
            ),

        "clusters_weak_similarity_F":
            cluster_counts.get(
                "F_weak_similarity",
                0
            ),

        "clusters_no_detected_similarity_G":
            cluster_counts.get(
                "G_no_detected_similarity",
                0
            ),

        "clusters_not_exact_or_near_identical_AB":
            n_not_AB,

        "clusters_divergent_from_final18_D_to_G":
            n_divergent_D_to_G,

        "clusters_without_close_final18_match_EFG":
            n_without_close,

        "rescue_evidence_level":
            rescue_level,

        "integrated_interpretation":
            impact,

        "changes_global_functional_repertoire":
            changes_global,
    })


# ============================================================
# PRODUCER-LEVEL COMPARISON
# ============================================================

producer_rows = []


for cap in CAPABILITIES:

    for producer in [
        "L",
        "M",
    ]:

        fc = final_counts(
            cap,
            producer
        )

        final_rep = representation_from_counts(
            fc
        )

        rr = [
            x
            for x in rescue_records
            if (
                x["capability"] == cap
                and x["producer"] == producer
            )
        ]

        n_high = sum(
            int(
                x[
                    "high_confidence_rescue"
                ]
            )
            for x in rr
        )

        n_strong = sum(
            x[
                "functional_strength"
            ] == "strong"
            for x in rr
        )

        if final_rep == "not_recovered":

            if n_high > 0:
                impact = (
                    "candidate_new_producer_level_"
                    "capability_from_rescue"
                )

            elif rr:
                impact = (
                    "limited_rescue_signal_in_"
                    "producer_without_final18_support"
                )

            else:
                impact = (
                    "not_recovered_in_final18_or_rescue"
                )

        elif final_rep == "represented_partial_only":

            if n_high > 0:
                impact = (
                    "rescue_strengthens_producer_level_"
                    "partial_final18_evidence"
                )

            elif rr:
                impact = (
                    "supplementary_rescue_evidence_"
                    "for_partial_final18_state"
                )

            else:
                impact = (
                    "partial_final18_only"
                )

        else:

            if rr:
                impact = (
                    "supplements_known_producer_"
                    "capability"
                )

            else:
                impact = (
                    "represented_in_final18_no_"
                    "additional_rescue_signal"
                )

        producer_rows.append({
            "producer":
                producer,

            "capability":
                cap,

            "final18_n_MAGs":
                fc["n_MAGs"],

            "final18_strong_MAGs":
                fc["strong"],

            "final18_partial_MAGs":
                fc["partial"],

            "final18_absent_MAGs":
                fc["absent"],

            "final18_producer_representation":
                final_rep,

            "rescue_records":
                len(rr),

            "rescue_strong_function_records":
                n_strong,

            "rescue_high_confidence_records":
                n_high,

            "rescue_unique_contigs":
                len({
                    x["contig_id"]
                    for x in rr
                    if x["contig_id"]
                }),

            "producer_level_interpretation":
                impact,
        })


# ============================================================
# CLUSTER × CAPABILITY TABLE
# ============================================================

cluster_rows = []


for cap in CAPABILITIES:

    for cid in sorted(
        cap_cluster_genes[
            cap
        ]
    ):

        genes = sorted(
            cap_cluster_genes[
                cap
            ][
                cid
            ]
        )

        match_class = (
            cap_cluster_class[
                cap
            ][
                cid
            ]
        )

        producers = sorted({
            rec["producer"]
            for rec in rescue_records
            if (
                rec["capability"] == cap
                and any(
                    g in genes
                    for g in split_ids(
                        rec["genes"]
                    )
                )
            )
        })

        cluster_rows.append({
            "capability":
                cap,

            "cluster_id":
                cid,

            "n_supporting_genes":
                len(genes),

            "supporting_genes":
                ";".join(genes),

            "producers":
                ";".join(producers),

            "best_final18_match_class":
                match_class,

            "near_identical_to_final18":
                int(
                    match_class
                    in {
                        "A_exact_full_length",
                        "B_near_identical_full_length",
                    }
                ),

            "without_close_final18_match":
                int(
                    match_class
                    in {
                        "E_homologous_or_domain_level",
                        "F_weak_similarity",
                        "G_no_detected_similarity",
                    }
                ),
        })


# ============================================================
# HIGH-CONFIDENCE RESCUE RECORDS
# ============================================================

high_conf_rows = [
    row
    for row in rescue_records
    if int(
        row[
            "high_confidence_rescue"
        ]
    ) == 1
]


# ============================================================
# WRITE OUTPUTS
# ============================================================

write_tsv(
    OUT
    / "90A_global_capability_rescue_impact.tsv",
    global_rows,
    list(
        global_rows[0].keys()
    )
)

write_tsv(
    OUT
    / "90B_producer_capability_rescue_impact.tsv",
    producer_rows,
    list(
        producer_rows[0].keys()
    )
)

write_tsv(
    OUT
    / "90C_rescue_sequence_clusters_by_capability.tsv",
    cluster_rows,
    [
        "capability",
        "cluster_id",
        "n_supporting_genes",
        "supporting_genes",
        "producers",
        "best_final18_match_class",
        "near_identical_to_final18",
        "without_close_final18_match",
    ]
)

write_tsv(
    OUT
    / "90D_high_confidence_rescue_evidence.tsv",
    high_conf_rows,
    [
        "capability",
        "source_table",
        "producer",
        "contig_id",
        "functional_strength",
        "context_support_class",
        "context_level",
        "prediction_integrity",
        "descriptor",
        "genes",
        "n_genes",
        "high_confidence_rescue",
    ]
)


# ============================================================
# PRIORITY INTERPRETATION TABLE
# ============================================================

priority_rows = []

for row in global_rows:

    if row[
        "capability"
    ] in {
        "histamine",
        "tyramine",
        "cadaverine",
        "putrescine",
        "CEP_like",
    }:

        priority_rows.append(
            row
        )

write_tsv(
    OUT
    / "90E_priority_CEP_BA_integrated_summary.tsv",
    priority_rows,
    list(
        priority_rows[0].keys()
    )
)


# ============================================================
# GLOBAL SUMMARY
# ============================================================

n_final_strong_global = sum(
    row[
        "final18_global_representation"
    ] == "represented_strong"
    for row in global_rows
)

n_final_partial_global = sum(
    row[
        "final18_global_representation"
    ] == "represented_partial_only"
    for row in global_rows
)

n_final_absent_global = sum(
    row[
        "final18_global_representation"
    ] == "not_recovered"
    for row in global_rows
)

n_rescue_high_caps = sum(
    int(
        row[
            "rescue_high_confidence_records"
        ]
    ) > 0
    for row in global_rows
)

n_candidate_global_changes = sum(
    row[
        "changes_global_functional_repertoire"
    ] == "candidate_change_requires_validation"
    for row in global_rows
)

producer_new = [
    row
    for row in producer_rows
    if row[
        "producer_level_interpretation"
    ] == (
        "candidate_new_producer_level_"
        "capability_from_rescue"
    )
]

producer_strengthened = [
    row
    for row in producer_rows
    if row[
        "producer_level_interpretation"
    ] == (
        "rescue_strengthens_producer_level_"
        "partial_final18_evidence"
    )
]

summary = [
    (
        "capabilities_assessed",
        len(CAPABILITIES)
    ),

    (
        "final18_globally_strong_capabilities",
        n_final_strong_global
    ),

    (
        "final18_globally_partial_only_capabilities",
        n_final_partial_global
    ),

    (
        "final18_globally_not_recovered_capabilities",
        n_final_absent_global
    ),

    (
        "rescue_capabilities_with_any_evidence",
        len(caps_with_rescue)
    ),

    (
        "rescue_capabilities_with_high_confidence_evidence",
        n_rescue_high_caps
    ),

    (
        "high_confidence_rescue_records",
        len(high_conf_rows)
    ),

    (
        "candidate_new_global_capabilities_from_rescue",
        n_candidate_global_changes
    ),

    (
        "candidate_new_producer_level_capabilities",
        len(producer_new)
    ),

    (
        "producer_level_partial_capabilities_strengthened",
        len(producer_strengthened)
    ),

    (
        "capability_sequence_cluster_links",
        len(cluster_rows)
    ),
]


with (
    OUT
    / "functional_rescue_synthesis_summary.tsv"
).open(
    "w",
    newline="",
    encoding="utf-8"
) as fh:

    w = csv.writer(
        fh,
        delimiter="\t",
        lineterminator="\n"
    )

    w.writerow(
        [
            "metric",
            "value"
        ]
    )

    w.writerows(
        summary
    )


# ============================================================
# README
# ============================================================

readme = """\
PASO 90 - SINTESIS DEL RESCATE FUNCIONAL DE BINS INCOMPLETOS
==============================================================

OBJETIVO
--------
Determinar si la evidencia funcional recuperada de bins <90%:

1. añade una capacidad funcional no representada en los 18 MAGs finales;
2. fortalece una capacidad que en final18 sólo tenía evidencia parcial;
3. añade diversidad de secuencia a una capacidad ya conocida; o
4. constituye solamente evidencia positiva suplementaria de gen/contig.

UNIDAD DE INTERPRETACION
------------------------
No se cuentan genes como funciones independientes.

Se integran:
- estados funcionales curados de final18;
- evidencia funcional del Paso 87;
- calidad del contexto del Paso 88;
- clusters de proteínas del Paso 89;
- clases de similitud CORREGIDAS del Paso 89c.

Las clases antiguas de match de 89D/89H/89I/89J NO se utilizan para
interpretación final. Se usan 89K y 89L corregidos.

EVIDENCIA RESCUE DE MAYOR CONFIANZA
------------------------------------
Se requiere simultáneamente:

- regla funcional clasificada como strong;
- predicción proteica completa o sistema sin CDS parcial;
- contexto B50-89.99;
- contaminación del mejor bin <=10%.

Esto NO convierte al bin en un MAG ni autoriza atribución de especie.

FINAL18
-------
Cada capacidad se clasifica por MAG como:

- strong
- partial
- absent

La representación global es:

represented_strong:
    >=1 MAG final con evidencia strong.

represented_partial_only:
    ningún MAG strong, pero >=1 partial.

not_recovered:
    ningún MAG strong ni partial.

INTERPRETACION DE SECUENCIA
---------------------------
A = exact full length
B = near-identical full length
C = near-identical query fragment
D = close homolog
E = homologous/domain-level
F = weak similarity
G = no detected similarity

A/B se consideran representación casi idéntica en final18.

C representa recuperación fragmentaria casi idéntica y NO se interpreta
como diversidad proteica adicional.

D representa un homólogo cercano pero secuencialmente divergente.

E/F/G se contabilizan por separado como clusters sin un match cercano
a final18, pero esto NO demuestra automáticamente una nueva función.

IMPORTANTE
----------
Una capacidad funcional ya representada globalmente en final18 no se
etiqueta como "nueva" aunque aparezcan proteínas divergentes en el rescate.

Los resultados por productor L/M son análisis de soporte funcional por
origen del coensamblaje/bin. No equivalen a asignación taxonómica.

La ausencia en bins incompletos nunca se interpreta biológicamente.
"""

with (
    OUT
    / "README_step90.txt"
).open(
    "w",
    encoding="utf-8"
) as fh:

    fh.write(
        readme
    )


# ============================================================
# CONSOLE
# ============================================================

print("=" * 78)
print("PASO 90 - SINTESIS FUNCIONAL COMPLETADA")
print("=" * 78)

for key, value in summary:
    print(
        f"{key:<55} {value}"
    )

print()
print("CANDIDATOS A CAMBIO GLOBAL:")

if n_candidate_global_changes == 0:
    print(
        "  Ninguna capacidad del panel."
    )
else:
    for row in global_rows:
        if row[
            "changes_global_functional_repertoire"
        ] == "candidate_change_requires_validation":
            print(
                "  "
                + row["capability"]
            )

print()
print(
    "CANDIDATOS A CAMBIO A NIVEL DE PRODUCTOR:"
)

if not producer_new:
    print("  Ninguno.")
else:
    for row in producer_new:
        print(
            f"  {row['producer']} : "
            f"{row['capability']}"
        )

print()
print(
    "CAPACIDADES PARCIALES FINAL18 FORTALECIDAS "
    "A NIVEL DE PRODUCTOR:"
)

if not producer_strengthened:
    print("  Ninguna.")
else:
    for row in producer_strengthened:
        print(
            f"  {row['producer']} : "
            f"{row['capability']}"
        )

print()
print(
    f"Salida: {OUT}"
)

print()
print(
    "PASO 90 FINALIZÓ CORRECTAMENTE; "
    "ES SEGURO SALIR."
)
