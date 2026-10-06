#!/usr/bin/env python3

import os
import csv
from collections import defaultdict, Counter
from itertools import product

USER = os.environ.get("USER", "user")

ROOT = f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"

IN79A = os.path.join(
    ROOT,
    "59_remaining_cheese_function_audit"
)

MASTER_FILE = os.path.join(
    IN79A,
    "remaining_cheese_function_matrix_18MAGs.tsv"
)

HITS_FILE = os.path.join(
    IN79A,
    "strict_remaining_function_marker_hits.tsv"
)

EPS_FILE = os.path.join(
    IN79A,
    "EPS_capsule_like_locus_candidates.tsv"
)

OUT = os.path.join(
    ROOT,
    "60_remaining_cheese_function_curated"
)

os.makedirs(OUT, exist_ok=True)


# ============================================================
# UTILIDADES
# ============================================================

def read_tsv(path):
    rows = []
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            clean = {}
            for k, v in row.items():
                if k is None:
                    continue
                kk = k.strip().replace("\r", "")
                vv = "" if v is None else str(v).strip().replace("\r", "")
                clean[kk] = vv
            rows.append(clean)
    return rows


def write_tsv(path, rows, fieldnames):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def to_int(x):
    try:
        return int(float(str(x)))
    except Exception:
        return 0


def lower(x):
    return str(x or "").lower()


def contains_ko(row, ko):
    return ko.lower() in lower(row.get("KEGG_ko", ""))


def contains_ec(row, ec):
    vals = lower(row.get("EC", ""))
    return ec.lower() in vals


def gene_index(row):
    return to_int(row.get("gene_index", 0))


def semicolon_join(values):
    vals = sorted(set(v for v in values if v))
    return ";".join(vals)


# ============================================================
# CARGA
# ============================================================

master = read_tsv(MASTER_FILE)
hits = read_tsv(HITS_FILE)
eps_loci = read_tsv(EPS_FILE)

if len(master) != 18:
    raise RuntimeError(
        f"Se esperaban 18 MAGs en 79a y se encontraron {len(master)}"
    )

mag_order = [r["MAG"] for r in master]
master_by_mag = {r["MAG"]: r for r in master}

hits_by_mag = defaultdict(list)
for r in hits:
    hits_by_mag[r["MAG"]].append(r)


# ============================================================
# 1. CURACIÓN DE LIPÓLISIS
# ============================================================

lipase_rows = []

for r in hits:

    if r.get("block") != "lipolysis":
        continue

    if r.get("marker") != "explicit_lipase_candidate":
        continue

    desc = lower(r.get("Description"))
    pref = lower(r.get("Preferred_name"))
    pfam = lower(r.get("PFAMs"))
    ko = lower(r.get("KEGG_ko"))
    ec = lower(r.get("EC"))

    tier = ""
    reason = ""

    # --------------------------------------------------------
    # D: accesorio, no enzima lipolítica por sí mismo
    # --------------------------------------------------------
    if (
        pref == "lifo"
        or "lipase_chap" in pfam
        or "folding of the extracellular lipase" in desc
    ):
        tier = "D_lipase_accessory_not_enzyme"
        reason = "lipase_chaperone_or_folding_accessory"

    # --------------------------------------------------------
    # A: anotación lipasa relativamente específica
    # --------------------------------------------------------
    elif (
        "lipase (class 3)" in desc
        or "k01046" in ko
        or "3.1.1.3" in ec
        or "secretory lipase" in desc
    ):
        tier = "A_specific_lipase_candidate"

        reasons = []

        if "lipase (class 3)" in desc:
            reasons.append("class_3_lipase_annotation")

        if "k01046" in ko:
            reasons.append("K01046")

        if "3.1.1.3" in ec:
            reasons.append("EC_3.1.1.3")

        if "secretory lipase" in desc:
            reasons.append("explicit_secretory_lipase")

        reason = ";".join(reasons)

    # --------------------------------------------------------
    # B: esterase/lipase específica o arquitectura secretada
    # --------------------------------------------------------
    elif (
        "k12686" in ko
        or (
            "autotransporter" in pfam
            and "lipase_gdsl" in pfam
        )
        or "cog0657 esterase lipase" in desc
        or "esterase lipase" in desc
        or "k03928" in ko
        or "k03929" in ko
        or "k14731" in ko
    ):
        tier = "B_esterase_lipase_candidate"

        reasons = []

        if "k12686" in ko:
            reasons.append("K12686")

        if "autotransporter" in pfam:
            reasons.append("autotransporter_architecture")

        if "cog0657 esterase lipase" in desc:
            reasons.append("COG0657")

        if "esterase lipase" in desc:
            reasons.append("esterase_lipase_annotation")

        if "k03928" in ko:
            reasons.append("K03928")

        if "k03929" in ko:
            reasons.append("K03929")

        if "k14731" in ko:
            reasons.append("K14731")

        reason = ";".join(reasons)

    # --------------------------------------------------------
    # C: gran familia GDSL, demasiado amplia para inferir
    # hidrólisis de grasa láctea
    # --------------------------------------------------------
    elif (
        "gdsl" in desc
        or "lipase_gdsl" in pfam
    ):
        tier = "C_broad_GDSL_acylhydrolase"
        reason = "broad_GDSL_lipase_acylhydrolase_family"

    else:
        tier = "E_unresolved_lipolysis_candidate"
        reason = "insufficient_specificity"

    nr = dict(r)
    nr["curated_lipolysis_tier"] = tier
    nr["curation_reason"] = reason

    lipase_rows.append(nr)


lip_counts = defaultdict(Counter)

for r in lipase_rows:
    lip_counts[r["MAG"]][r["curated_lipolysis_tier"]] += 1


# ============================================================
# 2. CURACIÓN DE AROMA DERIVADO DE AMINOÁCIDOS
# ============================================================

aroma_rows = []

for r in hits:

    if r.get("block") != "amino_acid_aroma":
        continue

    marker = r.get("marker", "")
    pref = lower(r.get("Preferred_name"))
    desc = lower(r.get("Description"))
    ko = lower(r.get("KEGG_ko"))
    ec = lower(r.get("EC"))

    tier = ""
    interpretation = ""
    discordance = "0"

    # --------------------------------------------------------
    # BCAA aminotransferase
    # --------------------------------------------------------
    if marker == "branched_chain_aminotransferase":

        if "k00826" in ko:
            tier = "A_BCAT_IlvE_supported"
            interpretation = "branched_chain_amino_acid_transamination"

        elif "k02619" in ko or pref == "pabc":
            tier = "X_excluded_PabC_like"
            interpretation = "not_counted_as_BCAT_aroma_marker"

        else:
            tier = "C_BCAT_annotation_only"
            interpretation = "branched_chain_transaminase_candidate"

    # --------------------------------------------------------
    # Aromatic amino acid aminotransferases
    # --------------------------------------------------------
    elif marker == "aromatic_amino_acid_aminotransferase":

        if pref.startswith("arat") or "k00841" in ko:
            tier = "A_AraT_supported"
            interpretation = "aromatic_amino_acid_transamination"

        elif pref.startswith("tyrb") or "k00832" in ko:
            tier = "B_TyrB_broad_aromatic_aminotransferase"
            interpretation = "broad_aromatic_amino_acid_transamination"

        else:
            tier = "C_aromatic_aminotransferase_candidate"
            interpretation = "aromatic_transamination_candidate"

    # --------------------------------------------------------
    # Metionina gamma-liasa
    # --------------------------------------------------------
    elif marker == "methionine_gamma_lyase":

        tier = "A_methionine_gamma_lyase_supported"
        interpretation = "direct_sulfur_volatile_precursor_candidate"

    # --------------------------------------------------------
    # Otras liasas de aminoácidos azufrados
    # --------------------------------------------------------
    elif marker == "sulfur_amino_acid_lyase_candidate":

        if "k01761" in ko or "4.4.1.11" in ec:

            tier = "A_K01761_methionine_gamma_lyase_like"
            interpretation = "direct_sulfur_volatile_precursor_candidate"

            if (
                "methionine gamma-lyase" not in desc
                and not pref.startswith("mgl")
            ):
                discordance = "1"

        elif (
            pref.startswith("mccb")
            or "k17217" in ko
            or "volatile sulfur compounds" in desc
        ):
            tier = "B_cheese_relevant_CBL_like"
            interpretation = "sulfur_aroma_candidate"

        elif pref.startswith("metc") or "k01760" in ko:
            tier = "C_MetC_cystathionine_beta_lyase"
            interpretation = "sulfur_amino_acid_metabolism_possible_aroma_role"

        elif pref.startswith("metb") or "k01739" in ko:
            tier = "D_MetB_sulfur_amino_acid_metabolism"
            interpretation = "sulfur_amino_acid_metabolism_not_direct_aroma_proof"

        else:
            tier = "E_other_sulfur_lyase_candidate"
            interpretation = "uncertain_sulfur_aroma_relevance"

    else:
        tier = "E_other_aroma_marker"
        interpretation = "uncurated_marker"

    nr = dict(r)
    nr["curated_aroma_tier"] = tier
    nr["curated_interpretation"] = interpretation
    nr["annotation_discordance"] = discordance

    aroma_rows.append(nr)


aroma_by_mag = defaultdict(list)

for r in aroma_rows:
    aroma_by_mag[r["MAG"]].append(r)


# ============================================================
# 3. CURACIÓN DE EPS / CÁPSULA
# ============================================================

curated_eps_rows = []
eps_by_mag = defaultdict(list)

for r in eps_loci:

    classes = set(
        x.strip()
        for x in r.get("marker_classes", "").split(";")
        if x.strip()
    )

    original_tier = r.get("curated_tier", "")

    has_eps = bool(
        {
            "explicit_exopolysaccharide",
            "eps_named"
        } & classes
    )

    has_capsule = bool(
        {
            "capsular_polysaccharide",
            "cps_named",
            "capsule_export_regulation"
        } & classes
    )

    has_export = bool(
        {
            "flippase_or_polymerase",
            "polysaccharide_export_polymerization"
        } & classes
    )

    n_gt = to_int(r.get("n_glycosyltransferase_genes"))

    if has_eps and has_capsule:
        interpretation = "mixed_EPS_capsule_like_locus"

    elif has_eps:
        interpretation = "EPS_biosynthesis_like_locus"

    elif has_capsule:
        interpretation = "capsule_biosynthesis_like_locus"

    elif has_export and n_gt >= 1:
        interpretation = "generic_surface_polysaccharide_locus"

    elif n_gt >= 2:
        interpretation = "glycosyltransferase_rich_context"

    else:
        interpretation = "single_or_weak_polysaccharide_context"

    nr = dict(r)

    nr["has_explicit_EPS_marker"] = int(has_eps)
    nr["has_capsule_marker"] = int(has_capsule)
    nr["has_export_polymerization_marker"] = int(has_export)
    nr["curated_interpretation_79b"] = interpretation

    curated_eps_rows.append(nr)
    eps_by_mag[r["MAG"]].append(nr)


# ============================================================
# 4. RECONSTRUCCIÓN DE SISTEMAS OSMÓTICOS
# ============================================================

osm_rows_by_mag = defaultdict(list)

for r in hits:
    if r.get("block") == "osmotic_salt_stress":
        osm_rows_by_mag[r["MAG"]].append(r)


def find_marker_cluster(rows, required_markers, max_span=15):
    """
    Busca todos los marcadores requeridos en el mismo contig y dentro
    de una ventana pequeña de ORFs.
    """

    by_contig = defaultdict(lambda: defaultdict(list))

    for r in rows:
        marker = r.get("marker", "")
        if marker in required_markers:
            by_contig[r.get("contig", "")][marker].append(r)

    best = None

    for contig, d in by_contig.items():

        if not all(x in d for x in required_markers):
            continue

        lists = [d[x] for x in required_markers]

        for combo in product(*lists):

            idx = [gene_index(x) for x in combo]

            span = max(idx) - min(idx)

            if span <= max_span:

                if best is None or span < best[0]:
                    best = (span, contig, combo)

    return best


def find_pref_cluster(rows, prefixes, max_span=15):
    """
    Reconstrucción por Preferred_name, útil para OpuC/OpuA.
    """

    by_contig = defaultdict(lambda: defaultdict(list))

    for r in rows:

        pref = lower(r.get("Preferred_name"))

        for label, prefix in prefixes.items():

            if pref.startswith(prefix):
                by_contig[r.get("contig", "")][label].append(r)

    best = None

    for contig, d in by_contig.items():

        if not all(x in d for x in prefixes):
            continue

        lists = [d[x] for x in prefixes]

        for combo in product(*lists):

            idx = [gene_index(x) for x in combo]

            span = max(idx) - min(idx)

            if span <= max_span:

                if best is None or span < best[0]:
                    best = (span, contig, combo)

    return best


osm_system_rows = []
osm_summary = {}


for mag in mag_order:

    rows = osm_rows_by_mag.get(mag, [])

    detected_complete = []

    # --------------------------------------------------------
    # BetAB
    # --------------------------------------------------------
    systems = [
        (
            "glycine_betaine_BetAB",
            ["BetA", "BetB"],
            10
        ),
        (
            "trehalose_OtsAB",
            ["OtsA", "OtsB"],
            10
        ),
        (
            "ectoine_EctABC",
            ["Ect_A", "Ect_B", "Ect_C"],
            10
        ),
        (
            "ProU_VWX",
            ["ProU_V", "ProU_W", "ProU_X"],
            15
        ),
    ]

    for system_name, req, span_limit in systems:

        hit = find_marker_cluster(
            rows,
            req,
            max_span=span_limit
        )

        if hit is not None:

            span, contig, combo = hit

            detected_complete.append(system_name)

            osm_system_rows.append({
                "MAG": mag,
                "system": system_name,
                "evidence_level": "complete_multigene_system",
                "contig": contig,
                "ORF_span": span,
                "genes": semicolon_join(
                    x.get("gene_id", "") for x in combo
                )
            })

    # --------------------------------------------------------
    # OpuC: CA / CB / CC / CD
    # --------------------------------------------------------
    opuc = find_pref_cluster(
        rows,
        {
            "A": "opuca",
            "B": "opucb",
            "C": "opucc",
            "D": "opucd"
        },
        max_span=10
    )

    if opuc is not None:

        span, contig, combo = opuc

        detected_complete.append(
            "OpuC_compatible_solute_transport"
        )

        osm_system_rows.append({
            "MAG": mag,
            "system": "OpuC_compatible_solute_transport",
            "evidence_level": "complete_multigene_system",
            "contig": contig,
            "ORF_span": span,
            "genes": semicolon_join(
                x.get("gene_id", "") for x in combo
            )
        })

    # --------------------------------------------------------
    # OpuA: AA / AB / AC
    # --------------------------------------------------------
    opua = find_pref_cluster(
        rows,
        {
            "A": "opuaa",
            "B": "opuab",
            "C": "opuac"
        },
        max_span=10
    )

    if opua is not None:

        span, contig, combo = opua

        detected_complete.append(
            "OpuA_compatible_solute_transport"
        )

        osm_system_rows.append({
            "MAG": mag,
            "system": "OpuA_compatible_solute_transport",
            "evidence_level": "complete_multigene_system",
            "contig": contig,
            "ORF_span": span,
            "genes": semicolon_join(
                x.get("gene_id", "") for x in combo
            )
        })

    # --------------------------------------------------------
    # Transportadores que pueden funcionar como unidades
    # individuales
    # --------------------------------------------------------
    standalone = []

    for r in rows:

        marker = r.get("marker", "")
        pref = lower(r.get("Preferred_name"))

        if marker == "ProP_osmolyte_transport":
            standalone.append("ProP")

        elif marker == "explicit_betaine_transport":
            standalone.append("explicit_betaine_transporter")

        elif pref.startswith("opud"):
            standalone.append("OpuD_BCCT")

    standalone = sorted(set(standalone))

    # BetT se conserva como transporte de precursor, no se usa solo
    # para declarar mecanismo completo de osmoprotección.
    bett_n = sum(
        1
        for r in rows
        if r.get("marker") == "BetT_choline_transport"
    )

    n_rows = len(rows)

    if detected_complete:

        evidence = "strong_multigene_osmoadaptation_system"

    elif standalone:

        evidence = "standalone_osmolyte_transporter_supported"

    elif n_rows > 0:

        evidence = "partial_or_single_component_only"

    else:

        evidence = "not_detected"

    osm_summary[mag] = {
        "n_complete_osmotic_systems": len(
            set(detected_complete)
        ),
        "complete_osmotic_systems": semicolon_join(
            detected_complete
        ),
        "standalone_osmolyte_transporters": semicolon_join(
            standalone
        ),
        "BetT_precursor_transport_n": bett_n,
        "osmotic_stress_evidence_curated": evidence
    }


# ============================================================
# 5. MATRIZ FINAL 79b
# ============================================================

final_rows = []

for mag in mag_order:

    m = master_by_mag[mag]

    # --------------------------------------------------------
    # Lipólisis
    # --------------------------------------------------------
    lc = lip_counts.get(mag, Counter())

    n_lip_A = lc.get(
        "A_specific_lipase_candidate",
        0
    )

    n_lip_B = lc.get(
        "B_esterase_lipase_candidate",
        0
    )

    n_lip_C = lc.get(
        "C_broad_GDSL_acylhydrolase",
        0
    )

    n_lip_D = lc.get(
        "D_lipase_accessory_not_enzyme",
        0
    )

    if n_lip_A > 0:
        lip_ev = "specific_lipase_candidate_recovered"

    elif n_lip_B > 0:
        lip_ev = "esterase_lipase_candidate_recovered"

    elif n_lip_C > 0:
        lip_ev = "broad_GDSL_hydrolase_only"

    elif n_lip_D > 0:
        lip_ev = "lipase_accessory_only"

    else:
        lip_ev = "no_curated_lipase_candidate"

    # --------------------------------------------------------
    # Aroma
    # --------------------------------------------------------
    ar = aroma_by_mag.get(mag, [])

    bcat = sum(
        r["curated_aroma_tier"] == "A_BCAT_IlvE_supported"
        for r in ar
    )

    arat = sum(
        r["curated_aroma_tier"] == "A_AraT_supported"
        for r in ar
    )

    tyrb = sum(
        r["curated_aroma_tier"]
        == "B_TyrB_broad_aromatic_aminotransferase"
        for r in ar
    )

    mgl = sum(
        r["curated_aroma_tier"] in {
            "A_methionine_gamma_lyase_supported",
            "A_K01761_methionine_gamma_lyase_like"
        }
        for r in ar
    )

    cheese_cbl = sum(
        r["curated_aroma_tier"]
        == "B_cheese_relevant_CBL_like"
        for r in ar
    )

    metc = sum(
        r["curated_aroma_tier"]
        == "C_MetC_cystathionine_beta_lyase"
        for r in ar
    )

    metb = sum(
        r["curated_aroma_tier"]
        == "D_MetB_sulfur_amino_acid_metabolism"
        for r in ar
    )

    sulfur_direct = mgl + cheese_cbl

    transamination = bcat + arat + tyrb

    if sulfur_direct > 0 and transamination > 0:
        aroma_ev = (
            "sulfur_aroma_candidate_plus_"
            "amino_acid_transamination"
        )

    elif sulfur_direct > 0:
        aroma_ev = "sulfur_aroma_candidate"

    elif arat > 0 or tyrb > 0:
        aroma_ev = "aromatic_amino_acid_transamination_potential"

    elif bcat > 0:
        aroma_ev = "branched_chain_transamination_only"

    elif metc > 0 or metb > 0:
        aroma_ev = "sulfur_amino_acid_metabolism_only"

    else:
        aroma_ev = "not_detected"

    ketoacid_n = to_int(
        m.get("ketoacid_decarboxylase_n", 0)
    )

    if transamination > 0 and ketoacid_n == 0:
        ehrlich = (
            "transamination_recovered_but_"
            "no_ketoacid_decarboxylase_marker"
        )
    elif transamination > 0 and ketoacid_n > 0:
        ehrlich = "multiple_Ehrlich_route_markers_recovered"
    else:
        ehrlich = "no_resolved_Ehrlich_route"

    # --------------------------------------------------------
    # EPS / cápsula
    # --------------------------------------------------------
    erows = eps_by_mag.get(mag, [])

    strong_probable = [
        r for r in erows
        if (
            r.get("curated_tier", "").startswith("A_")
            or r.get("curated_tier", "").startswith("B_")
        )
    ]

    eps_only = sum(
        r["curated_interpretation_79b"]
        == "EPS_biosynthesis_like_locus"
        for r in strong_probable
    )

    capsule_only = sum(
        r["curated_interpretation_79b"]
        == "capsule_biosynthesis_like_locus"
        for r in strong_probable
    )

    mixed = sum(
        r["curated_interpretation_79b"]
        == "mixed_EPS_capsule_like_locus"
        for r in strong_probable
    )

    generic = sum(
        r["curated_interpretation_79b"]
        in {
            "generic_surface_polysaccharide_locus",
            "glycosyltransferase_rich_context"
        }
        for r in strong_probable
    )

    if mixed > 0:
        eps_ev = "mixed_EPS_capsule_like_locus_recovered"

    elif eps_only > 0 and capsule_only > 0:
        eps_ev = "separate_EPS_like_and_capsule_like_loci"

    elif eps_only > 0:
        eps_ev = "EPS_biosynthesis_like_locus_recovered"

    elif capsule_only > 0:
        eps_ev = "capsule_biosynthesis_like_locus_recovered"

    elif generic > 0:
        eps_ev = "generic_surface_polysaccharide_context"

    elif to_int(m.get("single_marker_EPS_contexts", 0)) > 0:
        eps_ev = "single_marker_context_only"

    else:
        eps_ev = "no_resolved_EPS_capsule_locus"

    # --------------------------------------------------------
    # Osmótico
    # --------------------------------------------------------
    osm = osm_summary[mag]

    row = {
        "MAG": mag,
        "genus": m.get("genus", ""),
        "species": m.get("species", ""),
        "completeness": m.get("completeness", ""),

        "specific_lipase_candidate_n": n_lip_A,
        "esterase_lipase_candidate_n": n_lip_B,
        "broad_GDSL_acylhydrolase_n": n_lip_C,
        "lipase_accessory_n": n_lip_D,
        "lipolysis_evidence_curated": lip_ev,

        "BCAT_IlvE_n": bcat,
        "AraT_n": arat,
        "TyrB_n": tyrb,
        "MGL_K01761_like_n": mgl,
        "cheese_relevant_CBL_like_n": cheese_cbl,
        "MetC_n": metc,
        "MetB_n": metb,
        "amino_acid_aroma_evidence_curated": aroma_ev,
        "Ehrlich_route_evidence": ehrlich,

        "EPS_like_loci_n": eps_only,
        "capsule_like_loci_n": capsule_only,
        "mixed_EPS_capsule_like_loci_n": mixed,
        "generic_polysaccharide_loci_n": generic,
        "EPS_capsule_evidence_curated": eps_ev,

        "acid_stress_evidence": m.get(
            "acid_stress_evidence",
            ""
        ),
        "acid_stress_systems": m.get(
            "acid_stress_systems",
            ""
        ),

        "n_complete_osmotic_systems":
            osm["n_complete_osmotic_systems"],

        "complete_osmotic_systems":
            osm["complete_osmotic_systems"],

        "standalone_osmolyte_transporters":
            osm["standalone_osmolyte_transporters"],

        "BetT_precursor_transport_n":
            osm["BetT_precursor_transport_n"],

        "osmotic_stress_evidence_curated":
            osm["osmotic_stress_evidence_curated"],

        "oxidative_stress_evidence": m.get(
            "oxidative_stress_evidence",
            ""
        ),
        "oxidative_stress_mechanisms": m.get(
            "oxidative_stress_mechanisms",
            ""
        ),
    }

    final_rows.append(row)


# ============================================================
# 6. MAGs RELACIONADOS CON BACTERIOCINAS
# ============================================================

BACT_ROLES = {
    "L2__L2_maxbin2.004_sub":
        "source_MAG_ATTRLOC001",

    "L3__concoct_29":
        "source_MAG_ATTRLOC003_ATTRLOC004",

    "M2__M2_maxbin2.004":
        "source_MAG_ATTRLOC005_ATTRLOC006_ATTRLOC007",

    "L2__L2_maxbin2.011_sub":
        "competitive_shared_context_ATTRLOC002_not_independent",
}


focus_rows = []

for r in final_rows:

    if r["MAG"] not in BACT_ROLES:
        continue

    nr = dict(r)
    nr["bacteriocin_analysis_role"] = BACT_ROLES[r["MAG"]]

    focus_rows.append(nr)


# ============================================================
# 7. PREVALENCIA DE CATEGORÍAS CURADAS
# ============================================================

prevalence_rows = []

summary_fields = [
    ("lipolysis", "lipolysis_evidence_curated"),
    (
        "amino_acid_aroma",
        "amino_acid_aroma_evidence_curated"
    ),
    (
        "EPS_capsule",
        "EPS_capsule_evidence_curated"
    ),
    (
        "osmotic_salt_stress",
        "osmotic_stress_evidence_curated"
    ),
    (
        "acid_stress",
        "acid_stress_evidence"
    ),
    (
        "oxidative_stress",
        "oxidative_stress_evidence"
    )
]

for block, field in summary_fields:

    cats = defaultdict(list)

    for r in final_rows:
        cats[r.get(field, "")].append(r["MAG"])

    for category, mags in sorted(cats.items()):

        prevalence_rows.append({
            "block": block,
            "category": category,
            "n_MAGs": len(mags),
            "percent_MAGs": f"{100 * len(mags) / 18:.2f}",
            "MAGs": ";".join(mags)
        })


# ============================================================
# 8. ESCRITURA
# ============================================================

lip_fields = list(hits[0].keys()) + [
    "curated_lipolysis_tier",
    "curation_reason"
]

write_tsv(
    os.path.join(
        OUT,
        "curated_lipolysis_candidates.tsv"
    ),
    lipase_rows,
    lip_fields
)


aroma_fields = list(hits[0].keys()) + [
    "curated_aroma_tier",
    "curated_interpretation",
    "annotation_discordance"
]

write_tsv(
    os.path.join(
        OUT,
        "curated_amino_acid_aroma_markers.tsv"
    ),
    aroma_rows,
    aroma_fields
)


if curated_eps_rows:

    eps_fields = list(eps_loci[0].keys()) + [
        "has_explicit_EPS_marker",
        "has_capsule_marker",
        "has_export_polymerization_marker",
        "curated_interpretation_79b"
    ]

    write_tsv(
        os.path.join(
            OUT,
            "curated_EPS_capsule_loci.tsv"
        ),
        curated_eps_rows,
        eps_fields
    )


write_tsv(
    os.path.join(
        OUT,
        "curated_osmotic_systems.tsv"
    ),
    osm_system_rows,
    [
        "MAG",
        "system",
        "evidence_level",
        "contig",
        "ORF_span",
        "genes"
    ]
)


final_fields = list(final_rows[0].keys())

write_tsv(
    os.path.join(
        OUT,
        "remaining_cheese_functions_curated_matrix_18MAGs.tsv"
    ),
    final_rows,
    final_fields
)


focus_fields = [
    "MAG",
    "genus",
    "species",
    "completeness",
    "bacteriocin_analysis_role"
] + [
    x for x in final_fields
    if x not in {
        "MAG",
        "genus",
        "species",
        "completeness"
    }
]

write_tsv(
    os.path.join(
        OUT,
        "remaining_functions_bacteriocin_focus_curated.tsv"
    ),
    focus_rows,
    focus_fields
)


write_tsv(
    os.path.join(
        OUT,
        "curated_remaining_function_prevalence.tsv"
    ),
    prevalence_rows,
    [
        "block",
        "category",
        "n_MAGs",
        "percent_MAGs",
        "MAGs"
    ]
)


# ============================================================
# 9. README
# ============================================================

readme = f"""PASO 79b - CURACION DE FUNCIONES TECNOLOGICAS RESTANTES

MAGs analizados: 18

PRINCIPIOS DE CURACION

1. LIPOLISIS
A_specific_lipase_candidate:
    Lipasa clase 3, K01046, EC 3.1.1.3 o anotacion
    explicitamente secretory lipase.

B_esterase_lipase_candidate:
    Esterasas/lipasas o arquitecturas lipoliticas mas especificas,
    pero sin evidencia suficiente para inferir hidrolisis de grasa
    lactea.

C_broad_GDSL_acylhydrolase:
    Familias GDSL amplias. No deben interpretarse como evidencia
    directa de lipolisis de trigliceridos de leche.

D_lipase_accessory_not_enzyme:
    Proteinas accesorias de plegamiento/secrecion de lipasas.

La presencia genomica representa potencial, no actividad enzimatica.

2. AROMA DE AMINOACIDOS
IlvE/K00826 se conserva como evidencia de transaminacion de BCAA.
PabC/K02619 no se cuenta como BCAT tecnologicamente interpretable.
AraT y TyrB se separan.
K01761/EC 4.4.1.11 se trata como candidato tipo
methionine-gamma-lyase.
MccB/Cystathionine beta-lyase relacionado con compuestos sulfurados
se conserva como candidato de aroma.
MetB/MetC por si solos no demuestran produccion de aroma.

La ausencia de ketoacid decarboxylase en 79a impide declarar
reconstruida una ruta completa tipo Ehrlich.

3. EPS / CAPSULA
EPS-like y capsule-like se mantienen separados.
Un locus de polisacarido superficial no demuestra fenotipo de
produccion de EPS ni efecto de textura en queso.

4. OSMOADAPTACION
Solo se consideran sistemas multigenicos fuertes cuando los
componentes requeridos aparecen colocalizados:
BetAB, OtsAB, EctABC, ProU-VWX, OpuC u OpuA.

ProP/OpuD pueden conservarse como transportadores individuales
funcionales.

BetT aislado se registra como transporte de precursor y no eleva
por si solo al MAG a sistema osmoadaptativo completo.

5. INTERPRETACION
Todos los resultados son potencial genomico recuperado en MAGs.
No equivalen a expresion, actividad enzimatica ni fenotipo.

Salida:
{OUT}
"""

with open(
    os.path.join(
        OUT,
        "README_step79b.txt"
    ),
    "w",
    encoding="utf-8"
) as fh:
    fh.write(readme)


# ============================================================
# 10. RESUMEN CONSOLA
# ============================================================

print("=" * 60)
print("PASO 79b COMPLETADO")
print("=" * 60)
print(f"MAGs analizados: {len(final_rows)}")
print()

print("LIPOLISIS CURADA")
for category in [
    "A_specific_lipase_candidate",
    "B_esterase_lipase_candidate",
    "C_broad_GDSL_acylhydrolase",
    "D_lipase_accessory_not_enzyme"
]:
    genes = sum(
        1
        for r in lipase_rows
        if r["curated_lipolysis_tier"] == category
    )
    mags = len({
        r["MAG"]
        for r in lipase_rows
        if r["curated_lipolysis_tier"] == category
    })

    print(
        f"{category:<40} "
        f"{mags:>2}/18 MAGs  {genes:>3} genes"
    )

print()
print("SISTEMAS OSMOTICOS COMPLETOS")

system_counts = Counter(
    r["system"]
    for r in osm_system_rows
)

for system, n in sorted(system_counts.items()):
    print(f"{system:<40} {n:>2} MAGs")

print()
print("EPS/CAPSULA CURADO")

eps_counts = Counter(
    r["EPS_capsule_evidence_curated"]
    for r in final_rows
)

for k, v in sorted(eps_counts.items()):
    print(f"{k:<45} {v:>2}/18")

print()
print("AROMA AA CURADO")

aroma_counts = Counter(
    r["amino_acid_aroma_evidence_curated"]
    for r in final_rows
)

for k, v in sorted(aroma_counts.items()):
    print(f"{k:<55} {v:>2}/18")

print()
print(f"Salida: {OUT}")
print(
    "PASO 79b FINALIZÓ CORRECTAMENTE; "
    "ES SEGURO SALIR."
)
