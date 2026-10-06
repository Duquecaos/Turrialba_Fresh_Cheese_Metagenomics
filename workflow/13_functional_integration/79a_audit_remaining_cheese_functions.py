#!/usr/bin/env python3

import csv
import os
import re
from collections import defaultdict, Counter

USER = os.environ["USER"]
ROOT = f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"

MASTER = os.path.join(
    ROOT,
    "22_final_representative_mags",
    "representative_mags_master.tsv"
)

GENES = os.path.join(
    ROOT,
    "53_metabolic_annotation_consolidated",
    "gene_function_inventory_all_predicted.tsv"
)

OUTDIR = os.path.join(
    ROOT,
    "59_remaining_cheese_function_audit"
)

os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# UTILIDADES
# ============================================================

def clean(x):
    if x is None:
        return ""
    x = str(x).strip()
    if x in {"-", "NA", "N/A", "None", "nan"}:
        return ""
    return x


def lower(x):
    return clean(x).lower()


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fields):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore"
        )
        w.writeheader()
        for r in rows:
            w.writerow(r)


def parse_gene_id(gene_id):
    """
    Ejemplo:
    MAG___contig000024_400
    """
    gene_id = clean(gene_id)

    if "___" not in gene_id:
        return "", None

    suffix = gene_id.split("___", 1)[1]

    m = re.match(r"(.+)_([0-9]+)$", suffix)

    if not m:
        return suffix, None

    return m.group(1), int(m.group(2))


def same_contig_near(list_a, list_b, max_gap=10):
    for a in list_a:
        for b in list_b:
            if (
                a["contig"] and
                a["contig"] == b["contig"] and
                a["gene_index"] is not None and
                b["gene_index"] is not None and
                abs(a["gene_index"] - b["gene_index"]) <= max_gap
            ):
                return True
    return False


def component_cluster(marker_lists, components, min_components, max_span=15):
    """
    Detecta si >= min_components componentes diferentes
    están próximos en el mismo contig.
    """
    by_contig = defaultdict(list)

    for comp in components:
        for g in marker_lists.get(comp, []):
            if g["gene_index"] is not None:
                by_contig[g["contig"]].append(
                    (g["gene_index"], comp, g)
                )

    for contig, vals in by_contig.items():
        vals.sort()

        for i in range(len(vals)):
            start = vals[i][0]
            seen = set()

            for j in range(i, len(vals)):
                if vals[j][0] - start > max_span:
                    break

                seen.add(vals[j][1])

                if len(seen) >= min_components:
                    return True

    return False


# ============================================================
# MASTER DE MAGs
# ============================================================

master_rows = read_tsv(MASTER)

meta = {}
mag_order = []

for r in master_rows:
    mag = clean(r.get("representative_MAG"))

    if not mag:
        continue

    mag_order.append(mag)

    meta[mag] = {
        "genus": clean(r.get("genus")),
        "species": clean(r.get("species")),
        "completeness": clean(r.get("completeness")),
        "contamination": clean(r.get("contamination")),
    }

if len(mag_order) != 18:
    raise RuntimeError(
        f"Se esperaban 18 MAGs finales y se encontraron {len(mag_order)}"
    )


# ============================================================
# LECTURA DE INVENTARIO DE GENES
# ============================================================

raw_genes = read_tsv(GENES)

genes = []

for r in raw_genes:

    mag = clean(r.get("MAG"))

    gene_id = (
        clean(r.get("gene_id"))
        or clean(r.get("#query"))
        or clean(r.get("query"))
    )

    if not mag or mag not in meta or not gene_id:
        continue

    contig, gene_index = parse_gene_id(gene_id)

    g = {
        "MAG": mag,
        "gene_id": gene_id,
        "contig": contig,
        "gene_index": gene_index,

        "Preferred_name": clean(r.get("Preferred_name")),
        "Description": clean(r.get("Description")),
        "KEGG_ko": clean(r.get("KEGG_ko")),
        "EC": clean(r.get("EC")),
        "PFAMs": clean(r.get("PFAMs")),

        "_name": lower(r.get("Preferred_name")),
        "_desc": lower(r.get("Description")),
        "_ko": lower(r.get("KEGG_ko")),
        "_pfam": lower(r.get("PFAMs")),
    }

    genes.append(g)


genes_by_mag = defaultdict(list)
genes_by_mag_contig = defaultdict(list)

for g in genes:
    genes_by_mag[g["MAG"]].append(g)
    genes_by_mag_contig[(g["MAG"], g["contig"])].append(g)


# ============================================================
# DEFINICIÓN DE MARCADORES
# ============================================================

marker_hits = []
marker_seen = set()

stress_hits = defaultdict(lambda: defaultdict(list))


def add_hit(g, block, marker):
    key = (g["gene_id"], block, marker)

    if key in marker_seen:
        return

    marker_seen.add(key)

    marker_hits.append({
        "block": block,
        "marker": marker,
        "MAG": g["MAG"],
        "gene_id": g["gene_id"],
        "contig": g["contig"],
        "gene_index": (
            "" if g["gene_index"] is None else g["gene_index"]
        ),
        "Preferred_name": g["Preferred_name"],
        "Description": g["Description"],
        "KEGG_ko": g["KEGG_ko"],
        "EC": g["EC"],
        "PFAMs": g["PFAMs"],
    })


# ============================================================
# 1. LIPÓLISIS / ESTERASAS
# ============================================================

for g in genes:

    n = g["_name"]
    d = g["_desc"]
    p = g["_pfam"]

    lipoyl_related = any(
        x in d
        for x in (
            "lipoyl",
            "lipoate",
            "lipoic acid",
            "lipoylation",
        )
    )

    phospholipase = (
        "phospholipase" in d
        or re.match(r"^pl[acd][a-z0-9_-]*$", n or "") is not None
    )

    explicit_lipase = (
        (
            re.search(r"\blipase\b", d) is not None
            or "triacylglycerol lipase" in d
            or "lipase family" in d
            or re.search(r"\blipase\b", p) is not None
        )
        and not lipoyl_related
        and not phospholipase
    )

    explicit_esterase = (
        re.search(r"\besterase\b", d) is not None
        or "carboxylesterase" in d
    )

    if explicit_lipase:
        add_hit(g, "lipolysis", "explicit_lipase_candidate")

    if explicit_esterase:
        add_hit(g, "lipolysis", "explicit_esterase_candidate")

    if phospholipase:
        add_hit(g, "lipolysis", "phospholipase")


# ============================================================
# 2. AMINOÁCIDOS RELACIONADOS CON AROMA
# ============================================================

for g in genes:

    n = g["_name"]
    d = g["_desc"]

    # Branched-chain amino-acid aminotransferase
    if (
        n in {"ilve", "bcat", "bca_t", "bcaat"}
        or "branched-chain amino acid aminotransferase" in d
        or "branched chain amino acid aminotransferase" in d
    ):
        add_hit(
            g,
            "amino_acid_aroma",
            "branched_chain_aminotransferase"
        )

    # Aromatic amino-acid aminotransferase
    if (
        n in {"arat", "tyrb"}
        or "aromatic amino acid aminotransferase" in d
        or "aromatic-amino-acid aminotransferase" in d
        or "aromatic-amino-acid transaminase" in d
    ):
        add_hit(
            g,
            "amino_acid_aroma",
            "aromatic_amino_acid_aminotransferase"
        )

    # Methionine gamma-lyase
    if (
        n in {"mdea", "mgl"}
        or "methionine gamma-lyase" in d
        or "methionine gamma lyase" in d
    ):
        add_hit(
            g,
            "amino_acid_aroma",
            "methionine_gamma_lyase"
        )

    # Sulfur amino-acid lyases related but not equivalent to MGL
    if (
        "cystathionine beta-lyase" in d
        or "cystathionine gamma-lyase" in d
        or "cystathionine beta lyase" in d
        or "cystathionine gamma lyase" in d
    ):
        add_hit(
            g,
            "amino_acid_aroma",
            "sulfur_amino_acid_lyase_candidate"
        )

    # Keto-acid decarboxylation
    if (
        n == "kivd"
        or "alpha-keto acid decarboxylase" in d
        or "alpha keto acid decarboxylase" in d
        or "2-ketoacid decarboxylase" in d
    ):
        add_hit(
            g,
            "amino_acid_aroma",
            "ketoacid_decarboxylase"
        )


# ============================================================
# 3. EPS / POLISACÁRIDO CAPSULAR
# ============================================================

def eps_anchor_classes(g):

    n = g["_name"]
    d = g["_desc"]

    classes = set()

    if re.match(r"^eps[a-z0-9_-]+$", n or ""):
        classes.add("eps_named")

    if re.match(r"^cps[a-z0-9_-]+$", n or ""):
        classes.add("cps_named")

    if n in {"wzx", "wzy", "wzz"}:
        classes.add("flippase_or_polymerase")

    if n in {"wza", "wzb", "wzc"}:
        classes.add("capsule_export_regulation")

    if "exopolysaccharide" in d:
        classes.add("explicit_exopolysaccharide")

    if "capsular polysaccharide" in d:
        classes.add("capsular_polysaccharide")

    if (
        "polysaccharide export" in d
        or "polysaccharide polymerase" in d
        or "polysaccharide flippase" in d
    ):
        classes.add("polysaccharide_export_polymerization")

    return classes


def is_glycosyltransferase(g):

    d = g["_desc"]
    p = g["_pfam"]

    return (
        "glycosyltransferase" in d
        or "glycosyl transferase" in d
        or "glycosyltransf" in p
    )


eps_loci = []
eps_locus_counter = 0

eps_status_by_mag = defaultdict(
    lambda: {
        "A_strong": 0,
        "B_probable": 0,
        "C_single": 0
    }
)


for (mag, contig), contig_genes in genes_by_mag_contig.items():

    indexed = [
        g for g in contig_genes
        if g["gene_index"] is not None
    ]

    if not indexed:
        continue

    indexed.sort(key=lambda x: x["gene_index"])

    anchors = []

    for g in indexed:
        classes = eps_anchor_classes(g)

        if classes:
            anchors.append(
                (g["gene_index"], g, classes)
            )

            for c in classes:
                add_hit(
                    g,
                    "EPS_capsule",
                    c
                )

    if not anchors:
        continue

    # Agrupa anchors separados por <=20 ORFs
    anchor_groups = []
    current = [anchors[0]]

    for a in anchors[1:]:

        if a[0] - current[-1][0] <= 20:
            current.append(a)
        else:
            anchor_groups.append(current)
            current = [a]

    anchor_groups.append(current)

    for group in anchor_groups:

        min_anchor = min(x[0] for x in group)
        max_anchor = max(x[0] for x in group)

        start = max(1, min_anchor - 10)
        end = max_anchor + 10

        window = [
            g for g in indexed
            if start <= g["gene_index"] <= end
        ]

        classes = set()
        marker_gene_ids = []
        gt_gene_ids = []

        for g in window:

            cs = eps_anchor_classes(g)

            if cs:
                classes.update(cs)
                marker_gene_ids.append(g["gene_id"])

            if is_glycosyltransferase(g):
                classes.add("glycosyltransferase")
                gt_gene_ids.append(g["gene_id"])

        has_named = bool(
            classes.intersection({
                "eps_named",
                "cps_named",
                "explicit_exopolysaccharide",
                "capsular_polysaccharide",
            })
        )

        has_polymerization_export = bool(
            classes.intersection({
                "flippase_or_polymerase",
                "capsule_export_regulation",
                "polysaccharide_export_polymerization",
            })
        )

        has_gt = "glycosyltransferase" in classes

        n_anchor_genes = len(set(marker_gene_ids))

        if (
            has_named
            and has_polymerization_export
            and has_gt
        ):
            tier = "A_strong_EPS_capsule_like_locus"
            eps_status_by_mag[mag]["A_strong"] += 1

        elif (
            (has_named and has_gt)
            or
            (has_polymerization_export and has_gt)
            or
            n_anchor_genes >= 2
        ):
            tier = "B_probable_EPS_capsule_like_locus"
            eps_status_by_mag[mag]["B_probable"] += 1

        else:
            tier = "C_single_marker_context"
            eps_status_by_mag[mag]["C_single"] += 1

        eps_locus_counter += 1

        eps_loci.append({
            "locus_id": f"EPSLOC{eps_locus_counter:03d}",
            "MAG": mag,
            "genus": meta[mag]["genus"],
            "species": meta[mag]["species"],
            "contig": contig,
            "window_start_ORF": start,
            "window_end_ORF": end,
            "anchor_min_ORF": min_anchor,
            "anchor_max_ORF": max_anchor,
            "n_anchor_genes": n_anchor_genes,
            "n_glycosyltransferase_genes": len(
                set(gt_gene_ids)
            ),
            "marker_classes": ";".join(sorted(classes)),
            "marker_gene_ids": ";".join(
                sorted(set(marker_gene_ids))
            ),
            "glycosyltransferase_gene_ids": ";".join(
                sorted(set(gt_gene_ids))
            ),
            "curated_tier": tier,
        })


# ============================================================
# 4. ESTRÉS ÁCIDO
# ============================================================

for g in genes:

    n = g["_name"]
    d = g["_desc"]

    marker = None

    if n in {"gada", "gadb"} or "glutamate decarboxylase" in d:
        marker = "glutamate_decarboxylase"

    elif n == "gadc":
        marker = "glutamate_GABA_antiporter"

    elif n == "adia":
        marker = "arginine_decarboxylase_AR"

    elif n == "adic":
        marker = "arginine_agmatine_antiporter"

    elif (
        n == "arca"
        and "arginine deiminase" in d
    ):
        marker = "ADI_arcA"

    elif (
        n == "arcb"
        and (
            "ornithine carbamoyltransferase" in d
            or "ornithine transcarbamylase" in d
        )
    ):
        marker = "ADI_arcB"

    elif (
        n == "arcc"
        and "carbamate kinase" in d
    ):
        marker = "ADI_arcC"

    elif (
        n == "arcd"
        and (
            "arginine" in d
            and "antiporter" in d
        )
    ):
        marker = "ADI_arcD"

    elif n == "urea" and "urease" in d:
        marker = "urease_A"

    elif n == "ureb" and "urease" in d:
        marker = "urease_B"

    elif n == "urec" and "urease" in d:
        marker = "urease_C"

    if marker:
        add_hit(g, "acid_stress", marker)
        stress_hits[g["MAG"]][marker].append(g)


# ============================================================
# 5. ESTRÉS OSMÓTICO / SAL
# ============================================================

for g in genes:

    n = g["_name"]
    d = g["_desc"]

    markers = []

    if re.match(r"^opu[a-z0-9_-]+$", n or ""):
        markers.append("Opu_compatible_solute_transport")

    if n == "bett":
        markers.append("BetT_choline_transport")

    if n == "prop":
        markers.append("ProP_osmolyte_transport")

    if n in {"prov", "prow", "prox"}:
        markers.append(f"ProU_{n[-1].upper()}")

    if n == "beta":
        markers.append("BetA")

    if n == "betb":
        markers.append("BetB")

    if n in {"ecta", "ectb", "ectc", "ectd"}:
        markers.append(f"Ect_{n[-1].upper()}")

    if n == "otsa":
        markers.append("OtsA")

    if n == "otsb":
        markers.append("OtsB")

    # Descripciones explícitas como respaldo
    if (
        "glycine betaine transporter" in d
        or "betaine transporter" in d
    ):
        markers.append("explicit_betaine_transport")

    if "ectoine synthase" in d:
        markers.append("Ect_C")

    for marker in sorted(set(markers)):
        add_hit(g, "osmotic_salt_stress", marker)
        stress_hits[g["MAG"]][marker].append(g)


# ============================================================
# 6. ESTRÉS OXIDATIVO
# ============================================================

for g in genes:

    n = g["_name"]
    d = g["_desc"]

    markers = []

    if (
        n in {"kata", "kate", "katg"}
        or re.search(r"\bcatalase\b", d)
    ):
        markers.append("catalase")

    if (
        n in {"soda", "sodb", "sodc"}
        or "superoxide dismutase" in d
    ):
        markers.append("superoxide_dismutase")

    if (
        n in {"ahpc", "ahpf"}
        or "alkyl hydroperoxide reductase" in d
    ):
        markers.append("alkyl_hydroperoxide_reductase")

    if (
        n == "tpx"
        or "thiol peroxidase" in d
    ):
        markers.append("thiol_peroxidase")

    if (
        n == "bcp"
        or "bacterioferritin comigratory protein" in d
    ):
        markers.append("peroxiredoxin_Bcp")

    if n in {"msra", "msrb"}:
        markers.append("methionine_sulfoxide_repair")

    for marker in sorted(set(markers)):
        add_hit(g, "oxidative_stress", marker)
        stress_hits[g["MAG"]][marker].append(g)


# ============================================================
# EVALUACIÓN DE SISTEMAS DE ESTRÉS
# ============================================================

stress_summary = {}

for mag in mag_order:

    m = stress_hits[mag]

    # ----------------------------
    # ÁCIDO
    # ----------------------------

    acid_systems = []

    if same_contig_near(
        m.get("glutamate_decarboxylase", []),
        m.get("glutamate_GABA_antiporter", []),
        10
    ):
        acid_systems.append("glutamate_GAD_system")

    if same_contig_near(
        m.get("arginine_decarboxylase_AR", []),
        m.get("arginine_agmatine_antiporter", []),
        10
    ):
        acid_systems.append("arginine_decarboxylase_AR_system")

    if component_cluster(
        m,
        ["ADI_arcA", "ADI_arcB", "ADI_arcC", "ADI_arcD"],
        min_components=3,
        max_span=15
    ):
        acid_systems.append("arginine_deiminase_ADI_system")

    if component_cluster(
        m,
        ["urease_A", "urease_B", "urease_C"],
        min_components=3,
        max_span=10
    ):
        acid_systems.append("urease_core_system")

    acid_marker_n = sum(
        len(v)
        for k, v in m.items()
        if (
            k.startswith("ADI_")
            or k.startswith("urease_")
            or k in {
                "glutamate_decarboxylase",
                "glutamate_GABA_antiporter",
                "arginine_decarboxylase_AR",
                "arginine_agmatine_antiporter",
            }
        )
    )

    if acid_systems:
        acid_status = "system_supported"
    elif acid_marker_n > 0:
        acid_status = "marker_only"
    else:
        acid_status = "not_detected"

    # ----------------------------
    # OSMÓTICO
    # ----------------------------

    osmotic_mechanisms = []

    if m.get("BetT_choline_transport"):
        osmotic_mechanisms.append("BetT_transport")

    if m.get("ProP_osmolyte_transport"):
        osmotic_mechanisms.append("ProP_transport")

    if m.get("explicit_betaine_transport"):
        osmotic_mechanisms.append("explicit_betaine_transport")

    if component_cluster(
        m,
        ["ProU_V", "ProU_W", "ProU_X"],
        min_components=2,
        max_span=12
    ):
        osmotic_mechanisms.append("ProU_transport_system")

    opu_genes = m.get(
        "Opu_compatible_solute_transport",
        []
    )

    opu_supported = False

    by_contig_opu = defaultdict(list)

    for g in opu_genes:
        if g["gene_index"] is not None:
            by_contig_opu[g["contig"]].append(g["gene_index"])

    for vals in by_contig_opu.values():
        vals = sorted(vals)
        if len(vals) >= 2 and vals[-1] - vals[0] <= 20:
            opu_supported = True

    if opu_supported:
        osmotic_mechanisms.append("Opu_transport_cluster")
    elif opu_genes:
        osmotic_mechanisms.append("Opu_marker")

    if same_contig_near(
        m.get("BetA", []),
        m.get("BetB", []),
        10
    ):
        osmotic_mechanisms.append("glycine_betaine_BetAB")

    if component_cluster(
        m,
        ["Ect_A", "Ect_B", "Ect_C"],
        min_components=3,
        max_span=12
    ):
        osmotic_mechanisms.append("ectoine_EctABC")

    if same_contig_near(
        m.get("OtsA", []),
        m.get("OtsB", []),
        10
    ):
        osmotic_mechanisms.append("trehalose_OtsAB")

    osmotic_marker_n = sum(
        len(v)
        for k, v in m.items()
        if (
            k.startswith("Opu_")
            or k.startswith("ProU_")
            or k.startswith("Ect_")
            or k in {
                "BetT_choline_transport",
                "ProP_osmolyte_transport",
                "explicit_betaine_transport",
                "BetA",
                "BetB",
                "OtsA",
                "OtsB",
            }
        )
    )

    if osmotic_mechanisms:
        osmotic_status = "supported_mechanism"
    elif osmotic_marker_n > 0:
        osmotic_status = "marker_only"
    else:
        osmotic_status = "not_detected"

    # ----------------------------
    # OXIDATIVO
    # ----------------------------

    oxidative_markers = [
        "catalase",
        "superoxide_dismutase",
        "alkyl_hydroperoxide_reductase",
        "thiol_peroxidase",
        "peroxiredoxin_Bcp",
        "methionine_sulfoxide_repair",
    ]

    oxidative_present = [
        x for x in oxidative_markers
        if m.get(x)
    ]

    if len(oxidative_present) >= 3:
        oxidative_status = "broad_antioxidant_repertoire"

    elif len(oxidative_present) >= 1:
        oxidative_status = "supported"

    else:
        oxidative_status = "not_detected"

    stress_summary[mag] = {
        "acid_stress_evidence": acid_status,
        "acid_stress_systems": ";".join(
            sorted(set(acid_systems))
        ),
        "osmotic_salt_stress_evidence": osmotic_status,
        "osmotic_salt_mechanisms": ";".join(
            sorted(set(osmotic_mechanisms))
        ),
        "oxidative_stress_evidence": oxidative_status,
        "oxidative_stress_mechanisms": ";".join(
            sorted(set(oxidative_present))
        ),
    }


# ============================================================
# CONTEOS DE MARCADORES
# ============================================================

marker_count = defaultdict(Counter)

for h in marker_hits:
    marker_count[h["MAG"]][h["marker"]] += 1


# ============================================================
# MATRIZ REFINADA 18 MAGs
# ============================================================

matrix = []

for mag in mag_order:

    counts = marker_count[mag]

    lipase_n = counts["explicit_lipase_candidate"]
    esterase_n = counts["explicit_esterase_candidate"]
    phospholipase_n = counts["phospholipase"]

    if lipase_n > 0:
        lipolysis_status = "explicit_lipase_candidate_recovered"
    elif esterase_n > 0:
        lipolysis_status = "esterase_candidates_only"
    elif phospholipase_n > 0:
        lipolysis_status = "phospholipase_only"
    else:
        lipolysis_status = "no_specific_lipolytic_marker_recovered"

    aroma_markers = []

    aroma_marker_names = [
        "branched_chain_aminotransferase",
        "aromatic_amino_acid_aminotransferase",
        "methionine_gamma_lyase",
        "sulfur_amino_acid_lyase_candidate",
        "ketoacid_decarboxylase",
    ]

    for x in aroma_marker_names:
        if counts[x] > 0:
            aroma_markers.append(x)

    eps = eps_status_by_mag[mag]

    if eps["A_strong"] > 0:
        eps_status = "strong_EPS_capsule_like_locus"
    elif eps["B_probable"] > 0:
        eps_status = "probable_EPS_capsule_like_locus"
    elif eps["C_single"] > 0:
        eps_status = "single_marker_context_only"
    else:
        eps_status = "no_EPS_capsule_like_locus_recovered"

    row = {
        "MAG": mag,
        "genus": meta[mag]["genus"],
        "species": meta[mag]["species"],
        "completeness": meta[mag]["completeness"],

        "explicit_lipase_candidate_n": lipase_n,
        "explicit_esterase_candidate_n": esterase_n,
        "phospholipase_n": phospholipase_n,
        "lipolysis_evidence": lipolysis_status,

        "branched_chain_aminotransferase_n":
            counts["branched_chain_aminotransferase"],

        "aromatic_amino_acid_aminotransferase_n":
            counts["aromatic_amino_acid_aminotransferase"],

        "methionine_gamma_lyase_n":
            counts["methionine_gamma_lyase"],

        "sulfur_amino_acid_lyase_candidate_n":
            counts["sulfur_amino_acid_lyase_candidate"],

        "ketoacid_decarboxylase_n":
            counts["ketoacid_decarboxylase"],

        "n_aroma_marker_types": len(aroma_markers),
        "aroma_marker_types": ";".join(aroma_markers),

        "strong_EPS_capsule_like_loci":
            eps["A_strong"],

        "probable_EPS_capsule_like_loci":
            eps["B_probable"],

        "single_marker_EPS_contexts":
            eps["C_single"],

        "EPS_capsule_evidence": eps_status,

        **stress_summary[mag],
    }

    matrix.append(row)


# ============================================================
# PREVALENCIA DE MARCADORES
# ============================================================

prev = defaultdict(
    lambda: {
        "genes": 0,
        "mags": set()
    }
)

for h in marker_hits:
    key = (h["block"], h["marker"])

    prev[key]["genes"] += 1
    prev[key]["mags"].add(h["MAG"])


prevalence_rows = []

for (block, marker), v in sorted(prev.items()):

    mags = sorted(v["mags"])

    prevalence_rows.append({
        "block": block,
        "marker": marker,
        "n_MAGs": len(mags),
        "percent_MAGs": f"{100 * len(mags) / 18:.2f}",
        "n_genes": v["genes"],
        "MAGs": ";".join(mags),
    })


# ============================================================
# BACTERIOCIN FOCUS
# ============================================================

bacteriocin_roles = {
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

for row in matrix:

    mag = row["MAG"]

    if mag not in bacteriocin_roles:
        continue

    rr = dict(row)
    rr["bacteriocin_analysis_role"] = bacteriocin_roles[mag]

    focus_rows.append(rr)


# ============================================================
# ESCRITURA
# ============================================================

hit_fields = [
    "block",
    "marker",
    "MAG",
    "gene_id",
    "contig",
    "gene_index",
    "Preferred_name",
    "Description",
    "KEGG_ko",
    "EC",
    "PFAMs",
]

write_tsv(
    os.path.join(
        OUTDIR,
        "strict_remaining_function_marker_hits.tsv"
    ),
    sorted(
        marker_hits,
        key=lambda x: (
            x["block"],
            x["marker"],
            x["MAG"],
            x["contig"],
            str(x["gene_index"])
        )
    ),
    hit_fields
)


write_tsv(
    os.path.join(
        OUTDIR,
        "strict_remaining_function_marker_prevalence.tsv"
    ),
    prevalence_rows,
    [
        "block",
        "marker",
        "n_MAGs",
        "percent_MAGs",
        "n_genes",
        "MAGs",
    ]
)


matrix_fields = [
    "MAG",
    "genus",
    "species",
    "completeness",

    "explicit_lipase_candidate_n",
    "explicit_esterase_candidate_n",
    "phospholipase_n",
    "lipolysis_evidence",

    "branched_chain_aminotransferase_n",
    "aromatic_amino_acid_aminotransferase_n",
    "methionine_gamma_lyase_n",
    "sulfur_amino_acid_lyase_candidate_n",
    "ketoacid_decarboxylase_n",
    "n_aroma_marker_types",
    "aroma_marker_types",

    "strong_EPS_capsule_like_loci",
    "probable_EPS_capsule_like_loci",
    "single_marker_EPS_contexts",
    "EPS_capsule_evidence",

    "acid_stress_evidence",
    "acid_stress_systems",

    "osmotic_salt_stress_evidence",
    "osmotic_salt_mechanisms",

    "oxidative_stress_evidence",
    "oxidative_stress_mechanisms",
]


write_tsv(
    os.path.join(
        OUTDIR,
        "remaining_cheese_function_matrix_18MAGs.tsv"
    ),
    matrix,
    matrix_fields
)


write_tsv(
    os.path.join(
        OUTDIR,
        "EPS_capsule_like_locus_candidates.tsv"
    ),
    eps_loci,
    [
        "locus_id",
        "MAG",
        "genus",
        "species",
        "contig",
        "window_start_ORF",
        "window_end_ORF",
        "anchor_min_ORF",
        "anchor_max_ORF",
        "n_anchor_genes",
        "n_glycosyltransferase_genes",
        "marker_classes",
        "marker_gene_ids",
        "glycosyltransferase_gene_ids",
        "curated_tier",
    ]
)


focus_fields = [
    "MAG",
    "genus",
    "species",
    "completeness",
    "bacteriocin_analysis_role",
] + [
    x for x in matrix_fields
    if x not in {
        "MAG",
        "genus",
        "species",
        "completeness",
    }
]


write_tsv(
    os.path.join(
        OUTDIR,
        "remaining_functions_bacteriocin_focus_MAGs.tsv"
    ),
    focus_rows,
    focus_fields
)


# ============================================================
# README
# ============================================================

readme = os.path.join(
    OUTDIR,
    "README_step79a.txt"
)

with open(readme, "w", encoding="utf-8") as fh:

    fh.write(
"""PASO 79a — AUDITORÍA REFINADA DE FUNCIONES RESTANTES DEL QUESO

Objetivo
========
Reducir falsos positivos de los paneles amplios del Paso 76 para:
1. lipólisis/esterasas,
2. metabolismo de aminoácidos relacionado con aroma,
3. EPS/polisacárido capsular,
4. estrés ácido,
5. estrés osmótico/salino,
6. estrés oxidativo.

Principios de interpretación
============================
- Los resultados representan potencial genómico, no actividad.
- Una esterase no equivale automáticamente a lipólisis de grasa láctea.
- Una lipasa anotada se mantiene como candidato mientras no exista
  validación bioquímica/localización experimental.
- Los marcadores de aminoácidos indican capacidad potencial de generar
  intermediarios relacionados con aroma; no demuestran producción del
  compuesto volátil final.
- Los loci clasificados como EPS_capsule_like pueden corresponder a
  exopolisacárido, cápsula u otros polisacáridos de superficie.
  No deben denominarse EPS lácteo específico sin evidencia adicional.
- Los genes de estrés se interpretan como mecanismos recuperados en el MAG.
- Ausencia de marcador no equivale necesariamente a ausencia biológica.
- La completitud del MAG debe considerarse en toda inferencia negativa.

Este paso NO reanota proteínas ni ejecuta eggNOG nuevamente.
"""
    )


# ============================================================
# RESUMEN FINAL
# ============================================================

n_explicit_lipase = sum(
    1 for r in matrix
    if r["explicit_lipase_candidate_n"] > 0
)

n_aroma = sum(
    1 for r in matrix
    if r["n_aroma_marker_types"] > 0
)

n_eps_strong = sum(
    1 for r in matrix
    if r["strong_EPS_capsule_like_loci"] > 0
)

n_eps_probable = sum(
    1 for r in matrix
    if r["probable_EPS_capsule_like_loci"] > 0
)

n_acid_system = sum(
    1 for r in matrix
    if r["acid_stress_evidence"] == "system_supported"
)

n_osmo = sum(
    1 for r in matrix
    if r["osmotic_salt_stress_evidence"] == "supported_mechanism"
)

n_oxid_broad = sum(
    1 for r in matrix
    if r["oxidative_stress_evidence"] == "broad_antioxidant_repertoire"
)


print("=" * 60)
print("PASO 79a COMPLETADO")
print("=" * 60)

print(f"MAGs analizados:                    {len(matrix)}")
print(f"Hits estrictos totales:             {len(marker_hits)}")
print(f"Loci EPS/cápsula auditados:          {len(eps_loci)}")
print()

print(
    f"Lipasa explícita candidata:          "
    f"{n_explicit_lipase}/18 MAGs"
)

print(
    f"Marcadores de aroma AA:              "
    f"{n_aroma}/18 MAGs"
)

print(
    f"EPS/cápsula-like fuerte:              "
    f"{n_eps_strong}/18 MAGs"
)

print(
    f"EPS/cápsula-like probable:            "
    f"{n_eps_probable}/18 MAGs"
)

print(
    f"Sistema de resistencia ácida:         "
    f"{n_acid_system}/18 MAGs"
)

print(
    f"Mecanismo osmótico/salino:            "
    f"{n_osmo}/18 MAGs"
)

print(
    f"Repertorio antioxidante amplio:        "
    f"{n_oxid_broad}/18 MAGs"
)

print()
print(f"Salida: {OUTDIR}")
print(
    "PASO 79a FINALIZÓ CORRECTAMENTE; "
    "ES SEGURO SALIR."
)
