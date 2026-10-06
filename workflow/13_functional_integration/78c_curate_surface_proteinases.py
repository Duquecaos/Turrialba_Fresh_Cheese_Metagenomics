#!/usr/bin/env python3

import csv
import os
import re
from collections import Counter, defaultdict

USER = os.environ["USER"]
ROOT = f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"

INDIR = os.path.join(
    ROOT,
    "57_proteolytic_system_reconstruction"
)

AUDIT = os.path.join(
    INDIR,
    "surface_proteinase_candidate_audit.tsv"
)

MATRIX = os.path.join(
    INDIR,
    "proteolytic_capability_matrix_18MAGs.tsv"
)

OUTDIR = os.path.join(
    ROOT,
    "58_surface_proteinase_curated"
)

os.makedirs(OUTDIR, exist_ok=True)

OUT_CANDIDATES = os.path.join(
    OUTDIR,
    "curated_surface_proteinase_candidates.tsv"
)

OUT_SUMMARY = os.path.join(
    OUTDIR,
    "surface_proteinase_summary_by_MAG.tsv"
)

OUT_MATRIX = os.path.join(
    OUTDIR,
    "proteolytic_capability_matrix_18MAGs_78c.tsv"
)

OUT_FOCUS = os.path.join(
    OUTDIR,
    "surface_proteinase_bacteriocin_focus_MAGs.tsv"
)

OUT_README = os.path.join(
    OUTDIR,
    "README_78c.txt"
)


# ============================================================
# FUNCIONES AUXILIARES
# ============================================================

def read_tsv(path):
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader), reader.fieldnames


def to_int(value, default=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def yesno(value):
    return 1 if value else 0


def clean(value):
    if value is None:
        return ""
    return str(value).strip()


# ============================================================
# LEER INPUTS
# ============================================================

audit_rows, audit_fields = read_tsv(AUDIT)
matrix_rows, matrix_fields = read_tsv(MATRIX)

required_audit = {
    "MAG",
    "gene_id",
    "contig",
    "gene_index",
    "aa_length",
    "Preferred_name",
    "Description",
    "KEGG_ko",
    "PFAMs",
    "LPXTG_like_Cterminal",
}

missing = required_audit - set(audit_fields)

if missing:
    raise RuntimeError(
        "Faltan columnas requeridas en surface_proteinase_candidate_audit.tsv: "
        + ", ".join(sorted(missing))
    )

required_matrix = {
    "MAG",
    "genus",
    "species",
    "completeness",
    "peptide_utilization_evidence",
}

missing = required_matrix - set(matrix_fields)

if missing:
    raise RuntimeError(
        "Faltan columnas requeridas en proteolytic_capability_matrix_18MAGs.tsv: "
        + ", ".join(sorted(missing))
    )

if len(matrix_rows) != 18:
    raise RuntimeError(
        f"Se esperaban 18 MAGs en la matriz 78b y se encontraron {len(matrix_rows)}"
    )


# ============================================================
# REGLAS DE CURACIÓN
# ============================================================
#
# A = candidato CEP-like de alta confianza
#     - S8/S53
#     - >=1200 aa
#     - anotación explícita prt*
#     - anclaje LPXTG-like o dominio Gram_pos_anchor
#
# B = candidato CEP-like probable
#     - S8/S53
#     - >=1200 aa
#     - anotación explícita prt*
#     - arquitectura multidominio compatible
#     - sin requisito de anclaje recuperado
#
# C = subtilasa grande multidominio
#     - S8/S53
#     - >=1200 aa
#     - arquitectura accesoria
#     - pero sin anotación Prt suficientemente específica
#
# D = otras subtilasas S8/S53
#
# IMPORTANTE:
# Estas categorías representan potencial genómico recuperado.
# No demuestran secreción, expresión ni hidrólisis de caseína.
# ============================================================

explicit_prt_names = {
    "prtp",
    "prtb",
    "prth",
    "prts",
    "prtr",
}

accessory_patterns = {
    "PA": r"(^|,)PA($|,)",
    "SLAP": r"(^|,)SLAP($|,)",
    "fn3": r"fn3",
    "FIVAR": r"FIVAR",
    "Gram_pos_anchor": r"Gram_pos_anchor",
    "PKD": r"(^|,)PKD($|,)",
    "MAM": r"(^|,)MAM($|,)",
    "Big": r"(^|,)Big[_0-9]",
    "CARDB": r"CARDB",
    "Laminin_G": r"Laminin_G",
    "NHL": r"(^|,)NHL($|,)",
}

discordant_description_patterns = [
    r"gluconolactonase",
    r"xylan catabolic",
]

tier_order = {
    "A_high_confidence_CEP_like": 1,
    "B_probable_CEP_like": 2,
    "C_large_multidomain_subtilase_non_CEP": 3,
    "D_other_subtilase": 4,
}


def curate_candidate(row):

    aa = to_int(row["aa_length"])

    pref = clean(row["Preferred_name"])
    pref_l = pref.lower()

    desc = clean(row["Description"])
    desc_l = desc.lower()

    pfams = clean(row["PFAMs"])

    has_S8_S53 = (
        "Peptidase_S8" in pfams
        or "Peptidase_S53" in pfams
        or "subtilase" in desc_l
        or "subtilisin" in desc_l
        or "peptidase s8" in desc_l
        or "peptidase s53" in desc_l
    )

    explicit_prt = (
        pref_l in explicit_prt_names
        or bool(re.fullmatch(r"prt[a-z0-9]+", pref_l))
        or "cell envelope proteinase" in desc_l
        or "surface proteinase" in desc_l
    )

    lpxtg_flag = clean(row["LPXTG_like_Cterminal"]) == "1"

    gram_anchor = "Gram_pos_anchor" in pfams

    anchor_support = lpxtg_flag or gram_anchor

    accessory_hits = []

    for name, pattern in accessory_patterns.items():
        if re.search(pattern, pfams, flags=re.IGNORECASE):
            accessory_hits.append(name)

    discordant = any(
        re.search(pattern, desc_l)
        for pattern in discordant_description_patterns
    )

    reasons = []

    if has_S8_S53:
        reasons.append("S8_S53_supported")

    if aa >= 1200:
        reasons.append("large_>=1200aa")

    if explicit_prt:
        reasons.append("explicit_Prt_annotation")

    if lpxtg_flag:
        reasons.append("LPXTG_like_Cterminal")

    if gram_anchor:
        reasons.append("Gram_pos_anchor")

    if accessory_hits:
        reasons.append(
            "accessory_domains=" + ",".join(accessory_hits)
        )

    if discordant:
        reasons.append("discordant_primary_annotation")

    # --------------------------------------------------------
    # CLASIFICACIÓN CONSERVADORA
    # --------------------------------------------------------

    if (
        has_S8_S53
        and explicit_prt
        and aa >= 1200
        and anchor_support
        and not discordant
    ):
        tier = "A_high_confidence_CEP_like"

    elif (
        has_S8_S53
        and explicit_prt
        and aa >= 1200
        and len(accessory_hits) >= 2
        and not discordant
    ):
        tier = "B_probable_CEP_like"

    elif (
        has_S8_S53
        and aa >= 1200
        and len(accessory_hits) >= 1
    ):
        tier = "C_large_multidomain_subtilase_non_CEP"

    else:
        tier = "D_other_subtilase"

    return {
        **row,
        "has_S8_S53": yesno(has_S8_S53),
        "explicit_Prt_annotation": yesno(explicit_prt),
        "large_ge_1200aa": yesno(aa >= 1200),
        "LPXTG_like_support": yesno(lpxtg_flag),
        "Gram_pos_anchor_support": yesno(gram_anchor),
        "anchor_support_any": yesno(anchor_support),
        "n_accessory_architecture_domains": len(accessory_hits),
        "accessory_architecture_domains": ";".join(accessory_hits),
        "discordant_primary_annotation": yesno(discordant),
        "curated_tier_78c": tier,
        "curation_reason": ";".join(reasons),
    }


curated = [
    curate_candidate(row)
    for row in audit_rows
]

curated.sort(
    key=lambda r: (
        r["MAG"],
        tier_order[r["curated_tier_78c"]],
        -to_int(r["aa_length"]),
        r["gene_id"],
    )
)


# ============================================================
# GUARDAR CANDIDATOS CURADOS
# ============================================================

candidate_fields = audit_fields + [
    "has_S8_S53",
    "explicit_Prt_annotation",
    "large_ge_1200aa",
    "LPXTG_like_support",
    "Gram_pos_anchor_support",
    "anchor_support_any",
    "n_accessory_architecture_domains",
    "accessory_architecture_domains",
    "discordant_primary_annotation",
    "curated_tier_78c",
    "curation_reason",
]

with open(
    OUT_CANDIDATES,
    "w",
    encoding="utf-8",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        delimiter="\t",
        fieldnames=candidate_fields,
        extrasaction="ignore",
    )

    writer.writeheader()
    writer.writerows(curated)


# ============================================================
# RESUMEN POR MAG
# ============================================================

curated_by_mag = defaultdict(list)

for row in curated:
    curated_by_mag[row["MAG"]].append(row)

summary_rows = []

for mrow in matrix_rows:

    mag = mrow["MAG"]

    rows = curated_by_mag.get(mag, [])

    counts = Counter(
        r["curated_tier_78c"]
        for r in rows
    )

    best = None

    if rows:
        best = sorted(
            rows,
            key=lambda r: (
                tier_order[r["curated_tier_78c"]],
                -to_int(r["aa_length"]),
            )
        )[0]

    if counts["A_high_confidence_CEP_like"] > 0:
        evidence = "high_confidence_CEP_like_recovered"

    elif counts["B_probable_CEP_like"] > 0:
        evidence = "probable_CEP_like_recovered"

    elif counts["C_large_multidomain_subtilase_non_CEP"] > 0:
        evidence = "large_multidomain_subtilase_non_CEP_only"

    elif counts["D_other_subtilase"] > 0:
        evidence = "other_subtilase_only"

    else:
        evidence = "no_CEP_like_candidate_recovered"

    summary_rows.append({
        "MAG": mag,
        "genus": mrow.get("genus", ""),
        "species": mrow.get("species", ""),
        "completeness": mrow.get("completeness", ""),
        "A_high_confidence_CEP_like_n":
            counts["A_high_confidence_CEP_like"],
        "B_probable_CEP_like_n":
            counts["B_probable_CEP_like"],
        "C_large_multidomain_subtilase_non_CEP_n":
            counts["C_large_multidomain_subtilase_non_CEP"],
        "D_other_subtilase_n":
            counts["D_other_subtilase"],
        "best_surface_proteinase_evidence": evidence,
        "best_candidate_gene":
            best["gene_id"] if best else "",
        "best_candidate_name":
            best["Preferred_name"] if best else "",
        "best_candidate_aa_length":
            best["aa_length"] if best else "",
        "best_candidate_PFAMs":
            best["PFAMs"] if best else "",
        "best_candidate_tier":
            best["curated_tier_78c"] if best else "",
        "peptide_utilization_evidence":
            mrow.get("peptide_utilization_evidence", ""),
    })


summary_fields = [
    "MAG",
    "genus",
    "species",
    "completeness",
    "A_high_confidence_CEP_like_n",
    "B_probable_CEP_like_n",
    "C_large_multidomain_subtilase_non_CEP_n",
    "D_other_subtilase_n",
    "best_surface_proteinase_evidence",
    "best_candidate_gene",
    "best_candidate_name",
    "best_candidate_aa_length",
    "best_candidate_PFAMs",
    "best_candidate_tier",
    "peptide_utilization_evidence",
]

with open(
    OUT_SUMMARY,
    "w",
    encoding="utf-8",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        delimiter="\t",
        fieldnames=summary_fields,
    )

    writer.writeheader()
    writer.writerows(summary_rows)


# ============================================================
# MATRIZ 78b + CURACIÓN 78c
# ============================================================

summary_by_mag = {
    row["MAG"]: row
    for row in summary_rows
}

extra_fields = [
    "CEP_A_high_confidence_n",
    "CEP_B_probable_n",
    "large_multidomain_subtilase_non_CEP_n",
    "other_subtilase_n",
    "surface_proteinase_evidence_78c",
    "best_surface_proteinase_gene_78c",
    "best_surface_proteinase_tier_78c",
]

integrated_rows = []

for row in matrix_rows:

    mag = row["MAG"]
    s = summary_by_mag[mag]

    out = dict(row)

    out["CEP_A_high_confidence_n"] = \
        s["A_high_confidence_CEP_like_n"]

    out["CEP_B_probable_n"] = \
        s["B_probable_CEP_like_n"]

    out["large_multidomain_subtilase_non_CEP_n"] = \
        s["C_large_multidomain_subtilase_non_CEP_n"]

    out["other_subtilase_n"] = \
        s["D_other_subtilase_n"]

    out["surface_proteinase_evidence_78c"] = \
        s["best_surface_proteinase_evidence"]

    out["best_surface_proteinase_gene_78c"] = \
        s["best_candidate_gene"]

    out["best_surface_proteinase_tier_78c"] = \
        s["best_candidate_tier"]

    integrated_rows.append(out)

with open(
    OUT_MATRIX,
    "w",
    encoding="utf-8",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        delimiter="\t",
        fieldnames=matrix_fields + extra_fields,
        extrasaction="ignore",
    )

    writer.writeheader()
    writer.writerows(integrated_rows)


# ============================================================
# MAGs VINCULADOS AL ANÁLISIS DE BACTERIOCINAS
# ============================================================

focus_roles = {
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

for row in summary_rows:

    mag = row["MAG"]

    if mag not in focus_roles:
        continue

    x = dict(row)
    x["bacteriocin_analysis_role"] = focus_roles[mag]

    focus_rows.append(x)

focus_fields = [
    "MAG",
    "genus",
    "species",
    "completeness",
    "bacteriocin_analysis_role",
    "A_high_confidence_CEP_like_n",
    "B_probable_CEP_like_n",
    "C_large_multidomain_subtilase_non_CEP_n",
    "D_other_subtilase_n",
    "best_surface_proteinase_evidence",
    "best_candidate_gene",
    "best_candidate_name",
    "best_candidate_aa_length",
    "best_candidate_PFAMs",
    "best_candidate_tier",
    "peptide_utilization_evidence",
]

with open(
    OUT_FOCUS,
    "w",
    encoding="utf-8",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        delimiter="\t",
        fieldnames=focus_fields,
    )

    writer.writeheader()

    for mag in [
        "L2__L2_maxbin2.004_sub",
        "L3__concoct_29",
        "M2__M2_maxbin2.004",
        "L2__L2_maxbin2.011_sub",
    ]:
        for row in focus_rows:
            if row["MAG"] == mag:
                writer.writerow(row)


# ============================================================
# README
# ============================================================

with open(
    OUT_README,
    "w",
    encoding="utf-8"
) as fh:

    fh.write(
"""PASO 78c — Curación conservadora de proteinasas superficiales / CEP-like

OBJETIVO
-------
Separar candidatos compatibles con proteinasas de envoltura celular
de otras subtilasas S8/S53 recuperadas en los MAGs.

CATEGORÍAS
----------
A_high_confidence_CEP_like
  S8/S53 + >=1200 aa + anotación Prt explícita +
  evidencia de anclaje LPXTG-like o Gram_pos_anchor.

B_probable_CEP_like
  S8/S53 + >=1200 aa + anotación Prt explícita +
  arquitectura multidominio compatible, pero sin anclaje
  C-terminal inequívoco recuperado.

C_large_multidomain_subtilase_non_CEP
  S8/S53 grande y multidominio, pero sin evidencia suficiente
  para denominarla proteinasa CEP-like.

D_other_subtilase
  Otras subtilasas S8/S53 sin arquitectura suficiente para
  clasificación CEP-like.

LIMITACIONES
------------
La clasificación describe potencial genómico recuperado.
No demuestra expresión, secreción, actividad proteolítica ni
hidrólisis experimental de caseína.

Las proteínas con anotación primaria discordante, por ejemplo
'gluconolactonase activity', no se promueven a CEP-like sólo por
contener un dominio Peptidase_S8.

No se ejecutaron nuevas predicciones de señal peptídica,
transmembrana o localización subcelular en este paso.
"""
    )


# ============================================================
# RESUMEN EN TERMINAL
# ============================================================

global_counts = Counter(
    row["curated_tier_78c"]
    for row in curated
)

mag_evidence_counts = Counter(
    row["best_surface_proteinase_evidence"]
    for row in summary_rows
)

print("=" * 60)
print("PASO 78c COMPLETADO")
print("=" * 60)

print(f"MAGs analizados: {len(summary_rows)}")
print(f"Candidatos S8/S53 auditados: {len(curated)}")
print()

for tier in [
    "A_high_confidence_CEP_like",
    "B_probable_CEP_like",
    "C_large_multidomain_subtilase_non_CEP",
    "D_other_subtilase",
]:
    print(
        f"{tier:42s} "
        f"{global_counts[tier]:3d} genes"
    )

print()
print("EVIDENCIA POR MAG")

for status in [
    "high_confidence_CEP_like_recovered",
    "probable_CEP_like_recovered",
    "large_multidomain_subtilase_non_CEP_only",
    "other_subtilase_only",
    "no_CEP_like_candidate_recovered",
]:
    print(
        f"{status:45s} "
        f"{mag_evidence_counts[status]:2d}/18"
    )

print()
print(f"Salida: {OUTDIR}")
print(
    "PASO 78c FINALIZÓ CORRECTAMENTE; "
    "ES SEGURO SALIR."
)
