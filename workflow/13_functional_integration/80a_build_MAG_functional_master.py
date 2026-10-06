#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import csv
import sys

# ============================================================
# PASO 80a
# Integración funcional final por MAG
#
# Integra:
#   - calidad + taxonomía final dRep
#   - metabolismo refinado (77)
#   - aminas biógenas curadas (77b)
#   - sistema proteolítico (78b)
#   - proteinasa de envoltura / CEP-like (78c)
#   - lipólisis, aroma, EPS/cápsula y estrés (79b)
#   - rol de cada MAG en los loci candidatos de bacteriocinas
#
# IMPORTANTE:
#   * No calcula abundancia por muestra.
#   * No calcula expresión.
#   * No genera un "score tecnológico".
#   * ATTRLOC002 se mantiene como contexto competitivo compartido,
#     no como locus independiente.
# ============================================================

ROOT = Path("/scratch/global") / Path.home().name / "Shotgun_MAGs_Turrialba"

BASE = ROOT / "22_final_representative_mags" / "representative_mags_master.tsv"

METAB = (
    ROOT
    / "55_cheese_metabolic_capabilities"
    / "metabolic_capability_matrix_18MAGs.tsv"
)

BA = (
    ROOT
    / "55_cheese_metabolic_capabilities"
    / "curated_biogenic_amine_matrix_18MAGs.tsv"
)

PROT = (
    ROOT
    / "57_proteolytic_system_reconstruction"
    / "proteolytic_capability_matrix_18MAGs.tsv"
)

CEP = (
    ROOT
    / "58_surface_proteinase_curated"
    / "surface_proteinase_summary_by_MAG.tsv"
)

REM = (
    ROOT
    / "60_remaining_cheese_function_curated"
    / "remaining_cheese_functions_curated_matrix_18MAGs.tsv"
)

OUT = ROOT / "61_MAG_functional_master"
OUT.mkdir(parents=True, exist_ok=True)


# ============================================================
# Utilidades
# ============================================================

def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def read_tsv(path):
    if not path.exists():
        die(f"No existe archivo requerido: {path}")

    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        if reader.fieldnames is None:
            die(f"Archivo sin encabezado: {path}")

        rows = list(reader)

    return reader.fieldnames, rows


def index_rows(rows, key, source_name):
    out = {}

    for row in rows:
        value = row.get(key, "").strip()

        if not value:
            die(f"{source_name}: valor vacío en columna clave '{key}'")

        if value in out:
            die(f"{source_name}: MAG duplicado: {value}")

        out[value] = row

    return out


def require_columns(fieldnames, required, source_name):
    missing = [c for c in required if c not in fieldnames]

    if missing:
        die(
            f"{source_name}: faltan columnas requeridas: "
            + ", ".join(missing)
        )


def as_int(value):
    value = str(value).strip()
    if value in ("", "-", "NA", "nan", "None"):
        return 0
    try:
        return int(float(value))
    except ValueError:
        die(f"No se pudo convertir a entero: '{value}'")


def add_selected(master_row, source_row, fields):
    for field in fields:
        if field not in source_row:
            die(f"Columna ausente durante integración: {field}")

        if field in master_row:
            die(f"Columna duplicada en matriz maestra: {field}")

        master_row[field] = source_row[field]


def write_tsv(path, rows, fieldnames):
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=fieldnames,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


# ============================================================
# Lectura
# ============================================================

base_fields, base_rows = read_tsv(BASE)
metab_fields, metab_rows = read_tsv(METAB)
ba_fields, ba_rows = read_tsv(BA)
prot_fields, prot_rows = read_tsv(PROT)
cep_fields, cep_rows = read_tsv(CEP)
rem_fields, rem_rows = read_tsv(REM)


# ============================================================
# Validación de encabezados
# ============================================================

require_columns(
    base_fields,
    [
        "representative_MAG",
        "completeness",
        "contamination",
        "domain",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "classification",
    ],
    "representative_mags_master.tsv",
)

require_columns(
    metab_fields,
    [
        "MAG",
        "galactose_Leloir",
        "lactose_transport_betaGal",
        "lactose_PTS_LacEFG",
        "lactate_LDH",
        "acetate_PtaAckA",
        "formate_PFL",
        "citrate_fermentation",
        "acetoin_branch",
        "butanediol_branch",
        "n_strict_biogenic_amine_types",
        "strict_biogenic_amine_types",
    ],
    "metabolic_capability_matrix_18MAGs.tsv",
)

require_columns(
    ba_fields,
    [
        "MAG",
        "histamine",
        "tyramine",
        "cadaverine",
        "putrescine",
        "n_context_or_system_supported_BA_types",
        "n_enzyme_only_BA_types",
    ],
    "curated_biogenic_amine_matrix_18MAGs.tsv",
)

require_columns(
    prot_fields,
    [
        "MAG",
        "Opp_complete_clusters",
        "Opp_partial_clusters",
        "Dpp_complete_clusters",
        "Dpp_partial_clusters",
        "POT_Dtp_genes",
        "n_key_peptidase_types",
        "key_peptidase_types",
        "peptide_utilization_evidence",
    ],
    "proteolytic_capability_matrix_18MAGs.tsv",
)

require_columns(
    cep_fields,
    [
        "MAG",
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
    ],
    "surface_proteinase_summary_by_MAG.tsv",
)

require_columns(
    rem_fields,
    [
        "MAG",
        "lipolysis_evidence_curated",
        "amino_acid_aroma_evidence_curated",
        "Ehrlich_route_evidence",
        "EPS_capsule_evidence_curated",
        "acid_stress_evidence",
        "acid_stress_systems",
        "n_complete_osmotic_systems",
        "complete_osmotic_systems",
        "standalone_osmolyte_transporters",
        "BetT_precursor_transport_n",
        "osmotic_stress_evidence_curated",
        "oxidative_stress_evidence",
        "oxidative_stress_mechanisms",
    ],
    "remaining_cheese_functions_curated_matrix_18MAGs.tsv",
)


# ============================================================
# Índices
# ============================================================

base = index_rows(
    base_rows,
    "representative_MAG",
    "representative_mags_master.tsv",
)

metab = index_rows(
    metab_rows,
    "MAG",
    "metabolic_capability_matrix_18MAGs.tsv",
)

ba = index_rows(
    ba_rows,
    "MAG",
    "curated_biogenic_amine_matrix_18MAGs.tsv",
)

prot = index_rows(
    prot_rows,
    "MAG",
    "proteolytic_capability_matrix_18MAGs.tsv",
)

cep = index_rows(
    cep_rows,
    "MAG",
    "surface_proteinase_summary_by_MAG.tsv",
)

rem = index_rows(
    rem_rows,
    "MAG",
    "remaining_cheese_functions_curated_matrix_18MAGs.tsv",
)


# ============================================================
# Validar que los 6 archivos describen los mismos 18 MAGs
# ============================================================

master_mags = set(base)

if len(master_mags) != 18:
    die(
        f"Se esperaban 18 MAGs finales en master; "
        f"se encontraron {len(master_mags)}"
    )

for source_name, idx in [
    ("metabolismo 77", metab),
    ("aminas 77b", ba),
    ("proteólisis 78b", prot),
    ("CEP 78c", cep),
    ("funciones 79b", rem),
]:
    source_mags = set(idx)

    missing = sorted(master_mags - source_mags)
    extra = sorted(source_mags - master_mags)

    if missing or extra:
        die(
            f"{source_name}: conjunto de MAGs no coincide.\n"
            f"  Faltantes: {missing}\n"
            f"  Extra:     {extra}"
        )


# ============================================================
# Rol bacteriocina
#
# Seis loci biológicos independientes:
# ATTRLOC001
# ATTRLOC003
# ATTRLOC004
# ATTRLOC005
# ATTRLOC006
# ATTRLOC007
#
# ATTRLOC002 es contexto competitivo compartido y NO se suma.
# ============================================================

BACTERIOCIN_ROLE = {
    "L2__L2_maxbin2.004_sub": {
        "bacteriocin_analysis_role":
            "source_MAG_ATTRLOC001",
        "independent_bacteriocin_loci_n": "1",
        "independent_bacteriocin_loci":
            "ATTRLOC001",
        "shared_competitive_context_refs": "",
    },

    "L3__concoct_29": {
        "bacteriocin_analysis_role":
            "source_MAG_ATTRLOC003_ATTRLOC004",
        "independent_bacteriocin_loci_n": "2",
        "independent_bacteriocin_loci":
            "ATTRLOC003;ATTRLOC004",
        "shared_competitive_context_refs": "",
    },

    "M2__M2_maxbin2.004": {
        "bacteriocin_analysis_role":
            "source_MAG_ATTRLOC005_ATTRLOC006_ATTRLOC007",
        "independent_bacteriocin_loci_n": "3",
        "independent_bacteriocin_loci":
            "ATTRLOC005;ATTRLOC006;ATTRLOC007",
        "shared_competitive_context_refs": "",
    },

    "L2__L2_maxbin2.011_sub": {
        "bacteriocin_analysis_role":
            "competitive_shared_context_ATTRLOC002_not_independent",
        "independent_bacteriocin_loci_n": "0",
        "independent_bacteriocin_loci": "",
        "shared_competitive_context_refs":
            "ATTRLOC002",
    },
}


# ============================================================
# Columnas seleccionadas
# ============================================================

BASE_OUT_FIELDS = [
    "MAG",
    "completeness",
    "contamination",
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "classification",
    "representative_source",
    "drep_cluster",
    "n_cluster_members",
]

BACT_FIELDS = [
    "bacteriocin_analysis_role",
    "independent_bacteriocin_loci_n",
    "independent_bacteriocin_loci",
    "shared_competitive_context_refs",
]

METAB_FIELDS = [
    "galactose_Leloir",
    "lactose_transport_betaGal",
    "lactose_PTS_LacEFG",
    "lactate_LDH",
    "acetate_PtaAckA",
    "formate_PFL",
    "citrate_fermentation",
    "acetoin_branch",
    "butanediol_branch",
]

BA_FIELDS = [
    "histamine",
    "tyramine",
    "cadaverine",
    "putrescine",
    "n_context_or_system_supported_BA_types",
    "n_enzyme_only_BA_types",
]

PROT_FIELDS = [
    "Opp_complete_clusters",
    "Opp_partial_clusters",
    "Dpp_complete_clusters",
    "Dpp_partial_clusters",
    "POT_Dtp_genes",
    "n_key_peptidase_types",
    "key_peptidase_types",
    "peptide_utilization_evidence",
]

CEP_FIELDS = [
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
]

REM_FIELDS = [
    # Lipólisis
    "specific_lipase_candidate_n",
    "esterase_lipase_candidate_n",
    "broad_GDSL_acylhydrolase_n",
    "lipase_accessory_n",
    "lipolysis_evidence_curated",

    # Aroma de aminoácidos
    "BCAT_IlvE_n",
    "AraT_n",
    "TyrB_n",
    "MGL_K01761_like_n",
    "cheese_relevant_CBL_like_n",
    "MetC_n",
    "MetB_n",
    "amino_acid_aroma_evidence_curated",
    "Ehrlich_route_evidence",

    # EPS / cápsula
    "EPS_like_loci_n",
    "capsule_like_loci_n",
    "mixed_EPS_capsule_like_loci_n",
    "generic_polysaccharide_loci_n",
    "EPS_capsule_evidence_curated",

    # Estrés ácido
    "acid_stress_evidence",
    "acid_stress_systems",

    # Estrés osmótico / sal
    "n_complete_osmotic_systems",
    "complete_osmotic_systems",
    "standalone_osmolyte_transporters",
    "BetT_precursor_transport_n",
    "osmotic_stress_evidence_curated",

    # Estrés oxidativo
    "oxidative_stress_evidence",
    "oxidative_stress_mechanisms",
]


# ============================================================
# Construir matriz maestra
# ============================================================

master_rows = []

for mag in sorted(master_mags):

    b = base[mag]

    row = {
        "MAG": mag,
        "completeness": b.get("completeness", ""),
        "contamination": b.get("contamination", ""),
        "domain": b.get("domain", ""),
        "phylum": b.get("phylum", ""),
        "class": b.get("class", ""),
        "order": b.get("order", ""),
        "family": b.get("family", ""),
        "genus": b.get("genus", ""),
        "species": b.get("species", ""),
        "classification": b.get("classification", ""),
        "representative_source": b.get("representative_source", ""),
        "drep_cluster": b.get("drep_cluster", ""),
        "n_cluster_members": b.get("n_cluster_members", ""),
    }

    role = BACTERIOCIN_ROLE.get(
        mag,
        {
            "bacteriocin_analysis_role": "none",
            "independent_bacteriocin_loci_n": "0",
            "independent_bacteriocin_loci": "",
            "shared_competitive_context_refs": "",
        },
    )

    row.update(role)

    add_selected(row, metab[mag], METAB_FIELDS)
    add_selected(row, ba[mag], BA_FIELDS)

    # Resumen específico de aminas biógenas
    n_ctx = as_int(
        ba[mag]["n_context_or_system_supported_BA_types"]
    )
    n_enzyme = as_int(
        ba[mag]["n_enzyme_only_BA_types"]
    )

    if n_ctx > 0:
        row["biogenic_amine_evidence_summary"] = (
            "context_or_system_supported"
        )
    elif n_enzyme > 0:
        row["biogenic_amine_evidence_summary"] = (
            "enzyme_only_evidence"
        )
    else:
        row["biogenic_amine_evidence_summary"] = (
            "no_specific_BA_evidence"
        )

    add_selected(row, prot[mag], PROT_FIELDS)
    add_selected(row, cep[mag], CEP_FIELDS)
    add_selected(row, rem[mag], REM_FIELDS)

    master_rows.append(row)


# ============================================================
# Orden final de columnas
# ============================================================

FINAL_FIELDS = (
    BASE_OUT_FIELDS
    + BACT_FIELDS
    + METAB_FIELDS
    + BA_FIELDS
    + ["biogenic_amine_evidence_summary"]
    + PROT_FIELDS
    + CEP_FIELDS
    + REM_FIELDS
)

if len(FINAL_FIELDS) != len(set(FINAL_FIELDS)):
    dup = [
        x for x, n in Counter(FINAL_FIELDS).items()
        if n > 1
    ]
    die(f"Columnas duplicadas en FINAL_FIELDS: {dup}")


# ============================================================
# Validaciones finales
# ============================================================

if len(master_rows) != 18:
    die(
        f"Matriz maestra debería contener 18 filas; "
        f"contiene {len(master_rows)}"
    )

total_independent_loci = sum(
    as_int(r["independent_bacteriocin_loci_n"])
    for r in master_rows
)

if total_independent_loci != 6:
    die(
        f"Se esperaban 6 loci independientes de bacteriocinas; "
        f"se obtuvieron {total_independent_loci}"
    )

source_mags = [
    r["MAG"]
    for r in master_rows
    if as_int(r["independent_bacteriocin_loci_n"]) > 0
]

if len(source_mags) != 3:
    die(
        f"Se esperaban 3 MAGs fuente de loci independientes; "
        f"se obtuvieron {len(source_mags)}"
    )


# ============================================================
# 1. Matriz maestra
# ============================================================

MASTER_OUT = OUT / "MAG_functional_master_18MAGs.tsv"

write_tsv(
    MASTER_OUT,
    master_rows,
    FINAL_FIELDS,
)


# ============================================================
# 2. MAGs de interés para bacteriocinas
# ============================================================

focus_rows = [
    r for r in master_rows
    if r["bacteriocin_analysis_role"] != "none"
]

FOCUS_OUT = OUT / "bacteriocin_focus_MAG_functional_master.tsv"

write_tsv(
    FOCUS_OUT,
    focus_rows,
    FINAL_FIELDS,
)


# ============================================================
# 3. Tabla larga de estados funcionales
# ============================================================

FEATURE_DOMAINS = {
    # Fermentación / carbohidratos
    "galactose_Leloir": "carbohydrate_fermentation",
    "lactose_transport_betaGal": "carbohydrate_fermentation",
    "lactose_PTS_LacEFG": "carbohydrate_fermentation",
    "lactate_LDH": "carbohydrate_fermentation",
    "acetate_PtaAckA": "fermentation_end_products",
    "formate_PFL": "fermentation_end_products",
    "citrate_fermentation": "citrate_aroma",
    "acetoin_branch": "citrate_aroma",
    "butanediol_branch": "citrate_aroma",

    # Aminas
    "histamine": "biogenic_amines",
    "tyramine": "biogenic_amines",
    "cadaverine": "biogenic_amines",
    "putrescine": "biogenic_amines",
    "biogenic_amine_evidence_summary": "biogenic_amines",

    # Proteólisis
    "peptide_utilization_evidence": "proteolysis",
    "best_surface_proteinase_evidence": "proteolysis",

    # Lipólisis
    "lipolysis_evidence_curated": "lipolysis",

    # Aroma AA
    "amino_acid_aroma_evidence_curated": "amino_acid_aroma",
    "Ehrlich_route_evidence": "amino_acid_aroma",

    # EPS
    "EPS_capsule_evidence_curated": "surface_polysaccharides",

    # Estrés
    "acid_stress_evidence": "stress_adaptation",
    "osmotic_stress_evidence_curated": "stress_adaptation",
    "oxidative_stress_evidence": "stress_adaptation",

    # Bacteriocina
    "bacteriocin_analysis_role": "bacteriocin_context",
}

long_rows = []

for r in master_rows:
    for feature, domain_name in FEATURE_DOMAINS.items():
        long_rows.append(
            {
                "MAG": r["MAG"],
                "genus": r["genus"],
                "species": r["species"],
                "completeness": r["completeness"],
                "contamination": r["contamination"],
                "functional_domain": domain_name,
                "feature": feature,
                "state": r.get(feature, ""),
                "bacteriocin_analysis_role":
                    r["bacteriocin_analysis_role"],
            }
        )

LONG_OUT = OUT / "MAG_functional_master_long.tsv"

write_tsv(
    LONG_OUT,
    long_rows,
    [
        "MAG",
        "genus",
        "species",
        "completeness",
        "contamination",
        "functional_domain",
        "feature",
        "state",
        "bacteriocin_analysis_role",
    ],
)


# ============================================================
# 4. Prevalencia de estados funcionales
# ============================================================

prevalence_rows = []

for feature, domain_name in FEATURE_DOMAINS.items():

    by_state = {}

    for r in master_rows:
        state = r.get(feature, "").strip()

        if not state:
            state = "blank"

        by_state.setdefault(state, []).append(r["MAG"])

    for state in sorted(by_state):
        mags = sorted(by_state[state])

        prevalence_rows.append(
            {
                "functional_domain": domain_name,
                "feature": feature,
                "state": state,
                "n_MAGs": str(len(mags)),
                "percent_MAGs": f"{100.0 * len(mags) / 18:.2f}",
                "MAGs": ";".join(mags),
            }
        )

PREV_OUT = OUT / "MAG_functional_state_prevalence.tsv"

write_tsv(
    PREV_OUT,
    prevalence_rows,
    [
        "functional_domain",
        "feature",
        "state",
        "n_MAGs",
        "percent_MAGs",
        "MAGs",
    ],
)


# ============================================================
# 5. Resumen compacto de los tres MAGs fuente independientes
# ============================================================

SOURCE_SUMMARY_FIELDS = [
    "MAG",
    "genus",
    "species",
    "completeness",
    "contamination",
    "independent_bacteriocin_loci_n",
    "independent_bacteriocin_loci",

    "lactose_transport_betaGal",
    "lactose_PTS_LacEFG",
    "galactose_Leloir",
    "lactate_LDH",
    "citrate_fermentation",
    "acetoin_branch",
    "butanediol_branch",

    "peptide_utilization_evidence",
    "best_surface_proteinase_evidence",

    "lipolysis_evidence_curated",
    "amino_acid_aroma_evidence_curated",
    "Ehrlich_route_evidence",

    "EPS_capsule_evidence_curated",

    "acid_stress_evidence",
    "acid_stress_systems",

    "osmotic_stress_evidence_curated",
    "complete_osmotic_systems",

    "oxidative_stress_evidence",
    "oxidative_stress_mechanisms",

    "biogenic_amine_evidence_summary",
    "histamine",
    "tyramine",
    "cadaverine",
    "putrescine",
]

source_rows = [
    r for r in master_rows
    if as_int(r["independent_bacteriocin_loci_n"]) > 0
]

SOURCE_OUT = OUT / "independent_bacteriocin_source_MAG_profiles.tsv"

write_tsv(
    SOURCE_OUT,
    source_rows,
    SOURCE_SUMMARY_FIELDS,
)


# ============================================================
# 6. README
# ============================================================

README = OUT / "README_step80a.txt"

with README.open("w", encoding="utf-8") as fh:
    fh.write(
        "PASO 80a - MATRIZ MAESTRA FUNCIONAL POR MAG\n"
        "===========================================\n\n"

        "Este directorio integra las capacidades funcionales "
        "curadas de los 18 MAGs representativos finales.\n\n"

        "ARCHIVOS PRINCIPALES\n"
        "--------------------\n"
        "MAG_functional_master_18MAGs.tsv\n"
        "  Matriz principal, una fila por MAG.\n\n"

        "MAG_functional_master_long.tsv\n"
        "  Formato largo de estados funcionales, útil para "
        "visualización y análisis posterior.\n\n"

        "MAG_functional_state_prevalence.tsv\n"
        "  Prevalencia descriptiva de cada estado funcional "
        "entre los 18 MAGs.\n\n"

        "bacteriocin_focus_MAG_functional_master.tsv\n"
        "  Cuatro MAGs vinculados al análisis competitivo de "
        "loci de bacteriocinas.\n\n"

        "independent_bacteriocin_source_MAG_profiles.tsv\n"
        "  Los tres MAGs que contienen los seis loci biológicos "
        "independientes considerados en la integración final.\n\n"

        "INTERPRETACIÓN\n"
        "--------------\n"
        "1. Esta matriz representa POTENCIAL GENÓMICO recuperado.\n"
        "2. No demuestra expresión, actividad enzimática ni "
        "producción de metabolitos.\n"
        "3. Un estado ausente/no_detected significa que el marcador "
        "no fue recuperado en el MAG bajo las reglas utilizadas; "
        "no prueba ausencia biológica absoluta.\n"
        "4. No se calcula un puntaje tecnológico global porque "
        "implicaría ponderaciones arbitrarias entre funciones no "
        "equivalentes.\n"
        "5. ATTRLOC002 se conserva únicamente como contexto "
        "competitivo compartido asociado a Atlantibacter y no "
        "se cuenta entre los seis loci biológicos independientes.\n"
        "6. Los seis loci independientes se distribuyen en tres "
        "MAGs fuente: 1 en L. petauri, 2 en L. laudensis y "
        "3 en L. lactis.\n"
        "7. La proyección de estas capacidades a las 18 muestras "
        "se realizará en el Paso 80b usando evidencia de "
        "detección/abundancia de los MAGs.\n"
    )


# ============================================================
# Resumen final
# ============================================================

print("=" * 60)
print("PASO 80a COMPLETADO")
print("=" * 60)
print(f"MAGs integrados:                    {len(master_rows)}")
print(f"MAGs fuente bacteriocina:           {len(source_rows)}")
print(f"Loci independientes integrados:     {total_independent_loci}")
print(f"MAGs en foco competitivo total:     {len(focus_rows)}")
print(f"Variables funcionales matriz:       {len(FINAL_FIELDS)}")
print(f"Filas tabla larga:                  {len(long_rows)}")
print()
print(f"Salida: {OUT}")
print(
    "PASO 80a FINALIZÓ CORRECTAMENTE; "
    "ES SEGURO SALIR."
)
