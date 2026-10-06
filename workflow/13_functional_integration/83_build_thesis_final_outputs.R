#!/usr/bin/env Rscript

# ============================================================
# PASO 83
# Productos finales integrados para tesis
#
# NO recalcula:
# - MAGs
# - taxonomía
# - anotación funcional
# - abundancias
# - loci bacteriocina/RiPP
# - pruebas estadísticas
#
# Genera:
# Table83A: resumen compacto de 18 MAGs
# Table83B: seis loci candidatos independientes
# Table83C: tres MAGs fuente
# Table83D: inferencia longitudinal primaria
#
# Figure83A: composición A/B/C/D por locus
# Figure83B: abundancia longitudinal de MAGs fuente
# Figure83C: funciones tecnológicas de MAGs fuente
# Figure83D: distribuciones exactas de permutación
#
# Además:
# - tabla suplementaria funcional completa
# - manifiesto main/supplement
# - texto base para resultados
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

# ============================================================
# RUTAS
# ============================================================

USER <- Sys.getenv("USER")

ROOT <- file.path(
  "/scratch/global",
  USER,
  "Shotgun_MAGs_Turrialba"
)

MASTER <- file.path(
  ROOT,
  "61_MAG_functional_master",
  "MAG_functional_master_18MAGs.tsv"
)

LOCUS <- file.path(
  ROOT,
  "65_thesis_candidate_synthesis",
  "thesis_candidate_locus_synthesis_6.tsv"
)

SOURCE <- file.path(
  ROOT,
  "65_thesis_candidate_synthesis",
  "thesis_source_MAG_synthesis_3.tsv"
)

TRAJ <- file.path(
  ROOT,
  "64_bacteriocin_functional_integration",
  "source_MAG_locus_trajectories_54rows.tsv"
)

TESTS <- file.path(
  ROOT,
  "66_exact_longitudinal_inference",
  "exact_permutation_producer_tests.tsv"
)

EFFECT <- file.path(
  ROOT,
  "66_exact_longitudinal_inference",
  "exact_test_effect_interpretation.tsv"
)

NULL_FILE <- file.path(
  ROOT,
  "66_exact_longitudinal_inference",
  "exact_permutation_null_distributions.tsv"
)

SUBJECT <- file.path(
  ROOT,
  "66_exact_longitudinal_inference",
  "Llactis_subject_level_longitudinal_metrics.tsv"
)

OUT <- file.path(
  ROOT,
  "67_thesis_final_outputs"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# HELPERS
# ============================================================

required_files <- c(
  MASTER,
  LOCUS,
  SOURCE,
  TRAJ,
  TESTS,
  EFFECT,
  NULL_FILE,
  SUBJECT
)

missing_files <- required_files[
  !file.exists(
    required_files
  )
]

if (length(missing_files) > 0) {

  stop(
    paste(
      "Faltan archivos requeridos:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}

read_clean <- function(path) {

  read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE
  )
}

require_columns <- function(
  x,
  cols,
  label
) {

  missing <- setdiff(
    cols,
    names(x)
  )

  if (length(missing) > 0) {

    stop(
      paste(
        label,
        "faltan columnas:",
        paste(
          missing,
          collapse = ", "
        )
      )
    )
  }
}

# ============================================================
# CARGA
# ============================================================

master <- read_clean(
  MASTER
)

locus <- read_clean(
  LOCUS
)

source <- read_clean(
  SOURCE
)

traj <- read_clean(
  TRAJ
)

tests <- read_clean(
  TESTS
)

effect <- read_clean(
  EFFECT
)

null <- read_clean(
  NULL_FILE
)

subject <- read_clean(
  SUBJECT
)

# ============================================================
# QC
# ============================================================

if (nrow(master) != 18) {
  stop(
    "MASTER no contiene exactamente 18 MAGs."
  )
}

if (nrow(locus) != 6) {
  stop(
    "LOCUS no contiene exactamente 6 loci."
  )
}

if (nrow(source) != 3) {
  stop(
    "SOURCE no contiene exactamente 3 MAGs fuente."
  )
}

if (nrow(traj) != 54) {
  stop(
    "TRAJ no contiene exactamente 54 filas."
  )
}

require_columns(
  master,
  c(
    "MAG",
    "genus",
    "species",
    "completeness",
    "contamination",
    "bacteriocin_analysis_role",
    "galactose_Leloir",
    "lactose_transport_betaGal",
    "lactose_PTS_LacEFG",
    "lactate_LDH",
    "citrate_fermentation",
    "acetoin_branch",
    "butanediol_branch",
    "peptide_utilization_evidence",
    "best_surface_proteinase_evidence",
    "lipolysis_evidence_curated",
    "amino_acid_aroma_evidence_curated",
    "EPS_capsule_evidence_curated",
    "acid_stress_evidence",
    "osmotic_stress_evidence_curated",
    "oxidative_stress_evidence",
    "biogenic_amine_evidence_summary"
  ),
  "MASTER"
)

require_columns(
  locus,
  c(
    "locus_id",
    "source_MAG",
    "source_MAG_species",
    "candidate_names",
    "context_edge_status",
    "contains_interrupted_candidate",
    "primary_source_linked",
    "sensitivity_locus_signal",
    "A_source_linked",
    "B_locus_signal_source_unresolved",
    "C_source_MAG_without_locus_support",
    "D_no_support",
    "locus_MAG_linkage_pattern"
  ),
  "LOCUS"
)

require_columns(
  source,
  c(
    "source_MAG",
    "species",
    "completeness",
    "contamination",
    "n_independent_loci",
    "independent_loci",
    "n_primary_functional_features",
    "primary_functional_features",
    "biogenic_amine_detail",
    "L_source_MAG_median_RA_w0_pct",
    "L_source_MAG_median_RA_w1_pct",
    "L_source_MAG_median_RA_w2_pct",
    "M_source_MAG_median_RA_w0_pct",
    "M_source_MAG_median_RA_w1_pct",
    "M_source_MAG_median_RA_w2_pct"
  ),
  "SOURCE"
)

require_columns(
  traj,
  c(
    "sample",
    "producer",
    "producer_name",
    "biological_unit",
    "week",
    "source_MAG",
    "source_MAG_species",
    "source_MAG_RA_main_pct"
  ),
  "TRAJ"
)

require_columns(
  tests,
  c(
    "hypothesis",
    "analysis_role",
    "outcome",
    "mean_difference_M_minus_L",
    "median_difference_M_minus_L",
    "cliffs_delta_M_vs_L",
    "exact_two_sided_permutation_p"
  ),
  "TESTS"
)

require_columns(
  null,
  c(
    "hypothesis",
    "statistic_mean_M_minus_L"
  ),
  "NULL"
)

# ============================================================
# TABLA 83A
# RESUMEN COMPACTO 18 MAGs
# ============================================================

table83A <- master %>%
  transmute(

    MAG,

    taxon =
      if_else(
        !is.na(species) &
          species != "",
        species,
        genus
      ),

    completeness_pct =
      completeness,

    contamination_pct =
      contamination,

    bacteriocin_role =
      bacteriocin_analysis_role,

    galactose_Leloir,

    lactose_direct =
      lactose_transport_betaGal,

    lactose_PTS =
      lactose_PTS_LacEFG,

    lactate_LDH,

    citrate =
      citrate_fermentation,

    acetoin =
      acetoin_branch,

    butanediol =
      butanediol_branch,

    peptide_utilization =
      peptide_utilization_evidence,

    surface_proteinase =
      best_surface_proteinase_evidence,

    lipolysis =
      lipolysis_evidence_curated,

    amino_acid_aroma =
      amino_acid_aroma_evidence_curated,

    surface_polysaccharide =
      EPS_capsule_evidence_curated,

    acid_stress =
      acid_stress_evidence,

    osmotic_stress =
      osmotic_stress_evidence_curated,

    oxidative_stress =
      oxidative_stress_evidence,

    biogenic_amines =
      biogenic_amine_evidence_summary
  )

write_tsv(
  table83A,
  file.path(
    OUT,
    "Table83A_MAG_summary_18.tsv"
  )
)

# Tabla suplementaria sin reducir

write_tsv(
  master,
  file.path(
    OUT,
    "TableS83_full_MAG_functional_master.tsv"
  )
)

# ============================================================
# TABLA 83B
# SEIS LOCI
# ============================================================

table83B <- locus %>%
  transmute(

    locus_id,

    source_MAG,

    source_species =
      source_MAG_species,

    candidate =
      candidate_names,

    context_status =
      context_edge_status,

    interrupted_candidate =
      contains_interrupted_candidate,

    primary_source_linked_n18 =
      primary_source_linked,

    sensitivity_signal_n18 =
      sensitivity_locus_signal,

    A_source_linked,

    B_source_unresolved =
      B_locus_signal_source_unresolved,

    C_MAG_without_locus =
      C_source_MAG_without_locus_support,

    D_no_support,

    linkage_pattern =
      locus_MAG_linkage_pattern,

    fraction_source_MAG_with_locus =
      fraction_source_MAG_detections_with_locus_support,

    fraction_locus_signal_source_linked =
      fraction_primary_locus_signals_source_linked,

    L_trajectories =
      L_primary_trajectories,

    M_trajectories =
      M_primary_trajectories
  )

write_tsv(
  table83B,
  file.path(
    OUT,
    "Table83B_independent_candidate_loci_6.tsv"
  )
)

# ============================================================
# TABLA 83C
# TRES MAGs FUENTE
# ============================================================

table83C <- source %>%
  transmute(

    source_MAG,

    species,

    completeness_pct =
      completeness,

    contamination_pct =
      contamination,

    n_independent_loci,

    independent_loci,

    n_primary_functional_features,

    primary_functional_features,

    biogenic_amine_evidence =
      biogenic_amine_detail,

    L_RA_week0 =
      L_source_MAG_median_RA_w0_pct,

    L_RA_week1 =
      L_source_MAG_median_RA_w1_pct,

    L_RA_week2 =
      L_source_MAG_median_RA_w2_pct,

    M_RA_week0 =
      M_source_MAG_median_RA_w0_pct,

    M_RA_week1 =
      M_source_MAG_median_RA_w1_pct,

    M_RA_week2 =
      M_source_MAG_median_RA_w2_pct
  )

write_tsv(
  table83C,
  file.path(
    OUT,
    "Table83C_bacteriocin_source_MAGs_3.tsv"
  )
)

# ============================================================
# TABLA 83D
# INFERENCIA PRIMARIA 82b
# ============================================================

table83D <- tests %>%
  filter(
    analysis_role == "PRIMARY"
  ) %>%
  left_join(
    effect %>%
      select(
        hypothesis,
        cliff_magnitude,
        effect_direction
      ),
    by = "hypothesis"
  ) %>%
  select(
    hypothesis,
    outcome,
    mean_difference_M_minus_L,
    median_difference_M_minus_L,
    cliffs_delta_M_vs_L,
    cliff_magnitude,
    effect_direction,
    exact_two_sided_permutation_p,
    minimum_possible_two_sided_p
  )

write_tsv(
  table83D,
  file.path(
    OUT,
    "Table83D_primary_exact_longitudinal_tests.tsv"
  )
)

# ============================================================
# FIGURA 83A
# A/B/C/D POR LOCUS
# ============================================================

locus_long <- locus %>%
  select(
    locus_id,
    A_source_linked,
    B_locus_signal_source_unresolved,
    C_source_MAG_without_locus_support,
    D_no_support
  ) %>%
  pivot_longer(
    cols = -locus_id,
    names_to = "evidence_class",
    values_to = "n_samples"
  ) %>%
  mutate(
    evidence_class =
      factor(
        evidence_class,
        levels = c(
          "A_source_linked",
          "B_locus_signal_source_unresolved",
          "C_source_MAG_without_locus_support",
          "D_no_support"
        ),
        labels = c(
          "A: source-linked",
          "B: locus / source unresolved",
          "C: source MAG / locus unsupported",
          "D: no support"
        )
      )
  )

p83A <- ggplot(
  locus_long,
  aes(
    x = locus_id,
    y = n_samples,
    fill = evidence_class
  )
) +
  geom_col() +
  scale_y_continuous(
    breaks = seq(
      0,
      18,
      3
    ),
    limits = c(
      0,
      18
    )
  ) +
  labs(
    x = NULL,
    y = "Samples (n = 18)",
    fill = "Evidence",
    title =
      "Evidence structure of the six independent candidate loci",
    subtitle =
      "A-source-linked evidence is the primary criterion"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),
    legend.position = "bottom"
  )

ggsave(
  file.path(
    OUT,
    "Figure83A_candidate_locus_evidence.pdf"
  ),
  p83A,
  width = 8.5,
  height = 5.5
)

ggsave(
  file.path(
    OUT,
    "Figure83A_candidate_locus_evidence.png"
  ),
  p83A,
  width = 8.5,
  height = 5.5,
  dpi = 300
)

# ============================================================
# FIGURA 83B
# MAGs FUENTE POR PRODUCTOR / TIEMPO
# ============================================================

source_species_order <- c(
  "Lactococcus petauri",
  "Lactococcus_A laudensis",
  "Lactococcus lactis"
)

traj_plot <- traj %>%
  mutate(
    week =
      as.integer(
        week
      ),

    biological_unit =
      factor(
        biological_unit
      ),

    source_MAG_species =
      factor(
        source_MAG_species,
        levels =
          source_species_order
      )
  )

p83B <- ggplot(
  traj_plot,
  aes(
    x = week,
    y = source_MAG_RA_main_pct,
    group = biological_unit
  )
) +
  geom_line(
    linewidth = 0.75,
    alpha = 0.8
  ) +
  geom_point(
    size = 2.2
  ) +
  facet_grid(
    source_MAG_species ~ producer_name,
    scales = "free_y"
  ) +
  scale_x_continuous(
    breaks = c(
      0,
      1,
      2
    )
  ) +
  labs(
    x = "Week",
    y =
      "Relative abundance within recovered accepted MAGs (%)",
    title =
      "Longitudinal representation of bacteriocin-source MAGs",
    subtitle =
      "Each line represents one biological unit"
  ) +
  theme_bw(
    base_size = 10
  )

ggsave(
  file.path(
    OUT,
    "Figure83B_source_MAG_longitudinal_abundance.pdf"
  ),
  p83B,
  width = 9.5,
  height = 8
)

ggsave(
  file.path(
    OUT,
    "Figure83B_source_MAG_longitudinal_abundance.png"
  ),
  p83B,
  width = 9.5,
  height = 8,
  dpi = 300
)

# ============================================================
# FIGURA 83C
# FUNCIONES DE LOS 3 MAGs FUENTE
#
# Sólo features ya definidos como PRIMARY en 80a/81b.
# No convierte estados parciales en "presencia".
# ============================================================

key_features <- c(
  "galactose_Leloir_core",
  "lactose_transport_betaGal",
  "lactose_PTS_LacEFG",
  "lactate_LDH",
  "acetate_PtaAckA",
  "formate_PFL",
  "acetoin_branch",
  "butanediol_branch",
  "peptide_utilization_system",
  "CEP_like_surface_proteinase",
  "amino_acid_transamination_aroma_potential",
  "sulfur_aroma_candidate",
  "EPS_like_biosynthesis_context",
  "capsule_like_context",
  "acid_stress_system",
  "oxidative_stress_repertoire"
)

feature_labels <- c(
  galactose_Leloir_core =
    "Galactose Leloir",
  lactose_transport_betaGal =
    "Lactose transport + beta-gal",
  lactose_PTS_LacEFG =
    "Lactose PTS",
  lactate_LDH =
    "LDH",
  acetate_PtaAckA =
    "Pta-AckA",
  formate_PFL =
    "PFL",
  acetoin_branch =
    "Acetoin",
  butanediol_branch =
    "2,3-butanediol",
  peptide_utilization_system =
    "Peptide utilization",
  CEP_like_surface_proteinase =
    "CEP-like proteinase",
  amino_acid_transamination_aroma_potential =
    "AA transamination",
  sulfur_aroma_candidate =
    "Sulfur-aroma",
  EPS_like_biosynthesis_context =
    "EPS-like context",
  capsule_like_context =
    "Capsule-like context",
  acid_stress_system =
    "Acid-stress system",
  oxidative_stress_repertoire =
    "Oxidative-stress repertoire"
)

functional_long <- source %>%
  select(
    source_MAG,
    species,
    primary_functional_features
  ) %>%
  mutate(
    feature =
      str_split(
        primary_functional_features,
        ";"
      )
  ) %>%
  unnest(
    feature
  ) %>%
  filter(
    feature %in%
      key_features
  ) %>%
  mutate(
    primary_support = 1L
  ) %>%
  select(
    source_MAG,
    species,
    feature,
    primary_support
  ) %>%
  complete(
    source_MAG,
    species,
    feature =
      key_features,
    fill =
      list(
        primary_support = 0L
      )
  ) %>%
  mutate(
    feature =
      factor(
        feature,
        levels =
          rev(
            key_features
          ),
        labels =
          rev(
            feature_labels[
              key_features
            ]
          )
      ),

    species =
      factor(
        species,
        levels =
          source_species_order
      ),

    support_label =
      if_else(
        primary_support == 1,
        "Primary support",
        "Not primary-supported"
      )
  )

p83C <- ggplot(
  functional_long,
  aes(
    x = species,
    y = feature,
    fill = support_label
  )
) +
  geom_tile(
    linewidth = 0.3
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = NULL,
    title =
      "Primary technological-function evidence in bacteriocin-source MAGs",
    subtitle =
      "Only previously curated primary functional features are shown"
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 30,
        hjust = 1
      ),
    legend.position =
      "bottom"
  )

ggsave(
  file.path(
    OUT,
    "Figure83C_source_MAG_primary_functions.pdf"
  ),
  p83C,
  width = 8,
  height = 7
)

ggsave(
  file.path(
    OUT,
    "Figure83C_source_MAG_primary_functions.png"
  ),
  p83C,
  width = 8,
  height = 7,
  dpi = 300
)

# ============================================================
# FIGURA 83D
# DISTRIBUCIONES DE PERMUTACIÓN H1/H2
# ============================================================

primary_tests <- tests %>%
  filter(
    analysis_role == "PRIMARY"
  ) %>%
  select(
    hypothesis,
    mean_difference_M_minus_L,
    exact_two_sided_permutation_p
  )

null_primary <- null %>%
  filter(
    hypothesis %in%
      primary_tests$hypothesis
  ) %>%
  left_join(
    primary_tests,
    by = "hypothesis"
  ) %>%
  mutate(
    hypothesis_label =
      case_when(
        hypothesis ==
          "H1_Llactis_RA_change_week2_minus_week0" ~
          "H1: change in L. lactis MAG abundance",

        hypothesis ==
          "H2_direct_locus_count_change_week2_minus_week0" ~
          "H2: change in number of source-linked loci",

        TRUE ~
          hypothesis
      )
  )

p83D <- ggplot(
  null_primary,
  aes(
    x =
      statistic_mean_M_minus_L
  )
) +
  geom_histogram(
    bins = 10,
    boundary = 0
  ) +
  geom_vline(
    aes(
      xintercept =
        mean_difference_M_minus_L
    ),
    linewidth = 0.9,
    linetype = 2
  ) +
  facet_wrap(
    ~ hypothesis_label,
    scales = "free",
    ncol = 1
  ) +
  labs(
    x =
      "Mean difference: Marino - Lidieth",
    y =
      "Number of exact assignments",
    title =
      "Exact permutation distributions for the primary longitudinal hypotheses",
    subtitle =
      "All 20 possible 3-vs-3 producer-label assignments"
  ) +
  theme_bw(
    base_size = 10
  )

ggsave(
  file.path(
    OUT,
    "Figure83D_exact_permutation_primary_hypotheses.pdf"
  ),
  p83D,
  width = 8.5,
  height = 7
)

ggsave(
  file.path(
    OUT,
    "Figure83D_exact_permutation_primary_hypotheses.png"
  ),
  p83D,
  width = 8.5,
  height = 7,
  dpi = 300
)

# ============================================================
# MANIFIESTO: MAIN vs SUPPLEMENT
# ============================================================

manifest <- tribble(

  ~output,
  ~recommended_location,
  ~purpose,

  "Table83A_MAG_summary_18.tsv",
  "Supplement / Results reference",
  "Compact overview of the 18 representative MAGs",

  "Table83B_independent_candidate_loci_6.tsv",
  "MAIN TEXT",
  "Primary evidence for the six independent candidate bacteriocin/RiPP loci",

  "Table83C_bacteriocin_source_MAGs_3.tsv",
  "MAIN TEXT",
  "Integrated functional profile of the three independent source MAGs",

  "Table83D_primary_exact_longitudinal_tests.tsv",
  "MAIN TEXT or statistical supplement",
  "Exact exploratory producer comparison for the two predefined longitudinal hypotheses",

  "Figure83A_candidate_locus_evidence.pdf",
  "MAIN TEXT",
  "A/B/C/D evidence composition of candidate loci",

  "Figure83B_source_MAG_longitudinal_abundance.pdf",
  "MAIN TEXT",
  "Biological-unit trajectories of the three source MAGs",

  "Figure83C_source_MAG_primary_functions.pdf",
  "MAIN TEXT",
  "Integrated primary technological-function profiles",

  "Figure83D_exact_permutation_primary_hypotheses.pdf",
  "Statistical supplement or MAIN TEXT",
  "Exact null distributions and observed H1/H2 effects",

  "TableS83_full_MAG_functional_master.tsv",
  "SUPPLEMENT",
  "Full curated functional master matrix"
)

write_tsv(
  manifest,
  file.path(
    OUT,
    "thesis_output_manifest.tsv"
  )
)

# ============================================================
# TEXTO BASE PARA RESULTADOS
# ============================================================

H1 <- tests %>%
  filter(
    hypothesis ==
      "H1_Llactis_RA_change_week2_minus_week0"
  )

H2 <- tests %>%
  filter(
    hypothesis ==
      "H2_direct_locus_count_change_week2_minus_week0"
  )

LAC <- source %>%
  filter(
    source_MAG ==
      "M2__M2_maxbin2.004"
  )

result_text <- c(

  "# Texto base para resultados",

  "",

  "## Candidatos bacteriocina/RiPP",

  "",

  paste0(
    "Se integraron seis loci candidatos independientes asociados ",
    "con tres MAGs de Lactococcus. "
  ),

  paste0(
    "ATTRLOC003 y ATTRLOC004 presentaron concordancia exclusiva ",
    "con su MAG fuente, mientras que ATTRLOC005 mostró la mayor ",
    "prevalencia source-linked y fue recuperado en las 13 muestras ",
    "en las que su MAG fuente cumplió el criterio principal de ",
    "detección."
  ),

  paste0(
    "ATTRLOC001, ATTRLOC006 y ATTRLOC007 mostraron distintos grados ",
    "de discordancia entre la señal del locus y la detección del ",
    "MAG fuente, por lo que estas señales no se interpretaron como ",
    "atribuciones inequívocas a una misma población."
  ),

  "",

  "## Perfil funcional de los MAGs fuente",

  "",

  paste0(
    "El MAG M2__M2_maxbin2.004, clasificado como Lactococcus lactis, ",
    "presentó el repertorio tecnológico integrado más amplio entre ",
    "los tres MAGs fuente, con ",
    LAC$n_primary_functional_features,
    " funciones primarias curadas, tres loci candidatos independientes, ",
    "un sistema de utilización de péptidos y una proteinasa superficial ",
    "CEP-like de alta confianza."
  ),

  paste0(
    "No se recuperó evidencia específica curada de sistemas de ",
    "histamina, tiramina, cadaverina o putrescina en ninguno de los ",
    "tres MAGs fuente independientes."
  ),

  "",

  "## Dinámica longitudinal",

  "",

  paste0(
    "Para M2__M2_maxbin2.004, el cambio de abundancia relativa entre ",
    "las semanas 0 y 2 fue más negativo en Marino que en Lidieth ",
    "(diferencia media de cambios Marino-Lidieth = ",
    round(
      H1$mean_difference_M_minus_L,
      2
    ),
    " puntos porcentuales; Cliff's delta = ",
    H1$cliffs_delta_M_vs_L,
    "; p exacta bilateral = ",
    H1$exact_two_sided_permutation_p,
    ")."
  ),

  paste0(
    "El valor de p = ",
    H1$exact_two_sided_permutation_p,
    " corresponde al mínimo bilateral alcanzable con las 20 ",
    "asignaciones exactas posibles de tres unidades biológicas por ",
    "productor, por lo que la interpretación se centró en la magnitud ",
    "y dirección del efecto."
  ),

  paste0(
    "En contraste, el cambio en el número de loci ATTRLOC005-007 con ",
    "evidencia source-linked mostró un efecto menor ",
    "(Cliff's delta = ",
    round(
      H2$cliffs_delta_M_vs_L,
      3
    ),
    "; p exacta bilateral = ",
    H2$exact_two_sided_permutation_p,
    "), evidenciando que la dinámica de los loci no fue equivalente ",
    "a la dinámica de abundancia del MAG."
  ),

  "",

  "## Limitaciones interpretativas",

  "",

  paste0(
    "Los loci representan contextos genómicos candidatos y no demuestran ",
    "expresión, producción ni actividad antimicrobiana."
  ),

  paste0(
    "Las capacidades funcionales representan potencial genómico ",
    "recuperado y no actividad metabólica medida."
  ),

  paste0(
    "Las abundancias relativas corresponden al conjunto de MAGs ",
    "recuperados y aceptados y no deben interpretarse como porcentaje ",
    "de la microbiota total."
  ),

  paste0(
    "Los seis loci presentan truncamiento en al menos uno de sus ",
    "extremos, por lo que no se consideran BGC completos reconstruidos."
  )
)

writeLines(
  result_text,
  file.path(
    OUT,
    "thesis_results_base_text.md"
  )
)

# ============================================================
# CAPTIONS
# ============================================================

captions <- c(

  "Figure83A. Distribution of evidence classes for the six independent candidate bacteriocin/RiPP loci across the 18 samples. A indicates simultaneous coherent support for the locus and its source MAG; B indicates locus signal without resolved support for the source MAG; C indicates detection of the source MAG without primary locus support; and D indicates no primary support for either component.",

  "",

  "Figure83B. Longitudinal relative abundance of the three MAGs associated with independent candidate bacteriocin/RiPP loci. Each line represents one biological unit followed over weeks 0, 1 and 2. Relative abundance refers only to the set of recovered and accepted MAGs.",

  "",

  "Figure83C. Curated primary technological-function evidence in the three MAGs associated with the six independent candidate loci. Filled cells indicate functions previously classified as primary-supported; absence of a filled cell does not demonstrate biological absence from the original microbial population.",

  "",

  "Figure83D. Exact permutation distributions for the two predefined primary longitudinal hypotheses. The dashed line represents the observed difference between producers. With three biological units per producer, only 20 exact 3-vs-3 assignments are possible and the minimum attainable two-sided p-value is 0.10."
)

writeLines(
  captions,
  file.path(
    OUT,
    "thesis_figure_captions.txt"
  )
)

# ============================================================
# README
# ============================================================

readme <- c(

  "PASO 83 - PRODUCTOS FINALES PARA TESIS",
  "=======================================",
  "",
  "Este paso NO recalcula resultados bioinformaticos.",
  "",
  "Integra únicamente resultados previamente curados y validados.",
  "",
  "PRODUCTOS PRINCIPALES",
  "---------------------",
  "",
  "Table83B_independent_candidate_loci_6.tsv",
  "Table83C_bacteriocin_source_MAGs_3.tsv",
  "Table83D_primary_exact_longitudinal_tests.tsv",
  "",
  "Figure83A_candidate_locus_evidence.pdf/png",
  "Figure83B_source_MAG_longitudinal_abundance.pdf/png",
  "Figure83C_source_MAG_primary_functions.pdf/png",
  "Figure83D_exact_permutation_primary_hypotheses.pdf/png",
  "",
  "SUPLEMENTO",
  "----------",
  "",
  "Table83A_MAG_summary_18.tsv",
  "TableS83_full_MAG_functional_master.tsv",
  "",
  "TEXTO",
  "-----",
  "",
  "thesis_results_base_text.md",
  "thesis_figure_captions.txt",
  "",
  "CAUTELAS",
  "--------",
  "",
  "Los loci candidatos no demuestran expresión ni actividad.",
  "Las funciones representan potencial genómico.",
  "La abundancia MAG no equivale a abundancia de microbiota total.",
  "Los loci truncados no deben denominarse BGC completos.",
  "Los 18 tiempos no son 18 réplicas biológicas independientes."
)

writeLines(
  readme,
  file.path(
    OUT,
    "README_step83.txt"
  )
)

# ============================================================
# CONSOLA FINAL
# ============================================================

cat(
  "============================================================\n"
)

cat(
  "PASO 83 COMPLETADO\n"
)

cat(
  "============================================================\n"
)

cat(
  "MAGs en tabla compacta:              ",
  nrow(
    table83A
  ),
  "\n"
)

cat(
  "Loci candidatos independientes:      ",
  nrow(
    table83B
  ),
  "\n"
)

cat(
  "MAGs fuente independientes:          ",
  nrow(
    table83C
  ),
  "\n"
)

cat(
  "Hipótesis primarias integradas:      ",
  nrow(
    table83D
  ),
  "\n"
)

cat(
  "Figuras principales generadas:       4\n"
)

cat(
  "QC esperado loci A total:            ",
  sum(
    locus$A_source_linked
  ),
  "\n"
)

cat(
  "QC esperado A+B+C+D:                 ",
  sum(
    locus$A_source_linked +
      locus$B_locus_signal_source_unresolved +
      locus$C_source_MAG_without_locus_support +
      locus$D_no_support
  ),
  "\n"
)

cat(
  "\nSalida: ",
  OUT,
  "\n",
  sep = ""
)

cat(
  "PASO 83 FINALIZÓ CORRECTAMENTE; ES SEGURO SALIR.\n"
)
