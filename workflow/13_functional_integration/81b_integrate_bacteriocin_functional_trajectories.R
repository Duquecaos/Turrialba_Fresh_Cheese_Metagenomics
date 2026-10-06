#!/usr/bin/env Rscript

# ============================================================
# PASO 81b
# Integración longitudinal:
#
# locus candidato bacteriocina/RiPP
#        ↕
# MAG fuente
#        ↕
# capacidades funcionales curadas
#
# FUENTES
# -------
# Step 73:
#   evidencia directa locus-específica
#
# Step 74:
#   trayectoria longitudinal de los loci
#   (utilizada aquí como QC)
#
# Step 80a:
#   perfil funcional por MAG
#
# Step 80b:
#   proyección funcional por muestra
#   + abundancia de los MAGs vinculados a bacteriocinas
#
# PRINCIPIOS
# ----------
# 1. Detección de MAG != detección del locus.
# 2. Detección del locus != expresión ni producción de bacteriocina.
# 3. Potencial funcional != actividad metabólica.
# 4. Las abundancias son relativas dentro del conjunto de MAGs
#    recuperados/aceptados, no de toda la microbiota.
# 5. ATTRLOC002 NO forma parte de los seis loci biológicos
#    independientes.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

# ============================================================
# RUTAS
# ============================================================

user <- Sys.getenv("USER")

root <- file.path(
  "/scratch/global",
  user,
  "Shotgun_MAGs_Turrialba"
)

step73_file <- file.path(
  root,
  "51_bacteriocin_final_evidence",
  "final_independent_locus_evidence_18x6.tsv"
)

step74_file <- file.path(
  root,
  "52_bacteriocin_longitudinal_analysis",
  "longitudinal_trajectories.tsv"
)

master_file <- file.path(
  root,
  "61_MAG_functional_master",
  "MAG_functional_master_18MAGs.tsv"
)

rules_file <- file.path(
  root,
  "62_sample_functional_projection",
  "functional_projection_rules.tsv"
)

projection_file <- file.path(
  root,
  "62_sample_functional_projection",
  "sample_functional_projection_long.tsv"
)

focus_MAG_file <- file.path(
  root,
  "62_sample_functional_projection",
  "bacteriocin_focus_MAG_sample_projection.tsv"
)

outdir <- file.path(
  root,
  "64_bacteriocin_functional_integration"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# HELPERS
# ============================================================

check_file <- function(path) {
  if (!file.exists(path)) {
    stop(
      paste(
        "No existe archivo requerido:",
        path
      )
    )
  }
}

for (f in c(
  step73_file,
  step74_file,
  master_file,
  rules_file,
  projection_file,
  focus_MAG_file
)) {
  check_file(f)
}

write_clean <- function(x, filename) {
  write_tsv(
    x,
    file.path(
      outdir,
      filename
    ),
    na = ""
  )
}

split_semicolon <- function(x) {

  if (
    is.na(x) ||
    x == ""
  ) {
    return(
      character(0)
    )
  }

  str_split(
    x,
    ";",
    simplify = FALSE
  )[[1]]
}

# ============================================================
# CARGA
# ============================================================

loci <- read_tsv(
  step73_file,
  show_col_types = FALSE
)

traj74 <- read_tsv(
  step74_file,
  show_col_types = FALSE
)

master <- read_tsv(
  master_file,
  show_col_types = FALSE
)

rules <- read_tsv(
  rules_file,
  show_col_types = FALSE
)

proj <- read_tsv(
  projection_file,
  show_col_types = FALSE
)

focus <- read_tsv(
  focus_MAG_file,
  show_col_types = FALSE
)

# ============================================================
# QC DE COLUMNAS
# ============================================================

required_loci <- c(
  "sample",
  "producer",
  "producer_name",
  "biological_unit",
  "week",
  "subject",
  "locus_id",
  "source_MAG",
  "detection_class",
  "source_MAG_detection_class",
  "source_MAG_main_detection",
  "source_MAG_strict_detection",
  "evidence_tier_main",
  "evidence_tier_strict",
  "primary_source_linked_detection",
  "strict_source_linked_detection",
  "sensitivity_locus_signal",
  "coverage",
  "breadth",
  "breadth_expected",
  "breadth_to_expected_ratio",
  "filtered_read_pair_count",
  "coverage_per_million_input_pairs",
  "filtered_pairs_CPM",
  "primary_coverage_per_million",
  "primary_filtered_pairs_CPM",
  "count_as_independent_locus"
)

missing_loci <- setdiff(
  required_loci,
  names(loci)
)

if (length(missing_loci) > 0) {
  stop(
    paste(
      "Faltan columnas Step73:",
      paste(
        missing_loci,
        collapse = ", "
      )
    )
  )
}

required_focus <- c(
  "sample",
  "producer",
  "biological_unit",
  "week",
  "MAG",
  "genus",
  "species",
  "bacteriocin_analysis_role",
  "independent_bacteriocin_loci_n",
  "independent_bacteriocin_loci",
  "detection_class",
  "main_detection",
  "strict_detection",
  "relative_abundance_main_pct",
  "relative_abundance_strict_pct"
)

missing_focus <- setdiff(
  required_focus,
  names(focus)
)

if (length(missing_focus) > 0) {
  stop(
    paste(
      "Faltan columnas Step80b MAG-focus:",
      paste(
        missing_focus,
        collapse = ", "
      )
    )
  )
}

required_proj <- c(
  "sample",
  "functional_domain",
  "feature",
  "primary_supported_main",
  "primary_carriers_main",
  "primary_carrier_RA_main_pct",
  "primary_supported_strict",
  "primary_carrier_RA_strict_pct"
)

missing_proj <- setdiff(
  required_proj,
  names(proj)
)

if (length(missing_proj) > 0) {
  stop(
    paste(
      "Faltan columnas Step80b projection:",
      paste(
        missing_proj,
        collapse = ", "
      )
    )
  )
}

required_rules <- c(
  "feature",
  "domain",
  "source_field",
  "primary_states",
  "sensitivity_states"
)

missing_rules <- setdiff(
  required_rules,
  names(rules)
)

if (length(missing_rules) > 0) {
  stop(
    paste(
      "Faltan columnas de reglas:",
      paste(
        missing_rules,
        collapse = ", "
      )
    )
  )
}

# ============================================================
# TIPOS
# ============================================================

loci <- loci %>%
  mutate(
    biological_unit =
      as.integer(
        biological_unit
      ),
    week =
      as.integer(
        week
      ),
    primary_source_linked_detection =
      as.integer(
        primary_source_linked_detection
      ),
    strict_source_linked_detection =
      as.integer(
        strict_source_linked_detection
      ),
    sensitivity_locus_signal =
      as.integer(
        sensitivity_locus_signal
      ),
    source_MAG_main_detection =
      as.integer(
        source_MAG_main_detection
      ),
    source_MAG_strict_detection =
      as.integer(
        source_MAG_strict_detection
      ),
    count_as_independent_locus =
      as.integer(
        count_as_independent_locus
      )
  )

focus <- focus %>%
  mutate(
    biological_unit =
      as.integer(
        biological_unit
      ),
    week =
      as.integer(
        week
      ),
    main_detection =
      as.integer(
        main_detection
      ),
    strict_detection =
      as.integer(
        strict_detection
      ),
    relative_abundance_main_pct =
      as.numeric(
        relative_abundance_main_pct
      ),
    relative_abundance_strict_pct =
      as.numeric(
        relative_abundance_strict_pct
      )
  )

proj <- proj %>%
  mutate(
    week =
      as.integer(
        week
      ),
    biological_unit =
      as.integer(
        biological_unit
      ),
    primary_supported_main =
      as.integer(
        primary_supported_main
      ),
    primary_supported_strict =
      as.integer(
        primary_supported_strict
      ),
    primary_carrier_RA_main_pct =
      as.numeric(
        primary_carrier_RA_main_pct
      ),
    primary_carrier_RA_strict_pct =
      as.numeric(
        primary_carrier_RA_strict_pct
      )
  )

# ============================================================
# 1. QC ESTRUCTURAL STEP73
# ============================================================

if (nrow(loci) != 108) {
  stop(
    paste(
      "Se esperaban 108 pares muestra×locus;",
      "se encontraron",
      nrow(loci)
    )
  )
}

if (n_distinct(loci$sample) != 18) {
  stop(
    "Step73 no contiene exactamente 18 muestras."
  )
}

if (n_distinct(loci$locus_id) != 6) {
  stop(
    "Step73 no contiene exactamente 6 loci independientes."
  )
}

if (any(
  loci$count_as_independent_locus != 1
)) {
  stop(
    "La tabla Step73 contiene loci que no están marcados como independientes."
  )
}

expected_loci <- c(
  "ATTRLOC001",
  "ATTRLOC003",
  "ATTRLOC004",
  "ATTRLOC005",
  "ATTRLOC006",
  "ATTRLOC007"
)

if (!setequal(
  unique(
    loci$locus_id
  ),
  expected_loci
)) {
  stop(
    "El conjunto de loci independientes no coincide con Step73 final."
  )
}

# ============================================================
# 2. QC DE CONTEOS FINALES STEP73
# ============================================================

primary_total <- sum(
  loci$primary_source_linked_detection
)

strict_total <- sum(
  loci$strict_source_linked_detection
)

sensitivity_total <- sum(
  loci$sensitivity_locus_signal
)

if (primary_total != 45) {
  stop(
    paste(
      "Se esperaban 45 detecciones primary source-linked;",
      "se obtuvieron",
      primary_total
    )
  )
}

if (strict_total != 45) {
  stop(
    paste(
      "Se esperaban 45 detecciones strict source-linked;",
      "se obtuvieron",
      strict_total
    )
  )
}

if (sensitivity_total != 57) {
  stop(
    paste(
      "Se esperaban 57 señales locus sensitivity;",
      "se obtuvieron",
      sensitivity_total
    )
  )
}

main_strict_disagreement <- loci %>%
  filter(
    primary_source_linked_detection !=
      strict_source_linked_detection
  )

if (nrow(
  main_strict_disagreement
) != 0) {
  stop(
    "Step73 ya no presenta concordancia perfecta primary vs strict."
  )
}

# ============================================================
# 3. QC STEP73 vs STEP74
# Reconstruimos las 36 trayectorias subject×locus
# ============================================================

traj_reconstructed <- loci %>%
  select(
    producer,
    producer_name,
    biological_unit,
    subject,
    locus_id,
    week,
    primary_source_linked_detection
  ) %>%
  arrange(
    producer,
    biological_unit,
    locus_id,
    week
  ) %>%
  group_by(
    producer,
    producer_name,
    biological_unit,
    subject,
    locus_id
  ) %>%
  summarise(
    week0 =
      primary_source_linked_detection[
        week == 0
      ][1],
    week1 =
      primary_source_linked_detection[
        week == 1
      ][1],
    week2 =
      primary_source_linked_detection[
        week == 2
      ][1],
    trajectory =
      paste0(
        week0,
        week1,
        week2
      ),
    .groups = "drop"
  )

qc_traj <- traj_reconstructed %>%
  left_join(
    traj74,
    by = c(
      "producer",
      "producer_name",
      "biological_unit",
      "subject",
      "locus_id"
    ),
    suffix = c(
      "_reconstructed",
      "_step74"
    )
  ) %>%
  mutate(
    trajectory_match =
      trajectory_reconstructed ==
      trajectory_step74
  )

if (
  nrow(qc_traj) != 36 ||
  any(is.na(
    qc_traj$trajectory_match
  )) ||
  any(
    !qc_traj$trajectory_match
  )
) {
  stop(
    "QC Step73 vs Step74 falló."
  )
}

write_clean(
  qc_traj,
  "QC_step73_vs_step74_trajectories.tsv"
)

# ============================================================
# 4. MAGs FUENTE INDEPENDIENTES
# ============================================================

source_map <- loci %>%
  distinct(
    locus_id,
    source_MAG
  ) %>%
  arrange(
    source_MAG,
    locus_id
  )

if (n_distinct(
  source_map$source_MAG
) != 3) {
  stop(
    "Se esperaban exactamente 3 MAGs fuente independientes."
  )
}

source_summary <- source_map %>%
  group_by(
    source_MAG
  ) %>%
  summarise(
    n_independent_loci =
      n(),
    independent_loci =
      paste(
        locus_id,
        collapse = ";"
      ),
    .groups = "drop"
  )

# ============================================================
# 5. ABUNDANCIA DEL MAG FUENTE
# Excluir explícitamente Atlantibacter/ATTRLOC002 de esta
# integración independiente.
# ============================================================

focus_sources <- focus %>%
  filter(
    MAG %in%
      source_summary$source_MAG
  ) %>%
  rename(
    source_MAG = MAG,
    source_MAG_genus = genus,
    source_MAG_species = species,
    source_MAG_detection_class_focus =
      detection_class,
    source_MAG_main_detection_focus =
      main_detection,
    source_MAG_strict_detection_focus =
      strict_detection,
    source_MAG_RA_main_pct =
      relative_abundance_main_pct,
    source_MAG_RA_strict_pct =
      relative_abundance_strict_pct
  )

if (nrow(focus_sources) != 54) {
  stop(
    paste(
      "Se esperaban 54 filas sample×source_MAG;",
      "se encontraron",
      nrow(focus_sources)
    )
  )
}

# ============================================================
# 6. TABLA INTEGRADA 108 PARES LOCUS × MUESTRA
# ============================================================

locus_integrated <- loci %>%
  left_join(
    focus_sources %>%
      select(
        sample,
        source_MAG,
        source_MAG_genus,
        source_MAG_species,
        source_MAG_RA_main_pct,
        source_MAG_RA_strict_pct,
        source_MAG_main_detection_focus,
        source_MAG_strict_detection_focus,
        source_MAG_detection_class_focus
      ),
    by = c(
      "sample",
      "source_MAG"
    )
  ) %>%
  mutate(
    source_MAG_detection_match =
      source_MAG_main_detection ==
      source_MAG_main_detection_focus
  )

if (any(
  is.na(
    locus_integrated$source_MAG_RA_main_pct
  )
)) {
  stop(
    "Algún locus no pudo enlazarse con la abundancia de su MAG fuente."
  )
}

if (any(
  !locus_integrated$source_MAG_detection_match
)) {
  stop(
    "La detección del MAG fuente no coincide entre Step73 y Step80b."
  )
}

write_clean(
  locus_integrated,
  "locus_MAG_integrated_evidence_108rows.tsv"
)

# ============================================================
# 7. RESUMEN SAMPLE × SOURCE MAG
#
# 18 muestras × 3 MAGs = 54 filas
# ============================================================

source_sample <- locus_integrated %>%
  group_by(
    sample,
    producer,
    producer_name,
    biological_unit,
    week,
    subject,
    source_MAG,
    source_MAG_genus,
    source_MAG_species,
    source_MAG_RA_main_pct,
    source_MAG_RA_strict_pct,
    source_MAG_main_detection_focus,
    source_MAG_strict_detection_focus
  ) %>%
  summarise(
    n_independent_loci =
      n(),
    independent_loci =
      paste(
        locus_id,
        collapse = ";"
      ),
    n_primary_source_linked =
      sum(
        primary_source_linked_detection
      ),
    primary_detected_loci =
      paste(
        locus_id[
          primary_source_linked_detection == 1
        ],
        collapse = ";"
      ),
    n_strict_source_linked =
      sum(
        strict_source_linked_detection
      ),
    n_sensitivity_locus_signals =
      sum(
        sensitivity_locus_signal
      ),
    sensitivity_loci =
      paste(
        locus_id[
          sensitivity_locus_signal == 1
        ],
        collapse = ";"
      ),
    fraction_independent_loci_primary =
      n_primary_source_linked /
      n_independent_loci,
    fraction_independent_loci_sensitivity =
      n_sensitivity_locus_signals /
      n_independent_loci,
    all_independent_loci_primary =
      as.integer(
        n_primary_source_linked ==
          n_independent_loci
      ),
    any_independent_locus_primary =
      as.integer(
        n_primary_source_linked > 0
      ),
    .groups = "drop"
  )

if (nrow(source_sample) != 54) {
  stop(
    "La tabla source_MAG×sample no contiene 54 filas."
  )
}

write_clean(
  source_sample,
  "source_MAG_locus_trajectories_54rows.tsv"
)

# ============================================================
# 8. RESUMEN DE CONCORDANCIA POR LOCUS
# ============================================================

locus_summary <- locus_integrated %>%
  group_by(
    locus_id,
    source_MAG,
    source_MAG_genus,
    source_MAG_species
  ) %>%
  summarise(
    n_samples =
      n(),
    primary_source_linked =
      sum(
        primary_source_linked_detection
      ),
    strict_source_linked =
      sum(
        strict_source_linked_detection
      ),
    sensitivity_locus_signal =
      sum(
        sensitivity_locus_signal
      ),
    A_source_linked =
      sum(
        evidence_tier_main ==
          "A_source_linked"
      ),
    B_locus_signal_source_unresolved =
      sum(
        evidence_tier_main ==
          "B_locus_signal_source_unresolved"
      ),
    C_source_MAG_without_locus_support =
      sum(
        evidence_tier_main ==
          "C_source_MAG_without_locus_support"
      ),
    D_no_support =
      sum(
        evidence_tier_main ==
          "D_no_support"
      ),
    median_source_MAG_RA_when_primary =
      ifelse(
        any(
          primary_source_linked_detection == 1
        ),
        median(
          source_MAG_RA_main_pct[
            primary_source_linked_detection == 1
          ]
        ),
        NA_real_
      ),
    .groups = "drop"
  )

write_clean(
  locus_summary,
  "integrated_summary_by_locus.tsv"
)

# ============================================================
# 9. RECONSTRUIR FUNCIONES PRIMARY QUE CODIFICA CADA
# MAG FUENTE A PARTIR DE LAS REGLAS DE STEP80b
# ============================================================

source_function_rows <- list()
row_counter <- 1L

source_MAGs <- source_summary$source_MAG

for (mag in source_MAGs) {

  mrow <- master %>%
    filter(
      MAG == mag
    )

  if (nrow(mrow) != 1) {
    stop(
      paste(
        "No se encontró exactamente una fila del MAG:",
        mag
      )
    )
  }

  for (i in seq_len(
    nrow(rules)
  )) {

    source_field <-
      rules$source_field[i]

    state <-
      as.character(
        mrow[[source_field]][1]
      )

    primary_states <-
      split_semicolon(
        rules$primary_states[i]
      )

    sensitivity_states <-
      split_semicolon(
        rules$sensitivity_states[i]
      )

    source_function_rows[[row_counter]] <-
      tibble(
        source_MAG =
          mag,
        genus =
          mrow$genus[1],
        species =
          mrow$species[1],
        functional_domain =
          rules$domain[i],
        feature =
          rules$feature[i],
        source_field =
          source_field,
        source_MAG_functional_state =
          state,
        source_MAG_primary_carrier =
          as.integer(
            state %in%
              primary_states
          ),
        source_MAG_sensitivity_carrier =
          as.integer(
            state %in%
              sensitivity_states
          )
      )

    row_counter <-
      row_counter + 1L
  }
}

source_functions <- bind_rows(
  source_function_rows
)

source_functions_primary <- source_functions %>%
  filter(
    source_MAG_primary_carrier == 1
  )

write_clean(
  source_functions,
  "source_MAG_functional_rule_evaluation.tsv"
)

write_clean(
  source_functions_primary,
  "source_MAG_primary_functional_features.tsv"
)

# ============================================================
# 10. PERFIL FUNCIONAL COMPACTO DE CADA MAG FUENTE
# ============================================================

source_function_profile <- source_functions_primary %>%
  group_by(
    source_MAG,
    genus,
    species
  ) %>%
  summarise(
    n_primary_functional_features =
      n(),
    primary_functional_features =
      paste(
        feature,
        collapse = ";"
      ),
    .groups = "drop"
  ) %>%
  left_join(
    source_summary,
    by = "source_MAG"
  )

write_clean(
  source_function_profile,
  "source_MAG_primary_functional_profile_summary.tsv"
)

# ============================================================
# 11. CONTRIBUCIÓN DEL MAG FUENTE A CADA FUNCIÓN
#
# Solo para funciones que el MAG fuente porta como PRIMARY.
# ============================================================

source_function_contribution <- source_functions_primary %>%
  select(
    source_MAG,
    genus,
    species,
    functional_domain,
    feature,
    source_MAG_functional_state
  ) %>%
  crossing(
    sample =
      unique(
        proj$sample
      )
  ) %>%
  left_join(
    focus_sources %>%
      select(
        sample,
        source_MAG,
        producer,
        biological_unit,
        week,
        source_MAG_RA_main_pct,
        source_MAG_RA_strict_pct,
        source_MAG_main_detection_focus,
        source_MAG_strict_detection_focus
      ),
    by = c(
      "sample",
      "source_MAG"
    )
  ) %>%
  left_join(
    proj %>%
      select(
        sample,
        feature,
        primary_supported_main,
        primary_supported_strict,
        primary_carrier_RA_main_pct,
        primary_carrier_RA_strict_pct
      ),
    by = c(
      "sample",
      "feature"
    )
  ) %>%
  mutate(
    source_MAG_functional_contribution_RA_main_pct =
      if_else(
        source_MAG_main_detection_focus == 1,
        source_MAG_RA_main_pct,
        0
      ),
    source_MAG_functional_contribution_RA_strict_pct =
      if_else(
        source_MAG_strict_detection_focus == 1,
        source_MAG_RA_strict_pct,
        0
      ),
    source_MAG_share_of_functional_RA_main_pct =
      if_else(
        primary_carrier_RA_main_pct > 0,
        100 *
          source_MAG_functional_contribution_RA_main_pct /
          primary_carrier_RA_main_pct,
        NA_real_
      ),
    source_MAG_share_of_functional_RA_strict_pct =
      if_else(
        primary_carrier_RA_strict_pct > 0,
        100 *
          source_MAG_functional_contribution_RA_strict_pct /
          primary_carrier_RA_strict_pct,
        NA_real_
      )
  ) %>%
  left_join(
    source_sample %>%
      select(
        sample,
        source_MAG,
        n_independent_loci,
        n_primary_source_linked,
        n_sensitivity_locus_signals,
        fraction_independent_loci_primary,
        any_independent_locus_primary,
        all_independent_loci_primary
      ),
    by = c(
      "sample",
      "source_MAG"
    )
  )

# ============================================================
# QC:
# un MAG que es carrier primary no puede aportar más RA que
# toda la señal funcional comunitaria, salvo redondeo.
# ============================================================

bad_function_share <- source_function_contribution %>%
  filter(
    !is.na(
      source_MAG_share_of_functional_RA_main_pct
    ),
    source_MAG_share_of_functional_RA_main_pct >
      100.0001
  )

if (nrow(
  bad_function_share
) > 0) {
  stop(
    "Un MAG fuente aporta >100% de una señal funcional. Revisar integración."
  )
}

write_clean(
  source_function_contribution,
  "source_MAG_functional_contribution_by_sample.tsv"
)

# ============================================================
# 12. RESUMEN DE CUÁNTO DOMINA CADA MAG FUENTE SUS FUNCIONES
# ============================================================

source_function_dominance <- source_function_contribution %>%
  filter(
    source_MAG_main_detection_focus == 1,
    primary_carrier_RA_main_pct > 0
  ) %>%
  group_by(
    source_MAG,
    genus,
    species,
    functional_domain,
    feature
  ) %>%
  summarise(
    n_samples_source_MAG_detected =
      n(),
    median_source_MAG_RA_main_pct =
      median(
        source_MAG_RA_main_pct
      ),
    median_community_functional_RA_main_pct =
      median(
        primary_carrier_RA_main_pct
      ),
    median_source_MAG_share_of_functional_RA_pct =
      median(
        source_MAG_share_of_functional_RA_main_pct,
        na.rm = TRUE
      ),
    min_source_MAG_share_of_functional_RA_pct =
      min(
        source_MAG_share_of_functional_RA_main_pct,
        na.rm = TRUE
      ),
    max_source_MAG_share_of_functional_RA_pct =
      max(
        source_MAG_share_of_functional_RA_main_pct,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

write_clean(
  source_function_dominance,
  "source_MAG_functional_dominance_summary.tsv"
)

# ============================================================
# 13. RESUMEN LONGITUDINAL POR PRODUCTOR × SEMANA × MAG FUENTE
#
# Cada celda contiene 3 unidades biológicas.
# ============================================================

source_producer_week <- source_sample %>%
  group_by(
    producer,
    producer_name,
    week,
    source_MAG,
    source_MAG_genus,
    source_MAG_species,
    n_independent_loci
  ) %>%
  summarise(
    n_biological_units =
      n_distinct(
        biological_unit
      ),
    n_units_source_MAG_main_detected =
      sum(
        source_MAG_main_detection_focus
      ),
    prevalence_source_MAG_main =
      n_units_source_MAG_main_detected /
      n_biological_units,
    median_source_MAG_RA_main_pct =
      median(
        source_MAG_RA_main_pct
      ),
    min_source_MAG_RA_main_pct =
      min(
        source_MAG_RA_main_pct
      ),
    max_source_MAG_RA_main_pct =
      max(
        source_MAG_RA_main_pct
      ),
    total_primary_locus_detections =
      sum(
        n_primary_source_linked
      ),
    possible_locus_detections =
      n_biological_units *
      first(
        n_independent_loci
      ),
    fraction_possible_locus_detections =
      total_primary_locus_detections /
      possible_locus_detections,
    n_units_any_locus_primary =
      sum(
        any_independent_locus_primary
      ),
    n_units_all_loci_primary =
      sum(
        all_independent_loci_primary
      ),
    .groups = "drop"
  )

if (any(
  source_producer_week$n_biological_units != 3
)) {
  stop(
    "Alguna combinación productor×semana×MAG no contiene 3 unidades biológicas."
  )
}

write_clean(
  source_producer_week,
  "source_MAG_locus_summary_by_producer_week.tsv"
)

# ============================================================
# 14. MATRIZ SAMPLE × SOURCE MAG
# ABUNDANCIA + FRACCIÓN DE LOCI DIRECTAMENTE SOPORTADOS
# ============================================================

source_wide <- source_sample %>%
  select(
    sample,
    producer,
    biological_unit,
    week,
    source_MAG,
    source_MAG_RA_main_pct,
    fraction_independent_loci_primary
  ) %>%
  pivot_wider(
    names_from =
      source_MAG,
    values_from =
      c(
        source_MAG_RA_main_pct,
        fraction_independent_loci_primary
      ),
    names_sep =
      "__"
  )

write_clean(
  source_wide,
  "source_MAG_abundance_and_locus_fraction_matrix.tsv"
)

# ============================================================
# 15. FIGURA A
# Trayectorias de abundancia del MAG fuente.
# El tamaño del punto representa la fracción de sus loci
# con evidencia directa source-linked.
# ============================================================

plot_source <- source_sample %>%
  mutate(
    source_label =
      case_when(
        source_MAG ==
          "L2__L2_maxbin2.004_sub" ~
          "Lactococcus petauri | ATTRLOC001",
        source_MAG ==
          "L3__concoct_29" ~
          "Lactococcus_A laudensis | ATTRLOC003-004",
        source_MAG ==
          "M2__M2_maxbin2.004" ~
          "Lactococcus lactis | ATTRLOC005-007",
        TRUE ~
          source_MAG
      ),
    subject =
      factor(
        subject
      )
  )

pA <- ggplot(
  plot_source,
  aes(
    x = week,
    y =
      source_MAG_RA_main_pct,
    group =
      subject,
    color =
      producer
  )
) +
  geom_line(
    linewidth = 0.65,
    alpha = 0.75
  ) +
  geom_point(
    aes(
      size =
        fraction_independent_loci_primary
    ),
    alpha = 0.9
  ) +
  facet_wrap(
    ~ source_label,
    ncol = 1
  ) +
  scale_x_continuous(
    breaks = c(
      0,
      1,
      2
    )
  ) +
  scale_size_continuous(
    limits = c(
      0,
      1
    ),
    breaks = c(
      0,
      0.5,
      1
    ),
    name =
      "Fraction of source\nloci directly supported"
  ) +
  labs(
    x =
      "Week",
    y =
      "Source MAG relative abundance (%)",
    color =
      "Producer",
    title =
      "Source MAG trajectories and direct candidate-locus support",
    subtitle =
      "MAG relative abundance is within the recovered accepted MAG set; point size represents direct source-linked locus evidence"
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
    legend.position =
      "top"
  )

ggsave(
  file.path(
    outdir,
    "Figure81b_A_source_MAG_locus_trajectories.pdf"
  ),
  pA,
  width = 10,
  height = 11
)

ggsave(
  file.path(
    outdir,
    "Figure81b_A_source_MAG_locus_trajectories.png"
  ),
  pA,
  width = 10,
  height = 11,
  dpi = 300
)

# ============================================================
# 16. FIGURA B
# Qué proporción de la señal funcional comunitaria procede del
# MAG fuente de bacteriocina.
# ============================================================

feature_order <- rules$feature

pB_data <- source_function_contribution %>%
  filter(
    source_MAG_main_detection_focus == 1,
    !is.na(
      source_MAG_share_of_functional_RA_main_pct
    )
  ) %>%
  mutate(
    feature =
      factor(
        feature,
        levels =
          rev(
            feature_order
          )
      ),
    sample =
      factor(
        sample,
        levels =
          unique(
            proj %>%
              arrange(
                producer,
                biological_unit,
                week
              ) %>%
              pull(
                sample
              )
          )
      ),
    source_label =
      case_when(
        source_MAG ==
          "L2__L2_maxbin2.004_sub" ~
          "L. petauri",
        source_MAG ==
          "L3__concoct_29" ~
          "L. laudensis",
        source_MAG ==
          "M2__M2_maxbin2.004" ~
          "L. lactis",
        TRUE ~
          source_MAG
      )
  )

pB <- ggplot(
  pB_data,
  aes(
    x = sample,
    y = feature,
    fill =
      source_MAG_share_of_functional_RA_main_pct
  )
) +
  geom_tile() +
  facet_wrap(
    ~ source_label,
    scales = "free_y",
    ncol = 1
  ) +
  scale_fill_viridis_c(
    limits = c(
      0,
      100
    ),
    name =
      "Source MAG share\nof functional RA (%)"
  ) +
  labs(
    x = "Sample",
    y = NULL,
    title =
      "Contribution of bacteriocin-source MAGs to projected functional signals",
    subtitle =
      "Only primary functions encoded by each source MAG are shown"
  ) +
  theme_bw(
    base_size = 8
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),
    panel.grid =
      element_blank()
  )

ggsave(
  file.path(
    outdir,
    "Figure81b_B_source_MAG_functional_contribution.pdf"
  ),
  pB,
  width = 12,
  height = 14
)

ggsave(
  file.path(
    outdir,
    "Figure81b_B_source_MAG_functional_contribution.png"
  ),
  pB,
  width = 12,
  height = 14,
  dpi = 300
)

# ============================================================
# 17. FIGURA C
# Evidencia locus directa frente a abundancia del MAG fuente.
# DESCRIPTIVA: no se calcula p-value.
# ============================================================

pC <- ggplot(
  source_sample,
  aes(
    x =
      source_MAG_RA_main_pct,
    y =
      fraction_independent_loci_primary,
    shape =
      producer
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ source_MAG_genus,
    ncol = 1
  ) +
  scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = c(
      0,
      0.5,
      1
    )
  ) +
  labs(
    x =
      "Source MAG relative abundance (%)",
    y =
      "Fraction of independent loci with direct source-linked support",
    shape =
      "Producer",
    title =
      "Direct candidate-locus support versus source-MAG representation",
    subtitle =
      "Descriptive relationship only; repeated observations are not independent"
  ) +
  theme_bw(
    base_size = 10
  )

ggsave(
  file.path(
    outdir,
    "Figure81b_C_locus_support_vs_source_MAG_RA.pdf"
  ),
  pC,
  width = 9,
  height = 10
)

ggsave(
  file.path(
    outdir,
    "Figure81b_C_locus_support_vs_source_MAG_RA.png"
  ),
  pC,
  width = 9,
  height = 10,
  dpi = 300
)

# ============================================================
# 18. README
# ============================================================

writeLines(
  c(
    "PASO 81b - INTEGRACION BACTERIOCINA / MAG / FUNCION",
    "====================================================",
    "",
    "Este análisis mantiene separados tres niveles:",
    "",
    "1. Evidencia directa del locus candidato bacteriocina/RiPP.",
    "2. Detección y abundancia del MAG fuente.",
    "3. Potencial funcional codificado por el MAG.",
    "",
    "La detección del MAG NO sustituye la detección del locus.",
    "",
    "La detección directa del locus NO demuestra expresión,",
    "síntesis de péptido antimicrobiano ni actividad inhibitoria.",
    "",
    "La proyección funcional representa potencial genómico de",
    "MAGs recuperados y aceptados; no expresión ni metabolitos.",
    "",
    "ATTRLOC002 / Atlantibacter se excluye de los seis loci",
    "biológicos independientes y no entra en esta integración",
    "como fuente independiente.",
    "",
    "ARCHIVOS PRINCIPALES",
    "--------------------",
    "",
    "locus_MAG_integrated_evidence_108rows.tsv",
    "  Una fila por muestra x locus independiente.",
    "",
    "source_MAG_locus_trajectories_54rows.tsv",
    "  18 muestras x 3 MAGs fuente.",
    "",
    "source_MAG_primary_functional_features.tsv",
    "  Capacidades primary codificadas por cada MAG fuente.",
    "",
    "source_MAG_functional_contribution_by_sample.tsv",
    "  Contribución del MAG fuente a la señal funcional de la",
    "  comunidad recuperada en cada muestra.",
    "",
    "source_MAG_functional_dominance_summary.tsv",
    "  Resume cuánto de cada señal funcional procede del MAG",
    "  fuente cuando éste está detectado.",
    "",
    "source_MAG_locus_summary_by_producer_week.tsv",
    "  Resumen descriptivo longitudinal, n=3 unidades por",
    "  productor y semana.",
    "",
    "QC_step73_vs_step74_trajectories.tsv",
    "  Verifica identidad entre Step73 y Step74.",
    "",
    "No se realizan pruebas inferenciales en este paso."
  ),
  con =
    file.path(
      outdir,
      "README_step81b.txt"
    )
)

# ============================================================
# RESUMEN FINAL
# ============================================================

cat(
  strrep(
    "=",
    60
  ),
  "\n",
  sep = ""
)

cat(
  "PASO 81b COMPLETADO\n"
)

cat(
  strrep(
    "=",
    60
  ),
  "\n",
  sep = ""
)

cat(
  "Muestras:                         ",
  n_distinct(
    loci$sample
  ),
  "\n"
)

cat(
  "Loci independientes:              ",
  n_distinct(
    loci$locus_id
  ),
  "\n"
)

cat(
  "MAGs fuente independientes:       ",
  n_distinct(
    loci$source_MAG
  ),
  "\n"
)

cat(
  "Pares muestra-locus:              ",
  nrow(
    locus_integrated
  ),
  "\n"
)

cat(
  "Detecciones primary source-linked:",
  primary_total,
  "\n"
)

cat(
  "Detecciones strict source-linked: ",
  strict_total,
  "\n"
)

cat(
  "Señales sensitivity:              ",
  sensitivity_total,
  "\n"
)

cat(
  "Trayectorias Step73/74 idénticas: ",
  sum(
    qc_traj$trajectory_match
  ),
  "/",
  nrow(
    qc_traj
  ),
  "\n"
)

cat(
  "MAG×muestra integrados:           ",
  nrow(
    source_sample
  ),
  "\n"
)

cat(
  "\nSalida: ",
  outdir,
  "\n",
  sep = ""
)

cat(
  "PASO 81b FINALIZÓ CORRECTAMENTE; ES SEGURO SALIR.\n"
)
