#!/usr/bin/env Rscript

# ============================================================
# PASO 81a
# Dinámica longitudinal del potencial funcional recuperado
# y descomposición por MAG conductor.
#
# Diseño:
#   2 productores
#   3 unidades biológicas por productor
#   3 semanas repetidas: 0, 1, 2
#
# IMPORTANTE:
# - Las 18 muestras NO son 18 réplicas biológicas independientes.
# - La RA funcional es la SUMA de la abundancia relativa de MAGs
#   portadores dentro del conjunto de MAGs recuperados/aceptados.
# - No es abundancia de genes, expresión ni producción metabólica.
# - No se realizan pruebas inferenciales en este paso.
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

proj_file <- file.path(
  root,
  "62_sample_functional_projection",
  "sample_functional_projection_long.tsv"
)

abund_file <- file.path(
  root,
  "43_MAG_abundance_taxonomy",
  "MAG_abundance_taxonomy_long.tsv"
)

master_file <- file.path(
  root,
  "61_MAG_functional_master",
  "MAG_functional_master_18MAGs.tsv"
)

outdir <- file.path(
  root,
  "63_functional_longitudinal_descriptive"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# HELPERS
# ============================================================

stop_missing <- function(path) {
  if (!file.exists(path)) {
    stop(
      paste(
        "No existe archivo requerido:",
        path
      )
    )
  }
}

stop_missing(proj_file)
stop_missing(abund_file)
stop_missing(master_file)

write_tsv_clean <- function(x, filename) {
  write_tsv(
    x,
    file.path(
      outdir,
      filename
    ),
    na = ""
  )
}

# ============================================================
# CARGA
# ============================================================

proj <- read_tsv(
  proj_file,
  show_col_types = FALSE
)

abund <- read_tsv(
  abund_file,
  show_col_types = FALSE
)

master <- read_tsv(
  master_file,
  show_col_types = FALSE
)

# ============================================================
# QC BÁSICO
# ============================================================

required_proj <- c(
  "sample",
  "producer",
  "biological_unit",
  "week",
  "functional_domain",
  "feature",
  "primary_supported_main",
  "primary_carriers_main",
  "primary_carrier_RA_main_pct",
  "primary_supported_strict",
  "primary_carriers_strict",
  "primary_carrier_RA_strict_pct",
  "sensitivity_supported_main",
  "sensitivity_carrier_RA_main_pct"
)

missing_proj <- setdiff(
  required_proj,
  names(proj)
)

if (length(missing_proj) > 0) {
  stop(
    paste(
      "Faltan columnas en proyección:",
      paste(
        missing_proj,
        collapse = ", "
      )
    )
  )
}

required_abund <- c(
  "sample",
  "producer",
  "biological_unit",
  "week",
  "MAG",
  "detection_class",
  "main_detection",
  "strict_detection",
  "relative_abundance_main_pct",
  "relative_abundance_strict_pct"
)

missing_abund <- setdiff(
  required_abund,
  names(abund)
)

if (length(missing_abund) > 0) {
  stop(
    paste(
      "Faltan columnas en abundancia:",
      paste(
        missing_abund,
        collapse = ", "
      )
    )
  )
}

required_master <- c(
  "MAG",
  "genus",
  "species",
  "bacteriocin_analysis_role"
)

missing_master <- setdiff(
  required_master,
  names(master)
)

if (length(missing_master) > 0) {
  stop(
    paste(
      "Faltan columnas en master:",
      paste(
        missing_master,
        collapse = ", "
      )
    )
  )
}

if (nrow(proj) != 18 * 23) {
  stop(
    paste(
      "Se esperaban 414 filas muestra-función;",
      "se encontraron",
      nrow(proj)
    )
  )
}

if (n_distinct(proj$sample) != 18) {
  stop("La proyección no contiene exactamente 18 muestras.")
}

if (n_distinct(proj$feature) != 23) {
  stop("La proyección no contiene exactamente 23 funciones.")
}

# ============================================================
# TIPOS Y METADATA LONGITUDINAL
# ============================================================

proj <- proj %>%
  mutate(
    biological_unit =
      as.integer(
        biological_unit
      ),
    week =
      as.integer(
        week
      ),
    primary_supported_main =
      as.integer(
        primary_supported_main
      ),
    primary_supported_strict =
      as.integer(
        primary_supported_strict
      ),
    sensitivity_supported_main =
      as.integer(
        sensitivity_supported_main
      ),
    primary_carrier_RA_main_pct =
      as.numeric(
        primary_carrier_RA_main_pct
      ),
    primary_carrier_RA_strict_pct =
      as.numeric(
        primary_carrier_RA_strict_pct
      ),
    sensitivity_carrier_RA_main_pct =
      as.numeric(
        sensitivity_carrier_RA_main_pct
      ),
    subject_id =
      paste0(
        producer,
        biological_unit
      )
  )

abund <- abund %>%
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

# ============================================================
# ETIQUETAS PARA FIGURAS
#
# Se corrige SOLO la etiqueta visual del término de
# transaminación para no sobreinterpretarlo como aroma.
# El feature ID original se conserva en todos los TSV.
# ============================================================

feature_labels <- c(
  galactose_Leloir_core =
    "Leloir galactose core",
  lactose_transport_betaGal =
    "Lactose transport + beta-gal",
  lactose_PTS_LacEFG =
    "Lactose PTS LacEFG",
  lactate_LDH =
    "Lactate dehydrogenase",
  acetate_PtaAckA =
    "Acetate Pta-AckA",
  formate_PFL =
    "Pyruvate formate-lyase",
  citrate_fermentation =
    "Citrate fermentation",
  acetoin_branch =
    "Acetoin branch",
  butanediol_branch =
    "2,3-butanediol branch",
  histamine_potential =
    "Histamine potential",
  tyramine_potential =
    "Tyramine potential",
  cadaverine_potential =
    "Cadaverine potential",
  putrescine_potential =
    "Putrescine potential",
  peptide_utilization_system =
    "Peptide utilization",
  CEP_like_surface_proteinase =
    "CEP-like surface proteinase",
  lipolytic_enzyme_candidate =
    "Lipolytic enzyme candidate",
  amino_acid_transamination_aroma_potential =
    "Amino-acid transamination potential",
  sulfur_aroma_candidate =
    "Sulfur-aroma candidate",
  EPS_like_biosynthesis_context =
    "EPS-like biosynthesis context",
  capsule_like_context =
    "Capsule-like context",
  acid_stress_system =
    "Acid-stress system",
  osmoadaptation_multigene_system =
    "Multigene osmoadaptation",
  oxidative_stress_repertoire =
    "Oxidative-stress repertoire"
)

proj <- proj %>%
  mutate(
    feature_label =
      unname(
        feature_labels[
          feature
        ]
      ),
    feature_label =
      if_else(
        is.na(feature_label),
        feature,
        feature_label
      )
  )

# Mantener el orden original de Step 80b
feature_order <- proj %>%
  distinct(feature) %>%
  pull(feature)

label_order <- unname(
  feature_labels[
    feature_order
  ]
)

# ============================================================
# 1. ESTABILIDAD MAIN vs STRICT
# ============================================================

stability <- proj %>%
  group_by(
    functional_domain,
    feature,
    feature_label
  ) %>%
  summarise(
    n_samples_primary_main =
      sum(
        primary_supported_main
      ),
    n_samples_primary_strict =
      sum(
        primary_supported_strict
      ),
    n_presence_discordant =
      sum(
        primary_supported_main !=
          primary_supported_strict
      ),
    discordant_samples =
      paste(
        sample[
          primary_supported_main !=
            primary_supported_strict
        ],
        collapse = ";"
      ),
    mean_RA_main =
      mean(
        primary_carrier_RA_main_pct
      ),
    mean_RA_strict =
      mean(
        primary_carrier_RA_strict_pct
      ),
    median_abs_RA_main_strict_difference =
      median(
        abs(
          primary_carrier_RA_main_pct -
            primary_carrier_RA_strict_pct
        )
      ),
    max_abs_RA_main_strict_difference =
      max(
        abs(
          primary_carrier_RA_main_pct -
            primary_carrier_RA_strict_pct
        )
      ),
    .groups = "drop"
  ) %>%
  mutate(
    detection_robustness =
      if_else(
        n_presence_discordant == 0,
        "presence_concordant_main_vs_strict",
        "detection_threshold_sensitive"
      )
  )

write_tsv_clean(
  stability,
  "functional_main_vs_strict_stability.tsv"
)

# ============================================================
# 2. RESUMEN POR PRODUCTOR × SEMANA
#
# n = 3 unidades biológicas por celda.
# Se priorizan mediana e IQR.
# ============================================================

producer_week <- proj %>%
  group_by(
    functional_domain,
    feature,
    feature_label,
    producer,
    week
  ) %>%
  summarise(
    n_biological_units =
      n_distinct(
        subject_id
      ),
    n_units_primary_present =
      sum(
        primary_supported_main
      ),
    percent_units_primary_present =
      100 *
      n_units_primary_present /
      n_biological_units,
    median_carrier_RA_main_pct =
      median(
        primary_carrier_RA_main_pct
      ),
    Q1_carrier_RA_main_pct =
      quantile(
        primary_carrier_RA_main_pct,
        0.25,
        names = FALSE
      ),
    Q3_carrier_RA_main_pct =
      quantile(
        primary_carrier_RA_main_pct,
        0.75,
        names = FALSE
      ),
    min_carrier_RA_main_pct =
      min(
        primary_carrier_RA_main_pct
      ),
    max_carrier_RA_main_pct =
      max(
        primary_carrier_RA_main_pct
      ),
    median_carrier_RA_strict_pct =
      median(
        primary_carrier_RA_strict_pct
      ),
    .groups = "drop"
  )

if (any(
  producer_week$n_biological_units != 3
)) {
  stop(
    "Alguna combinación productor×semana no contiene 3 unidades biológicas."
  )
}

write_tsv_clean(
  producer_week,
  "functional_summary_by_producer_week.tsv"
)

# ============================================================
# 3. TRAYECTORIAS POR UNIDAD BIOLÓGICA
# ============================================================

subject_trajectories <- proj %>%
  select(
    producer,
    biological_unit,
    subject_id,
    sample,
    week,
    functional_domain,
    feature,
    feature_label,
    primary_supported_main,
    primary_carrier_RA_main_pct,
    primary_supported_strict,
    primary_carrier_RA_strict_pct
  ) %>%
  arrange(
    producer,
    biological_unit,
    feature,
    week
  )

write_tsv_clean(
  subject_trajectories,
  "functional_subject_trajectories_long.tsv"
)

# ============================================================
# 4. CAMBIOS INTRA-UNIDAD
# ============================================================

subject_changes <- proj %>%
  select(
    producer,
    biological_unit,
    subject_id,
    functional_domain,
    feature,
    feature_label,
    week,
    primary_carrier_RA_main_pct,
    primary_supported_main
  ) %>%
  pivot_wider(
    names_from = week,
    values_from = c(
      primary_carrier_RA_main_pct,
      primary_supported_main
    ),
    names_sep = "_week"
  ) %>%
  mutate(
    delta_RA_week1_minus_week0 =
      primary_carrier_RA_main_pct_week1 -
      primary_carrier_RA_main_pct_week0,
    delta_RA_week2_minus_week0 =
      primary_carrier_RA_main_pct_week2 -
      primary_carrier_RA_main_pct_week0,
    presence_trajectory =
      paste0(
        primary_supported_main_week0,
        primary_supported_main_week1,
        primary_supported_main_week2
      )
  ) %>%
  arrange(
    feature,
    producer,
    biological_unit
  )

write_tsv_clean(
  subject_changes,
  "functional_subject_changes.tsv"
)

# ============================================================
# 5. RESUMEN DE CAMBIO POR PRODUCTOR
#
# Descriptivo únicamente.
# ============================================================

producer_change <- subject_changes %>%
  group_by(
    functional_domain,
    feature,
    feature_label,
    producer
  ) %>%
  summarise(
    n_biological_units =
      n(),
    median_delta_RA_week1_minus_week0 =
      median(
        delta_RA_week1_minus_week0
      ),
    min_delta_RA_week1_minus_week0 =
      min(
        delta_RA_week1_minus_week0
      ),
    max_delta_RA_week1_minus_week0 =
      max(
        delta_RA_week1_minus_week0
      ),
    median_delta_RA_week2_minus_week0 =
      median(
        delta_RA_week2_minus_week0
      ),
    min_delta_RA_week2_minus_week0 =
      min(
        delta_RA_week2_minus_week0
      ),
    max_delta_RA_week2_minus_week0 =
      max(
        delta_RA_week2_minus_week0
      ),
    trajectories =
      paste(
        presence_trajectory,
        collapse = ";"
      ),
    .groups = "drop"
  )

write_tsv_clean(
  producer_change,
  "functional_change_summary_by_producer.tsv"
)

# ============================================================
# 6. DESCOMPONER CADA FUNCIÓN EN SUS MAGs PORTADORES
#
# Esto permite distinguir:
# "cambió una función"
# de
# "cambió la abundancia de un MAG que porta muchas funciones".
# ============================================================

carrier_rows <- proj %>%
  select(
    sample,
    producer,
    biological_unit,
    week,
    subject_id,
    functional_domain,
    feature,
    feature_label,
    primary_carriers_main,
    projected_carrier_RA_main_pct =
      primary_carrier_RA_main_pct
  ) %>%
  filter(
    !is.na(
      primary_carriers_main
    ),
    primary_carriers_main != ""
  ) %>%
  separate_rows(
    primary_carriers_main,
    sep = ";"
  ) %>%
  rename(
    MAG =
      primary_carriers_main
  ) %>%
  left_join(
    abund %>%
      select(
        sample,
        MAG,
        detection_class,
        main_detection,
        strict_detection,
        MAG_RA_main_pct =
          relative_abundance_main_pct,
        MAG_RA_strict_pct =
          relative_abundance_strict_pct
      ),
    by = c(
      "sample",
      "MAG"
    )
  ) %>%
  left_join(
    master %>%
      select(
        MAG,
        genus,
        species,
        bacteriocin_analysis_role
      ),
    by = "MAG"
  )

if (any(
  is.na(
    carrier_rows$MAG_RA_main_pct
  )
)) {
  stop(
    "Hay MAGs portadores que no pudieron enlazarse con la tabla de abundancia."
  )
}

if (any(
  carrier_rows$main_detection != 1
)) {
  stop(
    "Apareció un MAG portador primary_main sin main_detection=1."
  )
}

# ============================================================
# 7. QC: la suma de MAGs portadores debe reproducir Step 80b
# ============================================================

carrier_qc <- carrier_rows %>%
  group_by(
    sample,
    feature
  ) %>%
  summarise(
    recalculated_carrier_RA_main_pct =
      sum(
        MAG_RA_main_pct
      ),
    projected_carrier_RA_main_pct =
      first(
        projected_carrier_RA_main_pct
      ),
    absolute_difference =
      abs(
        recalculated_carrier_RA_main_pct -
          projected_carrier_RA_main_pct
      ),
    .groups = "drop"
  )

max_carrier_qc_difference <- max(
  carrier_qc$absolute_difference
)

if (
  max_carrier_qc_difference >
    1e-5
) {
  stop(
    paste(
      "QC portadores falló. Diferencia máxima:",
      max_carrier_qc_difference
    )
  )
}

write_tsv_clean(
  carrier_qc,
  "QC_carrier_RA_reconstruction.tsv"
)

# ============================================================
# 8. MAG CONDUCTOR POR MUESTRA × FUNCIÓN
# ============================================================

top_carriers <- carrier_rows %>%
  group_by(
    sample,
    producer,
    biological_unit,
    week,
    subject_id,
    functional_domain,
    feature,
    feature_label
  ) %>%
  arrange(
    desc(
      MAG_RA_main_pct
    ),
    .by_group = TRUE
  ) %>%
  mutate(
    rank_within_feature =
      row_number(),
    total_functional_carrier_RA_main_pct =
      sum(
        MAG_RA_main_pct
      )
  ) %>%
  filter(
    rank_within_feature == 1
  ) %>%
  mutate(
    top_carrier_fraction_of_functional_RA_pct =
      if_else(
        total_functional_carrier_RA_main_pct > 0,
        100 *
          MAG_RA_main_pct /
          total_functional_carrier_RA_main_pct,
        NA_real_
      )
  ) %>%
  ungroup() %>%
  select(
    sample,
    producer,
    biological_unit,
    week,
    subject_id,
    functional_domain,
    feature,
    feature_label,
    top_carrier_MAG = MAG,
    genus,
    species,
    bacteriocin_analysis_role,
    top_carrier_MAG_RA_main_pct =
      MAG_RA_main_pct,
    total_functional_carrier_RA_main_pct,
    top_carrier_fraction_of_functional_RA_pct
  )

write_tsv_clean(
  top_carriers,
  "functional_top_carrier_by_sample.tsv"
)

# ============================================================
# 9. DOMINANCIA DEL MAG CONDUCTOR POR FUNCIÓN
# ============================================================

driver_summary <- top_carriers %>%
  group_by(
    functional_domain,
    feature,
    feature_label
  ) %>%
  summarise(
    n_samples_with_primary_support =
      n(),
    median_top_carrier_fraction_pct =
      median(
        top_carrier_fraction_of_functional_RA_pct
      ),
    min_top_carrier_fraction_pct =
      min(
        top_carrier_fraction_of_functional_RA_pct
      ),
    max_top_carrier_fraction_pct =
      max(
        top_carrier_fraction_of_functional_RA_pct
      ),
    n_distinct_top_carrier_MAGs =
      n_distinct(
        top_carrier_MAG
      ),
    top_carrier_MAGs =
      paste(
        sort(
          unique(
            top_carrier_MAG
          )
        ),
        collapse = ";"
      ),
    .groups = "drop"
  )

write_tsv_clean(
  driver_summary,
  "functional_driver_summary.tsv"
)

# ============================================================
# 10. TABLA DE VARIABILIDAD FUNCIONAL
# ============================================================

variability <- proj %>%
  group_by(
    functional_domain,
    feature,
    feature_label
  ) %>%
  summarise(
    n_samples_present_main =
      sum(
        primary_supported_main
      ),
    prevalence_main_pct =
      100 *
      mean(
        primary_supported_main
      ),
    mean_carrier_RA_main_pct =
      mean(
        primary_carrier_RA_main_pct
      ),
    median_carrier_RA_main_pct =
      median(
        primary_carrier_RA_main_pct
      ),
    SD_carrier_RA_main_pct =
      sd(
        primary_carrier_RA_main_pct
      ),
    IQR_carrier_RA_main_pct =
      IQR(
        primary_carrier_RA_main_pct
      ),
    min_carrier_RA_main_pct =
      min(
        primary_carrier_RA_main_pct
      ),
    max_carrier_RA_main_pct =
      max(
        primary_carrier_RA_main_pct
      ),
    RA_range_pct =
      max_carrier_RA_main_pct -
      min_carrier_RA_main_pct,
    .groups = "drop"
  ) %>%
  left_join(
    stability %>%
      select(
        feature,
        n_presence_discordant,
        detection_robustness
      ),
    by = "feature"
  ) %>%
  left_join(
    driver_summary %>%
      select(
        feature,
        median_top_carrier_fraction_pct,
        n_distinct_top_carrier_MAGs,
        top_carrier_MAGs
      ),
    by = "feature"
  )

write_tsv_clean(
  variability,
  "functional_variability_and_driver_summary.tsv"
)

# ============================================================
# FIGURA 1
# HEATMAP DE RA DE MAGs PORTADORES
# ============================================================

sample_levels <- proj %>%
  distinct(
    producer,
    biological_unit,
    week,
    sample
  ) %>%
  arrange(
    producer,
    biological_unit,
    week
  ) %>%
  pull(
    sample
  )

heat <- proj %>%
  mutate(
    sample =
      factor(
        sample,
        levels =
          sample_levels
      ),
    feature_label =
      factor(
        feature_label,
        levels =
          rev(
            label_order
          )
      )
  )

p_heat <- ggplot(
  heat,
  aes(
    x = sample,
    y = feature_label,
    fill =
      primary_carrier_RA_main_pct
  )
) +
  geom_tile() +
  scale_fill_viridis_c(
    limits = c(
      0,
      100
    ),
    name =
      "Carrier MAG\nRA (%)"
  ) +
  labs(
    x = "Sample",
    y = NULL,
    title =
      "Projected functional potential among recovered MAGs",
    subtitle =
      "Sum of relative abundance of accepted MAG carriers; not whole-microbiome functional abundance"
  ) +
  theme_bw(
    base_size = 10
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
    "Figure81a_1_functional_RA_heatmap_main.png"
  ),
  p_heat,
  width = 12,
  height = 10,
  dpi = 300
)

ggsave(
  file.path(
    outdir,
    "Figure81a_1_functional_RA_heatmap_main.pdf"
  ),
  p_heat,
  width = 12,
  height = 10
)

# ============================================================
# FIGURA 2
# TRAYECTORIAS DE LAS 6 UNIDADES BIOLÓGICAS
# ============================================================

traj_plot_data <- proj %>%
  mutate(
    feature_label =
      factor(
        feature_label,
        levels =
          label_order
      )
  )

p_traj <- ggplot(
  traj_plot_data,
  aes(
    x = week,
    y =
      primary_carrier_RA_main_pct,
    group =
      subject_id,
    color =
      producer
  )
) +
  geom_line(
    linewidth = 0.55,
    alpha = 0.8
  ) +
  geom_point(
    size = 1.4
  ) +
  facet_wrap(
    ~ feature_label,
    ncol = 4,
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
      "Carrier MAG relative abundance (%)",
    color =
      "Producer",
    title =
      "Longitudinal trajectories of projected functional potential",
    subtitle =
      "Each line is one biological unit; values are relative to the recovered accepted MAG set"
  ) +
  theme_bw(
    base_size = 9
  ) +
  theme(
    legend.position =
      "top",
    strip.text =
      element_text(
        size = 7
      )
  )

ggsave(
  file.path(
    outdir,
    "Figure81a_2_functional_longitudinal_trajectories.pdf"
  ),
  p_traj,
  width = 14,
  height = 18
)

ggsave(
  file.path(
    outdir,
    "Figure81a_2_functional_longitudinal_trajectories.png"
  ),
  p_traj,
  width = 14,
  height = 18,
  dpi = 300
)

# ============================================================
# FIGURA 3
# CAMBIO SEMANA 2 - SEMANA 0 POR UNIDAD BIOLÓGICA
# ============================================================

change_plot_data <- subject_changes %>%
  mutate(
    feature_label =
      factor(
        feature_label,
        levels =
          rev(
            label_order
          )
      )
  )

p_change <- ggplot(
  change_plot_data,
  aes(
    x =
      delta_RA_week2_minus_week0,
    y =
      feature_label,
    shape =
      producer
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +
  geom_point(
    size = 2,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ producer,
    ncol = 2
  ) +
  labs(
    x =
      "Change in carrier MAG RA: week 2 - week 0 (percentage points)",
    y = NULL,
    shape =
      "Producer",
    title =
      "Within-biological-unit functional change"
  ) +
  theme_bw(
    base_size = 9
  )

ggsave(
  file.path(
    outdir,
    "Figure81a_3_week2_minus_week0_changes.pdf"
  ),
  p_change,
  width = 12,
  height = 10
)

ggsave(
  file.path(
    outdir,
    "Figure81a_3_week2_minus_week0_changes.png"
  ),
  p_change,
  width = 12,
  height = 10,
  dpi = 300
)

# ============================================================
# README
# ============================================================

readme <- file.path(
  outdir,
  "README_step81a.txt"
)

writeLines(
  c(
    "PASO 81a - DINAMICA FUNCIONAL LONGITUDINAL",
    "===========================================",
    "",
    "Este paso es DESCRIPTIVO.",
    "",
    "La unidad biologica independiente es la combinación",
    "productor x biological_unit (6 unidades en total).",
    "Las semanas 0, 1 y 2 son medidas repetidas.",
    "",
    "La variable continua principal es:",
    "primary_carrier_RA_main_pct",
    "",
    "Esta variable es la suma de la abundancia relativa de",
    "los MAGs aceptados que portan la capacidad funcional.",
    "",
    "NO significa:",
    "- abundancia de genes",
    "- expresión",
    "- actividad enzimática",
    "- concentración de metabolitos",
    "- porcentaje de toda la microbiota",
    "",
    "La abundancia está normalizada dentro del conjunto de",
    "MAGs recuperados y aceptados.",
    "",
    "ARCHIVOS CLAVE",
    "--------------",
    "functional_summary_by_producer_week.tsv",
    "  Resumen descriptivo de las 3 unidades por productor/semana.",
    "",
    "functional_subject_changes.tsv",
    "  Cambios intra-unidad entre semanas.",
    "",
    "functional_top_carrier_by_sample.tsv",
    "  MAG que más contribuye a cada señal funcional.",
    "",
    "functional_driver_summary.tsv",
    "  Evalúa si una función está dominada por uno o varios MAGs.",
    "",
    "functional_main_vs_strict_stability.tsv",
    "  Sensibilidad de cada capacidad al criterio de detección.",
    "",
    "functional_variability_and_driver_summary.tsv",
    "  Tabla compacta para decidir qué capacidades son",
    "  longitudinalmente informativas.",
    "",
    "No se realizan pruebas estadísticas inferenciales en 81a.",
    "La inferencia productor/tiempo, si se realiza, debe respetar",
    "n=3 unidades biológicas por productor y medidas repetidas."
  ),
  con = readme
)

# ============================================================
# RESUMEN EN CONSOLA
# ============================================================

cat(
  paste0(
    strrep("=", 60),
    "\n"
  )
)

cat(
  "PASO 81a COMPLETADO\n"
)

cat(
  paste0(
    strrep("=", 60),
    "\n"
  )
)

cat(
  "Muestras:                         ",
  n_distinct(
    proj$sample
  ),
  "\n"
)

cat(
  "Unidades biológicas:              ",
  n_distinct(
    proj$subject_id
  ),
  "\n"
)

cat(
  "Funciones:                        ",
  n_distinct(
    proj$feature
  ),
  "\n"
)

cat(
  "Funciones main/strict concordantes:",
  sum(
    stability$n_presence_discordant == 0
  ),
  "/",
  nrow(
    stability
  ),
  "\n"
)

cat(
  "Funciones sensibles al criterio:  ",
  sum(
    stability$n_presence_discordant > 0
  ),
  "/",
  nrow(
    stability
  ),
  "\n"
)

cat(
  "Diferencia máxima QC portadores:  ",
  format(
    max_carrier_qc_difference,
    scientific = TRUE
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
  "PASO 81a FINALIZÓ CORRECTAMENTE; ES SEGURO SALIR.\n"
)
