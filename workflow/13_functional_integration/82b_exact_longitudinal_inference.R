#!/usr/bin/env Rscript

# ============================================================
# PASO 82b
# Inferencia longitudinal exploratoria exacta
#
# Foco:
# M2__M2_maxbin2.004 = Lactococcus lactis
# ATTRLOC005 / ATTRLOC006 / ATTRLOC007
#
# UNIDAD INDEPENDIENTE:
# 6 unidades biológicas:
# L1 L2 L3 M1 M2 M3
#
# Cada unidad tiene 3 medidas repetidas:
# semana 0, 1 y 2
#
# IMPORTANTE:
# - NO tratar las 18 muestras como independientes.
# - NO usar t-test sobre las 18 muestras.
# - NO usar chi-cuadrado/Fisher sobre las 18 muestras.
# - NO interpretar p >= 0.10 como "ausencia de efecto".
#
# Con 3 vs 3 existen C(6,3)=20 asignaciones.
# El mínimo p exacto bilateral posible es 0.10.
#
# Las abundancias son relativas al conjunto de MAGs
# recuperados/aceptados, no a toda la microbiota.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
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

INPUT <- file.path(
  ROOT,
  "64_bacteriocin_functional_integration",
  "source_MAG_locus_trajectories_54rows.tsv"
)

OUT <- file.path(
  ROOT,
  "66_exact_longitudinal_inference"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(INPUT)) {
  stop(
    paste(
      "No existe:",
      INPUT
    )
  )
}

# ============================================================
# CARGA
# ============================================================

x <- read_tsv(
  INPUT,
  show_col_types = FALSE
)

required <- c(
  "sample",
  "producer",
  "producer_name",
  "biological_unit",
  "week",
  "subject",
  "source_MAG",
  "source_MAG_species",
  "source_MAG_RA_main_pct",
  "source_MAG_main_detection_focus",
  "n_independent_loci",
  "n_primary_source_linked",
  "n_sensitivity_locus_signals",
  "fraction_independent_loci_primary",
  "any_independent_locus_primary",
  "all_independent_loci_primary"
)

missing_cols <- setdiff(
  required,
  names(x)
)

if (length(missing_cols) > 0) {
  stop(
    paste(
      "Faltan columnas:",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  )
}

# ============================================================
# FILTRAR MAG DE L. LACTIS
# ============================================================

target_MAG <- "M2__M2_maxbin2.004"

dat <- x %>%
  filter(
    source_MAG == target_MAG
  ) %>%
  mutate(
    biological_unit =
      as.integer(
        biological_unit
      ),
    week =
      as.integer(
        week
      ),
    source_MAG_RA_main_pct =
      as.numeric(
        source_MAG_RA_main_pct
      ),
    source_MAG_main_detection_focus =
      as.integer(
        source_MAG_main_detection_focus
      ),
    n_independent_loci =
      as.integer(
        n_independent_loci
      ),
    n_primary_source_linked =
      as.integer(
        n_primary_source_linked
      ),
    n_sensitivity_locus_signals =
      as.integer(
        n_sensitivity_locus_signals
      )
  ) %>%
  arrange(
    producer,
    biological_unit,
    week
  )

# ============================================================
# QC ESTRUCTURAL
# ============================================================

if (nrow(dat) != 18) {
  stop(
    paste(
      "Se esperaban 18 filas para L. lactis;",
      "se encontraron",
      nrow(dat)
    )
  )
}

if (n_distinct(dat$subject) != 6) {
  stop(
    "Se esperaban exactamente 6 unidades biológicas."
  )
}

if (!setequal(
  unique(dat$week),
  c(0L, 1L, 2L)
)) {
  stop(
    "Las semanas no corresponden a 0,1,2."
  )
}

if (any(
  dat$n_independent_loci != 3
)) {
  stop(
    "L. lactis ya no tiene exactamente 3 loci independientes."
  )
}

qc_subject_weeks <- dat %>%
  count(
    subject
  )

if (any(
  qc_subject_weeks$n != 3
)) {
  stop(
    "Alguna unidad biológica no tiene exactamente 3 semanas."
  )
}

# ============================================================
# QC CONTRA RESULTADOS 81b YA CONGELADOS
# ============================================================

qc_pw <- dat %>%
  group_by(
    producer,
    week
  ) %>%
  summarise(
    total_primary_locus_detections =
      sum(
        n_primary_source_linked
      ),
    source_MAG_detected_units =
      sum(
        source_MAG_main_detection_focus
      ),
    .groups = "drop"
  )

expected_qc <- tibble(
  producer =
    c(
      "L","L","L",
      "M","M","M"
    ),
  week =
    c(
      0,1,2,
      0,1,2
    ),
  expected_locus_detections =
    c(
      5,6,4,
      8,6,5
    ),
  expected_MAG_detected_units =
    c(
      2,2,2,
      3,2,2
    )
)

qc_pw <- qc_pw %>%
  left_join(
    expected_qc,
    by = c(
      "producer",
      "week"
    )
  ) %>%
  mutate(
    locus_match =
      total_primary_locus_detections ==
      expected_locus_detections,
    MAG_match =
      source_MAG_detected_units ==
      expected_MAG_detected_units
  )

if (
  any(!qc_pw$locus_match) ||
  any(!qc_pw$MAG_match)
) {
  stop(
    "QC contra Step81b falló."
  )
}

write_tsv(
  qc_pw,
  file.path(
    OUT,
    "QC_step81b_Llactis_producer_week.tsv"
  )
)

# ============================================================
# TABLA A NIVEL DE UNIDAD BIOLÓGICA
# ============================================================

subject_metrics <- dat %>%
  group_by(
    producer,
    producer_name,
    biological_unit,
    subject,
    source_MAG,
    source_MAG_species
  ) %>%
  summarise(
    RA_w0 =
      source_MAG_RA_main_pct[
        week == 0
      ][1],

    RA_w1 =
      source_MAG_RA_main_pct[
        week == 1
      ][1],

    RA_w2 =
      source_MAG_RA_main_pct[
        week == 2
      ][1],

    loci_w0 =
      n_primary_source_linked[
        week == 0
      ][1],

    loci_w1 =
      n_primary_source_linked[
        week == 1
      ][1],

    loci_w2 =
      n_primary_source_linked[
        week == 2
      ][1],

    sensitivity_loci_w0 =
      n_sensitivity_locus_signals[
        week == 0
      ][1],

    sensitivity_loci_w1 =
      n_sensitivity_locus_signals[
        week == 1
      ][1],

    sensitivity_loci_w2 =
      n_sensitivity_locus_signals[
        week == 2
      ][1],

    n_weeks_MAG_main_detected =
      sum(
        source_MAG_main_detection_focus
      ),

    mean_RA_3weeks =
      mean(
        source_MAG_RA_main_pct
      ),

    median_RA_3weeks =
      median(
        source_MAG_RA_main_pct
      ),

    total_primary_locus_detections_3weeks =
      sum(
        n_primary_source_linked
      ),

    total_sensitivity_locus_signals_3weeks =
      sum(
        n_sensitivity_locus_signals
      ),

    .groups = "drop"
  ) %>%
  mutate(
    delta_RA_w2_minus_w0 =
      RA_w2 - RA_w0,

    delta_RA_w1_minus_w0 =
      RA_w1 - RA_w0,

    delta_loci_w2_minus_w0 =
      loci_w2 - loci_w0,

    delta_loci_w1_minus_w0 =
      loci_w1 - loci_w0,

    RA_trajectory =
      paste(
        sprintf(
          "%.6f",
          RA_w0
        ),
        sprintf(
          "%.6f",
          RA_w1
        ),
        sprintf(
          "%.6f",
          RA_w2
        ),
        sep = " -> "
      ),

    locus_trajectory =
      paste0(
        loci_w0,
        loci_w1,
        loci_w2
      )
  )

if (nrow(subject_metrics) != 6) {
  stop(
    "La tabla subject-level no contiene 6 filas."
  )
}

if (
  sum(
    subject_metrics$producer == "L"
  ) != 3 ||
  sum(
    subject_metrics$producer == "M"
  ) != 3
) {
  stop(
    "El diseño ya no es 3 unidades L vs 3 unidades M."
  )
}

write_tsv(
  subject_metrics,
  file.path(
    OUT,
    "Llactis_subject_level_longitudinal_metrics.tsv"
  )
)

# ============================================================
# DESCRIPTIVOS POR PRODUCTOR
# ============================================================

producer_summary <- subject_metrics %>%
  group_by(
    producer,
    producer_name
  ) %>%
  summarise(
    n_biological_units =
      n(),

    median_RA_w0 =
      median(
        RA_w0
      ),

    median_RA_w1 =
      median(
        RA_w1
      ),

    median_RA_w2 =
      median(
        RA_w2
      ),

    median_delta_RA_w2_minus_w0 =
      median(
        delta_RA_w2_minus_w0
      ),

    min_delta_RA_w2_minus_w0 =
      min(
        delta_RA_w2_minus_w0
      ),

    max_delta_RA_w2_minus_w0 =
      max(
        delta_RA_w2_minus_w0
      ),

    n_RA_decreased_w2_vs_w0 =
      sum(
        delta_RA_w2_minus_w0 < 0
      ),

    n_RA_increased_w2_vs_w0 =
      sum(
        delta_RA_w2_minus_w0 > 0
      ),

    median_loci_w0 =
      median(
        loci_w0
      ),

    median_loci_w1 =
      median(
        loci_w1
      ),

    median_loci_w2 =
      median(
        loci_w2
      ),

    median_delta_loci_w2_minus_w0 =
      median(
        delta_loci_w2_minus_w0
      ),

    min_delta_loci_w2_minus_w0 =
      min(
        delta_loci_w2_minus_w0
      ),

    max_delta_loci_w2_minus_w0 =
      max(
        delta_loci_w2_minus_w0
      ),

    n_locus_count_decreased =
      sum(
        delta_loci_w2_minus_w0 < 0
      ),

    n_locus_count_increased =
      sum(
        delta_loci_w2_minus_w0 > 0
      ),

    .groups = "drop"
  )

write_tsv(
  producer_summary,
  file.path(
    OUT,
    "Llactis_producer_descriptive_summary.tsv"
  )
)

# ============================================================
# EXACT PERMUTATION 3 vs 3
#
# H0:
# Las 6 unidades biológicas son intercambiables respecto
# al productor.
#
# Estadístico:
# media(M) - media(L)
#
# Exhaustivo:
# choose(6,3) = 20 asignaciones.
# ============================================================

cliffs_delta <- function(
  M,
  L
) {

  comparisons <- outer(
    M,
    L,
    FUN = "-"
  )

  (
    sum(
      comparisons > 0
    )
    -
    sum(
      comparisons < 0
    )
  ) /
    length(
      comparisons
    )
}


exact_3v3_test <- function(
  df,
  outcome,
  hypothesis,
  analysis_role
) {

  if (nrow(df) != 6) {
    stop(
      "exact_3v3_test requiere exactamente 6 sujetos."
    )
  }

  y <- df[[outcome]]

  names(y) <-
    df$subject

  actual_M <-
    which(
      df$producer == "M"
    )

  actual_L <-
    which(
      df$producer == "L"
    )

  observed <-
    mean(
      y[
        actual_M
      ]
    ) -
    mean(
      y[
        actual_L
      ]
    )

  observed_median_difference <-
    median(
      y[
        actual_M
      ]
    ) -
    median(
      y[
        actual_L
      ]
    )

  observed_cliff <-
    cliffs_delta(
      y[
        actual_M
      ],
      y[
        actual_L
      ]
    )

  permutations <-
    combn(
      seq_len(
        nrow(df)
      ),
      3,
      simplify = FALSE
    )

  null_rows <-
    lapply(
      seq_along(
        permutations
      ),
      function(i) {

        M_idx <-
          permutations[[i]]

        L_idx <-
          setdiff(
            seq_len(
              nrow(df)
            ),
            M_idx
          )

        stat <-
          mean(
            y[
              M_idx
            ]
          ) -
          mean(
            y[
              L_idx
            ]
          )

        tibble(
          hypothesis =
            hypothesis,
          outcome =
            outcome,
          permutation_id =
            i,
          permuted_M_subjects =
            paste(
              sort(
                df$subject[
                  M_idx
                ]
              ),
              collapse = ";"
            ),
          permuted_L_subjects =
            paste(
              sort(
                df$subject[
                  L_idx
                ]
              ),
              collapse = ";"
            ),
          statistic_mean_M_minus_L =
            stat
        )
      }
    ) %>%
    bind_rows()

  if (nrow(null_rows) != 20) {
    stop(
      "No se generaron exactamente 20 permutaciones."
    )
  }

  tolerance <- 1e-12

  p_two <-
    mean(
      abs(
        null_rows$
          statistic_mean_M_minus_L
      ) >=
        abs(
          observed
        ) -
        tolerance
    )

  result <-
    tibble(
      hypothesis =
        hypothesis,
      analysis_role =
        analysis_role,
      outcome =
        outcome,

      n_L =
        length(
          actual_L
        ),

      n_M =
        length(
          actual_M
        ),

      mean_L =
        mean(
          y[
            actual_L
          ]
        ),

      mean_M =
        mean(
          y[
            actual_M
          ]
        ),

      median_L =
        median(
          y[
            actual_L
          ]
        ),

      median_M =
        median(
          y[
            actual_M
          ]
        ),

      mean_difference_M_minus_L =
        observed,

      median_difference_M_minus_L =
        observed_median_difference,

      cliffs_delta_M_vs_L =
        observed_cliff,

      exact_two_sided_permutation_p =
        p_two,

      n_exact_assignments =
        nrow(
          null_rows
        ),

      minimum_possible_two_sided_p =
        0.10,

      interpretation_constraint =
        "exploratory_n3_vs_n3_effect_size_primary"
    )

  list(
    result =
      result,
    null =
      null_rows
  )
}

# ============================================================
# HIPÓTESIS PREDEFINIDAS
# ============================================================

tests <- list(

  exact_3v3_test(
    subject_metrics,
    "delta_RA_w2_minus_w0",
    "H1_Llactis_RA_change_week2_minus_week0",
    "PRIMARY"
  ),

  exact_3v3_test(
    subject_metrics,
    "delta_loci_w2_minus_w0",
    "H2_direct_locus_count_change_week2_minus_week0",
    "PRIMARY"
  ),

  exact_3v3_test(
    subject_metrics,
    "mean_RA_3weeks",
    "S1_Llactis_mean_RA_across_three_weeks",
    "SECONDARY"
  ),

  exact_3v3_test(
    subject_metrics,
    "total_primary_locus_detections_3weeks",
    "S2_total_direct_locus_detections_across_three_weeks",
    "SECONDARY"
  )
)

test_results <- bind_rows(
  lapply(
    tests,
    function(z)
      z$result
  )
)

null_distributions <- bind_rows(
  lapply(
    tests,
    function(z)
      z$null
  )
)

write_tsv(
  test_results,
  file.path(
    OUT,
    "exact_permutation_producer_tests.tsv"
  )
)

write_tsv(
  null_distributions,
  file.path(
    OUT,
    "exact_permutation_null_distributions.tsv"
  )
)

# ============================================================
# RANK / EFFECT INTERPRETATION
#
# NO cambia p-values.
# Sólo facilita lectura del tamaño/dirección del efecto.
# ============================================================

effect_interpretation <- test_results %>%
  mutate(
    cliff_abs =
      abs(
        cliffs_delta_M_vs_L
      ),

    cliff_magnitude =
      case_when(
        cliff_abs < 0.147 ~
          "negligible",
        cliff_abs < 0.33 ~
          "small",
        cliff_abs < 0.474 ~
          "medium",
        TRUE ~
          "large"
      ),

    effect_direction =
      case_when(
        mean_difference_M_minus_L < 0 ~
          "Marino_lower_or_more_negative_than_Lidieth",
        mean_difference_M_minus_L > 0 ~
          "Marino_higher_or_more_positive_than_Lidieth",
        TRUE ~
          "no_mean_difference"
      )
  ) %>%
  select(
    hypothesis,
    analysis_role,
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
  effect_interpretation,
  file.path(
    OUT,
    "exact_test_effect_interpretation.tsv"
  )
)

# ============================================================
# FIGURA 82B-A
# Trayectoria de abundancia relativa del MAG
# ============================================================

p_ra <- ggplot(
  dat,
  aes(
    x = week,
    y =
      source_MAG_RA_main_pct,
    group =
      subject
  )
) +
  geom_line(
    linewidth = 0.8,
    alpha = 0.8
  ) +
  geom_point(
    size = 2.6
  ) +
  facet_wrap(
    ~ producer_name,
    nrow = 1
  ) +
  scale_x_continuous(
    breaks =
      c(
        0,
        1,
        2
      )
  ) +
  labs(
    x =
      "Week",
    y =
      "L. lactis source-MAG relative abundance (%)",
    title =
      "Longitudinal representation of the Lactococcus lactis source MAG",
    subtitle =
      "Each line is one biological unit; abundance is relative to recovered accepted MAGs"
  ) +
  theme_bw(
    base_size = 11
  )

ggsave(
  file.path(
    OUT,
    "Figure82B_A_Llactis_RA_longitudinal.pdf"
  ),
  p_ra,
  width = 9,
  height = 5
)

ggsave(
  file.path(
    OUT,
    "Figure82B_A_Llactis_RA_longitudinal.png"
  ),
  p_ra,
  width = 9,
  height = 5,
  dpi = 300
)

# ============================================================
# FIGURA 82B-B
# Número de loci ATTRLOC005-007 directamente soportados
# ============================================================

p_loci <- ggplot(
  dat,
  aes(
    x = week,
    y =
      n_primary_source_linked,
    group =
      subject
  )
) +
  geom_line(
    linewidth = 0.8,
    alpha = 0.8
  ) +
  geom_point(
    size = 2.6
  ) +
  facet_wrap(
    ~ producer_name,
    nrow = 1
  ) +
  scale_x_continuous(
    breaks =
      c(
        0,
        1,
        2
      )
  ) +
  scale_y_continuous(
    breaks =
      0:3,
    limits =
      c(
        0,
        3
      )
  ) +
  labs(
    x =
      "Week",
    y =
      "Directly supported source-linked loci (0-3)",
    title =
      "Longitudinal support for ATTRLOC005-007",
    subtitle =
      "A-source-linked evidence; source-MAG support is part of the detection definition"
  ) +
  theme_bw(
    base_size = 11
  )

ggsave(
  file.path(
    OUT,
    "Figure82B_B_Llactis_locus_count_longitudinal.pdf"
  ),
  p_loci,
  width = 9,
  height = 5
)

ggsave(
  file.path(
    OUT,
    "Figure82B_B_Llactis_locus_count_longitudinal.png"
  ),
  p_loci,
  width = 9,
  height = 5,
  dpi = 300
)

# ============================================================
# FIGURA 82B-C
# Cambio individual w2-w0
# ============================================================

delta_long <- subject_metrics %>%
  select(
    producer,
    producer_name,
    subject,
    delta_RA_w2_minus_w0,
    delta_loci_w2_minus_w0
  ) %>%
  pivot_longer(
    cols =
      c(
        delta_RA_w2_minus_w0,
        delta_loci_w2_minus_w0
      ),
    names_to =
      "outcome",
    values_to =
      "delta"
  ) %>%
  mutate(
    outcome =
      recode(
        outcome,
        delta_RA_w2_minus_w0 =
          "MAG RA: week 2 - week 0",
        delta_loci_w2_minus_w0 =
          "Direct locus count: week 2 - week 0"
      )
  )

p_delta <- ggplot(
  delta_long,
  aes(
    x =
      producer_name,
    y =
      delta
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  geom_point(
    size = 3,
    position =
      position_jitter(
        width = 0.06,
        height = 0
      )
  ) +
  facet_wrap(
    ~ outcome,
    scales = "free_y"
  ) +
  labs(
    x =
      "Producer",
    y =
      "Within-unit change",
    title =
      "Subject-level longitudinal change from week 0 to week 2",
    subtitle =
      "Inference uses six biological units, not eighteen samples"
  ) +
  theme_bw(
    base_size = 11
  )

ggsave(
  file.path(
    OUT,
    "Figure82B_C_subject_level_change.pdf"
  ),
  p_delta,
  width = 9,
  height = 5
)

ggsave(
  file.path(
    OUT,
    "Figure82B_C_subject_level_change.png"
  ),
  p_delta,
  width = 9,
  height = 5,
  dpi = 300
)

# ============================================================
# README
# ============================================================

writeLines(
  c(
    "PASO 82b - INFERENCIA LONGITUDINAL EXPLORATORIA",
    "================================================",
    "",
    "Unidad independiente:",
    "6 unidades biologicas (3 Lidieth + 3 Marino).",
    "",
    "Medidas repetidas:",
    "semana 0, 1 y 2.",
    "",
    "HIPOTESIS PRIMARIAS",
    "-------------------",
    "",
    "H1:",
    "Diferencia entre productores en el cambio semana2-semana0",
    "de la abundancia relativa del MAG M2__M2_maxbin2.004.",
    "",
    "H2:",
    "Diferencia entre productores en el cambio semana2-semana0",
    "del numero de loci ATTRLOC005-007 con evidencia",
    "A_source_linked.",
    "",
    "ANALISIS SECUNDARIOS",
    "--------------------",
    "",
    "S1:",
    "Abundancia media del MAG durante las tres semanas.",
    "",
    "S2:",
    "Numero acumulado de detecciones directas de los tres loci",
    "durante las tres semanas.",
    "",
    "METODO",
    "------",
    "",
    "Permutacion exacta exhaustiva de las etiquetas de productor.",
    "",
    "Con seis unidades y grupos 3 vs 3 existen solamente",
    "C(6,3)=20 asignaciones.",
    "",
    "Por simetria, el minimo valor de p bilateral posible es 0.10.",
    "",
    "Por tanto, los resultados deben interpretarse principalmente",
    "mediante magnitud y direccion del efecto, no mediante un",
    "umbral p<0.05.",
    "",
    "Cliff's delta se informa como tamano de efecto ordinal.",
    "",
    "IMPORTANTE",
    "----------",
    "",
    "El numero de loci A_source_linked no es independiente de",
    "la deteccion del MAG fuente, porque el soporte del MAG forma",
    "parte de la definicion A_source_linked.",
    "",
    "No usar la relacion MAG/locus como prueba independiente de",
    "asociacion biologica.",
    "",
    "La deteccion de los loci no demuestra expresion ni produccion",
    "de bacteriocina.",
    "",
    "La abundancia relativa del MAG corresponde al conjunto de",
    "MAGs recuperados/aceptados, no a toda la microbiota."
  ),
  con =
    file.path(
      OUT,
      "README_step82b.txt"
    )
)

# ============================================================
# CONSOLA
# ============================================================

cat(
  "============================================================\n"
)

cat(
  "PASO 82b COMPLETADO\n"
)

cat(
  "============================================================\n"
)

cat(
  "Unidades biologicas:              ",
  nrow(
    subject_metrics
  ),
  "\n"
)

cat(
  "Lidieth:                          ",
  sum(
    subject_metrics$producer == "L"
  ),
  "\n"
)

cat(
  "Marino:                           ",
  sum(
    subject_metrics$producer == "M"
  ),
  "\n"
)

cat(
  "Asignaciones exactas por test:    20\n"
)

cat(
  "Hipotesis primarias:              2\n"
)

cat(
  "Analisis secundarios:             2\n"
)

cat(
  "Minimo p bilateral posible:       0.10\n"
)

cat(
  "\n===== RESULTADOS EXACTOS =====\n"
)

print(
  test_results,
  n = Inf,
  width = Inf
)

cat(
  "\nSalida: ",
  OUT,
  "\n",
  sep = ""
)

cat(
  "PASO 82b FINALIZÓ CORRECTAMENTE; ES SEGURO SALIR.\n"
)
