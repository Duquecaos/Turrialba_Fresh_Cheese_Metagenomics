#!/usr/bin/env Rscript

# ============================================================
# PASO 67
# Inferencia longitudinal sobre el conjunto de MAGs recuperados
#
# Diseño:
#   Productor: L / M
#   Unidad biológica: 1, 2, 3
#   Subject: L1-L3 / M1-M3
#   Semana: 0, 1, 2
#
# Main = análisis principal
# Strict = sensibilidad
# ============================================================


# ============================================================
# PAQUETES
# ============================================================

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "vegan",
  "permute"
)

missing <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing) > 0) {
  stop(
    "Faltan paquetes: ",
    paste(missing, collapse = ", ")
  )
}

library(readr)
library(dplyr)
library(tidyr)
library(vegan)
library(permute)


# ============================================================
# RUTAS
# ============================================================

user <- Sys.getenv("USER")

root <- file.path(
  "/scratch/global",
  user,
  "Shotgun_MAGs_Turrialba"
)

indir <- file.path(
  root,
  "43_MAG_abundance_taxonomy"
)

ecodir <- file.path(
  root,
  "44_MAG_ecology_descriptive"
)

outdir <- file.path(
  root,
  "45_MAG_longitudinal_inference"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


main_file <- file.path(
  indir,
  "relative_abundance_main_18x18.tsv"
)

strict_file <- file.path(
  indir,
  "relative_abundance_strict_18x18.tsv"
)

long_file <- file.path(
  indir,
  "MAG_abundance_taxonomy_long.tsv"
)

alpha_file <- file.path(
  ecodir,
  "MAGset_alpha_diversity.tsv"
)


# ============================================================
# CARGAR DATOS
# ============================================================

main_raw <- read_tsv(
  main_file,
  show_col_types = FALSE
)

strict_raw <- read_tsv(
  strict_file,
  show_col_types = FALSE
)

long <- read_tsv(
  long_file,
  show_col_types = FALSE
)

alpha <- read_tsv(
  alpha_file,
  show_col_types = FALSE
)


meta <- long %>%
  distinct(
    sample,
    producer,
    biological_unit,
    week
  ) %>%
  mutate(
    producer = factor(
      producer,
      levels = c("L", "M")
    ),

    biological_unit =
      factor(biological_unit),

    week_f = factor(
      week,
      levels = c(0, 1, 2)
    ),

    subject = factor(
      paste0(
        producer,
        biological_unit
      )
    )
  ) %>%
  arrange(
    producer,
    biological_unit,
    week_f
  )


stopifnot(
  nrow(meta) == 18,
  length(unique(meta$subject)) == 6
)


# ============================================================
# CONVERTIR MATRICES
# ============================================================

make_matrix <- function(x, sample_order) {

  x <- as.data.frame(x)

  rownames(x) <- x$sample

  x$sample <- NULL

  m <- as.matrix(x)

  storage.mode(m) <- "numeric"

  m <- m[
    sample_order,
    ,
    drop = FALSE
  ]

  return(m)
}


main_matrix <- make_matrix(
  main_raw,
  meta$sample
)

strict_matrix <- make_matrix(
  strict_raw,
  meta$sample
)


stopifnot(
  nrow(main_matrix) == 18,
  ncol(main_matrix) == 18,
  nrow(strict_matrix) == 18,
  ncol(strict_matrix) == 18
)


# ============================================================
# PERMUTACIONES PARA MEDIDAS REPETIDAS
#
# Solo se permutan semanas DENTRO de cada subject.
#
# Esto permite evaluar:
#   - semana
#   - productor:semana
#
# El término subject absorbe diferencias persistentes entre
# unidades biológicas.
#
# NO usamos el p de subject como prueba biológica.
# ============================================================

set.seed(20260914)

perm_repeated <- how(
  within = Within(
    type = "free"
  ),
  blocks = meta$subject,
  nperm = 9999
)


# ============================================================
# FUNCIÓN PERMANOVA LONGITUDINAL
# ============================================================

run_repeated_permanova <- function(
  matrix,
  label
) {

  result <- adonis2(
    matrix ~
      subject +
      week_f +
      producer:week_f,

    data = meta,

    permutations = perm_repeated,

    method = "bray",

    by = "terms"
  )

  tab <- as.data.frame(result)

  tab$term <- rownames(tab)

  rownames(tab) <- NULL

  tab <- tab %>%
    select(
      term,
      everything()
    ) %>%
    mutate(
      abundance_definition = label,

      interpretation = case_when(

        term == "subject" ~
          paste0(
            "Blocking/nuisance term; ",
            "p-value not interpreted"
          ),

        term == "week_f" ~
          paste0(
            "Within-subject temporal effect"
          ),

        term == "producer:week_f" ~
          paste0(
            "Differential temporal trajectory ",
            "between producers"
          ),

        TRUE ~ ""
      )
    )

  return(tab)
}


perm_main <- run_repeated_permanova(
  main_matrix,
  "main"
)

perm_strict <- run_repeated_permanova(
  strict_matrix,
  "strict"
)


write_tsv(
  perm_main,
  file.path(
    outdir,
    "PERMANOVA_repeated_main.tsv"
  )
)

write_tsv(
  perm_strict,
  file.path(
    outdir,
    "PERMANOVA_repeated_strict.tsv"
  )
)


# ============================================================
# PRODUCTOR:
# PROMEDIAR LAS 3 SEMANAS POR UNIDAD BIOLÓGICA
#
# Aquí sí tenemos 6 observaciones independientes.
# ============================================================

subject_meta <- meta %>%
  distinct(
    subject,
    producer
  ) %>%
  arrange(subject)


subject_mean_matrix <- function(matrix) {

  result <- lapply(
    levels(meta$subject),
    function(s) {

      idx <- which(
        meta$subject == s
      )

      colMeans(
        matrix[
          idx,
          ,
          drop = FALSE
        ]
      )
    }
  )

  result <- do.call(
    rbind,
    result
  )

  rownames(result) <-
    levels(meta$subject)

  return(result)
}


main_subject <- subject_mean_matrix(
  main_matrix
)

strict_subject <- subject_mean_matrix(
  strict_matrix
)


subject_meta <- subject_meta[
  match(
    rownames(main_subject),
    subject_meta$subject
  ),
  ,
  drop = FALSE
]


# ============================================================
# Con 6 subjects hay solamente C(6,3)=20 asignaciones únicas
# de 3 vs 3 productores.
#
# Usamos enumeración completa de las permutaciones posibles.
# ============================================================

perm_between <- how(
  within = Within(
    type = "free"
  ),
  complete = TRUE
)


run_producer_permanova <- function(
  matrix,
  label
) {

  res <- adonis2(
    matrix ~ producer,

    data = subject_meta,

    permutations = perm_between,

    method = "bray",

    by = "terms"
  )

  tab <- as.data.frame(res)

  tab$term <- rownames(tab)

  rownames(tab) <- NULL

  tab <- tab %>%
    select(
      term,
      everything()
    ) %>%
    mutate(
      abundance_definition = label,

      analysis_level =
        "Mean composition of each biological unit across weeks",

      independent_units = 6,

      producer_group_sizes =
        "3 L vs 3 M",

      caution =
        paste0(
          "Very low inferential resolution; ",
          "only 20 unique 3-vs-3 label allocations"
        )
    )

  return(tab)
}


producer_main <- run_producer_permanova(
  main_subject,
  "main"
)

producer_strict <- run_producer_permanova(
  strict_subject,
  "strict"
)


write_tsv(
  producer_main,
  file.path(
    outdir,
    "PERMANOVA_producer_subjectMeans_main.tsv"
  )
)

write_tsv(
  producer_strict,
  file.path(
    outdir,
    "PERMANOVA_producer_subjectMeans_strict.tsv"
  )
)


# ============================================================
# DISPERSIÓN ENTRE PRODUCTORES A NIVEL DE SUBJECT
#
# Evaluación auxiliar para interpretar el PERMANOVA de
# productor. Solo 3 unidades por grupo: interpretar con cautela.
# ============================================================

run_dispersion <- function(
  matrix,
  label
) {

  d <- vegdist(
    matrix,
    method = "bray"
  )

  bd <- betadisper(
    d,
    subject_meta$producer
  )

  pt <- permutest(
    bd,
    permutations = perm_between
  )

  tab <- as.data.frame(
    pt$tab
  )

  tab$term <- rownames(tab)

  rownames(tab) <- NULL

  tab <- tab %>%
    select(
      term,
      everything()
    ) %>%
    mutate(
      abundance_definition = label,
      independent_units = 6,
      caution =
        "Only 3 biological units per producer"
    )

  return(tab)
}


disp_main <- run_dispersion(
  main_subject,
  "main"
)

disp_strict <- run_dispersion(
  strict_subject,
  "strict"
)


write_tsv(
  disp_main,
  file.path(
    outdir,
    "PERMDISP_producer_subjectMeans_main.tsv"
  )
)

write_tsv(
  disp_strict,
  file.path(
    outdir,
    "PERMDISP_producer_subjectMeans_strict.tsv"
  )
)


# ============================================================
# SHANNON:
# ANOVA DE MEDIDAS REPETIDAS EXPLORATORIA
#
# Se usa únicamente para Shannon dentro del conjunto de MAGs.
# No representa alfa-diversidad del metagenoma completo.
# ============================================================

alpha2 <- alpha %>%
  mutate(
    producer = factor(
      producer,
      levels = c("L", "M")
    ),

    biological_unit =
      factor(biological_unit),

    week_f = factor(
      week,
      levels = c(0, 1, 2)
    ),

    subject = factor(subject)
  )


aov_main <- aov(
  Shannon_MAGset ~
    producer * week_f +
    Error(subject / week_f),

  data = alpha2
)


aov_strict <- aov(
  Shannon_MAGset_strict ~
    producer * week_f +
    Error(subject / week_f),

  data = alpha2
)


capture.output(
  summary(aov_main),
  file = file.path(
    outdir,
    "Shannon_repeated_ANOVA_main.txt"
  )
)

capture.output(
  summary(aov_strict),
  file = file.path(
    outdir,
    "Shannon_repeated_ANOVA_strict.txt"
  )
)


# ============================================================
# DESCRIPTIVOS DE SHANNON POR PRODUCTOR × SEMANA
# ============================================================

alpha_summary <- alpha2 %>%
  group_by(
    producer,
    week_f
  ) %>%
  summarise(
    n_subjects = n(),

    Shannon_main_mean =
      mean(
        Shannon_MAGset
      ),

    Shannon_main_sd =
      sd(
        Shannon_MAGset
      ),

    Shannon_strict_mean =
      mean(
        Shannon_MAGset_strict
      ),

    Shannon_strict_sd =
      sd(
        Shannon_MAGset_strict
      ),

    .groups = "drop"
  )


write_tsv(
  alpha_summary,
  file.path(
    outdir,
    "Shannon_producer_week_summary.tsv"
  )
)


# ============================================================
# MAGs DETECTADOS — SOLO DESCRIPTIVO
# ============================================================

richness_summary <- alpha2 %>%
  group_by(
    producer,
    week_f
  ) %>%
  summarise(
    n_subjects = n(),

    detected_MAGs_main_mean =
      mean(
        detected_MAGs
      ),

    detected_MAGs_main_range =
      paste0(
        min(detected_MAGs),
        "-",
        max(detected_MAGs)
      ),

    detected_MAGs_strict_mean =
      mean(
        detected_MAGs_strict
      ),

    detected_MAGs_strict_range =
      paste0(
        min(detected_MAGs_strict),
        "-",
        max(detected_MAGs_strict)
      ),

    .groups = "drop"
  )


write_tsv(
  richness_summary,
  file.path(
    outdir,
    "detected_MAGs_producer_week_summary.tsv"
  )
)


# ============================================================
# TRAYECTORIA:
# BRAY-CURTIS DE CADA SUBJECT CONTRA SU SEMANA 0
# ============================================================

baseline_distance <- function(
  matrix,
  label
) {

  out <- list()

  k <- 1

  for (s in levels(meta$subject)) {

    idx <- which(
      meta$subject == s
    )

    m <- meta[
      idx,
      ,
      drop = FALSE
    ]

    x <- matrix[
      idx,
      ,
      drop = FALSE
    ]

    base_idx <- which(
      m$week_f == "0"
    )

    for (j in seq_len(nrow(m))) {

      denominator <- sum(
        x[base_idx, ] +
        x[j, ]
      )

      bc <- if (
        denominator > 0
      ) {
        sum(
          abs(
            x[base_idx, ] -
              x[j, ]
          )
        ) /
          denominator
      } else {
        NA_real_
      }

      out[[k]] <- data.frame(
        abundance_definition = label,
        subject = s,
        producer =
          as.character(
            m$producer[j]
          ),
        week =
          as.character(
            m$week_f[j]
          ),
        BrayCurtis_from_week0 =
          bc
      )

      k <- k + 1
    }
  }

  bind_rows(out)
}


traj <- bind_rows(
  baseline_distance(
    main_matrix,
    "main"
  ),

  baseline_distance(
    strict_matrix,
    "strict"
  )
)


write_tsv(
  traj,
  file.path(
    outdir,
    "BrayCurtis_trajectory_from_week0.tsv"
  )
)


# ============================================================
# README METODOLÓGICO
# ============================================================

writeLines(
  c(
    "PASO 67 - Inferencia longitudinal del conjunto de MAGs",
    "",
    "IMPORTANT:",
    "",
    "1. The 18 samples are not treated as 18 independent biological replicates.",
    "   Six biological units (L1-L3, M1-M3) were repeatedly sampled at weeks 0, 1 and 2.",
    "",
    "2. Repeated-measures PERMANOVA:",
    "   Formula: subject + week + producer:week",
    "   Bray-Curtis distances.",
    "   Permutations are restricted within biological unit.",
    "   The subject term is a nuisance/blocking term and its permutation p-value is not interpreted.",
    "   Week tests the temporal effect.",
    "   Producer:week tests whether temporal trajectories differ between producers.",
    "",
    "3. Producer main effect:",
    "   The three weekly abundance profiles are averaged within each biological unit first.",
    "   PERMANOVA therefore uses six independent profiles: 3 L and 3 M.",
    "   Statistical resolution is intrinsically low because only 20 unique 3-vs-3 producer allocations exist.",
    "",
    "4. Shannon:",
    "   Repeated-measures ANOVA is exploratory because there are only three biological units per producer.",
    "   Shannon refers only to diversity within the recovered MAG set, not total metagenomic alpha diversity.",
    "",
    "5. Main is the primary abundance definition.",
    "   Strict is a sensitivity analysis."
  ),
  con = file.path(
    outdir,
    "README_longitudinal_inference.txt"
  )
)


# ============================================================
# SESSION INFO
# ============================================================

sink(
  file.path(
    outdir,
    "R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# FINAL
# ============================================================

cat(
  "============================================================\n"
)

cat(
  "PASO 67 COMPLETADO\n"
)

cat(
  "============================================================\n"
)

cat(
  "Unidades biológicas independientes: 6\n"
)

cat(
  "Mediciones por unidad: 3\n"
)

cat(
  "Muestras totales: 18\n"
)

cat(
  "\nPERMANOVA repetida MAIN:\n"
)

print(
  perm_main
)

cat(
  "\nPERMANOVA productor MAIN (promedios por unidad):\n"
)

print(
  producer_main
)

cat(
  "\nPERMDISP productor MAIN:\n"
)

print(
  disp_main
)

cat(
  "\nSalida: ",
  outdir,
  "\n",
  sep = ""
)
