#!/usr/bin/env Rscript

# ============================================================
# PASO 74
# Análisis longitudinal descriptivo de los seis loci
# bacteriocínicos/RiPP independientes
#
# NO realiza pruebas inferenciales productor × semana.
# n = 3 unidades biológicas por productor.
# ============================================================

options(stringsAsFactors = FALSE)

ROOT <- file.path(
  "/scratch/global",
  Sys.getenv("USER"),
  "Shotgun_MAGs_Turrialba"
)

INFILE <- file.path(
  ROOT,
  "51_bacteriocin_final_evidence",
  "final_independent_locus_evidence_18x6.tsv"
)

OUTDIR <- file.path(
  ROOT,
  "52_bacteriocin_longitudinal_analysis"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# DEPENDENCIAS
# ============================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop(
    "ERROR: ggplot2 no está instalado en este R."
  )
}

library(ggplot2)


# ============================================================
# LEER DATOS
# ============================================================

d <- read.delim(
  INFILE,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (nrow(d) != 108) {
  stop(
    paste0(
      "ERROR: se esperaban 108 filas; hay ",
      nrow(d)
    )
  )
}


# ============================================================
# VALIDACIONES
# ============================================================

expected_loci <- c(
  "ATTRLOC001",
  "ATTRLOC003",
  "ATTRLOC004",
  "ATTRLOC005",
  "ATTRLOC006",
  "ATTRLOC007"
)

expected_samples <- c(
  "L1_1","L1_2","L1_3",
  "L2_1","L2_2","L2_3",
  "L3_1","L3_2","L3_3",
  "M1_1","M1_2","M1_3",
  "M2_1","M2_2","M2_3",
  "M3_1","M3_2","M3_3"
)

if (!setequal(unique(d$locus_id), expected_loci)) {
  stop("ERROR: los seis ATTRLOC no coinciden.")
}

if (!setequal(unique(d$sample), expected_samples)) {
  stop("ERROR: las 18 muestras no coinciden.")
}

if (anyDuplicated(d[c("sample", "locus_id")])) {
  stop("ERROR: existen combinaciones sample × locus duplicadas.")
}


# ============================================================
# FACTORES
# ============================================================

d$locus_id <- factor(
  d$locus_id,
  levels = expected_loci
)

d$sample <- factor(
  d$sample,
  levels = rev(expected_samples)
)

d$week <- as.integer(d$week)

d$producer_name <- factor(
  d$producer_name,
  levels = c("Lidieth", "Marino")
)

d$biological_unit <- as.integer(
  d$biological_unit
)

d$subject <- paste0(
  d$producer,
  d$biological_unit
)

tier_levels <- c(
  "D_no_support",
  "C_source_MAG_without_locus_support",
  "B_locus_signal_source_unresolved",
  "A_source_linked"
)

d$evidence_tier_main <- factor(
  d$evidence_tier_main,
  levels = tier_levels
)


# ============================================================
# NUMÉRICOS
# ============================================================

num_fields <- c(
  "primary_source_linked_detection",
  "strict_source_linked_detection",
  "sensitivity_locus_signal",
  "coverage_per_million_input_pairs",
  "filtered_pairs_CPM",
  "primary_coverage_per_million",
  "primary_filtered_pairs_CPM",
  "coverage",
  "breadth",
  "reads_mean_PID"
)

for (x in num_fields) {
  d[[x]] <- suppressWarnings(
    as.numeric(d[[x]])
  )
}


# ============================================================
# 1. CONTROL MAIN VS STRICT
# ============================================================

discordant <- d[
  d$primary_source_linked_detection !=
    d$strict_source_linked_detection,
]

write.table(
  discordant,
  file.path(
    OUTDIR,
    "QC_main_vs_strict.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (nrow(discordant) != 0) {
  warning(
    "Hay discrepancias main vs strict."
  )
}


# ============================================================
# 2. RESUMEN POR LOCUS
# ============================================================

locus_summary <- do.call(
  rbind,
  lapply(
    expected_loci,
    function(loc) {

      x <- d[
        d$locus_id == loc,
      ]

      data.frame(
        locus_id = loc,
        source_MAG = x$source_MAG[1],

        n_samples = nrow(x),

        A_source_linked =
          sum(
            x$evidence_tier_main ==
              "A_source_linked",
            na.rm = TRUE
          ),

        B_source_unresolved =
          sum(
            x$evidence_tier_main ==
              "B_locus_signal_source_unresolved",
            na.rm = TRUE
          ),

        C_MAG_without_locus =
          sum(
            x$evidence_tier_main ==
              "C_source_MAG_without_locus_support",
            na.rm = TRUE
          ),

        D_no_support =
          sum(
            x$evidence_tier_main ==
              "D_no_support",
            na.rm = TRUE
          ),

        stringsAsFactors = FALSE
      )
    }
  )
)

write.table(
  locus_summary,
  file.path(
    OUTDIR,
    "summary_by_locus.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 3. RESUMEN POR MUESTRA
# ============================================================

sample_summary <- do.call(
  rbind,
  lapply(
    expected_samples,
    function(s) {

      x <- d[
        as.character(d$sample) == s,
      ]

      data.frame(
        sample = s,
        producer = x$producer[1],
        producer_name =
          as.character(x$producer_name[1]),
        biological_unit =
          x$biological_unit[1],
        week =
          x$week[1],

        A_source_linked =
          sum(
            x$primary_source_linked_detection,
            na.rm = TRUE
          ),

        locus_signal_sensitivity =
          sum(
            x$sensitivity_locus_signal,
            na.rm = TRUE
          ),

        unresolved_B =
          sum(
            x$evidence_tier_main ==
              "B_locus_signal_source_unresolved",
            na.rm = TRUE
          ),

        stringsAsFactors = FALSE
      )
    }
  )
)

write.table(
  sample_summary,
  file.path(
    OUTDIR,
    "summary_by_sample.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 4. PRODUCTOR × SEMANA
# SOLO DESCRIPTIVO
# ============================================================

pw <- aggregate(
  A_source_linked ~
    producer_name + week,
  data = sample_summary,
  FUN = function(x) {
    c(
      mean = mean(x),
      median = median(x),
      min = min(x),
      max = max(x)
    )
  }
)

pw_out <- data.frame(
  producer_name = pw$producer_name,
  week = pw$week,
  mean_loci = pw$A_source_linked[, "mean"],
  median_loci = pw$A_source_linked[, "median"],
  min_loci = pw$A_source_linked[, "min"],
  max_loci = pw$A_source_linked[, "max"]
)

write.table(
  pw_out,
  file.path(
    OUTDIR,
    "summary_producer_week.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 5. TRAYECTORIAS 0/1 DURANTE SEMANAS 0-1-2
#
# Una trayectoria por:
# productor × unidad biológica × locus
#
# Ejemplos:
# 111 = estable presente
# 000 = estable ausente
# 011 = ganancia tras semana 0
# 110 = pérdida en semana 2
# ============================================================

trajectory_rows <- list()
k <- 1

subjects <- unique(
  d[c(
    "producer",
    "producer_name",
    "biological_unit",
    "subject"
  )]
)

for (i in seq_len(nrow(subjects))) {

  ss <- subjects[i, ]

  for (loc in expected_loci) {

    x <- d[
      d$subject == ss$subject &
        d$locus_id == loc,
    ]

    x <- x[
      order(x$week),
    ]

    if (nrow(x) != 3) {
      stop(
        paste(
          "ERROR trayectoria:",
          ss$subject,
          loc,
          "n =",
          nrow(x)
        )
      )
    }

    pattern <- paste0(
      x$primary_source_linked_detection,
      collapse = ""
    )

    trajectory_rows[[k]] <- data.frame(
      producer =
        ss$producer,
      producer_name =
        as.character(ss$producer_name),
      biological_unit =
        ss$biological_unit,
      subject =
        ss$subject,
      locus_id =
        loc,
      week0 =
        x$primary_source_linked_detection[
          x$week == 0
        ],
      week1 =
        x$primary_source_linked_detection[
          x$week == 1
        ],
      week2 =
        x$primary_source_linked_detection[
          x$week == 2
        ],
      trajectory =
        pattern,
      stringsAsFactors = FALSE
    )

    k <- k + 1
  }
}

trajectories <- do.call(
  rbind,
  trajectory_rows
)

write.table(
  trajectories,
  file.path(
    OUTDIR,
    "longitudinal_trajectories.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


trajectory_summary <- as.data.frame(
  table(
    trajectories$trajectory
  ),
  stringsAsFactors = FALSE
)

names(trajectory_summary) <- c(
  "trajectory",
  "count"
)

trajectory_summary <- trajectory_summary[
  order(
    -trajectory_summary$count,
    trajectory_summary$trajectory
  ),
]

write.table(
  trajectory_summary,
  file.path(
    OUTDIR,
    "trajectory_pattern_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 6. PREVALENCIA PRODUCTOR × SEMANA × LOCUS
# ============================================================

prev <- aggregate(
  primary_source_linked_detection ~
    producer_name + week + locus_id,
  data = d,
  FUN = sum
)

names(prev)[
  names(prev) ==
    "primary_source_linked_detection"
] <- "detected"

prev$n <- 3
prev$prevalence <- prev$detected / prev$n

write.table(
  prev,
  file.path(
    OUTDIR,
    "primary_prevalence_producer_week_locus.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# FIGURA 1
# HEATMAP DE TIERS
# ============================================================

p1 <- ggplot(
  d,
  aes(
    x = locus_id,
    y = sample,
    fill = evidence_tier_main
  )
) +
  geom_tile(
    linewidth = 0.3
  ) +
  labs(
    title =
      "Evidencia de contextos bacteriocínicos/RiPP",
    subtitle =
      "Seis loci independientes; ATTRLOC002 excluido",
    x =
      "Contexto candidato",
    y =
      "Muestra",
    fill =
      "Evidencia"
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
    panel.grid =
      element_blank()
  )

ggsave(
  file.path(
    OUTDIR,
    "Figura74A_heatmap_evidence_tiers.pdf"
  ),
  p1,
  width = 9,
  height = 7
)

ggsave(
  file.path(
    OUTDIR,
    "Figura74A_heatmap_evidence_tiers.png"
  ),
  p1,
  width = 9,
  height = 7,
  dpi = 300
)


# ============================================================
# FIGURA 2
# NÚMERO DE LOCI A POR UNIDAD BIOLÓGICA
# ============================================================

sample_summary$week <- as.integer(
  sample_summary$week
)

sample_summary$subject <- paste0(
  sample_summary$producer,
  sample_summary$biological_unit
)

p2 <- ggplot(
  sample_summary,
  aes(
    x = week,
    y = A_source_linked,
    group = subject
  )
) +
  geom_line() +
  geom_point(
    size = 2.5
  ) +
  facet_wrap(
    ~ producer_name
  ) +
  scale_x_continuous(
    breaks = c(0, 1, 2)
  ) +
  scale_y_continuous(
    breaks = 0:6,
    limits = c(0, 6)
  ) +
  labs(
    title =
      "Trayectoria longitudinal de loci fuente-vinculados",
    x =
      "Semana",
    y =
      "Número de loci con evidencia A"
  ) +
  theme_bw(
    base_size = 11
  )

ggsave(
  file.path(
    OUTDIR,
    "Figura74B_loci_por_unidad_longitudinal.pdf"
  ),
  p2,
  width = 8,
  height = 5
)

ggsave(
  file.path(
    OUTDIR,
    "Figura74B_loci_por_unidad_longitudinal.png"
  ),
  p2,
  width = 8,
  height = 5,
  dpi = 300
)


# ============================================================
# FIGURA 3
# PREVALENCIA DESCRIPTIVA
# ============================================================

p3 <- ggplot(
  prev,
  aes(
    x = week,
    y = prevalence,
    group = locus_id
  )
) +
  geom_line() +
  geom_point(
    size = 2
  ) +
  facet_grid(
    locus_id ~ producer_name
  ) +
  scale_x_continuous(
    breaks = c(0, 1, 2)
  ) +
  scale_y_continuous(
    breaks = c(
      0,
      1/3,
      2/3,
      1
    ),
    limits = c(0, 1)
  ) +
  labs(
    title =
      "Prevalencia descriptiva de evidencia fuente-vinculada",
    subtitle =
      "Cada punto representa 0/3, 1/3, 2/3 o 3/3 unidades biológicas",
    x =
      "Semana",
    y =
      "Proporción de unidades biológicas"
  ) +
  theme_bw(
    base_size = 10
  )

ggsave(
  file.path(
    OUTDIR,
    "Figura74C_prevalencia_productor_semana_locus.pdf"
  ),
  p3,
  width = 9,
  height = 11
)

ggsave(
  file.path(
    OUTDIR,
    "Figura74C_prevalencia_productor_semana_locus.png"
  ),
  p3,
  width = 9,
  height = 11,
  dpi = 300
)


# ============================================================
# FIGURA 4
# SEÑAL NORMALIZADA SOLO PARA EVIDENCIA A
# log10(1 + coverage por millón de pares)
# ============================================================

d$signal_A <- ifelse(
  d$primary_source_linked_detection == 1,
  d$primary_coverage_per_million,
  0
)

d$log_signal_A <- log10(
  1 + d$signal_A
)

p4 <- ggplot(
  d,
  aes(
    x = locus_id,
    y = sample,
    fill = log_signal_A
  )
) +
  geom_tile(
    linewidth = 0.3
  ) +
  labs(
    title =
      "Señal normalizada de loci fuente-vinculados",
    subtitle =
      "log10(1 + cobertura por millón de fragmentos paired de entrada)",
    x =
      "Contexto candidato",
    y =
      "Muestra",
    fill =
      "Señal\nnormalizada"
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
    panel.grid =
      element_blank()
  )

ggsave(
  file.path(
    OUTDIR,
    "Figura74D_heatmap_normalized_signal.pdf"
  ),
  p4,
  width = 9,
  height = 7
)

ggsave(
  file.path(
    OUTDIR,
    "Figura74D_heatmap_normalized_signal.png"
  ),
  p4,
  width = 9,
  height = 7,
  dpi = 300
)


# ============================================================
# PDF MULTIPÁGINA
# ============================================================

pdf(
  file.path(
    OUTDIR,
    "Paso74_figuras_longitudinales.pdf"
  ),
  width = 9,
  height = 8
)

print(p1)
print(p2)
print(p3)
print(p4)

dev.off()


# ============================================================
# README / INTERPRETACIÓN
# ============================================================

sink(
  file.path(
    OUTDIR,
    "README_step74.txt"
  )
)

cat(
  "PASO 74 - ANALISIS LONGITUDINAL DESCRIPTIVO\n"
)

cat(
  "===========================================\n\n"
)

cat(
  "Muestras: 18\n"
)

cat(
  "Productores: 2\n"
)

cat(
  "Unidades biologicas por productor: 3\n"
)

cat(
  "Semanas: 0, 1, 2\n"
)

cat(
  "Loci independientes: 6\n\n"
)

cat(
  "El analisis principal utiliza evidencia A_source_linked:\n"
)

cat(
  "contexto del locus + MAG fuente respaldados simultaneamente.\n\n"
)

cat(
  "No se realizaron pruebas inferenciales productor x semana,\n"
)

cat(
  "debido al numero reducido de unidades biologicas independientes.\n"
)

cat(
  "Las prevalencias son descriptivas.\n\n"
)

cat(
  "La señal cuantitativa se expresa como cobertura del locus\n"
)

cat(
  "por millon de fragmentos paired de entrada y NO debe\n"
)

cat(
  "interpretarse como abundancia relativa de toda la microbiota\n"
)

cat(
  "ni como numero de copias del locus.\n"
)

sink()


# ============================================================
# CONSOLA
# ============================================================

cat(
  "============================================================\n"
)

cat(
  "PASO 74 COMPLETADO\n"
)

cat(
  "============================================================\n"
)

cat(
  "Filas analizadas:       ",
  nrow(d),
  "\n",
  sep = ""
)

cat(
  "Detecciones A:          ",
  sum(
    d$primary_source_linked_detection,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "Discordancias main/strict: ",
  nrow(discordant),
  "\n",
  sep = ""
)

cat(
  "Trayectorias evaluadas: ",
  nrow(trajectories),
  "\n",
  sep = ""
)

cat(
  "Salida: ",
  OUTDIR,
  "\n",
  sep = ""
)
