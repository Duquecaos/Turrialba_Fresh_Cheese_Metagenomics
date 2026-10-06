#!/usr/bin/env Rscript

# ============================================================
# PASO 66
# Análisis ecológico descriptivo de 18 MAGs representativos
# ============================================================

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "vegan",
  "scales"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Faltan paquetes de R: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(scales)


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

qcdir <- file.path(
  root,
  "41_MAG_abundance_matrices"
)

outdir <- file.path(
  root,
  "44_MAG_ecology_descriptive"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

long_file <- file.path(
  indir,
  "MAG_abundance_taxonomy_long.tsv"
)

main_matrix_file <- file.path(
  indir,
  "relative_abundance_main_18x18.tsv"
)

strict_matrix_file <- file.path(
  indir,
  "relative_abundance_strict_18x18.tsv"
)

qc_file <- file.path(
  qcdir,
  "sample_QC_summary.tsv"
)


# ============================================================
# CARGA
# ============================================================

long <- read_tsv(
  long_file,
  show_col_types = FALSE
)

main_mat_raw <- read_tsv(
  main_matrix_file,
  show_col_types = FALSE
)

strict_mat_raw <- read_tsv(
  strict_matrix_file,
  show_col_types = FALSE
)

qc <- read_tsv(
  qc_file,
  show_col_types = FALSE
)


# ============================================================
# VALIDACIONES
# ============================================================

stopifnot(
  nrow(long) == 324,
  nrow(main_mat_raw) == 18,
  nrow(strict_mat_raw) == 18,
  ncol(main_mat_raw) == 19,
  ncol(strict_mat_raw) == 19
)

samples <- long %>%
  distinct(
    sample,
    producer,
    biological_unit,
    week
  ) %>%
  mutate(
    week = as.integer(week),
    biological_unit = as.integer(biological_unit),
    subject = paste0(
      producer,
      biological_unit
    )
  ) %>%
  arrange(
    producer,
    biological_unit,
    week
  )

stopifnot(nrow(samples) == 18)


sample_order <- samples$sample


# ============================================================
# ETIQUETAS TAXONÓMICAS
# ============================================================

taxa <- long %>%
  distinct(
    MAG,
    phylum,
    class,
    order,
    family,
    genus,
    species,
    taxon_label,
    taxon_label_rank
  )


# Hacer un poco más informativo el placeholder GTDB GCA-*
taxa <- taxa %>%
  mutate(
    taxon_display = case_when(
      grepl("^GCA-", taxon_label) &
        !is.na(family) &
        family != "" ~
        paste(family, taxon_label),

      TRUE ~ taxon_label
    )
  )


# Si dos MAGs llegaran a tener la misma etiqueta,
# añadir el MAG para que no se fusionen accidentalmente.
duplicates <- taxa %>%
  count(taxon_display) %>%
  filter(n > 1) %>%
  pull(taxon_display)

taxa <- taxa %>%
  mutate(
    taxon_display = if_else(
      taxon_display %in% duplicates,
      paste0(
        taxon_display,
        " [",
        MAG,
        "]"
      ),
      taxon_display
    )
  )


long <- long %>%
  left_join(
    taxa %>%
      select(
        MAG,
        taxon_display
      ),
    by = "MAG"
  ) %>%
  mutate(
    week = as.integer(week),
    biological_unit =
      as.integer(biological_unit),

    subject = paste0(
      producer,
      biological_unit
    ),

    relative_abundance_main_pct =
      as.numeric(
        relative_abundance_main_pct
      ),

    relative_abundance_strict_pct =
      as.numeric(
        relative_abundance_strict_pct
      ),

    coverage =
      as.numeric(coverage),

    breadth =
      as.numeric(breadth)
  )


# ============================================================
# ORDEN TAXONÓMICO PARA FIGURAS
# ============================================================

taxon_order <- long %>%
  group_by(
    MAG,
    taxon_display
  ) %>%
  summarise(
    mean_RA =
      mean(
        relative_abundance_main_pct
      ),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_RA))

long <- long %>%
  mutate(
    taxon_display = factor(
      taxon_display,
      levels = taxon_order$taxon_display
    ),

    sample = factor(
      sample,
      levels = sample_order
    )
  )


# ============================================================
# 1. STACKED BARS LONGITUDINALES
# ============================================================

p_stack <- ggplot(
  long,
  aes(
    x = factor(week),
    y = relative_abundance_main_pct,
    fill = taxon_display
  )
) +
  geom_col(
    width = 0.82
  ) +
  facet_grid(
    producer ~ biological_unit,
    labeller = labeller(
      producer = c(
        L = "Lidieth",
        M = "Marino"
      )
    )
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  scale_fill_viridis_d(
    option = "turbo",
    end = 0.95
  ) +
  labs(
    x = "Semana",
    y = paste0(
      "Abundancia relativa dentro del ",
      "conjunto de MAGs (%)"
    ),
    fill = "MAG / taxón",
    title = paste0(
      "Composición longitudinal de los ",
      "MAGs recuperados"
    ),
    subtitle = paste0(
      "Criterio principal de reclutamiento; ",
      "no representa la microbiota total"
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    strip.background =
      element_rect(fill = "grey95"),
    panel.grid.major.x =
      element_blank()
  )


ggsave(
  file.path(
    outdir,
    "Figura1_stacked_MAGs_longitudinal.pdf"
  ),
  p_stack,
  width = 14,
  height = 8
)

ggsave(
  file.path(
    outdir,
    "Figura1_stacked_MAGs_longitudinal.png"
  ),
  p_stack,
  width = 14,
  height = 8,
  dpi = 300
)


# ============================================================
# 2. HEATMAP
# ============================================================

heat_data <- long %>%
  mutate(
    log_RA =
      log10(
        relative_abundance_main_pct +
          0.01
      )
  )


p_heat <- ggplot(
  heat_data,
  aes(
    x = sample,
    y = taxon_display,
    fill = log_RA
  )
) +
  geom_tile() +
  scale_fill_viridis_c(
    name = expression(
      log[10] * "(RA % + 0.01)"
    )
  ) +
  labs(
    x = "Muestra",
    y = "MAG / taxón",
    title = paste0(
      "Abundancia relativa de MAGs ",
      "recuperados"
    )
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid = element_blank()
  )


ggsave(
  file.path(
    outdir,
    "Figura2_heatmap_MAGs.pdf"
  ),
  p_heat,
  width = 12,
  height = 9
)

ggsave(
  file.path(
    outdir,
    "Figura2_heatmap_MAGs.png"
  ),
  p_heat,
  width = 12,
  height = 9,
  dpi = 300
)


# ============================================================
# 3. TOP 8 MAGs — TRAYECTORIAS
# ============================================================

top8_MAGs <- long %>%
  group_by(
    MAG,
    taxon_display
  ) %>%
  summarise(
    mean_RA =
      mean(
        relative_abundance_main_pct
      ),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_RA)) %>%
  slice_head(n = 8) %>%
  pull(MAG)


trajectory <- long %>%
  filter(
    MAG %in% top8_MAGs
  )


p_traj <- ggplot(
  trajectory,
  aes(
    x = week,
    y = relative_abundance_main_pct,
    group = subject,
    linetype = producer
  )
) +
  geom_line(
    linewidth = 0.7,
    alpha = 0.8
  ) +
  geom_point(
    size = 1.8
  ) +
  facet_wrap(
    ~ taxon_display,
    scales = "free_y",
    ncol = 2
  ) +
  scale_x_continuous(
    breaks = c(0, 1, 2)
  ) +
  labs(
    x = "Semana",
    y = "Abundancia relativa (%)",
    linetype = "Productor",
    title = paste0(
      "Trayectorias longitudinales ",
      "de los MAGs más abundantes"
    )
  ) +
  theme_bw(base_size = 10)


ggsave(
  file.path(
    outdir,
    "Figura3_trayectorias_top8_MAGs.pdf"
  ),
  p_traj,
  width = 11,
  height = 12
)


# ============================================================
# MATRICES NUMÉRICAS
# ============================================================

main_df <- as.data.frame(
  main_mat_raw
)

strict_df <- as.data.frame(
  strict_mat_raw
)

rownames(main_df) <-
  main_df$sample

rownames(strict_df) <-
  strict_df$sample

main_df$sample <- NULL
strict_df$sample <- NULL

main_matrix <- as.matrix(
  main_df
)

strict_matrix <- as.matrix(
  strict_df
)

storage.mode(main_matrix) <-
  "numeric"

storage.mode(strict_matrix) <-
  "numeric"


# Forzar mismo orden de muestras
main_matrix <- main_matrix[
  sample_order,
  ,
  drop = FALSE
]

strict_matrix <- strict_matrix[
  sample_order,
  ,
  drop = FALSE
]


# ============================================================
# 4. BRAY-CURTIS + PCoA
# ============================================================

bray_main <- vegdist(
  main_matrix,
  method = "bray"
)

bray_strict <- vegdist(
  strict_matrix,
  method = "bray"
)


pcoa_main <- cmdscale(
  bray_main,
  k = 2,
  eig = TRUE,
  add = TRUE
)

pcoa_strict <- cmdscale(
  bray_strict,
  k = 2,
  eig = TRUE,
  add = TRUE
)


variance_main <- 100 *
  pcoa_main$eig[
    1:2
  ] /
  sum(
    pcoa_main$eig[
      pcoa_main$eig > 0
    ]
  )


variance_strict <- 100 *
  pcoa_strict$eig[
    1:2
  ] /
  sum(
    pcoa_strict$eig[
      pcoa_strict$eig > 0
    ]
  )


pcoa_main_df <- data.frame(
  sample =
    rownames(pcoa_main$points),

  PCoA1 =
    pcoa_main$points[, 1],

  PCoA2 =
    pcoa_main$points[, 2]
) %>%
  left_join(
    samples,
    by = "sample"
  )


pcoa_strict_df <- data.frame(
  sample =
    rownames(pcoa_strict$points),

  PCoA1 =
    pcoa_strict$points[, 1],

  PCoA2 =
    pcoa_strict$points[, 2]
) %>%
  left_join(
    samples,
    by = "sample"
  )


p_pcoa_main <- ggplot(
  pcoa_main_df,
  aes(
    x = PCoA1,
    y = PCoA2,
    fill = producer,
    shape = factor(week),
    label = sample
  )
) +
  geom_point(
    size = 4,
    stroke = 0.8
  ) +
  scale_shape_manual(
    values = c(
      "0" = 21,
      "1" = 22,
      "2" = 24
    )
  ) +
  geom_text(
    nudge_y = 0.025,
    size = 3,
    check_overlap = TRUE
  ) +
  labs(
    x = sprintf(
      "PCoA1 (%.1f%%)",
      variance_main[1]
    ),
    y = sprintf(
      "PCoA2 (%.1f%%)",
      variance_main[2]
    ),
    fill = "Productor",
    shape = "Semana",
    title = "PCoA Bray-Curtis — criterio principal"
  ) +
  theme_bw(base_size = 11)


p_pcoa_strict <- ggplot(
  pcoa_strict_df,
  aes(
    x = PCoA1,
    y = PCoA2,
    fill = producer,
    shape = factor(week),
    label = sample
  )
) +
  geom_point(
    size = 4,
    stroke = 0.8
  ) +
  scale_shape_manual(
    values = c(
      "0" = 21,
      "1" = 22,
      "2" = 24
    )
  ) +
  geom_text(
    nudge_y = 0.025,
    size = 3,
    check_overlap = TRUE
  ) +
  labs(
    x = sprintf(
      "PCoA1 (%.1f%%)",
      variance_strict[1]
    ),
    y = sprintf(
      "PCoA2 (%.1f%%)",
      variance_strict[2]
    ),
    fill = "Productor",
    shape = "Semana",
    title = "PCoA Bray-Curtis — criterio estricto"
  ) +
  theme_bw(base_size = 11)


ggsave(
  file.path(
    outdir,
    "Figura4_PCoA_main.pdf"
  ),
  p_pcoa_main,
  width = 8,
  height = 6
)

ggsave(
  file.path(
    outdir,
    "Figura5_PCoA_strict.pdf"
  ),
  p_pcoa_strict,
  width = 8,
  height = 6
)


# ============================================================
# 5. DIVERSIDAD DENTRO DEL CONJUNTO DE MAGs
# ============================================================

alpha_main <- data.frame(
  sample = rownames(main_matrix),

  Shannon_MAGset =
    diversity(
      main_matrix / 100,
      index = "shannon"
    ),

  detected_MAGs =
    rowSums(
      main_matrix > 0
    )
)


alpha_strict <- data.frame(
  sample =
    rownames(strict_matrix),

  Shannon_MAGset_strict =
    diversity(
      strict_matrix / 100,
      index = "shannon"
    ),

  detected_MAGs_strict =
    rowSums(
      strict_matrix > 0
    )
)


alpha <- alpha_main %>%
  left_join(
    alpha_strict,
    by = "sample"
  ) %>%
  left_join(
    samples,
    by = "sample"
  )


write_tsv(
  alpha,
  file.path(
    outdir,
    "MAGset_alpha_diversity.tsv"
  )
)


# ============================================================
# 6. SENSIBILIDAD MAIN vs STRICT
# ============================================================

cell_spearman <- cor(
  as.vector(main_matrix),
  as.vector(strict_matrix),
  method = "spearman"
)


lower_main <- as.vector(
  as.dist(bray_main)
)

lower_strict <- as.vector(
  as.dist(bray_strict)
)


distance_spearman <- cor(
  lower_main,
  lower_strict,
  method = "spearman"
)


sample_sensitivity <- lapply(
  sample_order,
  function(s) {

    x <- main_matrix[s, ]
    y <- strict_matrix[s, ]

    # Bray-Curtis entre la composición main
    # y strict de la MISMA muestra.
    denom <- sum(x + y)

    bc <- if (
      denom > 0
    ) {
      sum(abs(x - y)) / denom
    } else {
      NA_real_
    }

    data.frame(
      sample = s,
      main_vs_strict_BrayCurtis = bc,
      main_MAGs =
        sum(x > 0),
      strict_MAGs =
        sum(y > 0)
    )
  }
) %>%
  bind_rows() %>%
  left_join(
    samples,
    by = "sample"
  )


write_tsv(
  sample_sensitivity,
  file.path(
    outdir,
    "main_vs_strict_sample_sensitivity.tsv"
  )
)


global_sensitivity <- data.frame(
  metric = c(
    "Spearman_all_324_abundances",
    "Spearman_pairwise_BrayCurtis"
  ),

  value = c(
    cell_spearman,
    distance_spearman
  )
)


write_tsv(
  global_sensitivity,
  file.path(
    outdir,
    "main_vs_strict_global_sensitivity.tsv"
  )
)


# ============================================================
# 7. QC + ALPHA
# ============================================================

qc2 <- qc %>%
  mutate(
    bowtie2_alignment_numeric =
      as.numeric(
        sub(
          "%",
          "",
          bowtie2_overall_alignment_rate
        )
      )
  ) %>%
  left_join(
    alpha,
    by = c(
      "sample",
      "producer",
      "biological_unit",
      "week"
    )
  )


write_tsv(
  qc2,
  file.path(
    outdir,
    "sample_QC_ecology.tsv"
  )
)


# ============================================================
# PDF MULTIPÁGINA
# ============================================================

pdf(
  file.path(
    outdir,
    "MAG_ecology_descriptive_report.pdf"
  ),
  width = 12,
  height = 8.5,
  onefile = TRUE
)

print(p_stack)
print(p_heat)
print(p_traj)
print(p_pcoa_main)
print(p_pcoa_strict)

dev.off()


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
# RESUMEN FINAL
# ============================================================

cat(
  "============================================================\n"
)

cat(
  "PASO 66 COMPLETADO\n"
)

cat(
  "============================================================\n"
)

cat(
  "Muestras: ",
  nrow(samples),
  "\n",
  sep = ""
)

cat(
  "MAGs: ",
  ncol(main_matrix),
  "\n",
  sep = ""
)

cat(
  sprintf(
    "Spearman abundancias main vs strict: %.4f\n",
    cell_spearman
  )
)

cat(
  sprintf(
    "Spearman distancias Bray-Curtis main vs strict: %.4f\n",
    distance_spearman
  )
)

cat(
  "Salida: ",
  outdir,
  "\n",
  sep = ""
)

cat(
  "No se realizaron pruebas inferenciales en este paso.\n"
)
