#!/usr/bin/env Rscript

# ============================================================
# PASO 68
# MAGs que explican los cambios temporales desde semana 0
# ============================================================

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2"
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
library(ggplot2)


# ============================================================
# RUTAS
# ============================================================

user <- Sys.getenv("USER")

root <- file.path(
  "/scratch/global",
  user,
  "Shotgun_MAGs_Turrialba"
)

infile <- file.path(
  root,
  "43_MAG_abundance_taxonomy",
  "MAG_abundance_taxonomy_long.tsv"
)

trajectory_file <- file.path(
  root,
  "45_MAG_longitudinal_inference",
  "BrayCurtis_trajectory_from_week0.tsv"
)

outdir <- file.path(
  root,
  "46_MAG_temporal_drivers"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# CARGA
# ============================================================

dat <- read_tsv(
  infile,
  show_col_types = FALSE
)

trajectory <- read_tsv(
  trajectory_file,
  show_col_types = FALSE
)


stopifnot(
  nrow(dat) == 324
)


# ============================================================
# PREPARAR TAXONOMÍA Y DISEÑO
# ============================================================

dat <- dat %>%
  mutate(
    week = as.integer(week),

    biological_unit =
      as.integer(biological_unit),

    subject = paste0(
      producer,
      biological_unit
    ),

    RA_main =
      as.numeric(
        relative_abundance_main_pct
      )
  )


# Etiqueta más informativa para GCA-* sin especie asignada
dat <- dat %>%
  mutate(
    taxon_display = case_when(

      grepl("^GCA-", taxon_label) &
        !is.na(family) &
        family != "" ~
        paste(
          family,
          taxon_label
        ),

      TRUE ~ taxon_label
    )
  )


# Evitar fusionar accidentalmente dos MAGs con la misma etiqueta
dup_labels <- dat %>%
  distinct(
    MAG,
    taxon_display
  ) %>%
  count(
    taxon_display
  ) %>%
  filter(n > 1) %>%
  pull(
    taxon_display
  )


dat <- dat %>%
  mutate(
    taxon_display = if_else(
      taxon_display %in% dup_labels,
      paste0(
        taxon_display,
        " [",
        MAG,
        "]"
      ),
      taxon_display
    )
  )


# ============================================================
# BASELINE SEMANA 0
# ============================================================

baseline <- dat %>%
  filter(
    week == 0
  ) %>%
  select(
    subject,
    MAG,
    RA_week0 = RA_main
  )


changes <- dat %>%
  filter(
    week %in% c(1, 2)
  ) %>%
  left_join(
    baseline,
    by = c(
      "subject",
      "MAG"
    )
  ) %>%
  mutate(
    delta_RA_pct =
      RA_main - RA_week0,

    abs_delta_RA_pct =
      abs(delta_RA_pct),

    # Porque ambas composiciones suman 100:
    # Bray-Curtis = sum(abs(delta)) / 200
    BC_contribution =
      abs_delta_RA_pct / 200
  )


stopifnot(
  nrow(changes) == 216
)


# ============================================================
# VALIDAR QUE LAS CONTRIBUCIONES RECONSTRUYAN BRAY-CURTIS
# ============================================================

bc_reconstructed <- changes %>%
  group_by(
    subject,
    producer,
    week
  ) %>%
  summarise(
    BrayCurtis_reconstructed =
      sum(
        BC_contribution
      ),
    .groups = "drop"
  )


bc_expected <- trajectory %>%
  filter(
    abundance_definition == "main",
    week %in% c(1, 2)
  ) %>%
  mutate(
    week = as.integer(week)
  ) %>%
  select(
    subject,
    producer,
    week,
    BrayCurtis_from_week0
  )


validation <- bc_reconstructed %>%
  left_join(
    bc_expected,
    by = c(
      "subject",
      "producer",
      "week"
    )
  ) %>%
  mutate(
    absolute_difference =
      abs(
        BrayCurtis_reconstructed -
          BrayCurtis_from_week0
      )
  )


if (
  any(
    validation$absolute_difference >
      1e-8
  )
) {
  stop(
    paste0(
      "ERROR: la descomposición no reproduce ",
      "Bray-Curtis."
    )
  )
}


write_tsv(
  validation,
  file.path(
    outdir,
    "BrayCurtis_contribution_validation.tsv"
  )
)


# ============================================================
# PORCENTAJE DE LA DISTANCIA EXPLICADO POR CADA MAG
# ============================================================

changes <- changes %>%
  group_by(
    subject,
    week
  ) %>%
  mutate(
    total_abs_change =
      sum(
        abs_delta_RA_pct
      ),

    percent_of_BC_change =
      if_else(
        total_abs_change > 0,
        100 *
          abs_delta_RA_pct /
          total_abs_change,
        0
      )
  ) %>%
  ungroup()


# ============================================================
# TABLA COMPLETA POR UNIDAD BIOLÓGICA
# ============================================================

driver_long <- changes %>%
  select(
    subject,
    producer,
    biological_unit,
    week,
    MAG,
    taxon_display,
    RA_week0,
    RA_main,
    delta_RA_pct,
    abs_delta_RA_pct,
    BC_contribution,
    percent_of_BC_change,
    detection_class,
    breadth,
    coverage
  ) %>%
  arrange(
    subject,
    week,
    desc(
      percent_of_BC_change
    )
  )


write_tsv(
  driver_long,
  file.path(
    outdir,
    "MAG_changes_from_week0_long.tsv"
  )
)


# ============================================================
# TOP DRIVERS POR UNIDAD BIOLÓGICA
# ============================================================

top_subject <- driver_long %>%
  group_by(
    subject,
    producer,
    week
  ) %>%
  slice_max(
    order_by =
      percent_of_BC_change,
    n = 5,
    with_ties = FALSE
  ) %>%
  ungroup()


write_tsv(
  top_subject,
  file.path(
    outdir,
    "top5_drivers_per_subject_week.tsv"
  )
)


# ============================================================
# RESUMEN POR PRODUCTOR × SEMANA × MAG
# ============================================================

producer_summary <- changes %>%
  group_by(
    producer,
    week,
    MAG,
    taxon_display
  ) %>%
  summarise(
    n_units = n(),

    mean_RA_week0 =
      mean(
        RA_week0
      ),

    mean_RA_week =
      mean(
        RA_main
      ),

    mean_delta_RA_pct =
      mean(
        delta_RA_pct
      ),

    median_delta_RA_pct =
      median(
        delta_RA_pct
      ),

    mean_abs_delta_RA_pct =
      mean(
        abs_delta_RA_pct
      ),

    mean_BC_contribution =
      mean(
        BC_contribution
      ),

    mean_percent_of_BC_change =
      mean(
        percent_of_BC_change
      ),

    n_increase =
      sum(
        delta_RA_pct > 0
      ),

    n_decrease =
      sum(
        delta_RA_pct < 0
      ),

    n_no_change =
      sum(
        abs(delta_RA_pct) <
          1e-10
      ),

    .groups = "drop"
  ) %>%
  mutate(
    dominant_direction = case_when(

      n_increase > n_decrease ~
        "increase",

      n_decrease > n_increase ~
        "decrease",

      TRUE ~
        "mixed/tie"
    ),

    direction_consistency =
      pmax(
        n_increase,
        n_decrease
      ) /
      n_units
  ) %>%
  arrange(
    producer,
    week,
    desc(
      mean_BC_contribution
    )
  )


write_tsv(
  producer_summary,
  file.path(
    outdir,
    "producer_week_MAG_driver_summary.tsv"
  )
)


# ============================================================
# TOP 8 DRIVERS POR PRODUCTOR × SEMANA
# ============================================================

top_producer <- producer_summary %>%
  group_by(
    producer,
    week
  ) %>%
  slice_max(
    order_by =
      mean_BC_contribution,
    n = 8,
    with_ties = FALSE
  ) %>%
  ungroup()


write_tsv(
  top_producer,
  file.path(
    outdir,
    "top8_drivers_by_producer_week.tsv"
  )
)


# ============================================================
# FIGURA 1
# Distancia respecto a semana 0 por unidad biológica
# ============================================================

traj_main <- trajectory %>%
  filter(
    abundance_definition == "main"
  ) %>%
  mutate(
    week = as.integer(week)
  )


p1 <- ggplot(
  traj_main,
  aes(
    x = week,
    y = BrayCurtis_from_week0,
    group = subject,
    linetype = producer
  )
) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_point(
    size = 2.5
  ) +
  facet_wrap(
    ~ producer,
    labeller = labeller(
      producer = c(
        L = "Lidieth",
        M = "Marino"
      )
    )
  ) +
  scale_x_continuous(
    breaks = c(0, 1, 2)
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  labs(
    x = "Semana",
    y = "Bray-Curtis respecto a semana 0",
    linetype = "Productor",
    title =
      "Magnitud del cambio composicional por unidad biológica"
  ) +
  theme_bw(
    base_size = 11
  )


# ============================================================
# FIGURA 2
# Heatmap de cambio porcentual desde semana 0
# ============================================================

heat_taxa <- changes %>%
  group_by(
    MAG,
    taxon_display
  ) %>%
  summarise(
    total_abs_change =
      sum(
        abs_delta_RA_pct
      ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(
      total_abs_change
    )
  )


changes_heat <- changes %>%
  mutate(
    contrast = paste0(
      subject,
      ": 0→",
      week
    ),

    taxon_display = factor(
      taxon_display,
      levels =
        rev(
          heat_taxa$taxon_display
        )
    )
  )


p2 <- ggplot(
  changes_heat,
  aes(
    x = contrast,
    y = taxon_display,
    fill = delta_RA_pct
  )
) +
  geom_tile() +
  scale_fill_gradient2(
    midpoint = 0,
    name = "Δ abundancia\n(puntos %)"
  ) +
  labs(
    x = "Unidad biológica y contraste",
    y = "MAG / taxón",
    title =
      "Cambios de abundancia relativa respecto a semana 0"
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


# ============================================================
# FIGURA 3
# Principales contribuyentes a Bray-Curtis por productor
# ============================================================

top_driver_taxa <- top_producer %>%
  distinct(
    MAG,
    taxon_display
  ) %>%
  pull(
    MAG
  )


driver_plot_data <- producer_summary %>%
  filter(
    MAG %in%
      top_driver_taxa
  ) %>%
  mutate(
    contrast =
      paste0(
        "0→",
        week
      )
  )


p3 <- ggplot(
  driver_plot_data,
  aes(
    x = mean_percent_of_BC_change,
    y = reorder(
      taxon_display,
      mean_percent_of_BC_change
    )
  )
) +
  geom_col() +
  facet_grid(
    producer ~ contrast,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(
      producer = c(
        L = "Lidieth",
        M = "Marino"
      )
    )
  ) +
  labs(
    x =
      "Contribución media al cambio Bray-Curtis (%)",
    y =
      "MAG / taxón",
    title =
      "MAGs que explican los cambios temporales"
  ) +
  theme_bw(
    base_size = 10
  )


# ============================================================
# GUARDAR FIGURAS
# ============================================================

ggsave(
  file.path(
    outdir,
    "Figura1_BrayCurtis_trajectory.pdf"
  ),
  p1,
  width = 9,
  height = 5
)

ggsave(
  file.path(
    outdir,
    "Figura1_BrayCurtis_trajectory.png"
  ),
  p1,
  width = 9,
  height = 5,
  dpi = 300
)


ggsave(
  file.path(
    outdir,
    "Figura2_delta_abundance_heatmap.pdf"
  ),
  p2,
  width = 12,
  height = 8
)


ggsave(
  file.path(
    outdir,
    "Figura3_temporal_drivers.pdf"
  ),
  p3,
  width = 12,
  height = 10
)


# ============================================================
# PDF MULTIPÁGINA
# ============================================================

pdf(
  file.path(
    outdir,
    "MAG_temporal_drivers_report.pdf"
  ),
  width = 12,
  height = 8.5,
  onefile = TRUE
)

print(p1)
print(p2)
print(p3)

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
# FINAL
# ============================================================

cat(
  "============================================================\n"
)

cat(
  "PASO 68 COMPLETADO\n"
)

cat(
  "============================================================\n"
)

cat(
  "Comparaciones individuales: ",
  nrow(validation),
  "\n",
  sep = ""
)

cat(
  "Descomposición Bray-Curtis validada: OK\n"
)

cat(
  "Salida: ",
  outdir,
  "\n",
  sep = ""
)
