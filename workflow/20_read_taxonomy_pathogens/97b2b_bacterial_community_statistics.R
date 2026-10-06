#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(width = 200)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
    stop("Usage: 97b2b_bacterial_community_statistics.R <ROOT>")
}

ROOT <- args[1]

IN <- file.path(
    ROOT,
    "98_read_taxonomy",
    "97B2A_matrix_QC"
)

OUT <- file.path(
    ROOT,
    "98_read_taxonomy",
    "97B2B_bacterial_community"
)

dir.create(
    OUT,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# 1. Dependencies
# ============================================================

required_packages <- c(
    "vegan",
    "permute",
    "ggplot2"
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
        paste(
            "Missing R packages:",
            paste(
                missing_packages,
                collapse = ", "
            )
        )
    )
}


suppressPackageStartupMessages(
    library(vegan)
)

suppressPackageStartupMessages(
    library(permute)
)

suppressPackageStartupMessages(
    library(ggplot2)
)


set.seed(20260929)


# ============================================================
# 2. Helpers
# ============================================================

write_tsv <- function(
    x,
    filename
) {

    write.table(
        x,
        file = file.path(
            OUT,
            filename
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
    )
}


read_tax_matrix <- function(
    filename
) {

    dat <- read.delim(
        file.path(
            IN,
            filename
        ),
        check.names = FALSE,
        quote = "",
        comment.char = ""
    )

    if (
        !all(
            c(
                "taxonomy_id",
                "name"
            )
            %in%
            colnames(dat)
        )
    ) {

        stop(
            paste(
                "Bad taxonomy matrix:",
                filename
            )
        )
    }


    info <- dat[
        ,
        c(
            "taxonomy_id",
            "name"
        )
    ]


    key <- make.unique(
        paste(
            info$taxonomy_id,
            info$name,
            sep = "|"
        )
    )

    mat <- as.matrix(
        dat[
            ,
            setdiff(
                colnames(dat),
                c(
                    "taxonomy_id",
                    "name"
                )
            )
        ]
    )

    storage.mode(
        mat
    ) <- "numeric"

    rownames(
        mat
    ) <- key


    info$key <- key


    list(
        info = info,
        matrix = mat
    )
}


normalize_columns <- function(
    m
) {

    sums <- colSums(
        m
    )

    if (
        any(
            sums <= 0
        )
    ) {

        stop(
            "A filtered matrix has sample(s) with sum <= 0"
        )
    }

    sweep(
        m,
        2,
        sums,
        "/"
    )
}


save_plot <- function(
    p,
    stem,
    width = 10,
    height = 6
) {

    pdf(
        file.path(
            OUT,
            paste0(
                stem,
                ".pdf"
            )
        ),
        width = width,
        height = height,
        useDingbats = FALSE
    )

    print(
        p
    )

    dev.off()


    svg(
        file.path(
            OUT,
            paste0(
                stem,
                ".svg"
            )
        ),
        width = width,
        height = height
    )

    print(
        p
    )

    dev.off()
}


adonis_table <- function(
    object,
    analysis_name
) {

    d <- as.data.frame(
        object
    )

    d$term <- rownames(
        d
    )

    rownames(
        d
    ) <- NULL

    d$analysis <- analysis_name

    d[
        ,
        c(
            "analysis",
            "term",
            setdiff(
                colnames(d),
                c(
                    "analysis",
                    "term"
                )
            )
        )
    ]
}


# ============================================================
# 3. Read bacterial matrix + metadata
# ============================================================

bac <- read_tax_matrix(
    "97B2A_genus_relative_bacteria.tsv"
)

M <- bac$matrix


metadata <- read.delim(
    file.path(
        IN,
        "97B2A_sample_metadata.tsv"
    ),
    check.names = FALSE
)


expected_samples <- colnames(
    M
)


if (
    !setequal(
        metadata$sample,
        expected_samples
    )
) {

    stop(
        "Metadata samples do not match bacterial matrix"
    )
}


metadata <- metadata[
    match(
        expected_samples,
        metadata$sample
    ),
]


metadata$producer_code <- factor(
    metadata$producer_code,
    levels = c(
        "L",
        "M"
    )
)


metadata$time_code <- factor(
    metadata$within_unit_time_code,
    levels = c(
        "1",
        "2",
        "3"
    )
)


metadata$subject_id <- factor(
    metadata$subject_id
)


metadata$sample <- factor(
    metadata$sample,
    levels = expected_samples
)


write_tsv(
    metadata,
    "97B2B_sample_metadata.tsv"
)


# ============================================================
# 4. Filter definitions
#
# PRIMARY:
# >=0.01% bacterial relative composition
# in >=2 samples.
#
# STRICT sensitivity:
# >=0.1% in >=2 samples.
# ============================================================

keep_unfiltered <- (
    rowSums(
        M > 0
    )
    >= 1
)


keep_primary <- (
    rowSums(
        M >= 0.0001
    )
    >= 2
)


keep_strict <- (
    rowSums(
        M >= 0.001
    )
    >= 2
)


filter_definitions <- data.frame(
    filter = c(
        "unfiltered",
        "primary",
        "strict"
    ),

    abundance_threshold_fraction = c(
        0,
        0.0001,
        0.001
    ),

    abundance_threshold_pct = c(
        0,
        0.01,
        0.1
    ),

    minimum_samples = c(
        1,
        2,
        2
    ),

    interpretation = c(
        "all_Bracken_retained_bacterial_genera",
        "PRIMARY_ge0.01pct_in_ge2_samples",
        "SENSITIVITY_ge0.1pct_in_ge2_samples"
    )
)


write_tsv(
    filter_definitions,
    "97B2B_filter_definitions.tsv"
)


filter_sets <- list(
    unfiltered = keep_unfiltered,
    primary = keep_primary,
    strict = keep_strict
)


filtered_matrices <- list()


filter_audit <- list()


for (
    filter_name
    in names(
        filter_sets
    )
) {

    keep <- filter_sets[[filter_name]]


    raw_filtered <- M[
        keep,
        ,
        drop = FALSE
    ]


    retained_before_renorm <- colSums(
        raw_filtered
    )


    normalized <- normalize_columns(
        raw_filtered
    )


    filtered_matrices[[filter_name]] <- normalized


    tmp <- data.frame(
        filter = filter_name,
        sample = colnames(M),
        genera_retained_global = sum(
            keep
        ),
        composition_retained_before_renormalization = retained_before_renorm,
        composition_removed_before_renormalization = 1 - retained_before_renorm
    )


    filter_audit[[filter_name]] <- tmp


    info <- bac$info[
        keep,
        ,
        drop = FALSE
    ]


    matrix_out <- data.frame(
        taxonomy_id = info$taxonomy_id,
        name = info$name,
        normalized,
        check.names = FALSE
    )


    write_tsv(
        matrix_out,
        paste0(
            "97B2B_bacterial_genus_matrix_",
            filter_name,
            ".tsv"
        )
    )
}


filter_audit_df <- do.call(
    rbind,
    filter_audit
)


write_tsv(
    filter_audit_df,
    "97B2B_filter_audit.tsv"
)


# ============================================================
# 5. Primary alpha diversity
#
# These are descriptive composition metrics.
# "Observed genera" = retained classified genera,
# not an ecological richness estimator.
# ============================================================

PRIMARY <- filtered_matrices[["primary"]]


X <- t(
    PRIMARY
)


alpha <- data.frame(
    sample = rownames(
        X
    ),

    observed_retained_genera = rowSums(
        X > 0
    ),

    Shannon = diversity(
        X,
        index = "shannon"
    ),

    Simpson = diversity(
        X,
        index = "simpson"
    ),

    inverse_Simpson = diversity(
        X,
        index = "invsimpson"
    ),

    effective_Shannon_genera = exp(
        diversity(
            X,
            index = "shannon"
        )
    )
)


alpha <- merge(
    metadata,
    alpha,
    by = "sample",
    sort = FALSE
)


alpha <- alpha[
    match(
        expected_samples,
        alpha$sample
    ),
]


write_tsv(
    alpha,
    "97B2B_alpha_diversity_primary.tsv"
)


# ============================================================
# 6. Bray-Curtis + PCoA
# ============================================================

bray_primary <- vegdist(
    X,
    method = "bray"
)


write.table(
    as.matrix(
        bray_primary
    ),
    file = file.path(
        OUT,
        "97B2B_BrayCurtis_primary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    col.names = NA
)


pcoa <- wcmdscale(
    bray_primary,
    k = 2,
    eig = TRUE,
    add = "lingoes"
)


coords <- as.data.frame(
    pcoa$points
)


colnames(
    coords
) <- c(
    "PCoA1",
    "PCoA2"
)


coords$sample <- rownames(
    coords
)


positive_eig <- pcoa$eig[
    pcoa$eig > 0
]


pct1 <- (
    100
    * pcoa$eig[1]
    / sum(
        positive_eig
    )
)


pct2 <- (
    100
    * pcoa$eig[2]
    / sum(
        positive_eig
    )
)


coords <- merge(
    metadata,
    coords,
    by = "sample",
    sort = FALSE
)


coords <- coords[
    match(
        expected_samples,
        coords$sample
    ),
]


write_tsv(
    coords,
    "97B2B_PCoA_primary_coordinates.tsv"
)


eig_df <- data.frame(
    axis = seq_along(
        pcoa$eig
    ),

    eigenvalue = pcoa$eig
)


write_tsv(
    eig_df,
    "97B2B_PCoA_primary_eigenvalues.tsv"
)


# ============================================================
# 7. Repeated-measures PERMANOVA: TIME
#
# Permutations restricted within biological unit.
# ============================================================

control_time <- how(
    nperm = 9999,
    blocks = metadata$subject_id
)


permanova_time <- adonis2(
    bray_primary ~ time_code,
    data = metadata,
    permutations = control_time
)


permanova_time_table <- adonis_table(
    permanova_time,
    "primary_global_time_repeated"
)


write_tsv(
    permanova_time_table,
    "97B2B_PERMANOVA_time_primary.tsv"
)


# ============================================================
# 8. Time within each producer
# ============================================================

producer_time_results <- list()


for (
    producer
    in c(
        "L",
        "M"
    )
) {

    idx <- (
        metadata$producer_code
        == producer
    )


    Xsub <- X[
        idx,
        ,
        drop = FALSE
    ]


    msub <- droplevels(
        metadata[
            idx,
            ,
            drop = FALSE
        ]
    )


    dsub <- vegdist(
        Xsub,
        method = "bray"
    )


    control_sub <- how(
        nperm = 9999,
        blocks = msub$subject_id
    )


    fit <- adonis2(
        dsub ~ time_code,
        data = msub,
        permutations = control_sub
    )


    producer_time_results[[producer]] <- adonis_table(
        fit,
        paste0(
            "primary_time_within_producer_",
            producer
        )
    )
}


producer_time_table <- do.call(
    rbind,
    producer_time_results
)


write_tsv(
    producer_time_table,
    "97B2B_PERMANOVA_time_within_producer.tsv"
)


# ============================================================
# 9. Producer comparison at independent subject level
#
# Average the three time-code compositions within each
# biological unit. This avoids treating the 18 repeated
# samples as independent producer replicates.
# ============================================================

subjects <- levels(
    metadata$subject_id
)


subject_profiles <- do.call(
    rbind,
    lapply(
        subjects,
        function(
            subject
        ) {

            idx <- (
                metadata$subject_id
                == subject
            )

            colMeans(
                X[
                    idx,
                    ,
                    drop = FALSE
                ]
            )
        }
    )
)


rownames(
    subject_profiles
) <- subjects


subject_profiles <- subject_profiles /
    rowSums(
        subject_profiles
    )


subject_metadata <- unique(
    metadata[
        ,
        c(
            "subject_id",
            "producer_code"
        )
    ]
)


subject_metadata <- subject_metadata[
    match(
        subjects,
        subject_metadata$subject_id
    ),
]


subject_metadata$producer_code <- droplevels(
    subject_metadata$producer_code
)


subject_bray <- vegdist(
    subject_profiles,
    method = "bray"
)


producer_fit <- adonis2(
    subject_bray ~ producer_code,
    data = subject_metadata,
    permutations = 9999
)


producer_table <- adonis_table(
    producer_fit,
    "primary_producer_subject_level_mean"
)


write_tsv(
    producer_table,
    "97B2B_PERMANOVA_producer_subject_level.tsv"
)


subject_profile_out <- data.frame(
    subject_id = rownames(
        subject_profiles
    ),
    subject_profiles,
    check.names = FALSE
)


write_tsv(
    subject_profile_out,
    "97B2B_subject_mean_bacterial_profiles.tsv"
)


# ============================================================
# 10. Producer dispersion diagnostic
# ============================================================

bd_producer <- betadisper(
    subject_bray,
    group = subject_metadata$producer_code
)


bd_producer_perm <- permutest(
    bd_producer,
    permutations = 9999
)


capture.output(
    bd_producer_perm,
    file = file.path(
        OUT,
        "97B2B_PERMDISP_producer_subject_level.txt"
    )
)


producer_dispersion <- data.frame(
    subject_id = names(
        bd_producer$distances
    ),

    distance_to_centroid = as.numeric(
        bd_producer$distances
    )
)


producer_dispersion <- merge(
    subject_metadata,
    producer_dispersion,
    by = "subject_id",
    sort = FALSE
)


write_tsv(
    producer_dispersion,
    "97B2B_producer_dispersion_distances.tsv"
)


# ============================================================
# 11. Sensitivity analysis:
# unfiltered vs primary vs strict
# ============================================================

sensitivity_results <- list()


for (
    filter_name
    in names(
        filtered_matrices
    )
) {

    MF <- filtered_matrices[[filter_name]]


    XF <- t(
        MF
    )


    DF <- vegdist(
        XF,
        method = "bray"
    )


    controlF <- how(
        nperm = 9999,
        blocks = metadata$subject_id
    )


    time_fit <- adonis2(
        DF ~ time_code,
        data = metadata,
        permutations = controlF
    )


    time_row <- as.data.frame(
        time_fit
    )[
        "time_code",
        ,
        drop = FALSE
    ]


    subjectF <- do.call(
        rbind,
        lapply(
            subjects,
            function(
                subject
            ) {

                idx <- (
                    metadata$subject_id
                    == subject
                )

                colMeans(
                    XF[
                        idx,
                        ,
                        drop = FALSE
                    ]
                )
            }
        )
    )


    subjectF <- subjectF /
        rowSums(
            subjectF
        )


    producer_distF <- vegdist(
        subjectF,
        method = "bray"
    )


    producer_fitF <- adonis2(
        producer_distF ~ producer_code,
        data = subject_metadata,
        permutations = 9999
    )


    producer_row <- as.data.frame(
        producer_fitF
    )[
        "producer_code",
        ,
        drop = FALSE
    ]


    sensitivity_results[[filter_name]] <- data.frame(
        filter = filter_name,

        genera = nrow(
            MF
        ),

        time_R2 = time_row$R2,

        time_F = time_row$F,

        time_p = time_row[["Pr(>F)"]],

        producer_subject_R2 = producer_row$R2,

        producer_subject_F = producer_row$F,

        producer_subject_p = producer_row[["Pr(>F)"]]
    )
}


sensitivity_table <- do.call(
    rbind,
    sensitivity_results
)


write_tsv(
    sensitivity_table,
    "97B2B_PERMANOVA_filter_sensitivity.tsv"
)


# ============================================================
# 12. Top 15 bacterial genera plot
# ============================================================

mean_abundance <- rowMeans(
    PRIMARY
)


top15_keys <- names(
    sort(
        mean_abundance,
        decreasing = TRUE
    )
)[
    seq_len(
        min(
            15,
            length(
                mean_abundance
            )
        )
    )
]


key_to_name <- setNames(
    bac$info$name,
    bac$info$key
)


composition_long <- list()


for (
    sample
    in colnames(
        PRIMARY
    )
) {

    vals <- PRIMARY[
        ,
        sample
    ]


    top_vals <- vals[
        top15_keys
    ]


    tmp <- data.frame(
        sample = sample,
        genus = key_to_name[
            top15_keys
        ],
        relative_abundance = as.numeric(
            top_vals
        )
    )


    other <- 1 - sum(
        top_vals
    )


    tmp <- rbind(
        tmp,
        data.frame(
            sample = sample,
            genus = "Other",
            relative_abundance = other
        )
    )


    composition_long[[sample]] <- tmp
}


composition_long <- do.call(
    rbind,
    composition_long
)


composition_long$sample <- factor(
    composition_long$sample,
    levels = expected_samples
)


write_tsv(
    composition_long,
    "97B2B_top15_bacterial_composition_plot_data.tsv"
)


p_comp <- ggplot(
    composition_long,
    aes(
        x = sample,
        y = relative_abundance,
        fill = genus
    )
) +
    geom_col(
        width = 0.9
    ) +
    scale_y_continuous(
        labels = function(x) {
            paste0(
                round(
                    x * 100,
                    1
                ),
                "%"
            )
        },
        expand = c(
            0,
            0
        )
    ) +
    labs(
        x = "Sample",
        y = "Relative composition within bacterial reads",
        fill = "Genus",
        title = "Shotgun read-based bacterial composition"
    ) +
    theme_minimal(
        base_size = 11
    ) +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1
        ),
        panel.grid.major.x = element_blank()
    )


save_plot(
    p_comp,
    "97B2B_Fig1_bacterial_composition_top15",
    width = 12,
    height = 7
)


# ============================================================
# 13. PCoA figure
# ============================================================

p_pcoa <- ggplot(
    coords,
    aes(
        x = PCoA1,
        y = PCoA2,
        group = subject_id
    )
) +
    geom_path(
        alpha = 0.45
    ) +
    geom_point(
        aes(
            shape = producer_code,
            fill = time_code
        ),
        size = 3
    ) +
    geom_text(
        aes(
            label = sample
        ),
        nudge_y = 0.01,
        size = 3,
        check_overlap = TRUE
    ) +
    labs(
        x = sprintf(
            "PCoA1 (%.1f%%)",
            pct1
        ),
        y = sprintf(
            "PCoA2 (%.1f%%)",
            pct2
        ),
        shape = "Producer code",
        fill = "Time code",
        title = "Bray-Curtis PCoA of bacterial genus composition"
    ) +
    theme_minimal(
        base_size = 11
    )


save_plot(
    p_pcoa,
    "97B2B_Fig2_bacterial_BrayCurtis_PCoA",
    width = 9,
    height = 7
)


# ============================================================
# 14. Shannon trajectory figure
# ============================================================

p_alpha <- ggplot(
    alpha,
    aes(
        x = time_code,
        y = Shannon,
        group = subject_id
    )
) +
    geom_line(
        alpha = 0.6
    ) +
    geom_point(
        aes(
            shape = producer_code
        ),
        size = 3
    ) +
    geom_text(
        aes(
            label = sample
        ),
        nudge_y = 0.03,
        size = 3,
        check_overlap = TRUE
    ) +
    labs(
        x = "Time code",
        y = "Shannon diversity",
        shape = "Producer code",
        title = "Bacterial genus diversity across sample codes"
    ) +
    theme_minimal(
        base_size = 11
    )


save_plot(
    p_alpha,
    "97B2B_Fig3_bacterial_Shannon",
    width = 8,
    height = 6
)


# ============================================================
# 15. Fungal signals as fraction of ALL input pairs
#
# Descriptive only.
# ============================================================

master <- read.delim(
    file.path(
        IN,
        "97B2A_genus_master.tsv"
    ),
    check.names = FALSE
)


fungi <- master[
    master$taxonomy_category
    == "Fungi",
    ,
    drop = FALSE
]


fungi$fraction_all_input_pairs <- as.numeric(
    fungi$fraction_all_input_pairs
)


fungal_global <- aggregate(
    fraction_all_input_pairs
    ~ name,
    data = fungi,
    FUN = sum
)


fungal_global <- fungal_global[
    order(
        fungal_global$fraction_all_input_pairs,
        decreasing = TRUE
    ),
]


top_fungi <- head(
    fungal_global$name,
    10
)


fungal_plot <- list()


for (
    sample
    in expected_samples
) {

    fs <- fungi[
        fungi$sample
        == sample,
        ,
        drop = FALSE
    ]


    vals <- setNames(
        fs$fraction_all_input_pairs,
        fs$name
    )


    top_vals <- sapply(
        top_fungi,
        function(
            taxon
        ) {

            if (
                taxon
                %in%
                names(
                    vals
                )
            ) {
                vals[[taxon]]
            } else {
                0
            }
        }
    )


    all_fungal <- sum(
        fs$fraction_all_input_pairs
    )


    other <- max(
        0,
        all_fungal
        - sum(
            top_vals
        )
    )


    fungal_plot[[sample]] <- rbind(
        data.frame(
            sample = sample,
            genus = top_fungi,
            fraction_all_input_pairs = as.numeric(
                top_vals
            )
        ),
        data.frame(
            sample = sample,
            genus = "Other fungi",
            fraction_all_input_pairs = other
        )
    )
}


fungal_plot_df <- do.call(
    rbind,
    fungal_plot
)


fungal_plot_df$sample <- factor(
    fungal_plot_df$sample,
    levels = expected_samples
)


write_tsv(
    fungal_plot_df,
    "97B2B_top10_fungi_fraction_all_pairs.tsv"
)


p_fungi <- ggplot(
    fungal_plot_df,
    aes(
        x = sample,
        y = fraction_all_input_pairs,
        fill = genus
    )
) +
    geom_col() +
    scale_y_continuous(
        labels = function(x) {
            paste0(
                round(
                    x * 100,
                    2
                ),
                "%"
            )
        }
    ) +
    labs(
        x = "Sample",
        y = "Estimated fungal reads / all input pairs",
        fill = "Fungal genus",
        title = "Shotgun read-based fungal signal"
    ) +
    theme_minimal(
        base_size = 11
    ) +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1
        ),
        panel.grid.major.x = element_blank()
    )


save_plot(
    p_fungi,
    "97B2B_Fig4_fungal_signal",
    width = 12,
    height = 6
)


# ============================================================
# 16. Classification QC
# ============================================================

qc <- read.delim(
    file.path(
        IN,
        "97B2A_sample_QC_summary.tsv"
    ),
    check.names = FALSE
)


qc_categories <- list()


for (
    i
    in seq_len(
        nrow(
            qc
        )
    )
) {

    r <- qc[
        i,
    ]


    vals <- c(
        Bacteria =
            as.numeric(
                r$bacteria_pct_all_pairs
            ),

        Fungi =
            as.numeric(
                r$fungi_pct_all_pairs
            ),

        Human =
            as.numeric(
                r$human_pct_all_pairs
            ),

        Unclassified =
            100
            - as.numeric(
                r$kraken_classified_pct
            )
    )


    other <- max(
        0,
        100
        - sum(
            vals
        )
    )


    vals <- c(
        vals,
        Other = other
    )


    qc_categories[[r$sample]] <- data.frame(
        sample = r$sample,
        category = names(
            vals
        ),
        pct_all_input_pairs = as.numeric(
            vals
        )
    )
}


qc_long <- do.call(
    rbind,
    qc_categories
)


qc_long$sample <- factor(
    qc_long$sample,
    levels = expected_samples
)


write_tsv(
    qc_long,
    "97B2B_classification_QC_plot_data.tsv"
)


p_qc <- ggplot(
    qc_long,
    aes(
        x = sample,
        y = pct_all_input_pairs,
        fill = category
    )
) +
    geom_col() +
    labs(
        x = "Sample",
        y = "% of all input pairs",
        fill = "Classification",
        title = "Read-based taxonomic classification QC"
    ) +
    theme_minimal(
        base_size = 11
    ) +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1
        ),
        panel.grid.major.x = element_blank()
    )


save_plot(
    p_qc,
    "97B2B_Fig5_classification_QC",
    width = 12,
    height = 6
)


# ============================================================
# 17. Software versions
# ============================================================

software <- data.frame(
    software = c(
        "R",
        "vegan",
        "permute",
        "ggplot2"
    ),

    version = c(
        R.version.string,
        as.character(
            packageVersion(
                "vegan"
            )
        ),
        as.character(
            packageVersion(
                "permute"
            )
        ),
        as.character(
            packageVersion(
                "ggplot2"
            )
        )
    )
)


write_tsv(
    software,
    "97B2B_software_versions.tsv"
)


# ============================================================
# 18. Global summary
# ============================================================

primary_audit <- filter_audit_df[
    filter_audit_df$filter
    == "primary",
]


summary_rows <- data.frame(
    metric = c(
        "samples",
        "independent_biological_units",
        "producers",
        "original_bacterial_genera",
        "primary_retained_genera",
        "strict_retained_genera",
        "primary_min_composition_retained_pct",
        "primary_max_composition_retained_pct",
        "primary_filter",
        "primary_beta_metric",
        "alpha_metrics",
        "time_test",
        "producer_test",
        "producer_independent_n",
        "formal_producer_time_interaction",
        "species_level_inferential_analysis",
        "fungal_inferential_analysis",
        "amplicon_data_used"
    ),

    value = c(
        nrow(
            metadata
        ),

        length(
            unique(
                metadata$subject_id
            )
        ),

        length(
            unique(
                metadata$producer_code
            )
        ),

        nrow(
            M
        ),

        sum(
            keep_primary
        ),

        sum(
            keep_strict
        ),

        sprintf(
            "%.6f",
            100
            * min(
                primary_audit$
                    composition_retained_before_renormalization
            )
        ),

        sprintf(
            "%.6f",
            100
            * max(
                primary_audit$
                    composition_retained_before_renormalization
            )
        ),

        "genus_ge0.01pct_bacterial_composition_in_ge2_samples",

        "Bray-Curtis",

        "observed_retained_genera;Shannon;Simpson;inverse_Simpson;effective_Shannon",

        "PERMANOVA_time_factor_with_permutations_restricted_within_biological_unit",

        "PERMANOVA_on_mean_profile_of_each_of_6_independent_biological_units",

        "6_total_3_per_producer",

        "NO",

        "NO_primary_genus_level_only",

        "NO_descriptive_fungal_signal_only",

        "NO"
    )
)


write_tsv(
    summary_rows,
    "97B2B_global_summary.tsv"
)


# ============================================================
# 19. Methodological scope / guardrails
# ============================================================

scope <- data.frame(
    field = c(
        "primary_community",
        "primary_rank",
        "primary_filter",
        "filter_rationale",
        "zero_interpretation",
        "time_labels",
        "repeated_measure_unit",
        "time_PERMANOVA",
        "producer_PERMANOVA",
        "producer_test_limitation",
        "interaction_test",
        "alpha_diversity_interpretation",
        "fungal_analysis",
        "archaea_analysis",
        "viral_analysis",
        "human_signal",
        "Salmonella_or_other_health_taxa",
        "relative_abundance_interpretation",
        "amplicon_data"
    ),

    value = c(
        "bacterial_read_based_genus_composition",

        "genus",

        "ge0.01pct_within_bacteria_in_at_least_2_samples",

        "removes_low_abundance_long_tail_while_preserving_nearly_all_composition",

        "evaluated_same_database_zero_estimated_reads_at_retained_rank_not_structural_NA",

        "codes_1_2_3_only_no_biological_week_relabeling_here",

        "L1_L2_L3_M1_M2_M3",

        "restricted_permutations_within_biological_unit",

        "mean_profile_per_independent_biological_unit_before_testing_producer",

        "only_3_independent_biological_units_per_producer_interpret_R2_and_P_cautiously",

        "not_formally_tested",

        "composition_diversity_not_species_richness_estimator",

        "descriptive_separate_read_based_signal",

        "descriptive_only_due_very_low_signal",

        "Kraken_virus_fraction_not_used_as_primary_virome_evidence_use_dedicated_geNomad_vOTU_analysis",

        "excluded_from_microbiota",

        "screening_signal_only_requires_targeted_validation_before_species_or_food_safety_claim",

        "Bracken_estimated_read_composition_not_cellular_absolute_abundance",

        "not_used"
    )
)


write_tsv(
    scope,
    "97B2B_methodological_scope.tsv"
)


# ============================================================
# 20. Final validation
# ============================================================

required <- c(
    "97B2B_bacterial_genus_matrix_primary.tsv",
    "97B2B_filter_audit.tsv",
    "97B2B_alpha_diversity_primary.tsv",
    "97B2B_BrayCurtis_primary.tsv",
    "97B2B_PCoA_primary_coordinates.tsv",
    "97B2B_PERMANOVA_time_primary.tsv",
    "97B2B_PERMANOVA_time_within_producer.tsv",
    "97B2B_PERMANOVA_producer_subject_level.tsv",
    "97B2B_PERMANOVA_filter_sensitivity.tsv",
    "97B2B_global_summary.tsv",
    "97B2B_methodological_scope.tsv",
    "97B2B_Fig1_bacterial_composition_top15.svg",
    "97B2B_Fig2_bacterial_BrayCurtis_PCoA.svg",
    "97B2B_Fig3_bacterial_Shannon.svg",
    "97B2B_Fig4_fungal_signal.svg",
    "97B2B_Fig5_classification_QC.svg"
)


for (
    f
    in required
) {

    p <- file.path(
        OUT,
        f
    )

    if (
        !file.exists(
            p
        )
        ||
        file.info(
            p
        )$size <= 0
    ) {

        stop(
            paste(
                "Missing/empty output:",
                f
            )
        )
    }
}


cat(
    "SAMPLES=",
    nrow(
        metadata
    ),
    "\n",
    sep = ""
)

cat(
    "BIOLOGICAL_UNITS=",
    length(
        unique(
            metadata$subject_id
        )
    ),
    "\n",
    sep = ""
)

cat(
    "ORIGINAL_BACTERIAL_GENERA=",
    nrow(
        M
    ),
    "\n",
    sep = ""
)

cat(
    "PRIMARY_RETAINED_GENERA=",
    sum(
        keep_primary
    ),
    "\n",
    sep = ""
)

cat(
    "STRICT_RETAINED_GENERA=",
    sum(
        keep_strict
    ),
    "\n",
    sep = ""
)

cat(
    "97B2B=PASS\n"
)
