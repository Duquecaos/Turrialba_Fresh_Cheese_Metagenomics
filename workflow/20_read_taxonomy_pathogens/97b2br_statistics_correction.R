#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(width = 220)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
    stop("Usage: 97b2br_statistics_correction.R <ROOT>")
}

ROOT <- args[1]

IN <- file.path(
    ROOT,
    "98_read_taxonomy",
    "97B2B_bacterial_community"
)

OUT <- file.path(
    ROOT,
    "98_read_taxonomy",
    "97B2BR_statistics_correction"
)

dir.create(
    OUT,
    recursive = TRUE,
    showWarnings = FALSE
)

suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(permute))

set.seed(20260929)


write_tsv <- function(x, filename) {

    write.table(
        x,
        file = file.path(OUT, filename),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
    )
}


read_matrix <- function(filename) {

    x <- read.delim(
        file.path(IN, filename),
        check.names = FALSE,
        quote = "",
        comment.char = ""
    )

    taxa <- x[, c("taxonomy_id", "name")]

    mat <- as.matrix(
        x[, setdiff(
            colnames(x),
            c("taxonomy_id", "name")
        )]
    )

    storage.mode(mat) <- "numeric"

    rownames(mat) <- make.unique(
        paste(
            taxa$taxonomy_id,
            taxa$name,
            sep = "|"
        )
    )

    mat
}


anova_table <- function(
    fit,
    analysis
) {

    d <- as.data.frame(fit)

    d$term <- rownames(d)
    rownames(d) <- NULL

    d$analysis <- analysis

    d[, c(
        "analysis",
        "term",
        setdiff(
            names(d),
            c("analysis", "term")
        )
    )]
}


extract_term <- function(
    fit,
    term
) {

    d <- as.data.frame(fit)

    rn <- rownames(d)

    idx <- which(
        rn == term
    )

    if (length(idx) != 1) {
        stop(
            paste(
                "Could not uniquely find term",
                term,
                "in:",
                paste(rn, collapse = ",")
            )
        )
    }

    pcol <- grep(
        "^Pr",
        names(d),
        value = TRUE
    )

    if (length(pcol) != 1) {
        stop("Could not identify P-value column")
    }

    data.frame(
        term = term,
        Df = d$Df[idx],
        SumOfSqs = d$SumOfSqs[idx],
        R2 = d$R2[idx],
        F = d$F[idx],
        P = d[[pcol]][idx]
    )
}


# ============================================================
# Metadata
# ============================================================

metadata <- read.delim(
    file.path(
        IN,
        "97B2B_sample_metadata.tsv"
    ),
    check.names = FALSE
)

metadata$producer_code <- factor(
    metadata$producer_code,
    levels = c("L", "M")
)

metadata$time_code <- factor(
    metadata$time_code,
    levels = c("1", "2", "3")
)

metadata$subject_id <- factor(
    metadata$subject_id
)

samples <- as.character(
    metadata$sample
)


# ============================================================
# Matrices
# ============================================================

matrices <- list(
    unfiltered = read_matrix(
        "97B2B_bacterial_genus_matrix_unfiltered.tsv"
    ),

    primary = read_matrix(
        "97B2B_bacterial_genus_matrix_primary.tsv"
    ),

    strict = read_matrix(
        "97B2B_bacterial_genus_matrix_strict.tsv"
    )
)


for (nm in names(matrices)) {

    M <- matrices[[nm]]

    if (!setequal(
        colnames(M),
        samples
    )) {
        stop(
            paste(
                "Sample mismatch:",
                nm
            )
        )
    }

    matrices[[nm]] <- M[, samples, drop = FALSE]
}


# ============================================================
# Primary model:
# subject + time
# restricted permutations within subject
# ============================================================

PRIMARY <- matrices[["primary"]]

X <- t(PRIMARY)

bray <- vegdist(
    X,
    method = "bray"
)


control_time <- how(
    nperm = 9999,
    blocks = metadata$subject_id
)


fit_time_adjusted <- adonis2(
    bray ~ subject_id + time_code,
    data = metadata,
    permutations = control_time,
    by = "terms"
)


write_tsv(
    anova_table(
        fit_time_adjusted,
        "primary_subject_adjusted_time"
    ),
    "97B2BR_PERMANOVA_subject_adjusted_time.tsv"
)


time_primary <- extract_term(
    fit_time_adjusted,
    "time_code"
)


# ============================================================
# Time within each producer
# subject + time
# ============================================================

within_results <- list()


for (producer in c("L", "M")) {

    idx <- metadata$producer_code == producer

    meta_sub <- droplevels(
        metadata[idx, , drop = FALSE]
    )

    Xsub <- X[
        idx,
        ,
        drop = FALSE
    ]

    dsub <- vegdist(
        Xsub,
        method = "bray"
    )

    control_sub <- how(
        nperm = 9999,
        blocks = meta_sub$subject_id
    )

    fit_sub <- adonis2(
        dsub ~ subject_id + time_code,
        data = meta_sub,
        permutations = control_sub,
        by = "terms"
    )

    tab <- anova_table(
        fit_sub,
        paste0(
            "time_within_producer_",
            producer,
            "_subject_adjusted"
        )
    )

    within_results[[producer]] <- tab
}


within_table <- do.call(
    rbind,
    within_results
)

write_tsv(
    within_table,
    "97B2BR_PERMANOVA_time_within_producer.tsv"
)


# ============================================================
# Producer: six independent biological units
# Mean composition over their three repeated samples
# ============================================================

subjects <- levels(
    metadata$subject_id
)


subject_profiles <- do.call(
    rbind,
    lapply(
        subjects,
        function(subject) {

            idx <- metadata$subject_id == subject

            colMeans(
                X[idx, , drop = FALSE]
            )
        }
    )
)


rownames(subject_profiles) <- subjects

subject_profiles <- subject_profiles /
    rowSums(subject_profiles)


subject_metadata <- unique(
    metadata[, c(
        "subject_id",
        "producer_code"
    )]
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


fit_producer <- adonis2(
    subject_bray ~ producer_code,
    data = subject_metadata,
    permutations = 9999,
    by = "terms"
)


write_tsv(
    anova_table(
        fit_producer,
        "producer_subject_level"
    ),
    "97B2BR_PERMANOVA_producer_subject_level.tsv"
)


producer_primary <- extract_term(
    fit_producer,
    "producer_code"
)


# ============================================================
# Producer PERMDISP
# Lingoes correction avoids the previous negative squared
# distance warning.
# ============================================================

disp <- betadisper(
    subject_bray,
    group = subject_metadata$producer_code,
    type = "median",
    bias.adjust = TRUE,
    add = "lingoes"
)


disp_perm <- permutest(
    disp,
    permutations = 9999
)


capture.output(
    disp_perm,
    file = file.path(
        OUT,
        "97B2BR_PERMDISP_producer.txt"
    )
)


disp_df <- data.frame(
    subject_id = names(
        disp$distances
    ),

    distance_to_group_median = as.numeric(
        disp$distances
    )
)


disp_df <- merge(
    subject_metadata,
    disp_df,
    by = "subject_id",
    sort = FALSE
)


write_tsv(
    disp_df,
    "97B2BR_producer_dispersion_by_subject.tsv"
)


disp_summary <- aggregate(
    distance_to_group_median
    ~ producer_code,
    data = disp_df,
    FUN = function(x) {
        c(
            mean = mean(x),
            median = median(x),
            sd = sd(x)
        )
    }
)


disp_summary <- data.frame(
    producer_code =
        disp_summary$producer_code,

    mean_distance =
        disp_summary$
            distance_to_group_median[, "mean"],

    median_distance =
        disp_summary$
            distance_to_group_median[, "median"],

    sd_distance =
        disp_summary$
            distance_to_group_median[, "sd"]
)


write_tsv(
    disp_summary,
    "97B2BR_producer_dispersion_summary.tsv"
)


# ============================================================
# Correct filter sensitivity
# ============================================================

sensitivity <- list()


for (filter_name in names(matrices)) {

    MF <- matrices[[filter_name]]

    XF <- t(MF)

    DF <- vegdist(
        XF,
        method = "bray"
    )


    # Time adjusted for subject
    controlF <- how(
        nperm = 9999,
        blocks = metadata$subject_id
    )

    fit_timeF <- adonis2(
        DF ~ subject_id + time_code,
        data = metadata,
        permutations = controlF,
        by = "terms"
    )

    tr <- extract_term(
        fit_timeF,
        "time_code"
    )


    # Producer based on six independent subject means
    subjectF <- do.call(
        rbind,
        lapply(
            subjects,
            function(subject) {

                idx <- metadata$subject_id == subject

                colMeans(
                    XF[idx, , drop = FALSE]
                )
            }
        )
    )


    subjectF <- subjectF /
        rowSums(subjectF)


    producer_distF <- vegdist(
        subjectF,
        method = "bray"
    )


    fit_prodF <- adonis2(
        producer_distF ~ producer_code,
        data = subject_metadata,
        permutations = 9999,
        by = "terms"
    )


    pr <- extract_term(
        fit_prodF,
        "producer_code"
    )


    sensitivity[[filter_name]] <- data.frame(
        filter = filter_name,
        genera = nrow(MF),

        time_subject_adjusted_R2 = tr$R2,
        time_subject_adjusted_F = tr$F,
        time_subject_adjusted_P = tr$P,

        producer_subject_R2 = pr$R2,
        producer_subject_F = pr$F,
        producer_subject_P = pr$P
    )
}


sensitivity_table <- do.call(
    rbind,
    sensitivity
)


write_tsv(
    sensitivity_table,
    "97B2BR_filter_sensitivity_corrected.tsv"
)


# ============================================================
# Alpha diversity: repeated time test
# ============================================================

alpha <- read.delim(
    file.path(
        IN,
        "97B2B_alpha_diversity_primary.tsv"
    ),
    check.names = FALSE
)


alpha$time_code <- factor(
    alpha$time_code,
    levels = c("1", "2", "3")
)

alpha$subject_id <- factor(
    alpha$subject_id
)

alpha$producer_code <- factor(
    alpha$producer_code,
    levels = c("L", "M")
)


alpha_metrics <- c(
    "Shannon",
    "Simpson",
    "inverse_Simpson",
    "effective_Shannon_genera"
)


friedman_results <- list()


for (metric in alpha_metrics) {

    form <- as.formula(
        paste(
            metric,
            "~ time_code | subject_id"
        )
    )


    ft <- friedman.test(
        form,
        data = alpha
    )


    friedman_results[[paste0(
        "ALL_",
        metric
    )]] <- data.frame(
        population = "ALL",
        metric = metric,
        n_subjects = length(
            unique(
                alpha$subject_id
            )
        ),
        statistic = unname(
            ft$statistic
        ),
        df = unname(
            ft$parameter
        ),
        P = ft$p.value
    )


    for (producer in c("L", "M")) {

        sub <- droplevels(
            alpha[
                alpha$producer_code
                == producer,
                ,
                drop = FALSE
            ]
        )


        ft_sub <- friedman.test(
            form,
            data = sub
        )


        friedman_results[[paste0(
            producer,
            "_",
            metric
        )]] <- data.frame(
            population = producer,
            metric = metric,
            n_subjects = length(
                unique(
                    sub$subject_id
                )
            ),
            statistic = unname(
                ft_sub$statistic
            ),
            df = unname(
                ft_sub$parameter
            ),
            P = ft_sub$p.value
        )
    }
}


friedman_table <- do.call(
    rbind,
    friedman_results
)


write_tsv(
    friedman_table,
    "97B2BR_alpha_Friedman_time.tsv"
)


# ============================================================
# PCoA axis variance
# ============================================================

eig <- read.delim(
    file.path(
        IN,
        "97B2B_PCoA_primary_eigenvalues.tsv"
    ),
    check.names = FALSE
)


positive <- eig$eigenvalue[
    eig$eigenvalue > 0
]


axis_summary <- data.frame(
    axis = c(
        "PCoA1",
        "PCoA2"
    ),

    pct_positive_eigenvalue_variation = c(
        100 * eig$eigenvalue[1] / sum(positive),
        100 * eig$eigenvalue[2] / sum(positive)
    )
)


write_tsv(
    axis_summary,
    "97B2BR_PCoA_axis_variance.tsv"
)


# ============================================================
# Thesis-safe summary
# ============================================================

within_L <- extract_term(
    adonis2(
        vegdist(
            X[
                metadata$producer_code == "L",
                ,
                drop = FALSE
            ],
            method = "bray"
        ) ~ subject_id + time_code,
        data = droplevels(
            metadata[
                metadata$producer_code == "L",
                ,
                drop = FALSE
            ]
        ),
        permutations = how(
            nperm = 9999,
            blocks = droplevels(
                metadata[
                    metadata$producer_code == "L",
                    "subject_id"
                ]
            )
        ),
        by = "terms"
    ),
    "time_code"
)


within_M <- extract_term(
    adonis2(
        vegdist(
            X[
                metadata$producer_code == "M",
                ,
                drop = FALSE
            ],
            method = "bray"
        ) ~ subject_id + time_code,
        data = droplevels(
            metadata[
                metadata$producer_code == "M",
                ,
                drop = FALSE
            ]
        ),
        permutations = how(
            nperm = 9999,
            blocks = droplevels(
                metadata[
                    metadata$producer_code == "M",
                    "subject_id"
                ]
            )
        ),
        by = "terms"
    ),
    "time_code"
)


summary <- data.frame(
    metric = c(
        "primary_genera",
        "time_subject_adjusted_R2",
        "time_subject_adjusted_F",
        "time_subject_adjusted_P",
        "time_within_L_R2",
        "time_within_L_P",
        "time_within_M_R2",
        "time_within_M_P",
        "producer_subject_level_R2",
        "producer_subject_level_F",
        "producer_subject_level_P",
        "producer_independent_units",
        "producer_label_combinations",
        "producer_dispersion_correction",
        "producer_dispersion_interpretation_required",
        "primary_filter",
        "amplicon_data_used"
    ),

    value = c(
        nrow(PRIMARY),

        sprintf(
            "%.8f",
            time_primary$R2
        ),

        sprintf(
            "%.8f",
            time_primary$F
        ),

        sprintf(
            "%.8g",
            time_primary$P
        ),

        sprintf(
            "%.8f",
            within_L$R2
        ),

        sprintf(
            "%.8g",
            within_L$P
        ),

        sprintf(
            "%.8f",
            within_M$R2
        ),

        sprintf(
            "%.8g",
            within_M$P
        ),

        sprintf(
            "%.8f",
            producer_primary$R2
        ),

        sprintf(
            "%.8f",
            producer_primary$F
        ),

        sprintf(
            "%.8g",
            producer_primary$P
        ),

        "6_total_3_per_producer",

        choose(6, 3),

        "Lingoes_plus_bias_adjusted_betadisper",

        "YES_before_centroid_conclusion",

        "genus_ge0.01pct_within_bacteria_in_ge2_samples",

        "NO"
    )
)


write_tsv(
    summary,
    "97B2BR_global_summary.tsv"
)


# ============================================================
# Scope
# ============================================================

scope <- data.frame(
    field = c(
        "reason_for_correction",
        "original_job_132974",
        "time_model",
        "time_permutations",
        "producer_model",
        "producer_limitation",
        "PERMDISP",
        "filter_sensitivity",
        "alpha_time_test",
        "formal_interaction",
        "interpretation"
    ),

    value = c(
        "subject_adjusted_repeated_measure_model_and_fix_NA_filter_sensitivity",

        "technically_complete_but_superseded_for_inferential_statistics_by_97B2BR",

        "BrayCurtis_subject_id_plus_time_code",

        "restricted_within_biological_unit",

        "mean_profile_per_6_independent_biological_units_then_producer_PERMANOVA",

        "n3_per_producer_coarse_permutation_resolution",

        "Lingoes_corrected_bias_adjusted",

        "unfiltered_primary_strict_recomputed_without_NA_extraction_bug",

        "Friedman_repeated_measure",

        "not_tested",

        "Bracken_read_composition_not_absolute_cell_abundance"
    )
)


write_tsv(
    scope,
    "97B2BR_methodological_scope.tsv"
)


# ============================================================
# Validation
# ============================================================

required <- c(
    "97B2BR_PERMANOVA_subject_adjusted_time.tsv",
    "97B2BR_PERMANOVA_time_within_producer.tsv",
    "97B2BR_PERMANOVA_producer_subject_level.tsv",
    "97B2BR_PERMDISP_producer.txt",
    "97B2BR_producer_dispersion_by_subject.tsv",
    "97B2BR_producer_dispersion_summary.tsv",
    "97B2BR_filter_sensitivity_corrected.tsv",
    "97B2BR_alpha_Friedman_time.tsv",
    "97B2BR_PCoA_axis_variance.tsv",
    "97B2BR_global_summary.tsv",
    "97B2BR_methodological_scope.tsv"
)


for (f in required) {

    path <- file.path(
        OUT,
        f
    )

    if (
        !file.exists(path)
        ||
        file.info(path)$size <= 0
    ) {

        stop(
            paste(
                "Missing/empty:",
                f
            )
        )
    }
}


cat("97B2BR=PASS\n")
cat(
    "PRIMARY_GENERA=",
    nrow(PRIMARY),
    "\n",
    sep = ""
)
cat(
    "TIME_ADJUSTED_R2=",
    time_primary$R2,
    "\n",
    sep = ""
)
cat(
    "TIME_ADJUSTED_P=",
    time_primary$P,
    "\n",
    sep = ""
)
cat(
    "PRODUCER_R2=",
    producer_primary$R2,
    "\n",
    sep = ""
)
cat(
    "PRODUCER_P=",
    producer_primary$P,
    "\n",
    sep = ""
)
