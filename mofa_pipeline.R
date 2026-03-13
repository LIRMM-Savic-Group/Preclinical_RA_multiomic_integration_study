# =============================================================================
# MOFA2 pipeline: multi-omics factor analysis
#
# Design:
#   - Primary model:   Bernoulli likelihood for somatic CHIP variants;
#                      Gaussian for all other views
#   - Companion model: all-Gaussian (used solely for variance-explained metrics,
#                      as R² is not defined under Bernoulli)
#
# Sections:
#   0.  Paths and analysis settings
#   1.  Metadata and cohort definition
#   2.  Utility functions
#   3.  Load modality matrices
#   4.  Align to cohort samples and QC
#   5.  Train MOFA2 (Bernoulli + all-Gaussian companion)
#   6.  Convergence QC and variance explained
#   7.  Cox models: factor associations with time-to-progression
#   8.  Kaplan-Meier curves (Factor 6 tertiles and median split)
#   9.  Factor 6 composition: per-view variance + top loadings
#   10. Factor 6 view-level visualisations
#   11. Palindromic-status logistic regression (all factors)
#   12. Per-factor variance explained (Factors 7, 8, 10)
#   13. Factor 7 analysis (RNA, IFN score, cytokines)
#   14. Factor 8 analysis (RNA, cytokines, HLA)
#   15. Factor 10 analysis (RNA, cytokines, IFN score, rare germline)
#
# Usage:
#   Set `base_dir` below to the project folder containing all input files,
#   then source the script or run section by section.
# =============================================================================

# -----------------------------------------------------------------------------
# NOTE: set this to your project directory before running
# -----------------------------------------------------------------------------
base_dir <- "."   # e.g. "/path/to/RA_cohort_analysis/mofa_model_v3"
setwd(base_dir)

# =============================================================================
# Libraries  (loaded once at the top)
# =============================================================================
suppressPackageStartupMessages({
  library(MOFA2)
  library(basilisk)
  library(readxl)
  library(writexl)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(survival)
  library(survminer)
  library(broom)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(msigdbr)
  library(clusterProfiler)
})

set.seed(1)

# =============================================================================
# 0. Paths and analysis settings
# =============================================================================

# -- Input files --------------------------------------------------------------
meta_file <- "mofa_cohort_matrix_v1.0.xlsx"
rna_file  <- "mofa_rnaseq2.csv"
germ_file <- "mofa_germline.csv"
prot_file <- "mofa_proteomics.csv"
asc_file  <- "mofa_asc.csv"
ifn_file  <- "mofa_ifn.csv"
chip_file <- "mofa_chip.csv"
hla_file  <- "mofa_hla.csv"

# -- Output prefixes ----------------------------------------------------------
out_prefix <- "MOFA_seed1"
out_hdf5   <- paste0(out_prefix, "_bernoulliCHIP.hdf5")
out_hdf5_g <- paste0(out_prefix, "_allGaussian_companion.hdf5")

# -- Analysis variables -------------------------------------------------------
time_var  <- "time_to_prog_years"
event_var <- "event"
id_var    <- "sample_id"
chip_view <- "Somatic CHIP variants"

# -- Colour palette (used throughout) ----------------------------------------
COL_RED  <- "red"
COL_BLUE <- "blue"

# =============================================================================
# 1. Metadata and cohort definition
# =============================================================================
meta <- readxl::read_xlsx(meta_file) %>%
  mutate(
    sample_id    = as.character(sample_id),
    ccp_2or3_bin = suppressWarnings(as.integer(ccp_2or3_bin))
  )

stopifnot(all(c("sample_id", "ccp_2or3_bin") %in% colnames(meta)))

meta_mofa <- meta %>%
  filter(ccp_2or3_bin == 1) %>%
  filter(!is.na(sample_id), sample_id != "") %>%
  distinct(sample_id, .keep_all = TRUE)

samples <- meta_mofa$sample_id
cat("MOFA cohort samples:", length(samples), "\n")

# =============================================================================
# 2. Utility functions
# =============================================================================

# -- Matrix QC and alignment --------------------------------------------------

#' Remove features with all-NA, non-finite, or zero-variance values
drop_bad_features <- function(mat, view_name = "view") {
  if (!is.matrix(mat)) mat <- as.matrix(mat)

  all_na      <- apply(mat, 1, function(x) all(is.na(x)))
  all_nonfinite <- apply(mat, 1, function(x) {
    x <- x[!is.na(x)]
    length(x) == 0 || all(!is.finite(x))
  })
  zero_var    <- apply(mat, 1, function(x) {
    x <- x[is.finite(x)]
    length(x) < 2 || isTRUE(all.equal(stats::sd(x), 0))
  })

  keep <- !(all_na | all_nonfinite | zero_var)
  out  <- mat[keep, , drop = FALSE]
  if (nrow(out) == 0) warning("After filtering, view '", view_name, "' has 0 features.")
  out
}

#' Align a features-x-samples matrix to a target sample vector (adds NA columns for missing samples)
reindex_samples <- function(mat, samples, view_name = "view") {
  if (is.null(mat))         stop("reindex_samples: ", view_name, " is NULL.")
  if (is.data.frame(mat))   mat <- as.matrix(mat)
  if (!is.matrix(mat))      stop("reindex_samples: ", view_name, " is not a matrix/data.frame.")
  if (is.null(colnames(mat))) stop("reindex_samples: ", view_name, " has no colnames (sample IDs).")

  missing  <- setdiff(samples, colnames(mat))
  if (length(missing) > 0) {
    na_block <- matrix(NA_real_, nrow = nrow(mat), ncol = length(missing),
                       dimnames = list(rownames(mat), missing))
    mat <- cbind(mat, na_block)
  }
  mat[, samples, drop = FALSE]
}

#' Read a modality CSV into a features-x-samples matrix.
#'   Supports two layouts:
#'     Layout A – rows = samples (detected by a `sample_id` column)
#'     Layout B – rows = features (first matching column used as feature ID)
read_mofa_csv <- function(path,
                          sample_col = "sample_id",
                          feature_col_candidates = c("feature_id", "feature", "gene",
                                                     "SYMBOL", "id", "ID"),
                          prefix    = NULL,
                          view_name = "view") {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  df <- df[, !is.na(colnames(df)) & colnames(df) != "", drop = FALSE]

  if (sample_col %in% colnames(df)) {
    # Layout A: samples x features  →  transpose to features x samples
    df[[sample_col]] <- as.character(df[[sample_col]])
    feat_cols <- setdiff(colnames(df), sample_col)
    df_num    <- df[, feat_cols, drop = FALSE]
    for (nm in colnames(df_num)) df_num[[nm]] <- suppressWarnings(as.numeric(df_num[[nm]]))
    mat_sf <- as.matrix(df_num)
    rownames(mat_sf) <- df[[sample_col]]
    mat <- t(mat_sf)
  } else {
    # Layout B: features x samples
    feature_col <- intersect(feature_col_candidates, colnames(df))
    feature_col <- if (length(feature_col) == 0) colnames(df)[1] else feature_col[1]
    df[[feature_col]] <- as.character(df[[feature_col]])
    sample_cols <- setdiff(colnames(df), feature_col)
    df_num <- df[, sample_cols, drop = FALSE]
    for (nm in colnames(df_num)) df_num[[nm]] <- suppressWarnings(as.numeric(df_num[[nm]]))
    mat <- as.matrix(df_num)
    rownames(mat) <- df[[feature_col]]
  }

  if (!is.null(prefix)) rownames(mat) <- paste0(prefix, rownames(mat))

  if (is.null(colnames(mat)) || any(colnames(mat) == ""))
    stop("No sample colnames detected in: ", path)
  if (is.null(rownames(mat)) || any(rownames(mat) == ""))
    stop("No feature rownames detected in: ", path)
  if (any(!is.finite(mat) & !is.na(mat))) {
    bad <- which(!is.finite(mat) & !is.na(mat), arr.ind = TRUE)
    stop("Non-finite values in view '", view_name, "' (", path, "). ",
         "Example: feature=", rownames(mat)[bad[1, 1]],
         " sample=", colnames(mat)[bad[1, 2]],
         " value=", mat[bad[1, 1], bad[1, 2]])
  }
  mat
}

#' Print a brief QC summary for a features-x-samples matrix
qc_view <- function(mat, name) {
  cat(sprintf("%-30s features=%-5d samples=%-5d NA=%-6d nonfinite(non-NA)=%d\n",
              name, nrow(mat), ncol(mat),
              sum(is.na(mat)),
              sum(!is.finite(mat) & !is.na(mat))))
}

#' Extract the features-x-samples data matrix for a single view from a trained MOFA object
get_view_matrix <- function(mofa_obj, view_name) {
  X_view <- MOFA2::get_data(mofa_obj, views = view_name)[[view_name]]
  X <- X_view[[1]]
  if (is.data.frame(X)) X <- as.matrix(X)
  stopifnot(is.matrix(X))
  X
}

# -- Shared plot helpers -------------------------------------------------------

#' Join palindromic status (0/1) from metadata and factorise for plotting
add_palindromic <- function(df_in, df_meta = df, id_col = "sample_id") {
  df_in %>%
    left_join(df_meta %>% select(sample_id, palindromic),
              by = setNames("sample_id", id_col)) %>%
    mutate(palindromic = factor(palindromic, levels = c(0, 1),
                                labels = c("Non-palindromic", "Palindromic"))) %>%
    filter(!is.na(palindromic))
}

#' Compute per-facet Spearman labels (rho + p) for annotating faceted scatter plots
facet_spearman_labels <- function(df_in, xvar, yvar, facet_var = "palindromic") {
  df_in %>%
    group_by(.data[[facet_var]]) %>%
    summarise(
      rho = unname(cor.test(.data[[xvar]], .data[[yvar]],
                            method = "spearman", exact = FALSE)$estimate),
      p   = cor.test(.data[[xvar]], .data[[yvar]],
                     method = "spearman", exact = FALSE)$p.value,
      .groups = "drop"
    ) %>%
    mutate(label = paste0("Spearman \u03c1 = ", round(rho, 2),
                          "\np = ", format.pval(p, digits = 2, eps = 1e-3)))
}

# -- Feature label formatters --------------------------------------------------

#' Format HLA feature names for display (e.g. "HLA_drb1_DRB1_04_01_carrier" → "DRB1*04:01 (carrier)")
format_hla_feature <- function(x) {
  ifelse(
    x %in% c("HLA_se_carrier", "HLA_se_dosage"),
    ifelse(x == "HLA_se_carrier", "Shared epitope (carrier)", "Shared epitope (dosage)"),
    x
  ) |>
    str_remove("^HLA_drb1_") |>
    str_replace("^([A-Za-z0-9]+)_([0-9]{2})_([0-9]{2})_(carrier|dosage)$", "\\1*\\2:\\3 (\\4)") |>
    str_replace("^([A-Za-z0-9]+)_([0-9]{2})_([0-9]{2})$", "\\1*\\2:\\3")
}

#' Format germline feature names for display (e.g. "GERM_germ_RA_n_damaging_variants" → "Number of damaging variants (RA panel)")
format_germ_feature <- function(x) {
  x2 <- x %>% str_replace_all("_+", " ") %>% str_squish()

  panel <- case_when(
    str_detect(x2, "\\bgerm\\b.*\\bRA\\b")   ~ "RA panel",
    str_detect(x2, "\\bgerm\\b.*\\bDISC\\b") ~ "Discovery panel",
    TRUE ~ NA_character_
  )
  measure <- case_when(
    str_detect(x2, "\\bn\\b.*\\bdamaging\\b.*\\bvariants\\b") ~ "Number of damaging variants",
    str_detect(x2, "\\bn\\b.*\\bdamaging\\b.*\\bgenes\\b")    ~ "Number of damaging genes",
    TRUE ~ x2
  )
  ifelse(!is.na(panel), paste0(measure, " (", panel, ")"), measure)
}

#' Strip standard prefixes from cytokine/proteomics feature names for display
format_cyto_feature <- function(x) {
  x %>%
    str_remove("^PROT[_ ]") %>%
    str_remove("^Protein[_ ]") %>%
    str_remove("^CYTO[_ ]") %>%
    str_replace_all("_", " ") %>%
    str_squish()
}

# -- Hallmark GSEA helper ------------------------------------------------------

#' Build a ranked gene list from MOFA RNA loadings and run Hallmark GSEA
run_hallmark_gsea <- function(rna_loadings, seed = 1) {
  gene_ranks <- rna_loadings %>%
    mutate(gene = toupper(gene)) %>%
    group_by(gene) %>%
    summarise(loading = mean(loading), .groups = "drop") %>%
    arrange(desc(loading)) %>%
    { setNames(.$loading, .$gene) }
  gene_ranks <- sort(gene_ranks, decreasing = TRUE)

  hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
    select(gs_name, gene_symbol) %>%
    distinct()

  set.seed(seed)
  gsea_result <- clusterProfiler::GSEA(
    geneList     = gene_ranks,
    TERM2GENE    = hallmark,
    pvalueCutoff = 0.05,
    verbose      = FALSE
  )
  gsea_result
}

#' Clean Hallmark GSEA result descriptions and filter to adjusted p < 0.05
tidy_gsea_df <- function(gsea_result) {
  as.data.frame(gsea_result@result) %>%
    mutate(
      Description = Description %>%
        str_remove("^HALLMARK_") %>%
        str_replace_all("_", " ") %>%
        str_to_sentence() %>%
        str_replace_all("\\bmtorc1\\b", "mTORC1") %>%
        str_replace_all("\\bil6\\b",    "IL6") %>%
        str_replace_all("\\bjak\\b",    "JAK") %>%
        str_replace_all("\\bstat3\\b",  "STAT3") %>%
        str_replace_all("\\bnfkb\\b",   "NFKB") %>%
        str_replace_all("\\btnfa\\b",   "TNF"),
      direction = factor(
        ifelse(NES > 0, "Positive enrichment", "Negative enrichment"),
        levels = c("Positive enrichment", "Negative enrichment")
      )
    ) %>%
    filter(p.adjust < 0.05)
}

#' Plot top Hallmark pathways per direction as a bubble plot
plot_gsea_bubbles <- function(gsea_df, factor_label, rna_view, n_top = 10,
                              col_pos = COL_RED, col_neg = COL_BLUE) {
  plot_df <- gsea_df %>%
    group_by(direction) %>%
    arrange(desc(abs(NES))) %>%
    slice_head(n = n_top) %>%
    ungroup()

  ggplot(plot_df, aes(x = NES, y = reorder(Description, NES),
                      size = -log10(p.adjust), colour = direction)) +
    geom_point(alpha = 0.9) +
    facet_wrap(~ direction, scales = "free_y") +
    scale_colour_manual(
      values = c("Positive enrichment" = col_pos,
                 "Negative enrichment" = col_neg),
      guide = "none"
    ) +
    labs(
      x = "Normalized enrichment score (NES)", y = NULL,
      size = expression(-log[10]("adj. p")),
      title    = paste0("Hallmark pathway enrichment \u2013 ", factor_label),
      subtitle = paste0("RNA view: ", rna_view, " (ranked by MOFA weights)")
    ) +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
}

# -- Variance explained helper -------------------------------------------------

#' Extract per-view R² for a single factor and return as a tidy data frame
factor_variance_df <- function(var_exp, factor_name) {
  r2_pf <- var_exp$r2_per_factor[[1]]
  as.data.frame(t(r2_pf[factor_name, , drop = FALSE])) %>%
    rownames_to_column("view") %>%
    setNames(c("view", "r2")) %>%
    arrange(desc(r2))
}

#' Bar chart of per-view variance explained for a single factor
plot_factor_variance <- function(r2_df, factor_label) {
  ggplot(r2_df, aes(x = reorder(view, r2), y = r2)) +
    geom_col(fill = "grey20") +
    coord_flip() +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = "Variance explained (R\u00b2)",
         title = paste0("Variance explained by ", factor_label))
}

# -- Factor 8 HLA carrier helper -----------------------------------------------

#' Faceted box-jitter plot of a factor score by carrier status, stratified by palindromic subgroup.
#'   Per-facet Wilcoxon p annotated in the bottom-right corner.
plot_carrier_facet <- function(dat, factor_col, carrier_col, title_text,
                               col_palette = c("Non-carrier" = COL_BLUE,
                                               "Carrier"     = COL_RED)) {
  dat2 <- dat %>%
    transmute(
      sample_id,
      Factor  = .data[[factor_col]],
      palindromic,
      carrier = .data[[carrier_col]]
    ) %>%
    filter(!is.na(carrier)) %>%
    mutate(carrier_status = factor(ifelse(carrier > 0, "Carrier", "Non-carrier"),
                                   levels = c("Non-carrier", "Carrier")))

  stat_df <- dat2 %>%
    group_by(palindromic) %>%
    summarise(
      p = wilcox.test(Factor ~ carrier_status, exact = FALSE)$p.value,
      .groups = "drop"
    ) %>%
    mutate(label = paste0("Wilcoxon p = ", format.pval(p, digits = 2, eps = 1e-3)))

  ggplot(dat2, aes(x = carrier_status, y = Factor)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, fill = "grey90", colour = "black") +
    geom_jitter(aes(colour = carrier_status), width = 0.15, alpha = 0.75, size = 1.8) +
    geom_text(data = stat_df,
              aes(x = Inf, y = -Inf, label = label),
              inherit.aes = FALSE,
              hjust = 1.05, vjust = -0.6, size = 4) +
    facet_wrap(~ palindromic, nrow = 1) +
    scale_colour_manual(values = col_palette, guide = "none") +
    theme_bw(base_size = 12) +
    theme(strip.text   = element_text(face = "bold", size = 13),
          axis.text.x  = element_text(size = 14)) +
    labs(x = NULL, y = paste0("MOFA ", factor_col, " score"), title = title_text)
}

# -- Cytokine factor plot helper -----------------------------------------------

#' Plot top cytokine loadings + module score scatter for a given factor;
#'   returns a named list of ggplot objects.
plot_factor_cytokines <- function(mofa_obj, cyto_mat, df_meta,
                                  cyto_view, factor_name, n = 10) {
  W   <- MOFA2::get_weights(mofa_obj, views = cyto_view, factors = factor_name)[[cyto_view]]
  w   <- as.numeric(W[, factor_name, drop = TRUE])
  names(w) <- rownames(W)

  loadings <- tibble(feature = names(w), loading = w) %>%
    arrange(desc(abs(loading))) %>%
    slice_head(n = n) %>%
    mutate(label = format_cyto_feature(feature),
           sign  = ifelse(loading >= 0, "Positive", "Negative"),
           label = fct_reorder(label, loading))

  p_load <- ggplot(loadings, aes(x = label, y = loading, fill = sign)) +
    geom_col(width = 0.75) +
    coord_flip() +
    scale_fill_manual(values = c("Positive" = COL_RED, "Negative" = COL_BLUE), guide = "none") +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = paste0("MOFA weight (", factor_name, ")"),
         title = paste0("Cytokine loadings \u2013 ", factor_name))

  list(loadings = p_load)
}

# =============================================================================
# 3. Load modality matrices (features x samples)
# =============================================================================
rna_mat  <- read_mofa_csv(rna_file,  sample_col = "sample_id", prefix = "RNA_",  view_name = "Gene expression")
germ_mat <- read_mofa_csv(germ_file, sample_col = "sample_id", prefix = "GERM_", view_name = "Rare Germline variants")
prot_mat <- read_mofa_csv(prot_file, sample_col = "sample_id", prefix = "PROT_", view_name = "Inflammatory cytokines")
asc_mat  <- read_mofa_csv(asc_file,  sample_col = "sample_id", prefix = "ASC_",  view_name = "Inflammasome ASC specks")
ifn_mat  <- read_mofa_csv(ifn_file,  sample_col = "sample_id", prefix = "IFN_",  view_name = "Interferon score")
som_mat  <- read_mofa_csv(chip_file, sample_col = "sample_id", prefix = "SOM_",  view_name = "Somatic CHIP variants")
hla_mat  <- read_mofa_csv(hla_file,  sample_col = "sample_id", prefix = "HLA_",  view_name = "HLA genotype")

# =============================================================================
# 4. Align to cohort samples and QC
# =============================================================================
data_list <- list(
  "Gene expression"         = drop_bad_features(reindex_samples(rna_mat,  samples, "Gene expression"),         "Gene expression"),
  "Rare Germline variants"  = drop_bad_features(reindex_samples(germ_mat, samples, "Rare Germline variants"),  "Rare Germline variants"),
  "Inflammatory cytokines"  = drop_bad_features(reindex_samples(prot_mat, samples, "Inflammatory cytokines"),  "Inflammatory cytokines"),
  "Inflammasome ASC specks" = drop_bad_features(reindex_samples(asc_mat,  samples, "Inflammasome ASC specks"), "Inflammasome ASC specks"),
  "Interferon score"        = drop_bad_features(reindex_samples(ifn_mat,  samples, "Interferon score"),        "Interferon score"),
  "Somatic CHIP variants"   = drop_bad_features(reindex_samples(som_mat,  samples, "Somatic CHIP variants"),   "Somatic CHIP variants"),
  "HLA genotype"            = drop_bad_features(reindex_samples(hla_mat,  samples, "HLA genotype"),            "HLA genotype")
)

cat("\n--- QC: per-view (after sample alignment and feature filtering) ---\n")
for (nm in names(data_list)) qc_view(data_list[[nm]], nm)

empty_views <- names(data_list)[vapply(data_list, function(x) nrow(x) == 0, logical(1))]
if (length(empty_views) > 0) {
  warning("Dropping empty views: ", paste(empty_views, collapse = ", "))
  data_list <- data_list[setdiff(names(data_list), empty_views)]
}

# Validate CHIP view: must be strictly 0/1 (Bernoulli-compatible)
stopifnot(chip_view %in% names(data_list))
x <- data_list[[chip_view]]
stopifnot(all(x[!is.na(x)] %in% c(0, 1)))
storage.mode(x) <- "numeric"
data_list[[chip_view]] <- x

# =============================================================================
# 5. Train MOFA2
# =============================================================================

# -- Sample metadata (MOFA2 expects a column named 'sample') ------------------
sample_md <- meta_mofa %>%
  mutate(sample = as.character(sample_id)) %>%
  select(sample, everything()) %>%
  as.data.frame()

# -- 5A. Primary model: Bernoulli for CHIP, Gaussian for all other views ------
mofa_b       <- create_mofa(data_list)
samples_metadata(mofa_b) <- sample_md

model_opts_b <- get_default_model_options(mofa_b)
model_opts_b$num_factors <- 10

views_b <- views_names(mofa_b)
stopifnot(chip_view %in% views_b)
model_opts_b$likelihoods                <- rep("gaussian", length(views_b))
names(model_opts_b$likelihoods)         <- views_b
model_opts_b$likelihoods[chip_view]     <- "bernoulli"

train_opts_b <- get_default_training_options(mofa_b)
train_opts_b$convergence_mode <- "medium"
train_opts_b$maxiter          <- 5000
train_opts_b$seed             <- 1

data_opts_b  <- get_default_data_options(mofa_b)
data_opts_b$scale_views <- FALSE

mofa_b <- prepare_mofa(mofa_b,
                       model_options   = model_opts_b,
                       training_options = train_opts_b,
                       data_options    = data_opts_b)

mofa_trained <- run_mofa(mofa_b, use_basilisk = TRUE, outfile = out_hdf5)

# -- 5B. Companion model: all-Gaussian (for variance-explained metrics only) --
mofa_g       <- create_mofa(data_list)
samples_metadata(mofa_g) <- sample_md

model_opts_g <- get_default_model_options(mofa_g)
model_opts_g$num_factors <- 10
views_g <- views_names(mofa_g)
model_opts_g$likelihoods <- setNames(rep("gaussian", length(views_g)), views_g)

train_opts_g <- get_default_training_options(mofa_g)
train_opts_g$convergence_mode <- "medium"
train_opts_g$maxiter          <- 5000
train_opts_g$seed             <- 1

data_opts_g  <- get_default_data_options(mofa_g)
data_opts_g$scale_views <- FALSE

mofa_g <- prepare_mofa(mofa_g,
                       model_options    = model_opts_g,
                       training_options = train_opts_g,
                       data_options     = data_opts_g)

mofa_g_trained <- run_mofa(mofa_g, use_basilisk = TRUE, outfile = out_hdf5_g)

# -- Sanity check: confirm Bernoulli applied to CHIP view ---------------------
print(mofa_trained@model_options$likelihoods)

# =============================================================================
# 6. Convergence QC and variance explained
# =============================================================================

# ELBO convergence (primary Bernoulli model)
elbo       <- mofa_trained@training_stats$elbo
elbo_valid <- elbo[!is.nan(elbo)]
plot(elbo_valid, type = "l", xlab = "ELBO index", ylab = "ELBO",
     main = "MOFA2 convergence (Bernoulli CHIP model)")

# Variance explained from all-Gaussian companion
var_exp <- calculate_variance_explained(mofa_g_trained)
print(var_exp$r2_total)

p_var_heat <- plot_variance_explained(mofa_g_trained, x = "view", y = "factor") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_var_tot  <- plot_variance_explained(mofa_g_trained, x = "view", y = "factor",
                                      plot_total = TRUE)[[2]] +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_var_heat)
print(p_var_tot)

# =============================================================================
# 7. Cox models: factor associations with time-to-progression
# =============================================================================

# -- Build analysis data frame ------------------------------------------------
Z          <- MOFA2::get_factors(mofa_trained, factors = "all")[[1]]  # samples x factors
df_factors <- as.data.frame(Z) %>% rownames_to_column("sample_id")

df <- meta %>%
  mutate(sample_id = as.character(sample_id)) %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  filter(ccp_2or3_bin == 1) %>%
  inner_join(df_factors, by = "sample_id")

stopifnot(all(c(time_var, event_var) %in% colnames(df)))
df[[event_var]] <- as.integer(df[[event_var]])
df[[time_var]]  <- as.numeric(df[[time_var]])
df <- df %>% filter(!is.na(.data[[time_var]]), !is.na(.data[[event_var]]))

write_xlsx(df, path = "MOFA_master_table.xlsx")

# -- Survival object ----------------------------------------------------------
surv_obj   <- Surv(time = df[[time_var]], event = df[[event_var]])
factor_cols <- grep("^Factor", colnames(df), value = TRUE)
stopifnot(length(factor_cols) > 0)

# -- Univariable Cox (all factors) --------------------------------------------
uni_res <- lapply(factor_cols, function(f) {
  fit <- coxph(as.formula(paste0("surv_obj ~ scale(", f, ")")), data = df)
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(factor = f)
}) %>%
  bind_rows() %>%
  select(factor, estimate, conf.low, conf.high, p.value) %>%
  arrange(p.value) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"))

# -- Adjusted Cox (all factors; covariates retained if present in df) ---------
covars <- c("age_at_BL", "gender_bin", "palindromic", "smoking_fact")
covars <- covars[covars %in% colnames(df)]

adj_res <- lapply(factor_cols, function(f) {
  rhs <- paste(c(paste0("scale(", f, ")"), covars), collapse = " + ")
  fit <- coxph(as.formula(paste0("surv_obj ~ ", rhs)), data = df)
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl(paste0("^scale\\(", f, "\\)$"), term)) %>%
    mutate(factor = f) %>%
    select(factor, estimate, conf.low, conf.high, p.value)
}) %>%
  bind_rows() %>%
  arrange(p.value) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"))

print(uni_res)
print(adj_res)

# -- Main-effects and interaction models for Factor 6 -------------------------
factor_of_interest <- "Factor6"
stopifnot(factor_of_interest %in% factor_cols)

# Ensure palindromic is numeric 0/1 if it was stored as character "0"/"1"
if (is.factor(df$palindromic) || is.character(df$palindromic)) {
  if (all(na.omit(unique(as.character(df$palindromic))) %in% c("0", "1"))) {
    df$palindromic <- as.integer(as.character(df$palindromic))
  }
}

covars_no_pal <- setdiff(covars, "palindromic")

rhs_main <- paste(c(paste0("scale(", factor_of_interest, ")"), "palindromic", covars_no_pal),
                  collapse = " + ")
fit_no_int <- coxph(as.formula(paste0("surv_obj ~ ", rhs_main)), data = df)

rhs_int <- paste(c(paste0("scale(", factor_of_interest, ") * palindromic"), covars_no_pal),
                 collapse = " + ")
fit_int <- coxph(as.formula(paste0("surv_obj ~ ", rhs_int)), data = df)

lrt <- anova(fit_no_int, fit_int, test = "LRT")
print(lrt)
print(summary(fit_no_int))
print(summary(fit_int))

# -- Proportional hazards assumption (Schoenfeld residuals) -------------------
zph_no_int <- cox.zph(fit_no_int, transform = "km")
zph_int    <- cox.zph(fit_int,    transform = "km")

cat("\n--- PH test: main effects model ---\n")
print(zph_no_int)
cat("\n--- PH test: interaction model ---\n")
print(zph_int)

op <- par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
plot(zph_no_int, main = paste0("PH check: main (", factor_of_interest, ")"))
abline(h = 0, lty = 2)
plot(zph_int, main = paste0("PH check: interaction (", factor_of_interest, "\u00d7palindromic)"))
abline(h = 0, lty = 2)
par(op)

# -- Stratum-specific HR for Factor 6 (from interaction model) ----------------
b         <- coef(fit_int)
V         <- vcov(fit_int)
term_f6   <- "scale(Factor6)"
term_int  <- "scale(Factor6):palindromic"
stopifnot(all(c(term_f6, term_int) %in% names(b)))

beta0 <- b[term_f6]
se0   <- sqrt(V[term_f6, term_f6])
beta1 <- b[term_f6] + b[term_int]
se1   <- sqrt(V[term_f6, term_f6] + V[term_int, term_int] + 2 * V[term_f6, term_int])

hr_tab <- tibble(
  palindromic = c(0, 1),
  HR      = exp(c(beta0, beta1)),
  CI_low  = exp(c(beta0 - 1.96 * se0, beta1 - 1.96 * se1)),
  CI_high = exp(c(beta0 + 1.96 * se0, beta1 + 1.96 * se1))
) %>%
  mutate(
    pal_group = factor(palindromic, levels = c(0, 1),
                       labels = c("Non-palindromic", "Palindromic")),
    HR_txt = sprintf("%.2f", HR),
    CI_txt = sprintf("%.2f\u2013%.2f", CI_low, CI_high),
    label  = paste0(HR_txt, " (", CI_txt, ")")
  )

print(hr_tab)

# =============================================================================
# 8. Kaplan-Meier curves (Factor 6)
# =============================================================================
stopifnot(all(c(time_var, event_var, "Factor6") %in% colnames(df)))

make_tertiles <- function(x) {
  z <- as.numeric(scale(x))
  q <- as.numeric(stats::quantile(z, probs = c(1/3, 2/3), na.rm = TRUE))
  if (!is.finite(q[1]) || !is.finite(q[2]) || q[1] == q[2])
    stop("Tertile cutpoints are not valid (possibly too many ties).")
  cut(z, breaks = c(-Inf, q[1], q[2], Inf),
      labels = c("Low", "Mid", "High"), right = TRUE)
}

make_median_split <- function(x) {
  m <- stats::median(x, na.rm = TRUE)
  factor(ifelse(x < m, "Low", "High"), levels = c("Low", "High"))
}

# 8A. Tertile KM
df_km1 <- df %>%
  mutate(f6_tertile = make_tertiles(Factor6)) %>%
  filter(!is.na(f6_tertile), !is.na(.data[[time_var]]), !is.na(.data[[event_var]]))

fit_km1 <- survfit(Surv(df_km1[[time_var]], df_km1[[event_var]]) ~ f6_tertile, data = df_km1)

p_km1 <- ggsurvplot(
  fit_km1, data = df_km1,
  risk.table       = TRUE,
  pval             = TRUE,
  conf.int         = TRUE,
  surv.median.line = "hv",
  xlab             = "Time (years)",
  ylab             = "Progression-free probability",
  title            = "Factor 6 tertiles vs time-to-progression",
  legend.title     = "Factor 6",
  legend.labs      = c("Low", "Mid", "High"),
  palette          = c("blue", "black", "red")
)
print(p_km1)

# 8B. Median-split KM
df_km2 <- df %>%
  mutate(f6_split = make_median_split(Factor6)) %>%
  filter(!is.na(f6_split), !is.na(.data[[time_var]]), !is.na(.data[[event_var]]))

fit_km2 <- survfit(Surv(df_km2[[time_var]], df_km2[[event_var]]) ~ f6_split, data = df_km2)

p_km2 <- ggsurvplot(
  fit_km2, data = df_km2,
  risk.table       = TRUE,
  pval             = TRUE,
  conf.int         = TRUE,
  conf.int.level   = 0.95,
  surv.median.line = "hv",
  xlab             = "Time (years)",
  ylab             = "Progression-free probability",
  title            = "Factor 6 (median split) vs time-to-progression",
  legend.title     = "Factor 6",
  legend.labs      = c("Low", "High"),
  palette          = c("blue", "red")
)

# Reduce CI ribbon alpha
is_ribbon <- vapply(p_km2$plot$layers, function(l) inherits(l$geom, "GeomRibbon"), logical(1))
p_km2$plot$layers[[which(is_ribbon)[1]]]$aes_params$alpha <- 0.15

print(p_km2)

# =============================================================================
# 9. Factor 6 composition: per-view variance and top loadings
# =============================================================================

# 9A. Per-view variance explained (Gaussian companion)
r2_f6_df <- factor_variance_df(var_exp, factor_of_interest)
print(r2_f6_df)
print(plot_factor_variance(r2_f6_df, factor_of_interest))

# 9B. Top loadings per view (from Bernoulli primary model)
get_top_loadings <- function(mofa_obj, view, factor_name, n = 25) {
  W <- MOFA2::get_weights(mofa_obj, views = view, factors = factor_name)[[view]]
  w <- as.numeric(W[, factor_name, drop = TRUE])
  names(w) <- rownames(W)
  dfw <- tibble(feature = names(w), weight = w, abs_weight = abs(w))
  list(
    positive = dfw %>% arrange(desc(weight))  %>% slice_head(n = n),
    negative = dfw %>% arrange(weight)         %>% slice_head(n = n)
  )
}

top_views <- r2_f6_df$view[seq_len(min(5, nrow(r2_f6_df)))]
stopifnot(all(top_views %in% MOFA2::views_names(mofa_trained)))
tops <- setNames(
  lapply(top_views, get_top_loadings, mofa_obj = mofa_trained,
         factor_name = factor_of_interest, n = 25),
  top_views
)

for (v in names(tops)) {
  cat("\n====", v, "====\n")
  cat("Top positive loadings:\n"); print(head(tops[[v]]$positive, 10))
  cat("Top negative loadings:\n"); print(head(tops[[v]]$negative, 10))
}

# =============================================================================
# 10. Factor 6 view-level visualisations
# =============================================================================

# Factor 6 scores (aligned to df sample IDs)
f6 <- setNames(df$Factor6, df$sample_id)
Z  <- MOFA2::get_factors(mofa_trained, factors = "all")[[1]]

# -- 10A. RNA: Hallmark GSEA and module score correlation ---------------------
rna_view <- "Gene expression"
stopifnot(rna_view %in% names(MOFA2::get_data(mofa_trained)))

W_rna  <- MOFA2::get_weights(mofa_trained, views = rna_view, factors = factor_of_interest)[[rna_view]]
w_rna  <- setNames(as.numeric(W_rna[, factor_of_interest, drop = TRUE]), rownames(W_rna))
rna_loadings <- tibble(feature = names(w_rna), loading = w_rna,
                       gene = sub("^RNA_", "", names(w_rna)))

gsea_f6 <- run_hallmark_gsea(rna_loadings)
print(plot_gsea_bubbles(tidy_gsea_df(gsea_f6), factor_of_interest, rna_view))

X_rna   <- MOFA2::get_data(mofa_trained)[[rna_view]][["group1"]]
f6_rna  <- f6[colnames(X_rna)]
N_module <- 100

pos_features <- rna_loadings %>% arrange(desc(loading)) %>% slice_head(n = N_module) %>%
  pull(feature) %>% intersect(rownames(X_rna))
neg_features <- rna_loadings %>% arrange(loading) %>% slice_head(n = N_module) %>%
  pull(feature) %>% intersect(rownames(X_rna))
stopifnot(length(pos_features) >= 10, length(neg_features) >= 10)

score_df <- tibble(
  sample_id     = colnames(X_rna),
  Factor6       = as.numeric(f6_rna),
  F6_pos_module = colMeans(X_rna[pos_features, , drop = FALSE], na.rm = TRUE),
  F6_neg_module = colMeans(X_rna[neg_features, , drop = FALSE], na.rm = TRUE)
) %>% filter(!is.na(Factor6))

cor_pos <- cor(score_df$Factor6, score_df$F6_pos_module, use = "complete.obs")
cor_neg <- cor(score_df$Factor6, score_df$F6_neg_module, use = "complete.obs")

p_rna_corr <- ggplot(score_df, aes(x = Factor6, y = F6_neg_module)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw(base_size = 12) +
  labs(x = paste0(factor_of_interest, " score"),
       y = paste0("Mean expression of top ", N_module, " negative-loading genes"),
       title = paste0("RNA programme validation \u2013 ", factor_of_interest),
       subtitle = paste0("Pearson r = ", round(cor_neg, 2),
                         " (positive-module r = ", round(cor_pos, 2), ")"))
print(p_rna_corr)

# -- 10B. ASC specks vs Factor 6 ----------------------------------------------
asc_view <- "Inflammasome ASC specks"
asc_X    <- as.data.frame(t(data_list[[asc_view]])) %>% rownames_to_column("sample_id")
df_asc   <- tibble(sample_id = names(f6), Factor6 = as.numeric(f6)) %>%
  inner_join(asc_X, by = "sample_id")
stopifnot(nrow(df_asc) > 0)

for (asc_feat in c("ASC_ASC_NLRP3_events", "ASC_ASC_events")) {
  if (!asc_feat %in% colnames(df_asc)) next
  ct  <- cor.test(df_asc$Factor6, df_asc[[asc_feat]], method = "spearman", exact = FALSE)
  lab <- paste0("Spearman \u03c1 = ", round(unname(ct$estimate), 2),
                "\np = ", format.pval(ct$p.value, digits = 2, eps = 1e-3))
  ylabel <- if (grepl("NLRP3", asc_feat)) "ASC-NLRP3 specks (events)" else "Non-specific ASC specks (events)"
  print(
    ggplot(df_asc, aes(x = Factor6, y = .data[[asc_feat]])) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE, colour = "#B2182B", fill = "#B2182B", alpha = 0.20) +
      annotate("text", x = Inf, y = Inf, label = lab, hjust = 1.1, vjust = 1.2, size = 4) +
      theme_bw(base_size = 12) +
      labs(x = "MOFA Factor 6 score", y = ylabel,
           title = paste0("Factor 6 vs ", ylabel))
  )
}

# -- 10C. CHIP status vs Factor 6 ---------------------------------------------
X_chip      <- MOFA2::get_data(mofa_trained)[[chip_view]][["group1"]]
common_samp <- intersect(colnames(X_chip), names(f6))
stopifnot(length(common_samp) >= 10)

chip_any <- colSums(X_chip[, common_samp, drop = FALSE] > 0, na.rm = TRUE) > 0
chip_df  <- tibble(
  sample_id   = common_samp,
  Factor6     = as.numeric(f6[common_samp]),
  CHIP_status = factor(ifelse(chip_any, "CHIP+", "CHIP-"), levels = c("CHIP-", "CHIP+"))
) %>%
  left_join(df %>% select(sample_id, palindromic), by = "sample_id") %>%
  mutate(palindromic = factor(palindromic, levels = c(0, 1),
                              labels = c("Non-palindromic", "Palindromic")))

wtest   <- wilcox.test(Factor6 ~ CHIP_status, data = chip_df)
p_lab   <- ifelse(wtest$p.value < 0.001, "p < 0.001",
                  paste0("p = ", formatC(wtest$p.value, format = "f", digits = 3)))
y_max   <- max(chip_df$Factor6, na.rm = TRUE)
y_range <- diff(range(chip_df$Factor6, na.rm = TRUE))
y_bar   <- y_max + 0.10 * y_range
y_text  <- y_max + 0.14 * y_range

print(
  ggplot(chip_df, aes(x = CHIP_status, y = Factor6)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, fill = "grey90", colour = "black") +
    geom_jitter(aes(colour = CHIP_status), width = 0.15, alpha = 0.7, size = 1.8) +
    geom_segment(aes(x = 1, xend = 2, y = y_bar, yend = y_bar), linewidth = 0.6) +
    geom_segment(aes(x = 1, xend = 1, y = y_bar - 0.03 * y_range, yend = y_bar), linewidth = 0.6) +
    geom_segment(aes(x = 2, xend = 2, y = y_bar - 0.03 * y_range, yend = y_bar), linewidth = 0.6) +
    annotate("text", x = 1.5, y = y_text, label = p_lab, size = 5) +
    scale_colour_manual(values = c("CHIP-" = "#2166AC", "CHIP+" = "#B2182B"), guide = "none") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(size = 14)) +
    labs(x = NULL, y = "MOFA Factor 6 score", title = "Factor 6 by CHIP status")
)

# -- 10D. HLA: loadings, DRB1*04:01 carrier, shared epitope dosage -----------
hla_view        <- names(MOFA2::get_data(mofa_trained))[str_detect(
  tolower(names(MOFA2::get_data(mofa_trained))), "hla")][1]
if (is.na(hla_view)) stop("No HLA view found. Check names(get_data(mofa_trained)).")
cat("Using HLA view:", hla_view, "\n")

W_hla       <- MOFA2::get_weights(mofa_trained, views = hla_view, factors = factor_of_interest)[[hla_view]]
w_hla       <- setNames(as.numeric(W_hla[, factor_of_interest, drop = TRUE]), rownames(W_hla))
hla_loadings <- tibble(feature = names(w_hla), loading = w_hla) %>%
  arrange(desc(abs(loading)))

N_hla <- 12
plot_df_hla <- hla_loadings %>%
  arrange(desc(abs(loading))) %>%
  slice_head(n = 2 * N_hla) %>%
  mutate(label = format_hla_feature(feature),
         sign  = ifelse(loading >= 0, "Positive", "Negative"),
         label = fct_reorder(label, loading))

print(
  ggplot(plot_df_hla, aes(x = label, y = loading, fill = sign)) +
    geom_col(width = 0.75) +
    coord_flip() +
    scale_fill_manual(values = c("Negative" = COL_BLUE, "Positive" = COL_RED), guide = "none") +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = paste0("MOFA weight (", factor_of_interest, ")"),
         title = "HLA feature loadings \u2013 Factor 6")
)

X_hla       <- MOFA2::get_data(mofa_trained)[[hla_view]][["group1"]]
common_samp <- intersect(colnames(X_hla), names(f6))
stopifnot(length(common_samp) >= 10)
X_hla2      <- X_hla[, common_samp, drop = FALSE]
f6_use      <- as.numeric(f6[common_samp])

feat_0401   <- "HLA_drb1_DRB1_04_01_carrier"
feat_se_dos <- "HLA_se_dosage"
stopifnot(feat_0401 %in% rownames(X_hla2), feat_se_dos %in% rownames(X_hla2))

val_df_hla <- tibble(
  sample_id         = common_samp,
  Factor6           = f6_use,
  DRB1_0401_carrier = as.numeric(X_hla2[feat_0401, common_samp]),
  SE_dosage         = as.numeric(X_hla2[feat_se_dos, common_samp])
) %>%
  mutate(carrier_status = factor(ifelse(DRB1_0401_carrier > 0, "Carrier", "Non-carrier"),
                                 levels = c("Non-carrier", "Carrier"))) %>%
  filter(!is.na(Factor6))

# DRB1*04:01 carrier box plot
df_0401 <- val_df_hla %>% filter(!is.na(carrier_status))
wtest_0401 <- wilcox.test(Factor6 ~ carrier_status, data = df_0401)
p_lab  <- ifelse(wtest_0401$p.value < 0.001, "p < 0.001",
                 paste0("p = ", sprintf("%.3f", wtest_0401$p.value)))
y_max  <- max(df_0401$Factor6, na.rm = TRUE)
y_rng  <- diff(range(df_0401$Factor6, na.rm = TRUE))
y_bar  <- y_max + 0.10 * y_rng
y_txt  <- y_max + 0.14 * y_rng

print(
  ggplot(df_0401, aes(x = carrier_status, y = Factor6)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, fill = "grey90", colour = "black") +
    geom_jitter(aes(colour = carrier_status), width = 0.15, alpha = 0.7, size = 1.8) +
    geom_segment(aes(x = 1, xend = 2, y = y_bar, yend = y_bar), linewidth = 0.6) +
    geom_segment(aes(x = 1, xend = 1, y = y_bar - 0.03 * y_rng, yend = y_bar), linewidth = 0.6) +
    geom_segment(aes(x = 2, xend = 2, y = y_bar - 0.03 * y_rng, yend = y_bar), linewidth = 0.6) +
    annotate("text", x = 1.5, y = y_txt, label = p_lab, size = 5) +
    scale_colour_manual(values = c("Non-carrier" = "#B2182B", "Carrier" = "#2166AC"), guide = "none") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(size = 14)) +
    labs(x = NULL, y = "MOFA Factor 6 score", title = "Factor 6 by HLA-DRB1*04:01 carrier status")
)

# Shared epitope dosage scatter
df_se  <- val_df_hla %>% filter(!is.na(SE_dosage))
ct_se  <- cor.test(df_se$Factor6, df_se$SE_dosage, method = "spearman", exact = FALSE)
sub_lab <- paste0("Spearman \u03c1 = ", round(unname(ct_se$estimate), 2),
                  ", p = ", format.pval(ct_se$p.value, digits = 2))

print(
  ggplot(df_se, aes(x = SE_dosage, y = Factor6)) +
    geom_jitter(width = 0.08, height = 0, alpha = 0.65, size = 1.6) +
    geom_smooth(method = "lm", se = TRUE, colour = "#2166AC", fill = "#2166AC", alpha = 0.20) +
    scale_x_continuous(breaks = seq(floor(min(df_se$SE_dosage, na.rm = TRUE)),
                                    ceiling(max(df_se$SE_dosage, na.rm = TRUE)), by = 1)) +
    theme_bw(base_size = 12) +
    labs(x = "Shared epitope dosage", y = "MOFA Factor 6 score",
         title = "Factor 6 vs shared epitope dosage", subtitle = sub_lab)
)

# -- 10E. Germline burdens vs Factor 6 ----------------------------------------
germ_view   <- "Rare Germline variants"
germ_mat_f6 <- data_list[[germ_view]]
common_samp <- intersect(colnames(germ_mat_f6), names(f6))

germ_df <- as.data.frame(t(germ_mat_f6[, common_samp, drop = FALSE])) %>%
  rownames_to_column("sample_id") %>%
  mutate(Factor6 = as.numeric(f6[common_samp]))

disc_feats <- c("GERM_germ_DISC_n_damaging_variants", "GERM_germ_DISC_n_damaging_genes")
ra_feats   <- c("GERM_germ_RA_n_damaging_variants",   "GERM_germ_RA_n_damaging_genes")
stopifnot(all(disc_feats %in% colnames(germ_df)), all(ra_feats %in% colnames(germ_df)))

germ_df <- germ_df %>%
  mutate(DISC_damaging_sum = rowSums(across(all_of(disc_feats)), na.rm = TRUE),
         RA_damaging_sum   = rowSums(across(all_of(ra_feats)),   na.rm = TRUE))

ct_disc <- cor.test(germ_df$Factor6, germ_df$DISC_damaging_sum, method = "spearman", exact = FALSE)
ct_ra   <- cor.test(germ_df$Factor6, germ_df$RA_damaging_sum,   method = "spearman", exact = FALSE)

p_disc <- ggplot(germ_df, aes(x = DISC_damaging_sum, y = Factor6)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, colour = "black") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Spearman \u03c1 = ", round(unname(ct_disc$estimate), 2),
                          "\np = ", format.pval(ct_disc$p.value, digits = 2, eps = 1e-3)),
           hjust = 1.1, vjust = 1.2, size = 4) +
  theme_bw(base_size = 12) +
  labs(x = "DISC damaging burden (variants + genes)", y = "MOFA Factor 6 score",
       title = "Factor 6 vs DISC damaging burden")

p_ra <- ggplot(germ_df, aes(x = RA_damaging_sum, y = Factor6)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, colour = "black") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Spearman \u03c1 = ", round(unname(ct_ra$estimate), 2),
                          "\np = ", format.pval(ct_ra$p.value, digits = 2, eps = 1e-3)),
           hjust = 1.1, vjust = 1.2, size = 4) +
  theme_bw(base_size = 12) +
  labs(x = "RA damaging burden (variants + genes)", y = "MOFA Factor 6 score",
       title = "Factor 6 vs RA damaging burden")

print(p_disc + p_ra)

# Export low-Factor 6 sample IDs
low_F6_ids <- df %>%
  filter(!is.na(Factor6), Factor6 <= median(Factor6, na.rm = TRUE)) %>%
  select(sample_id)
write_xlsx(low_F6_ids, path = "low_Factor6_sample_ids.xlsx")

# Germline loadings and carrier validation for Factor 6
W_germ_f6  <- MOFA2::get_weights(mofa_trained, views = germ_view, factors = factor_of_interest)[[germ_view]]
w_germ_f6  <- setNames(as.numeric(W_germ_f6[, factor_of_interest, drop = TRUE]), rownames(W_germ_f6))
germ_loadings_f6 <- tibble(feature = names(w_germ_f6), loading = w_germ_f6)

N_germ <- 20
plot_df_germ <- germ_loadings_f6 %>%
  arrange(desc(abs(loading))) %>%
  slice_head(n = N_germ) %>%
  mutate(label = format_germ_feature(feature),
         sign  = ifelse(loading >= 0, "Positive", "Negative"),
         label = fct_reorder(label, loading))

print(
  ggplot(plot_df_germ, aes(x = label, y = loading, fill = sign)) +
    geom_col(width = 0.75) +
    coord_flip() +
    scale_fill_manual(values = c("Negative" = COL_BLUE, "Positive" = COL_RED), guide = "none") +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = paste0("MOFA weight (", factor_of_interest, ")"),
         title = "Germline feature loadings \u2013 Factor 6")
)

X_germ_f6   <- MOFA2::get_data(mofa_trained)[[germ_view]][["group1"]]
common_samp <- intersect(colnames(X_germ_f6), names(f6))
stopifnot(length(common_samp) >= 10)
X_germ2     <- X_germ_f6[, common_samp, drop = FALSE]
f6_use      <- as.numeric(f6[common_samp])

ct_germ <- cor.test(f6_use, colSums(X_germ2, na.rm = TRUE), method = "spearman", exact = FALSE)
df_burden_f6 <- tibble(
  sample_id   = common_samp,
  Factor6     = f6_use,
  Germ_burden = colSums(X_germ2, na.rm = TRUE)
) %>% filter(!is.na(Factor6), !is.na(Germ_burden))

print(
  ggplot(df_burden_f6, aes(x = Germ_burden, y = Factor6)) +
    geom_point(alpha = 0.7, size = 1.7) +
    geom_smooth(method = "lm", se = TRUE, colour = COL_BLUE, fill = COL_BLUE, alpha = 0.20) +
    theme_bw(base_size = 12) +
    labs(x = "Germline variant burden (sum across features)",
         y = paste0(factor_of_interest, " score"),
         title = paste0(factor_of_interest, " vs germline burden"),
         subtitle = paste0("Spearman \u03c1 = ", round(unname(ct_germ$estimate), 2),
                           ", p = ", format.pval(ct_germ$p.value, digits = 2)))
)

# Carrier status of top germline feature vs Factor 6
top_feat <- plot_df_germ %>% arrange(desc(abs(loading))) %>% slice_head(n = 1) %>% pull(feature)
if (top_feat %in% rownames(X_germ2)) {
  df_top <- tibble(
    sample_id     = common_samp,
    Factor6       = f6_use,
    dosage        = as.numeric(X_germ2[top_feat, common_samp])
  ) %>%
    filter(!is.na(dosage), !is.na(Factor6)) %>%
    mutate(carrier_status = factor(ifelse(dosage > 0, "Carrier", "Non-carrier"),
                                   levels = c("Non-carrier", "Carrier")))

  wtest <- wilcox.test(Factor6 ~ carrier_status, data = df_top)
  p_lab <- ifelse(wtest$p.value < 0.001, "p < 0.001",
                  paste0("p = ", sprintf("%.3f", wtest$p.value)))
  y_max  <- max(df_top$Factor6, na.rm = TRUE)
  y_rng  <- diff(range(df_top$Factor6, na.rm = TRUE))

  print(
    ggplot(df_top, aes(x = carrier_status, y = Factor6)) +
      geom_boxplot(outlier.shape = NA, width = 0.6, fill = "grey90", colour = "black") +
      geom_jitter(aes(colour = carrier_status), width = 0.15, alpha = 0.7, size = 1.8) +
      geom_segment(aes(x = 1, xend = 2, y = y_max + 0.10 * y_rng, yend = y_max + 0.10 * y_rng), linewidth = 0.6) +
      geom_segment(aes(x = 1, xend = 1, y = y_max + 0.07 * y_rng, yend = y_max + 0.10 * y_rng), linewidth = 0.6) +
      geom_segment(aes(x = 2, xend = 2, y = y_max + 0.07 * y_rng, yend = y_max + 0.10 * y_rng), linewidth = 0.6) +
      annotate("text", x = 1.5, y = y_max + 0.14 * y_rng, label = p_lab, size = 5) +
      scale_colour_manual(values = c("Non-carrier" = COL_RED, "Carrier" = COL_BLUE), guide = "none") +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(size = 14)) +
      labs(x = NULL, y = paste0(factor_of_interest, " score"),
           title = "Factor 6 vs germline variant carrier status (Discovery panel)",
           subtitle = format_germ_feature(top_feat))
  )
}

# =============================================================================
# 11. Palindromic-status logistic regression (all factors)
# =============================================================================
df <- df %>% mutate(palindromic = as.integer(palindromic))
factor_test <- grep("^Factor", colnames(df), value = TRUE)

uni_pal <- lapply(factor_test, function(f) {
  fit <- glm(as.formula(paste0("palindromic ~ scale(", f, ")")),
             data = df, family = binomial())
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == paste0("scale(", f, ")")) %>%
    mutate(factor = f) %>%
    select(factor, estimate, conf.low, conf.high, p.value)
}) %>%
  bind_rows() %>%
  arrange(p.value) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"))

covars_pal <- c("age_at_BL", "gender_bin", "smoking_fact")
covars_pal <- covars_pal[covars_pal %in% colnames(df)]

adj_pal <- lapply(factor_test, function(f) {
  rhs <- paste(c(paste0("scale(", f, ")"), covars_pal), collapse = " + ")
  fit <- glm(as.formula(paste0("palindromic ~ ", rhs)), data = df, family = binomial())
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == paste0("scale(", f, ")")) %>%
    mutate(factor = f) %>%
    select(factor, estimate, conf.low, conf.high, p.value)
}) %>%
  bind_rows() %>%
  arrange(p.value) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"))

print(uni_pal)
print(adj_pal)

# Forest plot of FDR-significant factors (adjusted model)
make_forest_2panel <- function(plot_df, title = "Forest plot") {
  xmin <- min(plot_df$conf.low,  na.rm = TRUE) / 1.2
  xmax <- max(plot_df$conf.high, na.rm = TRUE) * 1.2

  p_forest <- ggplot(plot_df, aes(x = estimate, y = label)) +
    geom_vline(xintercept = 1, linetype = 2) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0) +
    geom_point(size = 2) +
    scale_x_log10(limits = c(xmin, xmax)) +
    labs(x = "Odds ratio (log scale)", y = NULL, title = title) +
    theme_bw(base_size = 12) +
    theme(plot.title  = element_text(face = "bold"),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 11),
          plot.margin = margin(5.5, 5.5, 5.5, 5.5))

  p_table <- ggplot(plot_df, aes(y = label)) +
    geom_text(aes(x = 0, label = annot), hjust = 0, size = 4) +
    xlim(0, 1) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = 12) +
    theme(plot.margin = margin(33, 5.5, 5.5, 5.5))

  cowplot::plot_grid(p_forest, p_table, ncol = 2, align = "h", rel_widths = c(0.62, 0.38))
}

plot_df_pal_sig <- adj_pal %>%
  filter(p_adj < 0.05) %>%
  arrange(p.value) %>%
  mutate(
    label   = paste0("MOFA Factor ", gsub("Factor", "", factor), " (per SD)"),
    p_txt   = if_else(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
    fdr_txt = if_else(p_adj   < 0.001, "<0.001", sprintf("%.3f", p_adj)),
    or_txt  = sprintf("%.2f", estimate),
    ci_txt  = sprintf("%.2f\u2013%.2f", conf.low, conf.high),
    annot   = paste0(or_txt, " (", ci_txt, "), p=", p_txt, ", FDR=", fdr_txt),
    label   = factor(label, levels = rev(unique(label)))
  )

print(make_forest_2panel(plot_df_pal_sig,
                         title = "Palindromic status: adjusted logistic regression (FDR < 0.05)"))

# =============================================================================
# 12. Per-factor variance explained (Factors 7, 8, 10)
# =============================================================================
for (fac in c("Factor7", "Factor8", "Factor10")) {
  r2_df <- factor_variance_df(var_exp, fac)
  print(r2_df)
  print(plot_factor_variance(r2_df, fac))
}

# =============================================================================
# 13. Factor 7 analysis
# =============================================================================
factor_of_interest <- "Factor7"
pal_cols <- c("Non-palindromic" = COL_BLUE, "Palindromic" = COL_RED)

f7 <- setNames(Z[, factor_of_interest], rownames(Z))

# -- 13A. IFN score vs Factor 7 (faceted by palindromic status) ---------------
ifn_view    <- "Interferon score"
ifn_X       <- as.data.frame(t(data_list[[ifn_view]])) %>% rownames_to_column("sample_id")
df_ifn_f7   <- tibble(sample_id = names(f7), Factor7 = as.numeric(f7)) %>%
  inner_join(ifn_X, by = "sample_id")
stopifnot(nrow(df_ifn_f7) > 0)

ifn_feature <- intersect(colnames(df_ifn_f7), "IFN_IFN_protein_score")[1]
stopifnot(!is.na(ifn_feature))

df_ifn_f7 <- df_ifn_f7 %>%
  left_join(df %>% select(sample_id, palindromic), by = "sample_id") %>%
  mutate(palindromic = factor(palindromic, levels = c(0, 1),
                              labels = c("Non-palindromic", "Palindromic"))) %>%
  filter(!is.na(palindromic))

stat_ifn_f7 <- df_ifn_f7 %>%
  group_by(palindromic) %>%
  summarise(
    rho = unname(cor.test(Factor7, .data[[ifn_feature]], method = "spearman", exact = FALSE)$estimate),
    p   = cor.test(Factor7, .data[[ifn_feature]], method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = paste0("Spearman \u03c1 = ", round(rho, 2),
                        "\np = ", format.pval(p, digits = 2, eps = 1e-3)))

print(
  ggplot(df_ifn_f7, aes(x = Factor7, y = .data[[ifn_feature]])) +
    geom_point(alpha = 0.7, size = 1.6) +
    geom_smooth(aes(colour = palindromic, fill = palindromic),
                method = "lm", se = TRUE, alpha = 0.20, linewidth = 0.9) +
    geom_text(data = stat_ifn_f7, aes(x = -Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = -0.05, vjust = 1.10, size = 4) +
    facet_wrap(~ palindromic, nrow = 1) +
    scale_colour_manual(values = pal_cols, guide = "none") +
    scale_fill_manual(values = pal_cols, guide = "none") +
    theme_bw(base_size = 12) +
    theme(strip.text  = element_text(face = "bold", size = 13),
          axis.title  = element_text(size = 13),
          plot.title  = element_text(face = "bold")) +
    labs(x = "MOFA Factor 7 score", y = "IFN protein score (z-scored)",
         title = "Factor 7 vs IFN protein score stratified by palindromic status")
)

# Interaction test
fit_ifn_int <- lm(as.formula(paste0(ifn_feature, " ~ Factor7 * palindromic")), data = df_ifn_f7)
p_int <- summary(fit_ifn_int)$coefficients["Factor7:palindromicPalindromic", "Pr(>|t|)"]
cat("Factor 7 x palindromic interaction p =", format.pval(p_int, digits = 2, eps = 1e-3), "\n")

# -- 13B. RNA: loadings, GSEA, module validation (Factor 7) -------------------
rna_view     <- "Gene expression"
W_rna_f7     <- MOFA2::get_weights(mofa_trained, views = rna_view, factors = factor_of_interest)[[rna_view]]
w_rna_f7     <- setNames(as.numeric(W_rna_f7[, factor_of_interest, drop = TRUE]), rownames(W_rna_f7))
rna_loadings <- tibble(feature = names(w_rna_f7), loading = w_rna_f7,
                       gene = sub("^RNA_", "", names(w_rna_f7)))

# Top loadings barplot
N_rna <- 10
plot_df_rna <- rna_loadings %>%
  arrange(desc(abs(loading))) %>%
  slice_head(n = 2 * N_rna) %>%
  mutate(direction = ifelse(loading >= 0, "Positive", "Negative"),
         gene = fct_reorder(gene, loading))

print(
  ggplot(plot_df_rna, aes(x = gene, y = loading, fill = direction)) +
    geom_col(width = 0.75) + coord_flip() +
    scale_fill_manual(values = c("Positive" = COL_RED, "Negative" = COL_BLUE), guide = "none") +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = paste0("MOFA weight (", factor_of_interest, ")"),
         title = paste0("Top RNA loadings \u2013 ", factor_of_interest), subtitle = rna_view)
)

# Hallmark GSEA
gsea_f7 <- run_hallmark_gsea(rna_loadings)
print(plot_gsea_bubbles(tidy_gsea_df(gsea_f7), factor_of_interest, rna_view))

# Module score validation
X_rna_f7    <- MOFA2::get_data(mofa_trained)[[rna_view]][["group1"]]
common_samp <- intersect(colnames(X_rna_f7), names(f7))
X_rna2      <- X_rna_f7[, common_samp, drop = FALSE]
f7_use      <- as.numeric(f7[common_samp])

pos_genes <- rna_loadings %>% arrange(desc(loading)) %>% slice_head(n = 100) %>%
  pull(feature) %>% intersect(rownames(X_rna2))
neg_genes <- rna_loadings %>% arrange(loading) %>% slice_head(n = 100) %>%
  pull(feature) %>% intersect(rownames(X_rna2))
stopifnot(length(pos_genes) >= 20, length(neg_genes) >= 20)

score_rna_f7 <- tibble(
  sample_id  = common_samp,
  Factor7    = f7_use,
  Pos_module = colMeans(X_rna2[pos_genes, , drop = FALSE], na.rm = TRUE),
  Neg_module = colMeans(X_rna2[neg_genes, , drop = FALSE], na.rm = TRUE)
) %>% filter(!is.na(Factor7))

ct_f7 <- cor.test(score_rna_f7$Factor7, score_rna_f7$Pos_module, method = "spearman", exact = FALSE)
print(
  ggplot(score_rna_f7, aes(x = Factor7, y = Pos_module)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, colour = COL_RED, fill = COL_RED, alpha = 0.20) +
    theme_bw(base_size = 12) +
    labs(x = "Factor 7 score", y = "Mean expression of top 100 positive-loading genes",
         title = paste0("RNA programme associated with Factor 7"),
         subtitle = paste0("Spearman \u03c1 = ", round(unname(ct_f7$estimate), 2),
                           ", p = ", format.pval(ct_f7$p.value, digits = 2)))
)

# -- 13C. Cytokines vs Factor 7 (loadings + palindromic-faceted validation) ---
views_available <- names(MOFA2::get_data(mofa_trained))
cyto_view <- views_available[str_detect(tolower(views_available), "cytok")][1]
if (is.na(cyto_view)) cyto_view <- views_available[str_detect(tolower(views_available), "protein")][1]
if (is.na(cyto_view)) cyto_view <- views_available[str_detect(tolower(views_available), "proteom")][1]
if (is.na(cyto_view)) stop("Could not auto-detect cytokine/protein view. Set cyto_view manually.")
cat("Using cytokine view:", cyto_view, "\n")

W_cyto_f7  <- MOFA2::get_weights(mofa_trained, views = cyto_view, factors = factor_of_interest)[[cyto_view]]
w_cyto_f7  <- setNames(as.numeric(W_cyto_f7[, factor_of_interest, drop = TRUE]), rownames(W_cyto_f7))
cyto_loadings_f7 <- tibble(feature = names(w_cyto_f7), loading = w_cyto_f7)

N_cyto <- 45
plot_df_cyto <- cyto_loadings_f7 %>%
  arrange(desc(abs(loading))) %>% slice_head(n = N_cyto) %>%
  mutate(label = format_cyto_feature(feature),
         sign  = ifelse(loading >= 0, "Positive", "Negative"),
         label = fct_reorder(label, loading))

print(
  ggplot(plot_df_cyto, aes(x = label, y = loading, fill = sign)) +
    geom_col(width = 0.75) + coord_flip() +
    scale_fill_manual(values = c("Positive" = COL_RED, "Negative" = COL_BLUE), guide = "none") +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = "MOFA weight (Factor 7)", title = "Cytokine/protein loadings \u2013 Factor 7")
)

X_cyto_f7   <- MOFA2::get_data(mofa_trained)[[cyto_view]][["group1"]]
common_samp <- intersect(colnames(X_cyto_f7), names(f7))
X_cyto2     <- X_cyto_f7[, common_samp, drop = FALSE]
f7_use      <- as.numeric(f7[common_samp])

pos_feats <- cyto_loadings_f7 %>% arrange(desc(loading)) %>% slice_head(n = N_cyto) %>%
  pull(feature) %>% intersect(rownames(X_cyto2))
stopifnot(length(pos_feats) >= 5)

score_cyto_f7 <- tibble(
  sample_id       = common_samp,
  Factor7         = f7_use,
  Cyto_pos_module = colMeans(X_cyto2[pos_feats, , drop = FALSE], na.rm = TRUE)
) %>%
  filter(!is.na(Factor7)) %>%
  add_palindromic()

stat_cyto_f7 <- facet_spearman_labels(score_cyto_f7, "Factor7", "Cyto_pos_module")

print(
  ggplot(score_cyto_f7, aes(x = Factor7, y = Cyto_pos_module)) +
    geom_point(alpha = 0.7, size = 1.6) +
    geom_smooth(aes(colour = palindromic, fill = palindromic),
                method = "lm", se = TRUE, alpha = 0.20, linewidth = 0.9) +
    geom_text(data = stat_cyto_f7, aes(x = -Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = -0.05, vjust = 1.10, size = 4) +
    facet_wrap(~ palindromic, nrow = 1) +
    scale_colour_manual(values = pal_cols, guide = "none") +
    scale_fill_manual(values = pal_cols, guide = "none") +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(face = "bold", size = 13)) +
    labs(x = "MOFA Factor 7 score",
         y = paste0("Mean abundance of top ", N_cyto, " positive-loading cytokines"),
         title = "Factor 7 vs cytokine module stratified by palindromic status")
)

# =============================================================================
# 14. Factor 8 analysis (RNA, cytokines, HLA)
# =============================================================================
factor_of_interest <- "Factor8"
pal_cols <- c("Non-palindromic" = COL_BLUE, "Palindromic" = COL_RED)
car_cols <- c("Non-carrier" = COL_BLUE, "Carrier" = COL_RED)

f8 <- setNames(Z[, factor_of_interest], rownames(Z))
stopifnot(all(c("sample_id", "palindromic") %in% colnames(df)))
df$sample_id <- as.character(df$sample_id)

# -- 14A. RNA: loadings, GSEA, module validation (palindromic facets) ---------
rna_view <- "Gene expression"
W_rna_f8 <- MOFA2::get_weights(mofa_trained, views = rna_view, factors = factor_of_interest)[[rna_view]]
w_rna_f8 <- setNames(as.numeric(W_rna_f8[, factor_of_interest, drop = TRUE]), rownames(W_rna_f8))
rna_loadings <- tibble(feature = names(w_rna_f8), loading = w_rna_f8,
                       gene = sub("^RNA_", "", names(w_rna_f8)))

plot_df_rna8 <- rna_loadings %>%
  arrange(desc(abs(loading))) %>% slice_head(n = 20) %>%
  mutate(direction = ifelse(loading >= 0, "Positive", "Negative"),
         gene = fct_reorder(gene, loading))

print(
  ggplot(plot_df_rna8, aes(x = gene, y = loading, fill = direction)) +
    geom_col(width = 0.75) + coord_flip() +
    scale_fill_manual(values = c("Positive" = COL_RED, "Negative" = COL_BLUE), guide = "none") +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = paste0("MOFA weight (", factor_of_interest, ")"),
         title = "Top RNA loadings \u2013 Factor 8", subtitle = rna_view)
)

gsea_f8 <- run_hallmark_gsea(rna_loadings)
print(plot_gsea_bubbles(tidy_gsea_df(gsea_f8), factor_of_interest, rna_view))

X_rna_f8    <- MOFA2::get_data(mofa_trained)[[rna_view]][["group1"]]
common_samp <- intersect(colnames(X_rna_f8), names(f8))
X_rna2      <- X_rna_f8[, common_samp, drop = FALSE]
f8_use      <- as.numeric(f8[common_samp])

pos_genes <- rna_loadings %>% arrange(desc(loading)) %>% slice_head(n = 100) %>%
  pull(feature) %>% intersect(rownames(X_rna2))
neg_genes <- rna_loadings %>% arrange(loading) %>% slice_head(n = 100) %>%
  pull(feature) %>% intersect(rownames(X_rna2))
stopifnot(length(pos_genes) >= 20, length(neg_genes) >= 20)

score_rna_f8 <- tibble(
  sample_id  = common_samp,
  Factor8    = f8_use,
  Pos_module = colMeans(X_rna2[pos_genes, , drop = FALSE], na.rm = TRUE),
  Neg_module = colMeans(X_rna2[neg_genes, , drop = FALSE], na.rm = TRUE)
) %>%
  filter(!is.na(Factor8)) %>%
  add_palindromic()

stat_rna_f8 <- facet_spearman_labels(score_rna_f8, "Factor8", "Pos_module")

print(
  ggplot(score_rna_f8, aes(x = Factor8, y = Neg_module)) +
    geom_point(alpha = 0.7, size = 1.6) +
    geom_smooth(aes(colour = palindromic, fill = palindromic),
                method = "lm", se = TRUE, alpha = 0.20, linewidth = 0.9) +
    geom_text(data = stat_rna_f8, aes(x = -Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = -0.05, vjust = 1.10, size = 4) +
    facet_wrap(~ palindromic, nrow = 1) +
    scale_colour_manual(values = pal_cols, guide = "none") +
    scale_fill_manual(values = pal_cols, guide = "none") +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(face = "bold", size = 13)) +
    labs(x = "MOFA Factor 8 score",
         y = "Mean expression of top 100 positive-loading genes",
         title = "RNA programme validation \u2013 Factor 8")
)

# -- 14B. Cytokines vs Factor 8 -----------------------------------------------
W_cyto_f8  <- MOFA2::get_weights(mofa_trained, views = cyto_view, factors = factor_of_interest)[[cyto_view]]
w_cyto_f8  <- setNames(as.numeric(W_cyto_f8[, factor_of_interest, drop = TRUE]), rownames(W_cyto_f8))
cyto_loadings_f8 <- tibble(feature = names(w_cyto_f8), loading = w_cyto_f8)

N_cyto <- 45
plot_df_cyto8 <- cyto_loadings_f8 %>%
  arrange(desc(abs(loading))) %>% slice_head(n = N_cyto) %>%
  mutate(label = format_cyto_feature(feature),
         sign  = ifelse(loading >= 0, "Positive", "Negative"),
         label = fct_reorder(label, loading))

print(
  ggplot(plot_df_cyto8, aes(x = label, y = loading, fill = sign)) +
    geom_col(width = 0.75) + coord_flip() +
    scale_fill_manual(values = c("Positive" = COL_RED, "Negative" = COL_BLUE), guide = "none") +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = "MOFA weight (Factor 8)", title = "Cytokine/protein loadings \u2013 Factor 8")
)

X_cyto_f8   <- MOFA2::get_data(mofa_trained)[[cyto_view]][["group1"]]
common_samp <- intersect(colnames(X_cyto_f8), names(f8))
X_cyto2     <- X_cyto_f8[, common_samp, drop = FALSE]
f8_use      <- as.numeric(f8[common_samp])

pos_feats <- cyto_loadings_f8 %>% arrange(desc(loading)) %>% slice_head(n = N_cyto) %>%
  pull(feature) %>% intersect(rownames(X_cyto2))
stopifnot(length(pos_feats) >= 5)

score_cyto_f8 <- tibble(
  sample_id       = common_samp,
  Factor8         = f8_use,
  Cyto_pos_module = colMeans(X_cyto2[pos_feats, , drop = FALSE], na.rm = TRUE)
) %>%
  filter(!is.na(Factor8)) %>%
  add_palindromic()

stat_cyto_f8 <- facet_spearman_labels(score_cyto_f8, "Factor8", "Cyto_pos_module")

print(
  ggplot(score_cyto_f8, aes(x = Factor8, y = Cyto_pos_module)) +
    geom_point(alpha = 0.7, size = 1.6) +
    geom_smooth(aes(colour = palindromic, fill = palindromic),
                method = "lm", se = TRUE, alpha = 0.20, linewidth = 0.9) +
    geom_text(data = stat_cyto_f8, aes(x = -Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = -0.05, vjust = 1.10, size = 4) +
    facet_wrap(~ palindromic, nrow = 1) +
    scale_colour_manual(values = pal_cols, guide = "none") +
    scale_fill_manual(values = pal_cols, guide = "none") +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(face = "bold", size = 13)) +
    labs(x = "MOFA Factor 8 score",
         y = paste0("Mean abundance of top ", N_cyto, " positive-loading cytokines"),
         title = "Cytokine/protein programme validation \u2013 Factor 8")
)

# -- 14C. HLA: loadings and DRB1 carrier validation (Factor 8) ----------------
hla_view_f8 <- names(MOFA2::get_data(mofa_trained))[str_detect(
  tolower(names(MOFA2::get_data(mofa_trained))), "hla")][1]
if (is.na(hla_view_f8)) stop("No HLA view found. Check names(get_data(mofa_trained)).")

W_hla_f8   <- MOFA2::get_weights(mofa_trained, views = hla_view_f8, factors = factor_of_interest)[[hla_view_f8]]
w_hla_f8   <- setNames(as.numeric(W_hla_f8[, factor_of_interest, drop = TRUE]), rownames(W_hla_f8))
hla_loadings_f8 <- tibble(feature = names(w_hla_f8), loading = w_hla_f8)

plot_df_hla8 <- hla_loadings_f8 %>%
  arrange(desc(abs(loading))) %>% slice_head(n = 24) %>%
  mutate(label = format_hla_feature(feature),
         sign  = ifelse(loading >= 0, "Positive", "Negative"),
         label = fct_reorder(label, loading))

print(
  ggplot(plot_df_hla8, aes(x = label, y = loading, fill = sign)) +
    geom_col(width = 0.75) + coord_flip() +
    scale_fill_manual(values = c("Negative" = COL_BLUE, "Positive" = COL_RED), guide = "none") +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = "MOFA weight (Factor 8)", title = "HLA feature loadings \u2013 Factor 8")
)

X_hla_f8    <- MOFA2::get_data(mofa_trained)[[hla_view_f8]][["group1"]]
common_samp <- intersect(colnames(X_hla_f8), names(f8))
stopifnot(length(common_samp) >= 10)
X_hla2      <- X_hla_f8[, common_samp, drop = FALSE]
f8_use      <- as.numeric(f8[common_samp])

feat_0301 <- "HLA_drb1_DRB1_03_01_carrier"
feat_1101 <- "HLA_drb1_DRB1_11_01_carrier"
stopifnot(feat_0301 %in% rownames(X_hla2), feat_1101 %in% rownames(X_hla2))

val_df_hla8 <- tibble(
  sample_id         = common_samp,
  Factor8           = f8_use,
  DRB1_0301_carrier = as.numeric(X_hla2[feat_0301, common_samp]),
  DRB1_1101_carrier = as.numeric(X_hla2[feat_1101, common_samp])
) %>%
  left_join(df %>% select(sample_id, palindromic), by = "sample_id") %>%
  mutate(palindromic = factor(palindromic, levels = c(0, 1),
                              labels = c("Non-palindromic", "Palindromic"))) %>%
  filter(!is.na(palindromic), !is.na(Factor8))

print(plot_carrier_facet(val_df_hla8, "Factor8", "DRB1_0301_carrier",
                         "Factor 8 by HLA-DRB1*03:01 carrier status"))
print(plot_carrier_facet(val_df_hla8, "Factor8", "DRB1_1101_carrier",
                         "Factor 8 by HLA-DRB1*11:01 carrier status"))

# =============================================================================
# 15. Factor 10 analysis (RNA, cytokines, IFN score, rare germline)
# =============================================================================
factor_of_interest <- "Factor10"
# Factor 10 convention: palindromic = blue (negative association)
pal_cols_f10 <- c("Non-palindromic" = COL_RED, "Palindromic" = COL_BLUE)

f10 <- setNames(Z[, factor_of_interest], rownames(Z))
stopifnot(all(c("sample_id", "palindromic") %in% colnames(df)))
df$sample_id <- as.character(df$sample_id)

# -- 15A. RNA: loadings, GSEA, module validation (palindromic facets) ---------
rna_view  <- "Gene expression"
W_rna_f10 <- MOFA2::get_weights(mofa_trained, views = rna_view, factors = factor_of_interest)[[rna_view]]
w_rna_f10 <- setNames(as.numeric(W_rna_f10[, factor_of_interest, drop = TRUE]), rownames(W_rna_f10))
rna_loadings <- tibble(feature = names(w_rna_f10), loading = w_rna_f10,
                       gene = sub("^RNA_", "", names(w_rna_f10)))

plot_df_rna10 <- rna_loadings %>%
  arrange(desc(abs(loading))) %>% slice_head(n = 20) %>%
  mutate(direction = ifelse(loading >= 0, "Positive", "Negative"),
         gene = fct_reorder(gene, loading))

print(
  ggplot(plot_df_rna10, aes(x = gene, y = loading, fill = direction)) +
    geom_col(width = 0.75) + coord_flip() +
    scale_fill_manual(values = c("Positive" = COL_RED, "Negative" = COL_BLUE), guide = "none") +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = paste0("MOFA weight (", factor_of_interest, ")"),
         title = paste0("Top RNA loadings \u2013 ", factor_of_interest), subtitle = rna_view)
)

gsea_f10 <- run_hallmark_gsea(rna_loadings)
print(plot_gsea_bubbles(tidy_gsea_df(gsea_f10), factor_of_interest, rna_view,
                        col_pos = COL_RED, col_neg = COL_BLUE))

X_rna_f10   <- MOFA2::get_data(mofa_trained)[[rna_view]][["group1"]]
common_samp <- intersect(colnames(X_rna_f10), names(f10))
X_rna2      <- X_rna_f10[, common_samp, drop = FALSE]
f10_use     <- as.numeric(f10[common_samp])

pos_genes <- rna_loadings %>% arrange(desc(loading)) %>% slice_head(n = 100) %>%
  pull(feature) %>% intersect(rownames(X_rna2))
neg_genes <- rna_loadings %>% arrange(loading) %>% slice_head(n = 100) %>%
  pull(feature) %>% intersect(rownames(X_rna2))
stopifnot(length(pos_genes) >= 20, length(neg_genes) >= 20)

score_rna_f10 <- tibble(
  sample_id  = common_samp,
  Factor10   = f10_use,
  Pos_module = colMeans(X_rna2[pos_genes, , drop = FALSE], na.rm = TRUE),
  Neg_module = colMeans(X_rna2[neg_genes, , drop = FALSE], na.rm = TRUE)
) %>%
  filter(!is.na(Factor10)) %>%
  add_palindromic()

stat_rna_f10 <- facet_spearman_labels(score_rna_f10, "Factor10", "Pos_module")

print(
  ggplot(score_rna_f10, aes(x = Factor10, y = Pos_module)) +
    geom_point(alpha = 0.7, size = 1.6) +
    geom_smooth(aes(colour = palindromic, fill = palindromic),
                method = "lm", se = TRUE, alpha = 0.20, linewidth = 0.9) +
    geom_text(data = stat_rna_f10, aes(x = -Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = -0.05, vjust = 1.10, size = 4) +
    facet_wrap(~ palindromic, nrow = 1) +
    scale_colour_manual(values = pal_cols_f10, guide = "none") +
    scale_fill_manual(values = pal_cols_f10, guide = "none") +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(face = "bold", size = 13)) +
    labs(x = "MOFA Factor 10 score",
         y = "Mean expression of top 100 positive-loading genes",
         title = "Gene expression validation \u2013 Factor 10")
)

# -- 15B. Cytokines vs Factor 10 ----------------------------------------------
W_cyto_f10  <- MOFA2::get_weights(mofa_trained, views = cyto_view, factors = factor_of_interest)[[cyto_view]]
w_cyto_f10  <- setNames(as.numeric(W_cyto_f10[, factor_of_interest, drop = TRUE]), rownames(W_cyto_f10))
cyto_loadings_f10 <- tibble(feature = names(w_cyto_f10), loading = w_cyto_f10)

N_cyto <- 45
plot_df_cyto10 <- cyto_loadings_f10 %>%
  arrange(desc(abs(loading))) %>% slice_head(n = N_cyto) %>%
  mutate(label = format_cyto_feature(feature),
         sign  = ifelse(loading >= 0, "Positive", "Negative"),
         label = fct_reorder(label, loading))

print(
  ggplot(plot_df_cyto10, aes(x = label, y = loading, fill = sign)) +
    geom_col(width = 0.75) + coord_flip() +
    scale_fill_manual(values = c("Positive" = COL_RED, "Negative" = COL_BLUE), guide = "none") +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = paste0("MOFA weight (", factor_of_interest, ")"),
         title = "Cytokine/protein loadings \u2013 Factor 10")
)

X_cyto_f10  <- MOFA2::get_data(mofa_trained)[[cyto_view]][["group1"]]
common_samp <- intersect(colnames(X_cyto_f10), names(f10))
X_cyto2     <- X_cyto_f10[, common_samp, drop = FALSE]
f10_use     <- as.numeric(f10[common_samp])

pos_feats <- cyto_loadings_f10 %>% arrange(desc(loading)) %>% slice_head(n = N_cyto) %>%
  pull(feature) %>% intersect(rownames(X_cyto2))
stopifnot(length(pos_feats) >= 5)

score_cyto_f10 <- tibble(
  sample_id       = common_samp,
  Factor10        = f10_use,
  Cyto_pos_module = colMeans(X_cyto2[pos_feats, , drop = FALSE], na.rm = TRUE)
) %>%
  filter(!is.na(Factor10)) %>%
  add_palindromic()

stat_cyto_f10 <- facet_spearman_labels(score_cyto_f10, "Factor10", "Cyto_pos_module")

print(
  ggplot(score_cyto_f10, aes(x = Factor10, y = Cyto_pos_module)) +
    geom_point(alpha = 0.7, size = 1.6) +
    geom_smooth(aes(colour = palindromic, fill = palindromic),
                method = "lm", se = TRUE, alpha = 0.20, linewidth = 0.9) +
    geom_text(data = stat_cyto_f10, aes(x = Inf, y = -Inf, label = label),
              inherit.aes = FALSE, hjust = 1.05, vjust = -0.6, size = 4) +
    facet_wrap(~ palindromic, nrow = 1) +
    scale_colour_manual(values = pal_cols_f10, guide = "none") +
    scale_fill_manual(values = pal_cols_f10, guide = "none") +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(face = "bold", size = 13)) +
    labs(x = "MOFA Factor 10 score",
         y = paste0("Mean abundance of top ", N_cyto, " positive-loading cytokines"),
         title = "Cytokine/protein validation \u2013 Factor 10")
)

# -- 15C. IFN score by palindromic status (Factor 10 context) -----------------
ifn_view <- "Interferon score"
stopifnot(ifn_view %in% names(MOFA2::get_data(mofa_trained)))
X_ifn     <- MOFA2::get_data(mofa_trained)[[ifn_view]][["group1"]]
stopifnot(nrow(X_ifn) == 1)

ifn_df <- tibble(
  sample_id = colnames(X_ifn),
  IFN_score = as.numeric(X_ifn[1, colnames(X_ifn)])
) %>%
  left_join(df %>% select(sample_id, palindromic), by = "sample_id") %>%
  mutate(palindromic = factor(palindromic, levels = c(0, 1),
                              labels = c("Non-palindromic", "Palindromic"))) %>%
  filter(!is.na(IFN_score), !is.na(palindromic))

stopifnot(nrow(ifn_df) >= 10)

wtest_ifn <- wilcox.test(IFN_score ~ palindromic, data = ifn_df, exact = FALSE)
y_max_ifn <- max(ifn_df$IFN_score, na.rm = TRUE)
y_rng_ifn <- diff(range(ifn_df$IFN_score, na.rm = TRUE))

print(
  ggplot(ifn_df, aes(x = palindromic, y = IFN_score)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, fill = "grey90", colour = "black") +
    geom_jitter(aes(colour = palindromic), width = 0.15, alpha = 0.7, size = 1.8) +
    annotate("text", x = 2, y = y_max_ifn + 0.08 * y_rng_ifn,
             label = paste0("Wilcoxon p = ", format.pval(wtest_ifn$p.value, digits = 2, eps = 1e-3)),
             hjust = 1, size = 4) +
    scale_colour_manual(values = pal_cols_f10, guide = "none") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(size = 14)) +
    labs(x = NULL, y = "Interferon score", title = "IFN score by palindromic status")
)

# -- 15D. Rare germline: loadings + burden validation (palindromic facets) ----
germ_view_f10 <- "Rare Germline variants"
stopifnot(germ_view_f10 %in% names(MOFA2::get_data(mofa_trained)))

W_germ_f10  <- MOFA2::get_weights(mofa_trained, views = germ_view_f10, factors = factor_of_interest)[[germ_view_f10]]
w_germ_f10  <- setNames(as.numeric(W_germ_f10[, factor_of_interest, drop = TRUE]), rownames(W_germ_f10))
germ_loadings_f10 <- tibble(feature = names(w_germ_f10), loading = w_germ_f10)

plot_df_germ10 <- germ_loadings_f10 %>%
  arrange(desc(abs(loading))) %>% slice_head(n = 20) %>%
  mutate(label = format_germ_feature(feature),
         sign  = ifelse(loading >= 0, "Positive", "Negative"),
         label = fct_reorder(label, loading))

print(
  ggplot(plot_df_germ10, aes(x = label, y = loading, fill = sign)) +
    geom_col(width = 0.75) + coord_flip() +
    scale_fill_manual(values = c("Negative" = COL_BLUE, "Positive" = COL_RED), guide = "none") +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = paste0("MOFA weight (", factor_of_interest, ")"),
         title = paste0("Rare germline feature loadings \u2013 Factor 10"))
)

# Germline view matrix
X_germ_f10  <- MOFA2::get_data(mofa_trained)[[germ_view_f10]][["group1"]]
common_samp <- intersect(colnames(X_germ_f10), names(f10))
stopifnot(length(common_samp) >= 10)
X_germ2     <- X_germ_f10[, common_samp, drop = FALSE]
f10_use     <- as.numeric(f10[common_samp])

germ_burden_f10 <- colSums(X_germ2, na.rm = TRUE)
df_burden_f10   <- tibble(
  sample_id   = common_samp,
  Factor10    = f10_use,
  Germ_burden = as.numeric(germ_burden_f10)
) %>%
  filter(!is.na(Factor10), !is.na(Germ_burden)) %>%
  add_palindromic()

stat_burden_f10 <- facet_spearman_labels(df_burden_f10, "Factor10", "Germ_burden")

print(
  ggplot(df_burden_f10, aes(x = Germ_burden, y = Factor10)) +
    geom_point(alpha = 0.7, size = 1.7) +
    geom_smooth(aes(colour = palindromic, fill = palindromic),
                method = "lm", se = TRUE, alpha = 0.20, linewidth = 0.9) +
    geom_text(data = stat_burden_f10, aes(x = Inf, y = -Inf, label = label),
              inherit.aes = FALSE, hjust = 1.05, vjust = -0.6, size = 4) +
    facet_wrap(~ palindromic, nrow = 1) +
    scale_colour_manual(values = pal_cols_f10, guide = "none") +
    scale_fill_manual(values = pal_cols_f10, guide = "none") +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(face = "bold", size = 13)) +
    labs(x = "Rare germline variant burden (sum across features)",
         y = paste0(factor_of_interest, " score"),
         title = paste0(factor_of_interest, " vs rare germline burden"))
)

# Top germline feature carrier vs Factor 10 (faceted)
top_feat_f10 <- plot_df_germ10 %>% arrange(desc(abs(loading))) %>% slice_head(n = 1) %>% pull(feature)
stopifnot(top_feat_f10 %in% rownames(X_germ2))

df_top_f10 <- tibble(
  sample_id = common_samp,
  Factor10  = f10_use,
  dosage    = as.numeric(X_germ2[top_feat_f10, common_samp])
) %>%
  filter(!is.na(dosage), !is.na(Factor10)) %>%
  mutate(carrier_status = factor(ifelse(dosage > 0, "Carrier", "Non-carrier"),
                                 levels = c("Non-carrier", "Carrier"))) %>%
  add_palindromic()

stat_carrier_f10 <- df_top_f10 %>%
  group_by(palindromic) %>%
  summarise(p = wilcox.test(Factor10 ~ carrier_status, exact = FALSE)$p.value, .groups = "drop") %>%
  mutate(label = paste0("Wilcoxon p = ", format.pval(p, digits = 2, eps = 1e-3)))

print(
  ggplot(df_top_f10, aes(x = carrier_status, y = Factor10)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, fill = "grey90", colour = "black") +
    geom_jitter(aes(colour = palindromic), width = 0.15, alpha = 0.7, size = 1.8) +
    geom_text(data = stat_carrier_f10, aes(x = Inf, y = -Inf, label = label),
              inherit.aes = FALSE, hjust = 1.05, vjust = -0.6, size = 4) +
    facet_wrap(~ palindromic, nrow = 1) +
    scale_colour_manual(values = pal_cols_f10, guide = "none") +
    theme_bw(base_size = 12) +
    theme(strip.text  = element_text(face = "bold", size = 13),
          axis.text.x = element_text(size = 14)) +
    labs(x = NULL, y = paste0(factor_of_interest, " score"),
         title = "Factor 10 by rare germline carrier status",
         subtitle = format_germ_feature(top_feat_f10))
)
