# =============================================================================
# RA RNAseq Analysis Pipeline
# =============================================================================
# Set base_dir to your working directory before running
base_dir <- "."

# =============================================================================
# 0. Libraries
# =============================================================================
suppressPackageStartupMessages({
  library(tximport)
  library(rtracklayer)
  library(DESeq2)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(scales)
  library(RColorBrewer)
  library(VennDiagram)
  library(ggVennDiagram)
  library(ComplexUpset)
  library(ComplexHeatmap)
  library(circlize)
  library(vsn)
  library(pheatmap)
  library(enrichR)
})

# =============================================================================
# 1. Load input data
# =============================================================================

# Salmon quantification files (one per sample, named by sample ID)
files <- list.files(file.path(base_dir, "quants"), pattern = ".sf",
                    full.names = TRUE, recursive = TRUE)
names(files) <- stringr::str_split(files, pattern = "/", simplify = TRUE)[, 5] %>%
  stringr::str_replace("_quant", "")

# Transcript-to-gene mapping from GTF
gtf    <- rtracklayer::import(file.path(base_dir, "gencode.v46.basic.annotation.gtf"))
gtf_df <- as.data.frame(gtf)

tx2gene       <- gtf_df %>% filter(type == "transcript") %>% select(transcript_id, gene_id)
gene_name_map <- gtf_df %>% filter(type == "gene")       %>% select(gene_id, gene_name)

# Import with tximport
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)

# Experimental metadata
meta_data <- read.csv("RA_meta.csv", row.names = 1)
meta_data$RA        <- as.factor(meta_data$RA)
meta_data$Phenotype <- as.factor(meta_data$Phenotype)
meta_data$Group     <- as.factor(meta_data$Group)
meta_data$Sex       <- as.factor(meta_data$Sex)

# Verify sample order matches
stopifnot(rownames(meta_data) == colnames(txi.salmon$counts))

# Protein-coding gene list (from BioMart)
protein_coding_genes <- read.csv("mart_export.csv")

# =============================================================================
# 2. DESeq2: Model 1 — ~ Sex + Group
# =============================================================================

meta_data$Group <- relevel(meta_data$Group, ref = "HC")
meta_data$Sex   <- relevel(meta_data$Sex,   ref = "M")

dds <- DESeqDataSetFromTximport(txi.salmon, colData = meta_data, design = ~ Sex + Group)

# Low-count filter: keep genes with >= 10 counts in at least 11 samples
smallestGroupSize <- 11
dds <- dds[rowSums(counts(dds) >= 10) >= smallestGroupSize, ]

dds <- DESeq(dds)
resultsNames(dds)

# Save normalised counts
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file = "DESeq2_normalized_counts.csv", quote = FALSE)

# Annotated normalised counts
normcounts_anno <- read.csv("DESeq2_normalized_counts.csv", row.names = 1) %>%
  tibble::rownames_to_column(var = "gene_id") %>%
  left_join(gene_name_map, by = "gene_id")
write.csv(as.data.frame(normcounts_anno), file = "DESeq2_normalized_counts_anno.csv")

# Cook's distance outlier overview
par(mar = c(8, 5, 2, 2))
boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2)

# =============================================================================
# 3. Extract pairwise results and write gene lists (Model 1)
# =============================================================================

# Add gene symbol to a DESeq2 results object
add_gene_name <- function(res) {
  as.data.frame(res) %>%
    tibble::rownames_to_column(var = "gene_id") %>%
    left_join(gene_name_map, by = "gene_id") %>%
    arrange(padj, abs(log2FoldChange))
}

# Filter to protein-coding genes, annotate, sort, and write CSVs.
# Returns a list(full, sig); also writes .FULL.csv and .0.05.csv files.
make_de_list <- function(res, label, protein_coding_genes, lfc_thresh = 0.5, padj_thresh = 0.05) {
  coding <- res[rownames(res) %in% protein_coding_genes$Gene_ID, ]
  coding <- add_gene_name(coding)
  coding <- coding[order(coding$padj), ]
  sig    <- subset(coding, padj < padj_thresh &
                     (log2FoldChange > lfc_thresh | log2FoldChange < -lfc_thresh))
  write.csv(as.data.frame(sig),    file = paste0(label, ".0.05.csv"))
  write.csv(as.data.frame(coding), file = paste0(label, ".FULL.csv"))
  invisible(list(full = coding, sig = sig))
}

# vs HC contrasts
AR_vs_HC      <- results(dds, name = "Group_AR_vs_HC")
PR_vs_HC      <- results(dds, name = "Group_PR_vs_HC")
PROG_AR_vs_HC <- results(dds, name = "Group_Progressor_AR_vs_HC")
PROG_PR_vs_HC <- results(dds, name = "Group_Progressor_PR_vs_HC")

summary(AR_vs_HC,      alpha = 0.05)
summary(PR_vs_HC,      alpha = 0.05)
summary(PROG_AR_vs_HC, alpha = 0.05)
summary(PROG_PR_vs_HC, alpha = 0.05)

de1 <- make_de_list(AR_vs_HC,      "AR_vs_HC_coding",      protein_coding_genes)
de2 <- make_de_list(PR_vs_HC,      "PR_vs_HC_coding",      protein_coding_genes)
de3 <- make_de_list(PROG_AR_vs_HC, "PROG_AR_vs_HC_coding", protein_coding_genes)
de4 <- make_de_list(PROG_PR_vs_HC, "PROG_PR_vs_HC_coding", protein_coding_genes)

AR_vs_HC_coding        <- de1$full; AR_vs_HC_coding.0.05        <- de1$sig
PR_vs_HC_coding        <- de2$full; PR_vs_HC_coding.0.05        <- de2$sig
PROG_AR_vs_HC_coding   <- de3$full; PROG_AR_vs_HC_coding.0.05   <- de3$sig
PROG_PR_vs_HC_coding   <- de4$full; PROG_PR_vs_HC_coding.0.05   <- de4$sig

# Progressor vs non-progressor contrasts
PROG_AR_vs_AR <- results(dds, contrast = c("Group", "Progressor_AR", "AR"))
PROG_PR_vs_PR <- results(dds, contrast = c("Group", "Progressor_PR", "PR"))

summary(PROG_AR_vs_AR, alpha = 0.05)
summary(PROG_PR_vs_PR, alpha = 0.05)

de5 <- make_de_list(PROG_AR_vs_AR, "PROG_AR_vs_AR_coding", protein_coding_genes)
de6 <- make_de_list(PROG_PR_vs_PR, "PROG_PR_vs_PR_coding", protein_coding_genes)

PROG_AR_vs_AR_coding <- de5$full; PROG_AR_vs_AR_coding.0.05 <- de5$sig
PROG_PR_vs_PR_coding <- de6$full; PROG_PR_vs_PR_coding.0.05 <- de6$sig

# =============================================================================
# 4. DESeq2: Model 2 — ~ Sex + AR_vs_PR (collapsed groups)
# =============================================================================

meta_data$AR_vs_PR <- case_when(
  meta_data$Group %in% c("AR", "Progressor_AR") ~ "AR_all",
  meta_data$Group %in% c("PR", "Progressor_PR") ~ "PR_all",
  TRUE ~ "HC"
)
meta_data$AR_vs_PR <- factor(meta_data$AR_vs_PR, levels = c("HC", "AR_all", "PR_all"))

dds2 <- DESeqDataSetFromTximport(txi.salmon, colData = meta_data, design = ~ Sex + AR_vs_PR)
dds2 <- dds2[rowSums(counts(dds2) >= 10) >= smallestGroupSize, ]
dds2 <- DESeq(dds2)
resultsNames(dds2)

AR_all_vs_PR_all <- results(dds2, contrast = c("AR_vs_PR", "AR_all", "PR_all"))
summary(AR_all_vs_PR_all, alpha = 0.05)

de7 <- make_de_list(AR_all_vs_PR_all, "AR_all_vs_PR_all_coding", protein_coding_genes)
AR_all_vs_PR_all_coding      <- de7$full
AR_all_vs_PR_all_coding.0.05 <- de7$sig

# =============================================================================
# 5. PCA (sex-corrected)
# =============================================================================

vsd  <- vst(dds, blind = TRUE)
mat  <- assay(vsd)
mat2 <- removeBatchEffect(mat, batch = colData(dds)$Sex)

pca2       <- prcomp(t(mat2), center = TRUE, scale. = FALSE)
percentVar <- pca2$sdev^2 / sum(pca2$sdev^2)

pca_df2           <- as.data.frame(pca2$x)
pca_df2$Group     <- colData(dds)$Group
pca_df2$Sex       <- colData(dds)$Sex
pca_df2$Sample    <- rownames(pca_df2)
pca_df2$RA        <- colData(dds)$RA
pca_df2$Phenotype <- colData(dds)$Phenotype

pca_df2 <- pca_df2 %>%
  mutate(Group_label = recode(Phenotype,
                              "AR" = "ACPA+ arthralgia",
                              "PR" = "ACPA+ PR",
                              "HC" = "Control"))

my_group_cols <- c(
  "Control"          = "red",
  "ACPA+ arthralgia" = "orange",
  "ACPA+ PR"         = "blue"
)

ggplot(pca_df2, aes(PC1, PC2)) +
  stat_ellipse(aes(fill = Group_label), geom = "polygon",
               level = 0.95, alpha = 0.3, colour = NA) +
  geom_point(aes(fill = Group_label), shape = 21, size = 2, stroke = 0.5) +
  scale_fill_manual(values = my_group_cols) +
  labs(x = sprintf("PC1: %.1f%% var", percentVar[1] * 100),
       y = sprintf("PC2: %.1f%% var", percentVar[2] * 100)) +
  theme_grey(base_size = 14) +
  theme(panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white", linewidth = 0.3))

# Mean-SD plot to verify VST stabilisation
meanSdPlot(assay(vsd))

# =============================================================================
# 6. Volcano plots
# =============================================================================

# Annotate DE direction and gene labels, return a plotting data frame
prep_volcano_df <- function(res_df, lfc_thresh = 0.5, padj_thresh = 0.05) {
  res_df$diffexpressed <- "NO"
  res_df$diffexpressed[res_df$log2FoldChange >  lfc_thresh & res_df$padj < padj_thresh] <- "UP"
  res_df$diffexpressed[res_df$log2FoldChange < -lfc_thresh & res_df$padj < padj_thresh] <- "DOWN"
  res_df$delabel <- NA
  res_df$delabel[res_df$diffexpressed != "NO"] <- res_df$gene_name[res_df$diffexpressed != "NO"]
  data.frame(
    log2FC        = res_df$log2FoldChange,
    logpv         = -log10(res_df$padj),
    diffexpressed = res_df$diffexpressed,
    delabel       = res_df$delabel
  )
}

# Standard volcano ggplot
plot_volcano <- function(df, title, xlim = c(-6, 6), ylim = c(0, 30)) {
  ggplot(df, aes(x = log2FC, y = logpv, col = diffexpressed, label = delabel)) +
    geom_point() +
    scale_color_manual(values = c("DOWN" = "blue", "NO" = "darkgrey", "UP" = "red")) +
    ggtitle(title) +
    xlab(expression("log2 fold change")) +
    xlim(xlim) +
    ylab(expression("-log10 adjusted p-value")) +
    ylim(ylim) +
    geom_hline(yintercept = 1.3,          col = "black", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = c(-0.5, 0.5), col = "black", linetype = "dashed", linewidth = 1) +
    geom_text_repel(max.overlaps = 10, size = 5) +
    theme_minimal() +
    theme(
      plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.line.x     = element_line(color = "black", size = 1),
      axis.line.y     = element_line(color = "black", size = 1),
      axis.title      = element_text(color = "black", size = 12, face = "bold"),
      axis.text       = element_text(color = "black", size = 10),
      legend.position = "none"
    )
}

# Helper: load a FULL CSV, annotate, and plot
make_volcano <- function(csv_file, title, xlim = c(-6, 6), ylim = c(0, 30)) {
  df <- read.csv(csv_file, row.names = 1)
  plot_volcano(prep_volcano_df(df), title, xlim, ylim)
}

make_volcano("PROG_AR_vs_AR_coding.FULL.csv",
             "ACPA+ arthralgia progressors vs non-progressors", xlim = c(-3, 3), ylim = c(0, 5))
make_volcano("PROG_PR_vs_PR_coding.FULL.csv",
             "ACPA+ PR progressors vs non-progressors",         xlim = c(-3, 3), ylim = c(0, 5))
make_volcano("AR_vs_HC_coding.FULL.csv",
             "ACPA+ arthralgia non-progressors vs HC")
make_volcano("PR_vs_HC_coding.FULL.csv",
             "ACPA+ PR non-progressors vs HC")
make_volcano("PROG_AR_vs_HC_coding.FULL.csv",
             "ACPA+ arthralgia progressors vs HC")
make_volcano("PROG_PR_vs_HC_coding.FULL.csv",
             "ACPA+ PR progressors vs HC")
make_volcano("AR_all_vs_PR_all_coding.FULL.csv",
             "ACPA+ arthralgia vs ACPA+ PR",                    xlim = c(-3, 3), ylim = c(0, 5))

# =============================================================================
# 7. Venn diagrams
# =============================================================================

venn_list_up <- list(
  "ACPA+ arthralgia non-progressors" = AR_vs_HC_coding.0.05$gene_name[AR_vs_HC_coding.0.05$log2FoldChange > 0.5],
  "ACPA+ PR non-progressors"         = PR_vs_HC_coding.0.05$gene_name[PR_vs_HC_coding.0.05$log2FoldChange > 0.5],
  "ACPA+ arthralgia progressors"     = PROG_AR_vs_HC_coding.0.05$gene_name[PROG_AR_vs_HC_coding.0.05$log2FoldChange > 0.5],
  "ACPA+ PR progressors"             = PROG_PR_vs_HC_coding.0.05$gene_name[PROG_PR_vs_HC_coding.0.05$log2FoldChange > 0.5]
)

venn_list_down <- list(
  AR      = AR_vs_HC_coding.0.05$gene_name[AR_vs_HC_coding.0.05$log2FoldChange < -0.5],
  PR      = PR_vs_HC_coding.0.05$gene_name[PR_vs_HC_coding.0.05$log2FoldChange < -0.5],
  AR_PROG = PROG_AR_vs_HC_coding.0.05$gene_name[PROG_AR_vs_HC_coding.0.05$log2FoldChange < -0.5],
  PR_PROG = PROG_PR_vs_HC_coding.0.05$gene_name[PROG_PR_vs_HC_coding.0.05$log2FoldChange < -0.5]
)

tiff("RA_RNAseq_Venn.tiff",      units = "in", width = 8, height = 5, res = 300)
ggVennDiagram(venn_list_up,   label_alpha = 0,
              category.names = c("AR", "PR", "AR Prog", "PR Prog")) +
  ggplot2::scale_fill_gradient(low = "white", high = "red")
dev.off()

tiff("RA_RNAseq_Venn_down.tiff", units = "in", width = 8, height = 5, res = 300)
ggVennDiagram(venn_list_down, label_alpha = 0,
              category.names = c("AR", "PR", "AR Prog", "PR Prog")) +
  ggplot2::scale_fill_gradient(low = "white", high = "cyan")
dev.off()

# Overlap gene sets
central_overlap_genes_up   <- Reduce(intersect, venn_list_up)
central_overlap_genes_down <- Reduce(intersect, venn_list_down)

AR_PROG_unique_up <- setdiff(
  venn_list_up[["ACPA+ arthralgia progressors"]],
  c(venn_list_up[["ACPA+ arthralgia non-progressors"]],
    venn_list_up[["ACPA+ PR non-progressors"]],
    venn_list_up[["ACPA+ PR progressors"]])
)

AR_PROG_PR_PROG_overlap_up <- setdiff(
  intersect(venn_list_up[["ACPA+ arthralgia progressors"]],
            venn_list_up[["ACPA+ PR progressors"]]),
  c(venn_list_up[["ACPA+ arthralgia non-progressors"]],
    venn_list_up[["ACPA+ PR non-progressors"]])
)

write.csv(central_overlap_genes_up,   "central_overlap_genes_up.csv",   row.names = FALSE)
write.csv(central_overlap_genes_down, "central_overlap_genes_down.csv", row.names = FALSE)
write.csv(AR_PROG_unique_up,          "AR_PROG_unique_up.csv",          row.names = FALSE)
write.csv(AR_PROG_PR_PROG_overlap_up, "AR_PROG_PR_PROG_overlap_up.csv", row.names = FALSE)

# =============================================================================
# 8. UpSet plot (up- and down-regulated genes across comparisons)
# =============================================================================

up_list <- list(
  AR_vs_HC      = AR_vs_HC_coding.0.05      %>% filter(log2FoldChange > 0) %>% pull(gene_name),
  PR_vs_HC      = PR_vs_HC_coding.0.05      %>% filter(log2FoldChange > 0) %>% pull(gene_name),
  AR_PROG_vs_HC = PROG_AR_vs_HC_coding.0.05 %>% filter(log2FoldChange > 0) %>% pull(gene_name),
  PR_PROG_vs_HC = PROG_PR_vs_HC_coding.0.05 %>% filter(log2FoldChange > 0) %>% pull(gene_name)
)

down_list <- list(
  AR_vs_HC      = AR_vs_HC_coding.0.05      %>% filter(log2FoldChange < 0) %>% pull(gene_name),
  PR_vs_HC      = PR_vs_HC_coding.0.05      %>% filter(log2FoldChange < 0) %>% pull(gene_name),
  AR_PROG_vs_HC = PROG_AR_vs_HC_coding.0.05 %>% filter(log2FoldChange < 0) %>% pull(gene_name),
  PR_PROG_vs_HC = PROG_PR_vs_HC_coding.0.05 %>% filter(log2FoldChange < 0) %>% pull(gene_name)
)

df_wide_reg <- bind_rows(
  enframe(up_list,   name = "comparison", value = "gene") %>% unnest(gene) %>% mutate(Regulation = "Up"),
  enframe(down_list, name = "comparison", value = "gene") %>% unnest(gene) %>% mutate(Regulation = "Down")
) %>%
  mutate(Regulation = factor(Regulation, levels = c("Up", "Down"))) %>%
  distinct(gene, Regulation, comparison) %>%
  mutate(present = TRUE) %>%
  pivot_wider(id_cols = c(gene, Regulation), names_from = comparison,
              values_from = present, values_fill = list(present = FALSE))

upset(
  df_wide_reg,
  intersect        = names(up_list),
  name             = "comparison",
  base_annotations = list(
    'Intersection size' = intersection_size(
      mapping = aes(fill = Regulation),
      text    = list(aes(label = after_stat(count)))
    ) +
      scale_fill_manual(name = "Regulation", values = c(Up = "red", Down = "blue")) +
      theme(legend.position = "right", legend.title = element_text(size = 10))
  )
) + theme_minimal()

# =============================================================================
# 9. Pathway enrichment (enrichR)
# =============================================================================

# Run enrichR on a gene vector, tidy results, and plot a dotplot of top n terms
run_enrichr_dotplot <- function(genes, dbs, db_name, top_n = 10,
                                title = NULL, point_color_high = "red") {
  enriched <- enrichr(genes, dbs)
  res <- enriched[[db_name]] %>%
    separate(Overlap, into = c("OverlapCount", "TermSize"), sep = "/", convert = TRUE) %>%
    mutate(GeneRatio = OverlapCount / TermSize,
           logAdjP   = -log10(Adjusted.P.value)) %>%
    arrange(Adjusted.P.value) %>%
    slice_head(n = top_n)

  p <- ggplot(res, aes(x = logAdjP, y = reorder(Term, logAdjP),
                       size = OverlapCount, color = Combined.Score)) +
    geom_point() +
    scale_color_gradient(low = "grey", high = point_color_high) +
    labs(x = expression(-log[10] ~ "(adj. P-value)"), y = NULL,
         size = "Gene Count", color = "Combined Score", title = title) +
    theme_minimal() +
    theme(axis.text.y  = element_text(size = 11, face = "bold"),
          axis.title   = element_text(size = 12, face = "bold"),
          plot.title   = element_text(size = 14, face = "bold"))
  print(p)
  invisible(enriched)
}

enrichr_dbs <- c("GO_Biological_Process_2023", "GO_Molecular_Function_2023",
                 "GO_Cellular_Component_2023", "MSigDB_Hallmark_2020",
                 "Reactome_Pathways_2024")

upgenes_shared   <- read.csv("central_overlap_genes_up.csv")
run_enrichr_dotplot(upgenes_shared,   enrichr_dbs, "MSigDB_Hallmark_2020",
                    title = "Shared upregulated — Hallmark",   point_color_high = "red")

downgenes_shared <- read.csv("central_overlap_genes_down.csv")
run_enrichr_dotplot(downgenes_shared, enrichr_dbs, "MSigDB_Hallmark_2020",
                    title = "Shared downregulated — Hallmark", point_color_high = "blue")

prog_shared_up   <- read.csv("AR_PROG_PR_PROG_overlap_up.csv")
run_enrichr_dotplot(prog_shared_up,   enrichr_dbs, "Reactome_Pathways_2024",
                    title = "Progressor shared up — Reactome", point_color_high = "red")

prog_shared_down <- read.csv("AR_PROG_PR_PROG_overlap_down.csv")
run_enrichr_dotplot(prog_shared_down, enrichr_dbs, "Reactome_Pathways_2024",
                    title = "Progressor shared down — Reactome", point_color_high = "blue")

shared_Prog_up      <- read.csv("prog_shared_up.csv")
run_enrichr_dotplot(shared_Prog_up,      enrichr_dbs, "MSigDB_Hallmark_2020",
                    title = "AR_PROG unique up — Hallmark",    point_color_high = "red")

enriched_overlap_up <- read.csv("AR_PROG_PR_PROG_overlap_up.csv")
run_enrichr_dotplot(enriched_overlap_up, enrichr_dbs, "MSigDB_Hallmark_2020",
                    title = "AR_PROG + PR_PROG overlap up — Hallmark", point_color_high = "red")

# =============================================================================
# 10. Gene-level count plots (individual genes of interest)
# =============================================================================

# plotCounts boxplot for a single gene
plot_counts_box <- function(dds, ensg_id, gene_label,
                            intgroup = "Group",
                            group_colors = c(HC = "darkgrey", AR = "yellow", PR = "blue",
                                             Progressor_AR = "darkgreen", Progressor_PR = "purple")) {
  dat <- plotCounts(dds, gene = ensg_id, intgroup = intgroup, returnData = TRUE)
  ggplot(dat, aes(.data[[intgroup]], count, fill = .data[[intgroup]])) +
    geom_boxplot() +
    scale_fill_manual(name = intgroup, values = group_colors) +
    ggtitle(gene_label) +
    ylab("Normalised Counts") +
    theme_minimal() +
    theme(axis.title.x    = element_blank(),
          axis.text.x     = element_blank(),
          axis.ticks.x    = element_blank(),
          legend.position = "right")
}

genes_of_interest <- c(
  IFITM1 = "ENSG00000185885.17",
  IFITM2 = "ENSG00000185201.18",
  IFITM3 = "ENSG00000142089.17",
  IRF1   = "ENSG00000125347.15",
  IRF3   = "ENSG00000126456.16",
  IRF9   = "ENSG00000213928.10",
  IFIT1B = "ENSG00000204010.3"
)

for (gene_label in names(genes_of_interest)) {
  print(plot_counts_box(dds, genes_of_interest[[gene_label]], gene_label))
}

# =============================================================================
# 11. Heatmaps (ComplexHeatmap)
# =============================================================================

# Z-score normalisation (row-wise helper used inside vst_heatmap)
z_score_normalisation <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

# Build a ComplexHeatmap from a gene list CSV (columns: Gene_stable_ID_version, Gene_name).
# Pass outlier = "SAMPLE_ID" to drop one sample, or outlier = NULL to keep all.
vst_heatmap <- function(vsd, meta_data, gene_list_csv, out_tiff,
                        width = 8, height = 6, fontsize_rows = 10, outlier = NULL) {
  gene_list   <- read.csv(gene_list_csv, header = TRUE)
  ensembl_ids <- gene_list$Gene_stable_ID_version
  gene_names  <- gene_list$Gene_name
  selected    <- ensembl_ids[ensembl_ids %in% rownames(assay(vsd))]
  mat         <- assay(vsd)[selected, ]

  if (!is.null(outlier)) {
    mat       <- mat[, colnames(mat) != outlier, drop = FALSE]
    meta_data <- meta_data[rownames(meta_data) != outlier, , drop = FALSE]
  }

  mat           <- t(apply(mat, 1, z_score_normalisation))
  rownames(mat) <- make.unique(gene_names[match(selected, ensembl_ids)])

  ordered_samples          <- rownames(meta_data)[order(meta_data$Group)]
  mat                      <- mat[, ordered_samples]
  annotation_col           <- data.frame(Group = meta_data[ordered_samples, "Group"])
  rownames(annotation_col) <- ordered_samples

  ann_colors <- list(Group = c(HC = "red", AR = "yellow", Progressor_AR = "blue",
                                PR = "cyan", Progressor_PR = "purple"))
  ha      <- HeatmapAnnotation(df = annotation_col, col = ann_colors)
  col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

  tiff(out_tiff, units = "in", width = width, height = height, res = 300)
  draw(Heatmap(mat,
               name                 = "Expression (Z-score)",
               col                  = col_fun,
               top_annotation       = ha,
               cluster_rows         = TRUE,
               cluster_columns      = TRUE,
               show_row_names       = TRUE,
               show_column_names    = FALSE,
               row_names_gp         = gpar(fontsize = fontsize_rows),
               heatmap_legend_param = list(title = "Z-score"),
               row_dend_width       = unit(2, "cm"),
               column_dend_height   = unit(2, "cm")))
  dev.off()
}

vsd <- vst(dds, blind = FALSE)

# Heatmap 1: IFN overlap genes (with outlier removed; set outlier to sample ID string or NULL)
vst_heatmap(vsd, meta_data,
            gene_list_csv = "enrichment_overlap_genes_IFN_Comp.csv",
            out_tiff      = "overlap_heatmap.tiff",
            width = 8, height = 6, fontsize_rows = 10,
            outlier = NULL)

# Heatmap 2: IFN overlap genes (no outlier removal, smaller font for denser matrix)
vst_heatmap(vsd, meta_data,
            gene_list_csv = "enrichment_overlap_genes_IFN_Comp.csv",
            out_tiff      = "overlap_heatmap_IFN_Comp.tiff",
            width = 8, height = 6, fontsize_rows = 8,
            outlier = NULL)

# =============================================================================
# 12. Sample-to-sample distance heatmap
# =============================================================================

vsd       <- vst(dds, blind = FALSE)
res_de    <- results(dds, alpha = 0.05)
sig_genes <- rownames(subset(res_de, padj < 0.05))

vsd_mat_de          <- assay(vsd)[sig_genes, ]
sampleDists_de      <- dist(t(vsd_mat_de))
sampleDistMatrix_de <- as.matrix(sampleDists_de)
rownames(sampleDistMatrix_de) <- paste(vsd$Group, vsd$type, sep = "-")

pheatmap(sampleDistMatrix_de,
         clustering_distance_rows = sampleDists_de,
         clustering_distance_cols = sampleDists_de,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255))

# =============================================================================
# 13. IFN score and filtered gene count exports
# =============================================================================

genecounts <- read.csv("salmon_counts_genes.csv", stringsAsFactors = FALSE)

IFN_scoring_genes        <- read.csv("IFNScoreGenes.csv", stringsAsFactors = FALSE, header = TRUE)[[1]]
overlap_and_unique_genes <- read.csv("overlap_vs_HC.csv", stringsAsFactors = FALSE, header = TRUE)[[1]]
tnfsig_genes             <- read.csv("tnfsig.csv",        stringsAsFactors = FALSE, header = TRUE)[[1]]

write.csv(genecounts %>% filter(gene_name %in% IFN_scoring_genes),        "IFNScoringGenes_output.csv",          row.names = FALSE)
write.csv(genecounts %>% filter(gene_name %in% overlap_and_unique_genes), "overlap_and_unique_genes_output.csv", row.names = FALSE)
write.csv(genecounts %>% filter(gene_name %in% tnfsig_genes),             "tnfsig_genes_output.csv",             row.names = FALSE)

# =============================================================================
# 14. MOFA export: VST z-scored top-5000 variable genes (gene symbols)
# =============================================================================

vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)

map <- gene_name_map %>%
  distinct(gene_id, gene_name) %>%
  mutate(gene_id = as.character(gene_id), gene_name = as.character(gene_name))

symbols <- map$gene_name[match(rownames(mat), map$gene_id)]
keep    <- !is.na(symbols) & symbols != ""
mat     <- mat[keep, , drop = FALSE]
symbols <- symbols[keep]

# For duplicate symbols, keep the ENSG row with highest variance
gene_var <- apply(mat, 1, var, na.rm = TRUE)
best_idx <- tapply(seq_along(symbols), symbols, function(ix) ix[which.max(gene_var[ix])])
mat_sym  <- mat[unlist(best_idx), , drop = FALSE]
rownames(mat_sym) <- names(best_idx)

# Top 5,000 most variable symbols
top_n    <- 5000
sym_var  <- apply(mat_sym, 1, var, na.rm = TRUE)
keep_sym <- names(sort(sym_var, decreasing = TRUE))[1:min(top_n, length(sym_var))]
mat_sym  <- mat_sym[keep_sym, , drop = FALSE]

# Z-score per gene
mat_z             <- t(scale(t(mat_sym)))
mat_z[is.na(mat_z)] <- 0

write.csv(mat_z, file = "MOFA_RNAseq_VST_z_5000_SYMBOL_geneXsample.csv", quote = FALSE)
write.csv(data.frame(gene_symbol = rownames(mat_z)),
          "MOFA_RNAseq_5000_SYMBOL_list.csv", row.names = FALSE, quote = FALSE)
