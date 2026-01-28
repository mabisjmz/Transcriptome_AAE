############################################################
# Title: RNA-seq analysis pipeline for Aedes aegypti
# Populations: New Orleans & San Nicolás
# Treatments: Control, Broflanilide, Imidacloprid
#
# Main tasks:
#   1. Import featureCounts outputs and build count matrix
#   2. Create metadata table and DESeq2 object
#   3. Perform QC (PCA and sample distance heatmap)
#   4. Identify top variable genes and plot heatmap
#   5. Differential expression for 7 biologically relevant contrasts
#   6. Gene Ontology enrichment using NCBI gene2go + clusterProfiler
#   7. Generate main figures (Figs 1–8) and Supplementary Data S1–S3
############################################################

###############################
# 0. Load packages and set working directory
###############################

# Core differential expression and data handling
library(DESeq2)
library(tidyverse)    # dplyr, readr, tibble, stringr, purrr, ggplot2, forcats

# QC and visualization
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(matrixStats)
library(patchwork)

# GO enrichment
library(data.table)
library(clusterProfiler)

# Export to Excel for supplementary material
library(writexl)

# Working directory (change if necessary)
setwd("/media/lfme/Datos/trans_manus")

# Create output folders (if they do not exist)
dir.create("counts_matriz", showWarnings = FALSE)
dir.create("qc",            showWarnings = FALSE)
dir.create("heatmaps",      showWarnings = FALSE)
dir.create("degs",          showWarnings = FALSE)
dir.create("results/GO",    recursive    = TRUE, showWarnings = FALSE)
dir.create("figures_final", showWarnings = FALSE)
dir.create("supplementary", showWarnings = FALSE)

############################################################
# 1. Import featureCounts data and build count matrix
############################################################

# List all featureCounts files with pattern *_counts.txt
files <- list.files(pattern = "_counts\\.txt$")
stopifnot(length(files) > 0)  # basic check that files exist

# Sample names are the filenames without suffix "_counts.txt"
samples <- stringr::str_remove(files, "_counts\\.txt$")

# Helper function: read a single featureCounts file
# Assumes featureCounts format: first column = gene ID, last column = raw counts
read_count_file <- function(file, sample_name) {
  df <- readr::read_tsv(file, comment = "#", show_col_types = FALSE)
  df <- df %>%
    dplyr::select(
      GeneID = 1,                 # first column: gene identifier
      count  = dplyr::last_col()  # last column: counts for this sample
    )
  colnames(df)[colnames(df) == "count"] <- sample_name
  df
}

# Read all count files into a list
counts_list <- purrr::map2(files, samples, read_count_file)

# Merge all by GeneID (full join) and sort by GeneID
counts_tbl <- purrr::reduce(counts_list, dplyr::full_join, by = "GeneID") %>%
  dplyr::arrange(GeneID)

# Replace any missing values (genes absent in some samples) with 0
counts_tbl[is.na(counts_tbl)] <- 0

# Save full count table to disk
readr::write_tsv(counts_tbl,
                 file.path("counts_matriz", "matriz_counts_todos.txt"))

# Convert to numeric matrix (genes x samples) for DESeq2
counts_mat <- counts_tbl %>%
  tibble::column_to_rownames("GeneID") %>%
  as.matrix()

mode(counts_mat) <- "numeric"

cat("Count matrix dimensions (genes x samples):\n")
print(dim(counts_mat))

############################################################
# 2. Build sample metadata (sample_info / colData)
############################################################

# Experimental design:
#   CNO = New Orleans, Control
#   NOB = New Orleans, Broflanilide
#   NOI = New Orleans, Imidacloprid
#   CSN = San Nicolás, Control
#   SNB = San Nicolás, Broflanilide
#   SNI = San Nicolás, Imidacloprid

sample_info <- tibble::tibble(
  sample    = colnames(counts_mat),
  condition = stringr::str_extract(sample, "^[A-Z]{3}"),
  replicate = as.integer(stringr::str_extract(sample, "[0-9]+$"))
) %>%
  dplyr::mutate(
    population = dplyr::case_when(
      condition %in% c("CNO", "NOB", "NOI") ~ "NewOrleans",
      condition %in% c("CSN", "SNB", "SNI") ~ "SanNicolas",
      TRUE ~ NA_character_
    ),
    treatment = dplyr::case_when(
      condition %in% c("CNO", "CSN") ~ "Control",
      condition %in% c("NOB", "SNB") ~ "Broflanilide",
      condition %in% c("NOI", "SNI") ~ "Imidacloprid",
      TRUE ~ NA_character_
    ),
    # Set factor levels to control order in plots and tables
    condition  = factor(condition,
                        levels = c("CNO", "NOB", "NOI", "CSN", "SNB", "SNI")),
    population = factor(population,
                        levels = c("NewOrleans", "SanNicolas")),
    treatment  = factor(treatment,
                        levels = c("Control", "Broflanilide", "Imidacloprid"))
  )

# Convert to data.frame for DESeq2 and set rownames
coldata <- as.data.frame(sample_info)
rownames(coldata) <- coldata$sample

# Save metadata table used for DESeq2
readr::write_tsv(tibble::as_tibble(coldata),
                 file.path("counts_matriz", "sample_info.txt"))

cat("Sample metadata (first rows):\n")
print(head(coldata))

############################################################
# 3. Create DESeq2 object, filter low counts, run DESeq
############################################################

# Create DESeqDataSet with design ~ condition
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = coldata,
  design    = ~ condition
)

# Filter out genes with very low total counts across all samples
dds <- dds[rowSums(DESeq2::counts(dds)) > 10, ]

# Run DESeq2 pipeline (size factors, dispersion, Wald tests)
dds <- DESeq2::DESeq(dds)

############################################################
# 4. Variance-stabilizing transformation and PCA
############################################################

# Variance-stabilizing transformation (used for distance and heatmaps)
vsd <- DESeq2::vst(dds, blind = FALSE)

# Order samples: population (NO first), then treatment, then replicate
sample_info_ordered <- sample_info %>%
  dplyr::arrange(population, treatment, replicate)

sample_order <- sample_info_ordered$sample

# Reorder transformed matrix and colData accordingly
vsd <- vsd[, sample_order]
coldata_ordered <- coldata[sample_order, , drop = FALSE]

# PCA by condition (saved as QC figure, not main figure)
pca_plot <- DESeq2::plotPCA(vsd, intgroup = "condition") +
  ggplot2::ggtitle("PCA by condition")

ggplot2::ggsave(
  filename = file.path("qc", "PCA_condition.png"),
  plot     = pca_plot,
  width    = 6,
  height   = 6,
  dpi      = 300
)

############################################################
# 5. Sample-to-sample distance heatmap (Figure 1)
############################################################

# Euclidean distance between samples on VST-transformed counts
sample_dist <- dist(t(SummarizedExperiment::assay(vsd)))
sample_dist_mat <- as.matrix(sample_dist)

# Reorder rows and columns by sample_order
sample_dist_mat <- sample_dist_mat[sample_order, sample_order]

# Create "nice" labels for controls:
#   CNO1–3 -> NOC1–3 (New Orleans Control)
#   CSN1–3 -> SNC1–3 (San Nicolás Control)
# Other samples (NOB*, NOI*, SNB*, SNI*) keep original sample ID
sample_info_ordered <- sample_info_ordered %>%
  dplyr::mutate(
    sample_label = dplyr::case_when(
      condition == "CNO" ~ paste0("NOC", replicate),
      condition == "CSN" ~ paste0("SNC", replicate),
      TRUE               ~ sample
    )
  )

labels_vec <- sample_info_ordered$sample_label

# Set row/column names of distance matrix to the pretty labels
rownames(sample_dist_mat) <- labels_vec
colnames(sample_dist_mat) <- labels_vec

# Build annotation data frame for heatmap (rows and columns share info)
ann_df <- sample_info_ordered %>%
  dplyr::select(sample_label, treatment, population, condition) %>%
  tibble::column_to_rownames("sample_label")

# Colors for annotations
pop_colors <- c(
  NewOrleans = "#0072B2",
  SanNicolas = "#E69F00"
)

treat_colors <- c(
  Control      = "#7F7F7F",
  Broflanilide = "#009E73",
  Imidacloprid = "#D55E00"
)

cond_colors <- c(
  CNO = "#B3CDE3",
  NOB = "#8DD3C7",
  NOI = "#BC80BD",
  CSN = "#FED9A6",
  SNB = "#CCEBC5",
  SNI = "#FB8072"
)

annotation_colors <- list(
  population = pop_colors,
  treatment  = treat_colors,
  condition  = cond_colors
)

# Gap index between populations (after last New Orleans sample)
n_NO <- sum(sample_info_ordered$population == "NewOrleans")
gap_index <- n_NO

# Legend breaks for colorbar (distance scale)
dist_range    <- range(sample_dist_mat, na.rm = TRUE)
legend_breaks <- pretty(dist_range, n = 4)

# Generate and save Figure 1
pheatmap::pheatmap(
  sample_dist_mat,
  annotation_col    = ann_df,
  annotation_row    = ann_df,
  annotation_colors = annotation_colors,
  display_numbers   = TRUE,
  number_format     = "%.2f",
  number_color      = "black",
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  gaps_row          = gap_index,
  gaps_col          = gap_index,
  main              = "Sample-to-sample Euclidean distance (vst)",
  fontsize          = 9,
  fontsize_row      = 8,
  fontsize_col      = 8,
  legend_breaks     = legend_breaks,
  legend_labels     = legend_breaks,
  filename          = file.path("figures_final",
                                "Fig1_heatmap_sample_distances.png"),
  width             = 14,
  height            = 10,
  dpi               = 300
)

############################################################
# 6. Heatmap of top 10 most variable genes (Figure 3)
############################################################

# Identify top 10 genes with highest variance across samples (VST counts)
top10_genes <- head(order(matrixStats::rowVars(SummarizedExperiment::assay(vsd)),
                          decreasing = TRUE), 10)

# Extract VST values for these genes and reorder columns
mat_top10 <- SummarizedExperiment::assay(vsd)[top10_genes, sample_order]

# Center each row (gene) by subtracting its mean (z-score like)
mat_top10_centered <- mat_top10 - rowMeans(mat_top10)

# Use the same labels as in the distance heatmap
colnames(mat_top10_centered) <- labels_vec

# Color palette for expression heatmap
heat_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100)

# Legend breaks for centered expression scale
val_range     <- range(mat_top10_centered, na.rm = TRUE)
legend_breaks <- pretty(val_range, n = 5)

# Generate and save Figure 3
pheatmap::pheatmap(
  mat_top10_centered,
  color             = heat_colors,
  annotation_col    = ann_df,
  annotation_colors = annotation_colors,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  gaps_col          = gap_index,   # visual break between populations
  display_numbers   = TRUE,
  number_format     = "%.1f",
  number_color      = "black",
  fontsize          = 9,
  fontsize_row      = 8,
  fontsize_col      = 8,
  main              = "Top 10 most variable genes (vst, row-centered log2 expression)",
  legend_breaks     = legend_breaks,
  legend_labels     = legend_breaks,
  filename          = file.path("figures_final",
                                "Fig3_top10_variable_genes_heatmap.png"),
  width             = 10,
  height            = 6,
  dpi               = 300
)

############################################################
# 7. Differential expression contrasts (DESeq2)
############################################################

# Table to collect summary statistics for each contrast
summary_degs <- tibble::tibble()

# Helper function:
#   - Run DESeq2 contrast
#   - Export full and significant results
#   - Export list of significant genes (for GO)
#   - Update global summary table
run_contrast <- function(dds,
                         group1,
                         group2,
                         label,
                         alpha = 0.05,
                         lfc   = 1,
                         outdir = "degs") {
  
  # DESeq2 contrast: group1 vs group2 (log2FC > 0 = up in group1)
  res <- DESeq2::results(dds, contrast = c("condition", group1, group2))
  
  # Convert to data frame and add GeneID column
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("GeneID")
  
  # Full results (all genes)
  full_file <- file.path(outdir, paste0("DESeq2_", label, "_all_genes.csv"))
  readr::write_csv(res_df, full_file)
  
  # Significant DEGs (padj < alpha & |log2FC| > lfc)
  res_sig <- res_df %>%
    dplyr::filter(
      !is.na(padj),
      padj < alpha,
      abs(log2FoldChange) > lfc
    ) %>%
    dplyr::arrange(padj)
  
  sig_file <- file.path(outdir,
                        paste0("DESeq2_", label, "_DEGs_padj0.05_log2FC1.csv"))
  readr::write_csv(res_sig, sig_file)
  
  # Gene list for enrichment (one gene per line)
  gene_list_file <- file.path(outdir,
                              paste0("DEGs_", label, "_sig_gene_list.txt"))
  res_sig %>%
    dplyr::select(GeneID) %>%
    readr::write_tsv(gene_list_file)
  
  # Update global summary table (number of significant genes)
  tested <- sum(!is.na(res$padj))
  n_sig  <- nrow(res_sig)
  n_up   <- sum(res_sig$log2FoldChange > 0)
  n_down <- sum(res_sig$log2FoldChange < 0)
  
  group1_pop <- unique(as.character(coldata$population[coldata$condition == group1]))
  group2_pop <- unique(as.character(coldata$population[coldata$condition == group2]))
  group1_trt <- unique(as.character(coldata$treatment[coldata$condition == group1]))
  group2_trt <- unique(as.character(coldata$treatment[coldata$condition == group2]))
  
  summary_degs <<- dplyr::bind_rows(
    summary_degs,
    tibble::tibble(
      contrast                       = label,
      group1_condition               = group1,
      group2_condition               = group2,
      group1_population              = group1_pop,
      group2_population              = group2_pop,
      group1_treatment               = group1_trt,
      group2_treatment               = group2_trt,
      n_genes_tested                 = tested,
      n_DEGs_padj0.05_log2FC1        = n_sig,
      n_DEGs_padj0.05_log2FC1_up     = n_up,
      n_DEGs_padj0.05_log2FC1_down   = n_down,
      alpha                          = alpha,
      log2FC_threshold               = lfc
    )
  )
  
  invisible(res)
}

# 7.1 Insecticide vs control within each population
# New Orleans: NOB vs CNO (Broflanilide), NOI vs CNO (Imidacloprid)
res_NOB_vs_CNO <- run_contrast(dds, group1 = "NOB", group2 = "CNO",
                               label = "NOB_vs_CNO")
res_NOI_vs_CNO <- run_contrast(dds, group1 = "NOI", group2 = "CNO",
                               label = "NOI_vs_CNO")

# San Nicolás: SNB vs CSN (Broflanilide), SNI vs CSN (Imidacloprid)
res_SNB_vs_CSN <- run_contrast(dds, group1 = "SNB", group2 = "CSN",
                               label = "SNB_vs_CSN")
res_SNI_vs_CSN <- run_contrast(dds, group1 = "SNI", group2 = "CSN",
                               label = "SNI_vs_CSN")

# 7.2 Population differences under same treatment
res_CSN_vs_CNO <- run_contrast(dds, group1 = "CSN", group2 = "CNO",
                               label = "CSN_vs_CNO") # baseline controls
res_SNB_vs_NOB <- run_contrast(dds, group1 = "SNB", group2 = "NOB",
                               label = "SNB_vs_NOB") # broflanilide
res_SNI_vs_NOI <- run_contrast(dds, group1 = "SNI", group2 = "NOI",
                               label = "SNI_vs_NOI") # imidacloprid

# 7.3 Save summary table of DEGs per contrast
summary_file <- file.path("degs", "DEGs_summary_counts.csv")
readr::write_csv(summary_degs, summary_file)

cat("Summary of significant DEGs saved to:", summary_file, "\n")

############################################################
# 8. Multi-panel volcano figure (Figure 2)
############################################################

# Helper function: build a volcano panel from DESeq2 results
make_volcano_panel <- function(res_obj, title, subtitle) {
  EnhancedVolcano::EnhancedVolcano(
    res_obj,
    lab = rownames(res_obj),
    x   = "log2FoldChange",
    y   = "padj",
    xlab = bquote(Log[2]~"fold change"),
    ylab = bquote(-Log[10]~"adjusted p-value"),
    title     = title,
    subtitle  = subtitle,
    pCutoff   = 0.05,
    FCcutoff  = 1,
    pointSize = 1.2,
    labSize   = 2.5,
    col       = c("grey80", "#4daf4a", "#377eb8", "#e41a1c"),
    colAlpha  = 0.9,
    legendLabels = c(
      "NS",
      "|log2FC| > 1",
      "padj < 0.05",
      "padj < 0.05 & |log2FC| > 1"
    ),
    legendPosition   = "right",
    gridlines.major  = FALSE,
    gridlines.minor  = FALSE,
    drawConnectors   = TRUE,
    maxoverlapsConnectors = 10,
    widthConnectors  = 0.3
  ) +
    ggplot2::theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 11),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      axis.title    = element_text(size = 9),
      axis.text     = element_text(size = 8),
      legend.title  = element_text(size = 8),
      legend.text   = element_text(size = 7)
    )
}

# Individual volcano panels for each comparison
p_NOB_CNO <- make_volcano_panel(
  res_obj  = res_NOB_vs_CNO,
  title    = "New Orleans – Broflanilide vs control",
  subtitle = "NOB vs CNO"
)

p_NOI_CNO <- make_volcano_panel(
  res_obj  = res_NOI_vs_CNO,
  title    = "New Orleans – Imidacloprid vs control",
  subtitle = "NOI vs CNO"
)

p_SNB_CSN <- make_volcano_panel(
  res_obj  = res_SNB_vs_CSN,
  title    = "San Nicolás – Broflanilide vs control",
  subtitle = "SNB vs CSN"
)

p_SNI_CSN <- make_volcano_panel(
  res_obj  = res_SNI_vs_CSN,
  title    = "San Nicolás – Imidacloprid vs control",
  subtitle = "SNI vs CSN"
)

p_CSN_CNO <- make_volcano_panel(
  res_obj  = res_CSN_vs_CNO,
  title    = "Baseline population difference",
  subtitle = "CSN vs CNO"
)

p_SNB_NOB <- make_volcano_panel(
  res_obj  = res_SNB_vs_NOB,
  title    = "Population difference under broflanilide",
  subtitle = "SNB vs NOB"
)

p_SNI_NOI <- make_volcano_panel(
  res_obj  = res_SNI_vs_NOI,
  title    = "Population difference under imidacloprid",
  subtitle = "SNI vs NOI"
)

# Assemble Figure 2:
# Row 1: NO insecticide vs control
# Row 2: SN insecticide vs control
# Row 3: population differences under control, broflanilide, imidacloprid
fig2_volcano <- (p_NOB_CNO + p_NOI_CNO) /
  (p_SNB_CSN + p_SNI_CSN) /
  (p_CSN_CNO + p_SNB_NOB + p_SNI_NOI) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title      = "Differential gene expression across populations and insecticides"
  ) &
  theme(
    plot.tag   = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )

ggplot2::ggsave(
  filename = file.path("figures_final", "Fig2_volcano_all_comparisons.png"),
  plot     = fig2_volcano,
  width    = 14,
  height   = 12,
  dpi      = 300
)

############################################################
# 9. Gene Ontology (GO) enrichment with NCBI gene2go
############################################################

# 9.1 Locate and read gene2go file (local or download)
go_files_local <- list.files(pattern = "^gene2go")

if (length(go_files_local) == 0) {
  message("No local 'gene2go' file found. Downloading gene2go.gz from NCBI...")
  download.file(
    url      = "https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz",
    destfile = "gene2go.gz"
  )
  gene2go <- data.table::fread("gene2go.gz")
} else {
  message("Using local gene2go file: ", go_files_local[1])
  gene2go <- data.table::fread(go_files_local[1])
}

# In NCBI gene2go, the first column is '#tax_id'; remove the '#' prefix
colnames(gene2go) <- sub("^#", "", colnames(gene2go))

# 9.2 Filter entries for Aedes aegypti (tax_id 7159) and split by ontology
gene2go_aaeg <- gene2go[tax_id == 7159]

gene2go_bp <- gene2go_aaeg[Category == "Process",
                           .(GO_ID, GeneID, GO_term)]
gene2go_mf <- gene2go_aaeg[Category == "Function",
                           .(GO_ID, GeneID, GO_term)]
gene2go_cc <- gene2go_aaeg[Category == "Component",
                           .(GO_ID, GeneID, GO_term)]

TERM2GENE_BP <- gene2go_bp[, .(GO_ID, GeneID)]
TERM2GENE_MF <- gene2go_mf[, .(GO_ID, GeneID)]
TERM2GENE_CC <- gene2go_cc[, .(GO_ID, GeneID)]

TERM2NAME_BP <- unique(gene2go_bp[, .(GO_ID, GO_term)])
TERM2NAME_MF <- unique(gene2go_mf[, .(GO_ID, GO_term)])
TERM2NAME_CC <- unique(gene2go_cc[, .(GO_ID, GO_term)])

# 9.3 Helper to clean GeneIDs (remove 'LOC' prefix if present)
clean_geneid <- function(x) {
  x <- as.character(x)
  gsub("^LOC", "", x)
}

# Gene universe: all genes tested in DESeq2 that have GO annotation
universe_ids <- clean_geneid(rownames(dds))
universe_ids <- universe_ids[universe_ids %in% gene2go_aaeg$GeneID]

cat("Universe size (genes with GO annotation):", length(universe_ids), "\n")

# 9.4 Collect DESeq2 results for the four insecticide-vs-control contrasts
res_list_go <- list(
  NOB_vs_CNO = res_NOB_vs_CNO,  # Broflanilide – New Orleans
  NOI_vs_CNO = res_NOI_vs_CNO,  # Imidacloprid – New Orleans
  SNB_vs_CSN = res_SNB_vs_CSN,  # Broflanilide – San Nicolás
  SNI_vs_CSN = res_SNI_vs_CSN   # Imidacloprid – San Nicolás
)

# 9.5 Function to run GO enrichment (BP/MF/CC) for one contrast
run_go_enrichment <- function(res_obj,
                              contrast_name,
                              padj_cutoff = 0.05,
                              lfc_cutoff  = 1) {
  
  res_df <- as.data.frame(res_obj)
  res_df$GeneID <- clean_geneid(rownames(res_df))
  
  # Set of significant genes for GO (padj + log2FC threshold)
  sig_genes <- res_df %>%
    dplyr::filter(!is.na(padj),
                  padj < padj_cutoff,
                  abs(log2FoldChange) > lfc_cutoff) %>%
    dplyr::pull(GeneID) %>%
    unique()
  
  # Keep only genes that have GO annotation
  sig_genes <- sig_genes[sig_genes %in% gene2go_aaeg$GeneID]
  
  message("Contrast ", contrast_name, ": ",
          length(sig_genes), " significant genes with GO annotation")
  
  # Require at least 10 annotated genes to perform enrichment
  if (length(sig_genes) < 10) {
    warning("Too few annotated genes for GO enrichment in ", contrast_name)
    return(NULL)
  }
  
  # Biological Process
  ego_bp <- clusterProfiler::enricher(
    gene          = sig_genes,
    universe      = universe_ids,
    TERM2GENE     = TERM2GENE_BP,
    TERM2NAME     = TERM2NAME_BP,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05
  )
  
  # Molecular Function
  ego_mf <- clusterProfiler::enricher(
    gene          = sig_genes,
    universe      = universe_ids,
    TERM2GENE     = TERM2GENE_MF,
    TERM2NAME     = TERM2NAME_MF,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05
  )
  
  # Cellular Component
  ego_cc <- clusterProfiler::enricher(
    gene          = sig_genes,
    universe      = universe_ids,
    TERM2GENE     = TERM2GENE_CC,
    TERM2NAME     = TERM2NAME_CC,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05
  )
  
  # Export enrichment tables as TSV files (one per ontology)
  if (!is.null(ego_bp) && nrow(ego_bp@result) > 0) {
    readr::write_tsv(
      ego_bp@result,
      file.path("results/GO", paste0("GO_BP_", contrast_name, ".tsv"))
    )
  }
  
  if (!is.null(ego_mf) && nrow(ego_mf@result) > 0) {
    readr::write_tsv(
      ego_mf@result,
      file.path("results/GO", paste0("GO_MF_", contrast_name, ".tsv"))
    )
  }
  
  if (!is.null(ego_cc) && nrow(ego_cc@result) > 0) {
    readr::write_tsv(
      ego_cc@result,
      file.path("results/GO", paste0("GO_CC_", contrast_name, ".tsv"))
    )
  }
  
  # Optional: barplot of top 15 BP terms for this contrast
  if (!is.null(ego_bp) && nrow(ego_bp@result) > 0) {
    top_bp <- ego_bp@result %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = 15) %>%
      dplyr::mutate(minus_log10_p = -log10(p.adjust))
    
    bp_plot <- ggplot2::ggplot(
      top_bp,
      ggplot2::aes(x = reorder(Description, minus_log10_p),
                   y = minus_log10_p,
                   fill = minus_log10_p)
    ) +
      ggplot2::geom_col(color = "black") +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_gradient(
        name = "-log10(adj p-value)",
        low  = "#c6dbef",
        high = "#08306b"
      ) +
      ggplot2::xlab("GO Biological Process") +
      ggplot2::ylab("-log10(adj p-value)") +
      ggplot2::ggtitle(paste("Top GO BP terms –", contrast_name)) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title   = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.text.y  = element_text(size = 8),
        axis.text.x  = element_text(size = 8),
        axis.title   = element_text(size = 10),
        legend.position = "right",
        plot.margin  = margin(t = 10, r = 20, b = 10, l = 10)
      )
    
    ggplot2::ggsave(
      filename = file.path("results/GO",
                           paste0("GO_BP_barplot_", contrast_name, ".png")),
      plot   = bp_plot,
      width  = 12,
      height = 9,
      dpi    = 300
    )
  }
  
  invisible(list(BP = ego_bp, MF = ego_mf, CC = ego_cc))
}

# Run GO enrichment for the four insecticide vs control contrasts
go_results <- lapply(names(res_list_go), function(cname) {
  run_go_enrichment(res_list_go[[cname]], contrast_name = cname)
})
names(go_results) <- names(res_list_go)

message("GO enrichment finished. Tables saved in 'results/GO/'")

############################################################
# 10. GO dotplots across comparisons (Figures 4–6)
############################################################

# Helper: build data frame for one ontology (BP, CC, MF)
# using GO_*_*.tsv files in results/GO and global top_n terms
build_go_df_from_files <- function(ontology = c("BP", "CC", "MF"),
                                   top_n = 15) {
  ontology <- match.arg(ontology)
  
  go_dir  <- "results/GO"
  pattern <- paste0("^GO_", ontology, "_.*\\.tsv$")
  
  files_go <- list.files(go_dir, pattern = pattern, full.names = TRUE)
  if (length(files_go) == 0) {
    warning("No GO ", ontology, " files found in ", go_dir)
    return(NULL)
  }
  
  message("Reading ", length(files_go), " GO_", ontology,
          " files from ", go_dir, " ...")
  
  go_all <- purrr::map_dfr(files_go, function(f) {
    df <- readr::read_tsv(f, show_col_types = FALSE)
    fname <- basename(f)
    contrast <- sub(paste0("^GO_", ontology, "_"), "", fname)
    contrast <- sub("\\.tsv$", "", contrast)
    df$Contrast <- contrast
    df
  })
  
  # Keep only significant terms (p.adjust < 0.05)
  go_sig <- go_all %>%
    dplyr::filter(!is.na(p.adjust), p.adjust < 0.05) %>%
    dplyr::mutate(
      minus_log10_p = -log10(p.adjust),
      Description   = stringr::str_trunc(Description, 60)
    )
  
  if (nrow(go_sig) == 0) {
    warning("No significant GO ", ontology, " terms with p.adjust < 0.05")
    return(NULL)
  }
  
  # Select global top_n terms by minimum p.adjust across comparisons
  top_terms <- go_sig %>%
    dplyr::group_by(Description) %>%
    dplyr::summarise(min_p = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(min_p) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(Description)
  
  go_top <- go_sig %>%
    dplyr::filter(Description %in% top_terms)
  
  # Order contrasts on x-axis
  go_top <- go_top %>%
    dplyr::mutate(
      Contrast = factor(
        Contrast,
        levels = c("NOB_vs_CNO",
                   "NOI_vs_CNO",
                   "SNB_vs_CSN",
                   "SNI_vs_CSN")
      )
    )
  
  # Order terms by maximum significance across contrasts
  go_top <- go_top %>%
    dplyr::mutate(
      Description = forcats::fct_reorder(Description, minus_log10_p,
                                         .fun = max)
    )
  
  go_top
}

# Helper: build and save a dotplot for one ontology
make_go_dotplot <- function(go_df, ont_title, file_name, gradient_colors) {
  if (is.null(go_df) || nrow(go_df) == 0) {
    warning("No GO terms to plot for ", ont_title)
    return(NULL)
  }
  
  p <- ggplot2::ggplot(go_df, aes(x = Contrast, y = Description)) +
    ggplot2::geom_point(
      aes(size = Count, colour = minus_log10_p),
      alpha = 1
    ) +
    ggplot2::geom_text(
      aes(label = Count),
      vjust = -1.5,
      size  = 2.7
    ) +
    ggplot2::scale_colour_gradientn(
      colours = gradient_colors,
      name    = expression(-log[10]~"(adj p-value)")
    ) +
    ggplot2::scale_size_continuous(
      name  = "Gene count",
      range = c(2.5, 7)
    ) +
    ggplot2::labs(
      x = "Comparison",
      y = NULL,
      title = ont_title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.y      = element_text(size = 8),
      axis.text.x      = element_text(size = 9, angle = 45, hjust = 1),
      axis.title.x     = element_text(size = 10),
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 8),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank()
    )
  
  ggplot2::ggsave(
    filename = file.path("figures_final", file_name),
    plot     = p,
    width    = 10,
    height   = 12,
    dpi      = 300
  )
  
  p
}

# Figure 4 – BP (blue palette), 11 global terms
go_bp_df <- build_go_df_from_files("BP", top_n = 11)
p_bp <- make_go_dotplot(
  go_bp_df,
  ont_title       = "GO Biological Process enrichment across comparisons",
  file_name       = "Fig4_GO_BP_dotplot.png",
  gradient_colors = c("#deebf7", "#9ecae1", "#3182bd", "#08519c")
)

# Figure 5 – CC (green palette), 8 global terms
go_cc_df <- build_go_df_from_files("CC", top_n = 8)
p_cc <- make_go_dotplot(
  go_cc_df,
  ont_title       = "GO Cellular Component enrichment across comparisons",
  file_name       = "Fig5_GO_CC_dotplot.png",
  gradient_colors = c("#e5f5e0", "#a1d99b", "#31a354", "#006d2c")
)

# Figure 6 – MF (red palette), 16 global terms
go_mf_df <- build_go_df_from_files("MF", top_n = 16)
p_mf <- make_go_dotplot(
  go_mf_df,
  ont_title       = "GO Molecular Function enrichment across comparisons",
  file_name       = "Fig6_GO_MF_dotplot.png",
  gradient_colors = c("#fee0d2", "#fc9272", "#fb6a4a", "#a50f15")
)

############################################################
# 11. BP barplots by population and insecticide (Figures 7–8)
############################################################

# Helper: read top_n BP terms for a specific contrast
read_go_bp_top <- function(contrast_name, top_n = 15) {
  fname <- file.path("results", "GO", paste0("GO_BP_", contrast_name, ".tsv"))
  df <- readr::read_tsv(fname, show_col_types = FALSE) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      Description   = stringr::str_trunc(Description, 80),
      minus_log10_p = -log10(p.adjust),
      Description   = forcats::fct_reorder(Description, minus_log10_p)
    )
  df
}

# Helper: build a single BP barplot
make_go_bp_barplot <- function(go_df, panel_title, palette_colors) {
  ggplot2::ggplot(go_df, aes(x = minus_log10_p, y = Description)) +
    ggplot2::geom_col(aes(fill = minus_log10_p), color = "black") +
    ggplot2::scale_fill_gradientn(
      colours = palette_colors,
      name    = expression(-log[10]~"(adj p-value)")
    ) +
    ggplot2::labs(
      x     = expression(-log[10]~"(adj p-value)"),
      y     = NULL,
      title = panel_title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text.y  = element_text(size = 8),
      axis.text.x  = element_text(size = 9),
      axis.title.x = element_text(size = 10),
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8)
    )
}

# Figure 7 – New Orleans (NOI vs CNO, NOB vs CNO)
palette_NO <- c("#f3e5f5", "#ce93d8", "#ab47bc", "#6a1b9a")

bp_NOI <- read_go_bp_top("NOI_vs_CNO", top_n = 15)
bp_NOB <- read_go_bp_top("NOB_vs_CNO", top_n = 15)

p_NOI <- make_go_bp_barplot(
  bp_NOI,
  panel_title    = "New Orleans – Imidacloprid vs control (NOI vs CNO)",
  palette_colors = palette_NO
)

p_NOB <- make_go_bp_barplot(
  bp_NOB,
  panel_title    = "New Orleans – Broflanilide vs control (NOB vs CNO)",
  palette_colors = palette_NO
)

fig7_NO <- p_NOI / p_NOB +
  patchwork::plot_annotation(
    tag_levels = "A",
    title      = "Top GO Biological Process terms – New Orleans"
  ) &
  theme(
    plot.tag   = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
  )

ggplot2::ggsave(
  filename = file.path("figures_final", "Fig7_GO_BP_barplots_NewOrleans.png"),
  plot     = fig7_NO,
  width    = 10,
  height   = 10,
  dpi      = 300
)

# Figure 8 – San Nicolás (SNI vs CSN, SNB vs CSN)
palette_SN <- c("#fff3e0", "#ffcc80", "#fb8c00", "#e65100")

bp_SNI <- read_go_bp_top("SNI_vs_CSN", top_n = 15)
bp_SNB <- read_go_bp_top("SNB_vs_CSN", top_n = 15)

p_SNI <- make_go_bp_barplot(
  bp_SNI,
  panel_title    = "San Nicolás – Imidacloprid vs control (SNI vs CSN)",
  palette_colors = palette_SN
)

p_SNB <- make_go_bp_barplot(
  bp_SNB,
  panel_title    = "San Nicolás – Broflanilide vs control (SNB vs CSN)",
  palette_colors = palette_SN
)

fig8_SN <- p_SNI / p_SNB +
  patchwork::plot_annotation(
    tag_levels = "A",
    title      = "Top GO Biological Process terms – San Nicolás"
  ) &
  theme(
    plot.tag   = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
  )

ggplot2::ggsave(
  filename = file.path("figures_final", "Fig8_GO_BP_barplots_SanNicolas.png"),
  plot     = fig8_SN,
  width    = 10,
  height   = 10,
  dpi      = 300
)

############################################################
# 12. Supplementary Data S1–S3 (Excel workbooks)
############################################################

# S1: Sample metadata + raw counts + filtered counts
supp_S1_metadata <- sample_info %>%
  dplyr::select(sample, population, treatment, condition, replicate)

supp_S1_raw_counts <- counts_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("GeneID")

supp_S1_filtered_counts <- DESeq2::counts(dds, normalized = FALSE) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("GeneID")

supp_S1_list <- list(
  Sample_metadata = supp_S1_metadata,
  Raw_counts      = supp_S1_raw_counts,
  Filtered_counts = supp_S1_filtered_counts
)

writexl::write_xlsx(
  supp_S1_list,
  path = file.path(
    "supplementary",
    "Supplementary_Data_S1_RNAseq_counts_and_metadata.xlsx"
  )
)

# S2: Full DESeq2 results for all 7 contrasts + DEG summary

res_list_s2 <- list(
  NOB_vs_CNO = res_NOB_vs_CNO,
  NOI_vs_CNO = res_NOI_vs_CNO,
  SNB_vs_CSN = res_SNB_vs_CSN,
  SNI_vs_CSN = res_SNI_vs_CSN,
  SNB_vs_NOB = res_SNB_vs_NOB,
  SNI_vs_NOI = res_SNI_vs_NOI,
  CSN_vs_CNO = res_CSN_vs_CNO
)

res_list_s2_df <- lapply(res_list_s2, function(res_obj) {
  as.data.frame(res_obj) %>%
    tibble::rownames_to_column("GeneID")
})

# Recreate a compact DEG summary (using same thresholds as summary_degs)
deg_summary_s2 <- purrr::map_dfr(names(res_list_s2_df), function(name) {
  df <- res_list_s2_df[[name]]
  tested <- sum(!is.na(df$padj))
  sig_fc <- df %>%
    dplyr::filter(
      !is.na(padj),
      padj < 0.05,
      abs(log2FoldChange) > 1
    )
  up   <- sum(sig_fc$log2FoldChange > 0)
  down <- sum(sig_fc$log2FoldChange < 0)
  
  tibble::tibble(
    contrast                       = name,
    n_genes_tested                 = tested,
    n_DEGs_padj0.05_log2FC1        = nrow(sig_fc),
    n_DEGs_padj0.05_log2FC1_up     = up,
    n_DEGs_padj0.05_log2FC1_down   = down
  )
})

supp_S2_list <- lapply(names(res_list_s2_df), function(name) {
  res_list_s2_df[[name]]
})
names(supp_S2_list) <- paste0(names(res_list_s2_df), "_all_genes")
supp_S2_list[["DEG_summary"]] <- deg_summary_s2

writexl::write_xlsx(
  supp_S2_list,
  path = file.path(
    "supplementary",
    "Supplementary_Data_S2_DESeq2_results_all_comparisons.xlsx"
  )
)

# S3: Lists of significant DEGs + GO enrichment tables

deg_lists_sig <- lapply(res_list_s2_df, function(df) {
  df %>%
    dplyr::filter(
      !is.na(padj),
      padj < 0.05,
      abs(log2FoldChange) > 1
    ) %>%
    dplyr::mutate(direction = if_else(log2FoldChange > 0, "up", "down")) %>%
    dplyr::arrange(padj)
})
names(deg_lists_sig) <- paste0("DEGs_", names(res_list_s2_df), "_sig")

go_dir   <- "results/GO"
go_files <- list.files(go_dir, pattern = "^GO_.*\\.tsv$", full.names = TRUE)

if (length(go_files) > 0) {
  go_tables <- lapply(go_files, function(f) {
    readr::read_tsv(f, show_col_types = FALSE)
  })
  names(go_tables) <- basename(go_files) %>%
    sub("\\.tsv$", "", .)  # e.g., "GO_BP_NOI_vs_CNO"
  
  supp_S3_list <- c(deg_lists_sig, go_tables)
} else {
  warning("No GO_*.tsv files found in 'results/GO'. S3 will contain only DEG lists.")
  supp_S3_list <- deg_lists_sig
}

writexl::write_xlsx(
  supp_S3_list,
  path = file.path(
    "supplementary",
    "Supplementary_Data_S3_DEG_lists_and_GO_enrichment.xlsx"
  )
)

cat("\nPipeline finished.\n",
    "Main figures saved in 'figures_final/'.\n",
    "GO tables and barplots in 'results/GO/'.\n",
    "Supplementary files saved in 'supplementary/'.\n")
