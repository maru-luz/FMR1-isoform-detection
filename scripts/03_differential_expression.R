#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT 3: DIFFERENTIAL EXPRESSION ANALYSIS (limma-voom)
# ==============================================================================
# Description:
# This script performs Differential Transcript Expression (DTE) analysis on
# FMR1 isoforms using limma-voom with duplicateCorrelation to account for
# the repeated-measures structure of the design (multiple samples per woman).
#
# Since QC-based exclusion of samples from women A and D resulted in an
# incomplete block design, woman was modeled as a random effect using
# duplicateCorrelation (Ritchie et al., NAR 2015; Smyth et al., limma
# User's Guide, Bioconductor, section 9.7) rather than as a fixed effect.
# Technical replicates from woman E (bc14/bc16 for B2; bc15/bc17 for GC)
# were averaged prior to analysis to prevent their high within-replicate
# correlation from dominating the intra-woman correlation estimate.
# Both voom and duplicateCorrelation were run in two iterative rounds
# as recommended in the limma documentation.
#
# Three pairwise comparisons are performed:
#   - B1 vs GC
#   - B2 vs GC
#   - B1 vs B2
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(edgeR)
})

# ------------------------------------------------------------------------------
# 1. DATA LOADING
# ------------------------------------------------------------------------------
counts_file <- "filtered_isoforms.tsv"
counts_matrix <- read.table(counts_file,
                            header      = TRUE,
                            row.names   = 1,
                            sep         = "\t",
                            check.names = FALSE) %>%
  as.matrix()

# ------------------------------------------------------------------------------
# 2. AVERAGE TECHNICAL REPLICATES FROM DONOR E
# ------------------------------------------------------------------------------
# Woman E contributed technical replicates for B2 (bc14, bc16) and GC (bc15,
# bc17). These are averaged to a single value per tissue type per woman to
# avoid inflating the intra-woman correlation estimate in duplicateCorrelation.

E_B2 <- round(rowMeans(counts_matrix[, c("bc14", "bc16")]))
E_GC <- round(rowMeans(counts_matrix[, c("bc15", "bc17")]))

counts_matrix <- counts_matrix[, !colnames(counts_matrix) %in%
                                   c("bc14", "bc15", "bc16", "bc17")]
counts_matrix <- cbind(counts_matrix, E_B2 = E_B2, E_GC = E_GC)

cat("Samples after averaging woman E replicates:", colnames(counts_matrix), "\n\n")

# ------------------------------------------------------------------------------
# 3. METADATA CREATION
# ------------------------------------------------------------------------------
# Tissue and woman assignments for each sample (post-averaging)
col_data <- data.frame(
  barcode = colnames(counts_matrix),
  tissue  = c("B1", "B1", "B2", "GC",   # bc01(A), bc04(B), bc05(B), bc06(B)
              "B1", "B2",               # bc07(D), bc08(D)
              "B1", "B2", "GC",         # bc10(C), bc11(C), bc12(C)
              "B1", "B2", "GC"),        # bc13(E), E_B2, E_GC
  woman   = c("A", "B", "B", "B",
              "D", "D",
              "C", "C", "C",
              "E", "E", "E"),
  row.names = NULL
) %>% column_to_rownames("barcode")

col_data$tissue <- factor(col_data$tissue, levels = c("GC", "B1", "B2"))
col_data$woman  <- factor(col_data$woman)

# Verify sample alignment
stopifnot(all(colnames(counts_matrix) == rownames(col_data)))

# ------------------------------------------------------------------------------
# 4. NORMALIZATION
# ------------------------------------------------------------------------------
dge <- DGEList(counts = counts_matrix)
dge <- calcNormFactors(dge)

# ------------------------------------------------------------------------------
# 5. DESIGN MATRIX
# ------------------------------------------------------------------------------
# Woman is NOT included in the design formula, its contribution is accounted
# for via duplicateCorrelation (random effect approximation)
design <- model.matrix(~ tissue, data = col_data)

# ------------------------------------------------------------------------------
# 6. VOOM + DUPLICATECORRELATION (TWO ITERATIVE ROUNDS)
# ------------------------------------------------------------------------------
# Round 1: estimate intra-woman correlation with initial precision weights
cat("Round 1: voom + duplicateCorrelation...\n")
v1      <- voom(dge, design, plot = FALSE)
corfit1 <- duplicateCorrelation(v1, design, block = col_data$woman)
cat(sprintf("  Intra-woman correlation (round 1): %.3f\n", corfit1$consensus))

# Round 2: re-estimate precision weights incorporating the correlation,
# then re-estimate the correlation with improved weights
cat("Round 2: voom + duplicateCorrelation...\n")
v2      <- voom(dge, design,
                block       = col_data$woman,
                correlation = corfit1$consensus,
                plot        = FALSE)
corfit2 <- duplicateCorrelation(v2, design, block = col_data$woman)
cat(sprintf("  Intra-woman correlation (round 2): %.3f\n\n", corfit2$consensus))

# ------------------------------------------------------------------------------
# 7. MODEL FITTING
# ------------------------------------------------------------------------------
fit <- lmFit(v2, design,
             block       = col_data$woman,
             correlation = corfit2$consensus)

# Define pairwise contrasts
cont_matrix <- makeContrasts(
  B1_vs_GC = tissueB1,
  B2_vs_GC = tissueB2,
  B1_vs_B2 = tissueB1 - tissueB2,
  levels   = design
)

fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2)

# ------------------------------------------------------------------------------
# 8. EXTRACT AND SAVE RESULTS
# ------------------------------------------------------------------------------
dir.create("para_paper", showWarnings = FALSE)

comparisons <- c("B1_vs_GC", "B2_vs_GC", "B1_vs_B2")

for(comp_name in comparisons){
  res <- topTable(fit2, coef = comp_name, number = Inf, sort.by = "P") %>%
    as.data.frame()

  n_sig <- sum(res$adj.P.Val < 0.05 & abs(res$logFC) > 1, na.rm = TRUE)
  cat(sprintf("%s: %d significant isoforms (adj.P.Val < 0.05, |logFC| > 1)\n",
              comp_name, n_sig))

  output_file <- paste0("para_paper/DTE_limma_", comp_name, "_FMR1.csv")
  write.csv(res, file = output_file, row.names = TRUE)
}

cat("\n✔ Analysis complete. Results saved in para_paper/DTE_limma_*_FMR1.csv\n")
