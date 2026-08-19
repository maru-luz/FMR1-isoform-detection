#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT 3: DIFFERENTIAL EXPRESSION ANALYSIS (DESeq2)
# ==============================================================================
# Description:
# This script performs Differential Transcript Expression (DTE) analysis on 
# FMR1 isoforms using DESeq2. It takes the filtered raw count matrix 
# (excluding outlier samples bc02, bc03, bc09), assigns biological groups 
# (B1, B2, GC), and calculates pairwise comparisons.
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
})

# ------------------------------------------------------------------------------
# 1. DATA LOADING
# ------------------------------------------------------------------------------
counts_file <- "filtered_isoforms.tsv"
counts_matrix <- read.table(counts_file, 
                            header = TRUE, 
                            row.names = 1, 
                            sep = "\t", 
                            check.names = FALSE) %>% 
  as.matrix()

# ------------------------------------------------------------------------------
# 2. METADATA CREATION (colData)
# ------------------------------------------------------------------------------
# Define barcode lists according to study groups (excluding outliers)
b1_list <- c("bc01","bc04","bc07","bc10","bc13")
b2_list <- c("bc05","bc08","bc11","bc14","bc16")
gc_list <- c("bc06","bc12","bc15","bc17")

# Create metadata dataframe
col_data <- data.frame(barcode = colnames(counts_matrix)) %>%
  mutate(tejido = case_when(
    barcode %in% b1_list ~ "B1",
    barcode %in% b2_list ~ "B2",
    barcode %in% gc_list ~ "GC",
    TRUE ~ "Other"
  )) %>%
  column_to_rownames("barcode")

# Force factor conversion and set GC (Control Group) as the reference level
col_data$tejido <- factor(col_data$tejido, levels = c("GC", "B1", "B2"))

# ------------------------------------------------------------------------------
# 3. DESeq2 OBJECT CONSTRUCTION
# ------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = round(counts_matrix),
                              colData = col_data,
                              design = ~ tejido)

# ------------------------------------------------------------------------------
# 4. RUN ANALYSIS
# ------------------------------------------------------------------------------
dds <- DESeq(dds)

# ------------------------------------------------------------------------------
# 5. EXTRACT AND SAVE RESULTS
# ------------------------------------------------------------------------------
# Define the pairwise comparisons of interest
comparisons <- list(
  B1_vs_GC = c("tejido", "B1", "GC"),
  B2_vs_GC = c("tejido", "B2", "GC"),
  B1_vs_B2 = c("tejido", "B1", "B2")
)

# Loop to process, sort (by adjusted p-value), and save each table
for(comp_name in names(comparisons)){
  res <- results(dds, contrast = comparisons[[comp_name]])
  
  res_tab <- as.data.frame(res) %>%
    arrange(padj)
  
  output_file <- paste0("para_paper/DTE_filtr_final_", comp_name, "_FMR1.csv")
  write.csv(res_tab, file = output_file, row.names = TRUE)
}

cat("✔ Analysis complete. FMR1 differential expression result tables have been generated.\n")
