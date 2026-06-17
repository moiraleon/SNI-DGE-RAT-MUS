setwd("~/Desktop/Professional/Development/Masters Data Analysis Scripts/SNI-DGE-RAT-MUS")

# =========================================================
# User-defined settings
# =========================================================

ANALYSIS_NAME <- "SNI_SC_MUS_2_MALE_GSE306455_VS_GSE202166"

DESEQ2_RESULTS_FILE <- paste0(
  "final_matrix/deseq2_results_",
  ANALYSIS_NAME,
  ".csv"
)

SIGNIFICANT_GENES_FILE <- paste0(
  "final_matrix/significant_genes_",
  ANALYSIS_NAME,
  "_padj_0.05.csv"
)

PADJ_THRESHOLD <- 0.05

cat("\n==============================\n")
cat("Creating Significant Genes Table\n")
cat("==============================\n\n")

library(readr)
library(dplyr)

# =========================================================
# 1. Load DESeq2 results
# =========================================================

cat("Loading DESeq2 results...\n")

res_df <- read_csv(DESEQ2_RESULTS_FILE)

cat("Rows loaded:", nrow(res_df), "\n\n")

# =========================================================
# 2. Filter significant genes
# =========================================================

cat("Filtering significant genes...\n")

sig_genes <- res_df %>%
  filter(
    !is.na(padj),
    padj < PADJ_THRESHOLD
  ) %>%
  arrange(padj)

cat("Significant genes found:", nrow(sig_genes), "\n\n")

# =========================================================
# 3. Save significant genes
# =========================================================

write_csv(
  sig_genes,
  SIGNIFICANT_GENES_FILE
)

cat("Saved significant genes file:\n")
cat(SIGNIFICANT_GENES_FILE, "\n\n")

cat("==============================\n")
cat("PIPELINE COMPLETE\n")
cat("==============================\n")