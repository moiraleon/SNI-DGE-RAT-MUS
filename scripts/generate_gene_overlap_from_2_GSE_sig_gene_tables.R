setwd("~/Desktop/Professional/Development/Masters Data Analysis Scripts/SNI-DGE-RAT-MUS")

# =========================================================
# User-defined settings
# =========================================================

GSE_A <- "SNI_SC_MUS_DPI_7_FEMALE_GSE241361"
GSE_B <- "SNI_SC_MUS_DPI_7_MALE_GSE202166"

SIG_GENES_FILE_A <- paste0(
  "final_matrix/significant_genes_",
  GSE_A,
  "_SNI_vs_CTRL_padj_0.05.csv"
)

SIG_GENES_FILE_B <- paste0(
  "final_matrix/significant_genes_",
  GSE_B,
  "_SNI_vs_CTRL_padj_0.05.csv"
)

SHARED_GENES_OUTPUT_FILE <- paste0(
  "final_matrix/shared_significant_genes_",
  GSE_A,
  "_vs_",
  GSE_B,
  "_by_ensembl_gene_id.csv"
)

OVERLAP_SUMMARY_FILE <- paste0(
  "final_matrix/shared_significant_genes_",
  GSE_A,
  "_vs_",
  GSE_B,
  "_summary.csv"
)

cat("\n==============================\n")
cat("Starting Significant Gene Overlap Analysis\n")
cat("==============================\n\n")

library(readr)
library(dplyr)

# =========================================================
# 1. Load significant gene tables
# =========================================================

cat("STEP 1: Loading significant gene tables...\n")

sig_a <- read_csv(SIG_GENES_FILE_A)
sig_b <- read_csv(SIG_GENES_FILE_B)

cat(GSE_A, "significant genes:", nrow(sig_a), "\n")
cat(GSE_B, "significant genes:", nrow(sig_b), "\n\n")

# =========================================================
# 2. Check required column
# =========================================================

if (!"ensembl_gene_id" %in% colnames(sig_a)) {
  stop("ERROR: First significant gene file must contain ensembl_gene_id column.")
}

if (!"ensembl_gene_id" %in% colnames(sig_b)) {
  stop("ERROR: Second significant gene file must contain ensembl_gene_id column.")
}

# =========================================================
# 3. Compare by Ensembl gene ID
# =========================================================

cat("STEP 2: Finding shared significant genes by Ensembl gene ID...\n")

shared_genes <- sig_a %>%
  inner_join(
    sig_b,
    by = "ensembl_gene_id",
    suffix = c(paste0("_", GSE_A), paste0("_", GSE_B))
  )

cat("Shared significant genes found:", nrow(shared_genes), "\n\n")

# =========================================================
# 4. Create summary table
# =========================================================

summary_table <- tibble(
  metric = c(
    "gse_a",
    "gse_b",
    "significant_genes_gse_a",
    "significant_genes_gse_b",
    "shared_significant_genes",
    "gse_a_only_significant_genes",
    "gse_b_only_significant_genes"
  ),
  value = c(
    GSE_A,
    GSE_B,
    nrow(sig_a),
    nrow(sig_b),
    nrow(shared_genes),
    nrow(sig_a) - nrow(shared_genes),
    nrow(sig_b) - nrow(shared_genes)
  )
)

# =========================================================
# 5. Save shared genes
# =========================================================

write_csv(
  shared_genes,
  SHARED_GENES_OUTPUT_FILE
)

write_csv(
  summary_table,
  OVERLAP_SUMMARY_FILE
)

cat("Saved shared genes to:\n")
cat(SHARED_GENES_OUTPUT_FILE, "\n\n")

cat("==============================\n")
cat("OVERLAP ANALYSIS COMPLETE\n")
cat("==============================\n")