setwd("~/Desktop/Professional/Development/Masters Data Analysis Scripts/SNI-DGE-RAT-MUS")

# =========================================================
# User-defined settings
# =========================================================

METADATA_FILE <- "metadata/SNI_SC_MUS_DPI_7_MALE_GSE202166_METADATA.CSV"
COUNTS_FILE <- "final_matrix/gene_level_counts_matrix.csv"

FILE_NAME <- "SNI_SC_MUS_DPI_7_MALE_GSE202166"
TARGET_GSE <- "GSE202166"

DESEQ2_RESULTS_FILE <- paste0("final_matrix/deseq2_results_", FILE_NAME, "_SNI_vs_CTRL.csv")
DESEQ2_SUMMARY_FILE <- paste0("final_matrix/deseq2_summary_", FILE_NAME, "_SNI_vs_CTRL.csv")

cat("\n==============================\n")
cat("Starting One-GSE DESeq2 Differential Expression Pipeline\n")
cat("==============================\n\n")

library(readr)
library(dplyr)
library(tibble)
library(DESeq2)

# =========================================================
# 1. Load counts
# =========================================================

counts <- read_csv(COUNTS_FILE)

gene_info <- counts %>%
  select(ensembl_gene_id, gene_symbol)

count_matrix <- counts %>%
  select(-gene_symbol) %>%
  column_to_rownames("ensembl_gene_id")

# =========================================================
# 2. Load metadata and filter to one GSE
# =========================================================

metadata <- read_csv(METADATA_FILE)

if (!"sample_id" %in% colnames(metadata)) {
  stop("ERROR: metadata file must contain a sample_id column.")
}

if (!"condition" %in% colnames(metadata)) {
  stop("ERROR: metadata file must contain a condition column.")
}

if (!"gse" %in% colnames(metadata)) {
  stop("ERROR: metadata file must contain a gse column.")
}

metadata <- metadata %>%
  filter(gse == TARGET_GSE)

cat("Filtering analysis to:", TARGET_GSE, "\n")
cat("Samples retained:", nrow(metadata), "\n")
cat("Condition table:\n")
print(table(metadata$condition))
cat("\n")

metadata <- metadata %>%
  column_to_rownames("sample_id")

# =========================================================
# 3. Match counts to metadata
# =========================================================

missing_in_counts <- setdiff(rownames(metadata), colnames(count_matrix))

if (length(missing_in_counts) > 0) {
  stop(
    "ERROR: These samples are in metadata but missing from count matrix:\n",
    paste(missing_in_counts, collapse = "\n")
  )
}

count_matrix <- count_matrix[, rownames(metadata)]

stopifnot(all(colnames(count_matrix) == rownames(metadata)))

# =========================================================
# 4. Clean counts
# =========================================================

count_matrix <- as.matrix(count_matrix)
mode(count_matrix) <- "numeric"

if (any(is.na(count_matrix))) {
  stop("ERROR: Count matrix contains NA values.")
}

count_matrix <- round(count_matrix)

# =========================================================
# 5. Create DESeq2 dataset
# =========================================================

metadata$condition <- relevel(factor(metadata$condition), ref = "CTRL")

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ condition
)

cat("DESeq2 dataset created successfully.\n")
cat("Genes before filtering:", nrow(dds), "\n\n")

# =========================================================
# 6. Filter low-count genes
# =========================================================

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

cat("Genes retained after filtering:", nrow(dds), "\n")
cat("Genes removed:", sum(!keep), "\n\n")

# =========================================================
# 7. Run DESeq2
# =========================================================

dds <- DESeq(dds)

# =========================================================
# 8. Extract results
# =========================================================

res <- results(dds, contrast = c("condition", "SNI", "CTRL"))

res_df <- as.data.frame(res) %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(gene_info, by = "ensembl_gene_id") %>%
  relocate(gene_symbol, .after = ensembl_gene_id)

# =========================================================
# 9. Save results
# =========================================================

write_csv(res_df, DESEQ2_RESULTS_FILE)

summary_table <- tibble(
  metric = c(
    "target_gse",
    "total_genes_tested",
    "significant_genes_padj_0.05",
    "upregulated_log2FC_gt_1",
    "downregulated_log2FC_lt_minus_1"
  ),
  value = c(
    TARGET_GSE,
    nrow(res_df),
    sum(res_df$padj < 0.05, na.rm = TRUE),
    sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1, na.rm = TRUE),
    sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm = TRUE)
  )
)

write_csv(summary_table, DESEQ2_SUMMARY_FILE)

cat("Saved DESeq2 results to:", DESEQ2_RESULTS_FILE, "\n")
cat("Saved DESeq2 summary to:", DESEQ2_SUMMARY_FILE, "\n\n")

print(summary_table)

cat("\n==============================\n")
cat("ONE-GSE DESeq2 PIPELINE COMPLETE\n")
cat("==============================\n")