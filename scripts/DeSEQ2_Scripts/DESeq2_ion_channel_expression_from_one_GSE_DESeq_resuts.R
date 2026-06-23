setwd("~/Desktop/Professional/Development/Masters Data Analysis Scripts/SNI-DGE-RAT-MUS")

library(readr)
library(dplyr)

# =========================================================
# User-defined settings
# =========================================================

FILE_NAME <- "SNI_SC_MUS_DPI_7_FEMALE_GSE241361"

DESEQ2_RESULTS_FILE <- paste0(
  "final_matrix/deseq2_results_",
  FILE_NAME,
  "_SNI_vs_CTRL.csv"
)

ION_CHANNEL_RESULTS_FILE <- paste0(
  "final_matrix/ion_channels_ranked_",
  FILE_NAME,
  "_SNI_vs_CTRL.csv"
)

UPREGULATED_ION_CHANNELS_FILE <- paste0(
  "final_matrix/upregulated_ion_channels_",
  FILE_NAME,
  "_SNI_vs_CTRL.csv"
)

# =========================================================
# 1. Load DESeq2 results
# =========================================================

res_df <- read_csv(DESEQ2_RESULTS_FILE)

# =========================================================
# 2. Filter for ion channel genes
# =========================================================

ion_channels_df <- res_df %>%
  filter(
    grepl("^P2rx", gene_symbol) |
      grepl("^Scn", gene_symbol) |
      grepl("^Kcn", gene_symbol) |
      grepl("^Cacn", gene_symbol) |
      grepl("^Trp", gene_symbol) |
      grepl("^Clcn", gene_symbol) |
      grepl("^Piezo", gene_symbol)
  ) %>%
  arrange(desc(log2FoldChange))

# =========================================================
# 3. Filter for significantly upregulated ion channels
# =========================================================

#removing p value based on PMID: 34244727 sorting by log fold changes
# upregulated_ion_channels <- ion_channels_df %>%
#   filter(
#     !is.na(padj),
#     padj < 0.05,
#     log2FoldChange > 0
#   ) %>%
#   arrange(desc(log2FoldChange))

upregulated_ion_channels <- ion_channels_df %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange))

# =========================================================
# 4. Save outputs
# =========================================================

write_csv(ion_channels_df, ION_CHANNEL_RESULTS_FILE)
write_csv(upregulated_ion_channels, UPREGULATED_ION_CHANNELS_FILE)

# =========================================================
# 5. Print results
# =========================================================

cat("\nIon channel genes found:", nrow(ion_channels_df), "\n")
cat("Significantly upregulated ion channels:", nrow(upregulated_ion_channels), "\n\n")

cat("Top upregulated ion channels:\n")
print(
  upregulated_ion_channels %>%
    select(gene_symbol, log2FoldChange, pvalue, padj) %>%
    head(20)
)

cat("\nP2rx3 result:\n")
print(
  res_df %>%
    filter(gene_symbol == "P2rx3") %>%
    select(gene_symbol, log2FoldChange, pvalue, padj)
)

cat("\nSaved ranked ion channels to:", ION_CHANNEL_RESULTS_FILE, "\n")
cat("Saved upregulated ion channels to:", UPREGULATED_ION_CHANNELS_FILE, "\n")