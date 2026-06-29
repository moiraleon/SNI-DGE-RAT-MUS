setwd("~/Desktop/Professional/Development/Masters Data Analysis Scripts/SNI-DGE-RAT-MUS")

library(readr)
library(dplyr)
library(grid)
library(gridExtra)


# =========================================================
# User-defined settings
# =========================================================

FILE_NAME <-  "SNI_SC_MUS_2_MALE_GSE306455_VS_GSE202166"
#"SNI_SC_MUS_2_MALE_GSE306455_VS_GSE202166"

DESEQ2_RESULTS_FILE <- paste0(
  "final_matrix/deseq2_results_",
  FILE_NAME,
  ".csv"
)

ION_CHANNEL_RESULTS_FILE <- paste0(
  "final_matrix/v3_ion_channels_ranked_",
  FILE_NAME,
 ".csv"
)

UPREGULATED_ION_CHANNELS_FILE <- paste0(
  "final_matrix/v3_upregulated_ion_channels_",
  FILE_NAME,
  #"_SNI_vs_CTRL.csv"
  ".csv"
)
OUTPUT_DIR <- paste0(
"VIS_ion_channel_analysis/",
FILE_NAME
)

TOP_N <- 10

if (!dir.exists(OUTPUT_DIR)) {
dir.create(OUTPUT_DIR, recursive = TRUE)
}


# =========================================================
# 1. Load DESeq2 results
# =========================================================

res_df <- read_csv(DESEQ2_RESULTS_FILE)

# =========================================================
# 2. Filter for ion-channel-related genes
# =========================================================

ion_channel_pattern <- paste(
  c(
    # Purinergic ion channels
    "^P2rx[1-7]$",

    # Voltage-gated sodium channels and ENaC subunits
    "^Scn[0-9]+[a-z]+$",
    "^Scnn[0-9]+[a-z]+$",

    # Potassium-channel genes, subunits, and regulators
    "^Kcn[a-z0-9]+$",

    # Calcium-channel genes and associated subunits
    "^Cacn[a-z0-9]+$",

    # True TRP-channel families
    "^Trpc[1-7]$",
    "^Trpm[1-8]$",
    "^Trpa1$",
    "^Trpv[1-6]$",
    "^Mcoln[1-3]$",

    # Chloride channels
    "^Clcn[a-z0-9]+$",

    # Mechanosensitive channels
    "^Piezo[12]$",

    # Additional ion-channel families
    "^Hcn[1-4]$",
    "^Asic[1-5]$",
    "^Orai[1-3]$",
    "^Ano[12]$",

    # Ionotropic glutamate receptors
    "^Gria[1-4]$",
    "^Grik[1-5]$",
    "^Grin[1-3][a-z]?$",

    # Ionotropic GABA receptors
    "^Gabra[1-6]$",
    "^Gabrb[1-3]$",
    "^Gabrd$",
    "^Gabre$",
    "^Gabrg[1-3]$",
    "^Gabrr[1-3]$",

    # Glycine receptors
    "^Glyra[1-4]$",
    "^Glyrb$",

    # Nicotinic acetylcholine receptors
    "^Chrna[0-9]+$",
    "^Chrnb[1-4]$",
    "^Chrnd$",
    "^Chrne$",
    "^Chrng$",

    # Serotonin 5-HT3 receptor subunits
    "^Htr3[ab]$"
  ),
  collapse = "|"
)

ion_channels_df <- res_df %>%
  filter(
    !is.na(.data$gene_symbol),
    grepl(ion_channel_pattern, .data$gene_symbol)
  ) %>%
  arrange(desc(.data$log2FoldChange))

# =========================================================
# 3. Filter for significantly upregulated ion channel genes
# =========================================================

upregulated_ion_channels <- ion_channels_df %>%
  filter(
    !is.na(.data$padj),
    .data$padj < 0.05,
    .data$log2FoldChange > 1
  ) %>%
  arrange(desc(.data$log2FoldChange))

# =========================================================
# 4. Save outputs
# =========================================================

write_csv(ion_channels_df, ION_CHANNEL_RESULTS_FILE)
write_csv(upregulated_ion_channels, UPREGULATED_ION_CHANNELS_FILE)
# =========================================================
# 5. Create table of top-ranked ion-channel-related genes
# =========================================================

top_10_ion_channels <- ion_channels_df %>%
  slice_head(n = TOP_N) %>%
  mutate(
    classification = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange > 1 ~
        "Significant upregulated",
      !is.na(padj) & padj < 0.05 ~
        "Significant, smaller effect",
      TRUE ~
        "Not significant"
    )
  )

cat("\nTop ", TOP_N, " ranked ion-channel-related genes:\n", sep = "")

print(
  top_10_ion_channels %>%
    select(
      ensembl_gene_id,
      gene_symbol,
      baseMean,
      log2FoldChange,
      pvalue,
      padj,
      classification
    )
)

write_csv(
  top_10_ion_channels,
  file.path(
    OUTPUT_DIR,
    "top_10_ranked_ion_channel_genes.csv"
  )
)

top_10_display <- top_10_ion_channels %>%
  mutate(
    Rank = row_number(),
    baseMean = round(baseMean, 2),
    log2FoldChange = round(log2FoldChange, 2),
    pvalue = format(pvalue, scientific = TRUE, digits = 3),
    padj = if_else(
      is.na(padj),
      "NA",
      format(padj, scientific = TRUE, digits = 3)
    )
  ) %>%
  transmute(
    Rank,
    `Gene symbol` = gene_symbol,
    `Ensembl ID` = ensembl_gene_id,
    `Base mean` = baseMean,
    `log2 fold change` = log2FoldChange,
    `P-value` = pvalue,
    `Adjusted p-value` = padj,
    Classification = classification
  )

table_grob <- tableGrob(
  top_10_display,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 9,
    core = list(
      fg_params = list(hjust = 0.5, x = 0.5)
    ),
    colhead = list(
      fg_params = list(
        fontface = "bold",
        hjust = 0.5,
        x = 0.5
      )
    )
  )
)

table_png_file <- file.path(
  OUTPUT_DIR,
  "top_10_ranked_ion_channel_genes.png"
)

png(
  filename = table_png_file,
  width = 15,
  height = 5,
  units = "in",
  res = 300
)

grid.newpage()

grid.text(
  "Top 10 Ion-Channel-Related Genes Ranked by log2 Fold Change",
  x = 0.5,
  y = 0.97,
  gp = gpar(
    fontsize = 14,
    fontface = "bold"
  )
)

grid.draw(table_grob)

dev.off()

cat("\nSaved top 10 ion-channel table to:\n")
cat(file.path(OUTPUT_DIR, "top_10_ranked_ion_channel_genes.csv"), "\n")
cat(table_png_file, "\n")

# =========================================================
# 6. Print results
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