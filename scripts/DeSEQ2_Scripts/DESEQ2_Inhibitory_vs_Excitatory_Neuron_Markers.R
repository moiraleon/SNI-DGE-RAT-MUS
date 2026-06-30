setwd("~/Desktop/Professional/Development/Masters Data Analysis Scripts/SNI-DGE-RAT-MUS")

library(readr)
library(dplyr)
library(grid)
library(gridExtra)

# =========================================================
# User-defined settings
# =========================================================

#FILE_NAME <- "SNI_SC_MUS_2_MALE_GSE306455_VS_GSE202166"
# FILE_NAME <- "SNI_SC_MUS_DPI_7_FEMALE_GSE241361_SNI_vs_CTRL"
# FILE_NAME <- "SNI_SC_MUS_DPI_7_MALE_GSE202166_SNI_vs_CTRL"

DESEQ2_RESULTS_FILE <- paste0(
  "final_matrix/deseq2_results_",
  FILE_NAME,
  ".csv"
)

OUTPUT_DIR <- paste0(
  "VIS_interneuron_marker_analysis/",
  FILE_NAME
)

ALL_MARKERS_FILE <- file.path(
  OUTPUT_DIR,
  "excitatory_inhibitory_interneuron_marker_results.csv"
)

STRICT_MARKERS_FILE <- file.path(
  OUTPUT_DIR,
  "significant_excitatory_inhibitory_markers.csv"
)

CORE_SUMMARY_FILE <- file.path(
  OUTPUT_DIR,
  "core_interneuron_marker_summary.csv"
)

TABLE_PNG_FILE <- file.path(
  OUTPUT_DIR,
  "excitatory_inhibitory_interneuron_marker_table.png"
)

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# =========================================================
# 1. Load DESeq2 results
# =========================================================

res_df <- read_csv(DESEQ2_RESULTS_FILE)

required_columns <- c(
  "ensembl_gene_id",
  "gene_symbol",
  "baseMean",
  "log2FoldChange",
  "pvalue",
  "padj"
)

missing_columns <- setdiff(required_columns, names(res_df))

if (length(missing_columns) > 0) {
  stop(
    "The DESeq2 results file is missing these required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# =========================================================
# 2. Define dorsal-horn interneuron marker panel
# =========================================================
#
# Core identity markers:
# - Inhibitory: Gad1, Gad2, Slc32a1, Slc6a5, Pax2
# - Excitatory: Slc17a6
#
# Subtype-associated markers are exploratory and should not
# be interpreted individually as proof of excitatory or
# inhibitory circuit shifts.
# =========================================================

marker_panel <- tibble::tibble(
  gene_symbol = c(
    # Inhibitory core identity markers
    "Gad1",
    "Gad2",
    "Slc32a1",
    "Slc6a5",
    "Pax2",

    # Inhibitory subtype-associated markers
    "Pdyn",
    "Gal",
    "Npy",
    "Nos1",

    # Excitatory core identity marker
    "Slc17a6",

    # Excitatory subtype-associated markers
    "Sst",
    "Tac1",
    "Cck",
    "Grp",
    "Rorb",
    "Crh",
    "Prkcg"
  ),

  marker_group = c(
    rep("Inhibitory", 9),
    rep("Excitatory", 8)
  ),

  marker_level = c(
    rep("Core identity", 5),
    rep("Subtype-associated", 4),
    "Core identity",
    rep("Subtype-associated", 7)
  ),

  interpretation = c(
    "GABA synthesis",
    "GABA synthesis",
    "Vesicular inhibitory transmitter transporter",
    "Glycinergic inhibitory marker",
    "Inhibitory interneuron-associated transcription factor",

    "Inhibitory subtype-associated marker",
    "Inhibitory subtype-associated marker",
    "Inhibitory subtype-associated marker",
    "Inhibitory subtype-associated marker",

    "Vesicular glutamate transporter; broad excitatory marker",

    "Excitatory subtype-associated marker",
    "Excitatory subtype-associated marker",
    "Excitatory subtype-associated marker",
    "Excitatory subtype-associated marker",
    "Excitatory subtype-associated marker",
    "Excitatory subtype-associated marker",
    "Excitatory subtype-associated marker"
  )
)

# =========================================================
# 3. Match marker panel to DESeq2 results
# =========================================================

res_df_unique <- res_df %>%
  filter(!is.na(.data$gene_symbol)) %>%
  distinct(.data$gene_symbol, .keep_all = TRUE)

marker_results <- marker_panel %>%
  left_join(
    res_df_unique %>%
      select(
        gene_symbol,
        ensembl_gene_id,
        baseMean,
        log2FoldChange,
        pvalue,
        padj
      ),
    by = "gene_symbol"
  ) %>%
  mutate(
    classification = case_when(
      is.na(.data$log2FoldChange) ~
        "Not found in DESeq2 results",

      !is.na(.data$padj) &
        .data$padj < 0.05 &
        .data$log2FoldChange > 1 ~
        "Significant increased",

      !is.na(.data$padj) &
        .data$padj < 0.05 &
        .data$log2FoldChange < -1 ~
        "Significant decreased",

      !is.na(.data$padj) &
        .data$padj < 0.05 ~
        "Significant, smaller effect",

      TRUE ~
        "Not significant"
    ),

    marker_group = factor(
      .data$marker_group,
      levels = c("Inhibitory", "Excitatory")
    ),

    marker_level = factor(
      .data$marker_level,
      levels = c("Core identity", "Subtype-associated")
    )
  ) %>%
  arrange(
    .data$marker_group,
    .data$marker_level,
    desc(.data$log2FoldChange)
  )

# =========================================================
# 4. Create strict significant-marker subset
# =========================================================

strict_markers <- marker_results %>%
  filter(
    .data$classification %in% c(
      "Significant increased",
      "Significant decreased"
    )
  ) %>%
  arrange(
    .data$marker_group,
    desc(abs(.data$log2FoldChange))
  )

# =========================================================
# 5. Summarize core identity markers
# =========================================================

core_marker_summary <- marker_results %>%
  filter(.data$marker_level == "Core identity") %>%
  group_by(.data$marker_group) %>%
  summarise(
    markers_in_panel = n(),
    markers_present_in_results = sum(!is.na(.data$log2FoldChange)),
    significantly_increased = sum(
      .data$classification == "Significant increased"
    ),
    significantly_decreased = sum(
      .data$classification == "Significant decreased"
    ),
    significant_smaller_effect = sum(
      .data$classification == "Significant, smaller effect"
    ),
    .groups = "drop"
  )

# =========================================================
# 6. Save CSV outputs
# =========================================================

write_csv(marker_results, ALL_MARKERS_FILE)
write_csv(strict_markers, STRICT_MARKERS_FILE)
write_csv(core_marker_summary, CORE_SUMMARY_FILE)

# =========================================================
# 7. Create thesis-ready marker table
# =========================================================

marker_display <- marker_results %>%
  mutate(
    baseMean_display = if_else(
      is.na(.data$baseMean),
      "NA",
      format(round(.data$baseMean, 2), nsmall = 2)
    ),

    log2FC_display = if_else(
      is.na(.data$log2FoldChange),
      "NA",
      format(round(.data$log2FoldChange, 2), nsmall = 2)
    ),

    pvalue_display = if_else(
      is.na(.data$pvalue),
      "NA",
      format(.data$pvalue, scientific = TRUE, digits = 3)
    ),

    padj_display = if_else(
      is.na(.data$padj),
      "NA",
      format(.data$padj, scientific = TRUE, digits = 3)
    )
  ) %>%
  transmute(
    `Marker group` = marker_group,
    `Marker level` = marker_level,
    `Gene symbol` = gene_symbol,
    `Ensembl ID` = ensembl_gene_id,
    `Base mean` = baseMean_display,
    `log2 fold change` = log2FC_display,
    `P-value` = pvalue_display,
    `Adjusted p-value` = padj_display,
    Classification = classification
  )

marker_table_grob <- tableGrob(
  marker_display,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 8,
    core = list(
      fg_params = list(
        hjust = 0.5,
        x = 0.5
      )
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

png(
  filename = TABLE_PNG_FILE,
  width = 18,
  height = 9,
  units = "in",
  res = 300
)

grid.newpage()

grid.text(
  paste(
    "Excitatory and Inhibitory Interneuron Marker Results",
    FILE_NAME,
    sep = "\n"
  ),
  x = 0.5,
  y = 0.97,
  gp = gpar(
    fontsize = 15,
    fontface = "bold"
  )
)

grid.draw(marker_table_grob)

dev.off()

# =========================================================
# 8. Print results
# =========================================================

cat("\n==============================================\n")
cat("EXCITATORY / INHIBITORY INTERNEURON ANALYSIS\n")
cat("==============================================\n")

cat("\nCore marker summary:\n")
print(core_marker_summary)

cat("\nStrict significant interneuron markers:\n")

if (nrow(strict_markers) == 0) {
  cat(
    "No interneuron markers met both criteria: ",
    "adjusted p < 0.05 and |log2FC| > 1.\n",
    sep = ""
  )
} else {
  print(
    strict_markers %>%
      select(
        marker_group,
        marker_level,
        gene_symbol,
        baseMean,
        log2FoldChange,
        pvalue,
        padj,
        classification
      )
  )
}

cat("\nAll marker results:\n")

print(
  marker_results %>%
    select(
      marker_group,
      marker_level,
      gene_symbol,
      log2FoldChange,
      pvalue,
      padj,
      classification
    )
)

cat("\nSaved files:\n")
cat(ALL_MARKERS_FILE, "\n")
cat(STRICT_MARKERS_FILE, "\n")
cat(CORE_SUMMARY_FILE, "\n")
cat(TABLE_PNG_FILE, "\n")