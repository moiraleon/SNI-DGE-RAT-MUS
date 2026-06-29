setwd("~/Desktop/Professional/Development/Masters Data Analysis Scripts/SNI-DGE-RAT-MUS")

cat("\n============================================\n")
cat("Starting 14 DPI Visualization Pipeline\n")
cat("============================================\n\n")

# =========================================================

# Install packages once only, if needed

# Keep these lines commented after installation.

# =========================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages(c(

"readr",

"dplyr",

"tibble",

"ggplot2",

"ggrepel",

"gridExtra"

))

# =========================================================

# Load libraries

# =========================================================

library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)

cat("Libraries loaded successfully.\n\n")

# =========================================================

# User-defined settings

# =========================================================

DESEQ2_RESULTS_FILE <- "final_matrix/PHASE_2/COMB_2/deseq2_results_SNI_SC_MUS_2_MALE_GSE306455_VS_GSE202166.csv"

OUTPUT_DIR <- "VIS_temporal_analysis/COMB_2_14DPI_Male_GSE306455_GSE202166"

PADJ_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 1

TISSUE <- "Spinal cord"
SEX <- "Male"
DPI <- 14
DATASETS <- "GSE306455 and GSE202166"

# =========================================================

# Create output directory

# =========================================================

if (!dir.exists(OUTPUT_DIR)) {
dir.create(OUTPUT_DIR, recursive = TRUE)
}

if (!dir.exists(OUTPUT_DIR)) {
stop("ERROR: Failed to create output directory: ", OUTPUT_DIR)
}

cat("Output directory ready:\n")
cat(OUTPUT_DIR, "\n\n")

# =========================================================

# 1. Load DESeq2 results

# =========================================================

cat("STEP 1: Loading DESeq2 results...\n")

res <- read_csv(DESEQ2_RESULTS_FILE, show_col_types = FALSE)

required_columns <- c(
"ensembl_gene_id",
"gene_symbol",
"log2FoldChange",
"padj"
)

missing_columns <- setdiff(required_columns, colnames(res))

if (length(missing_columns) > 0) {
stop(
"ERROR: Missing required columns:\n",
paste(missing_columns, collapse = ", ")
)
}

cat("Results loaded:", nrow(res), "genes\n\n")

# =========================================================

# 2. Clean gene IDs and classify genes

# =========================================================

cat("STEP 2: Cleaning and classifying genes...\n")

res <- res %>%
mutate(
ensembl_clean = sub("[.].*$", "", ensembl_gene_id),

gene_label = if_else(
  is.na(gene_symbol) | gene_symbol == "",
  ensembl_clean,
  gene_symbol
),

neg_log10_padj = -log10(pmax(padj, 1e-300)),

significance = case_when(
  padj < PADJ_CUTOFF & log2FoldChange > LOG2FC_CUTOFF ~ "Upregulated",
  padj < PADJ_CUTOFF & log2FoldChange < -LOG2FC_CUTOFF ~ "Downregulated",
  padj < PADJ_CUTOFF ~ "Significant, smaller effect",
  TRUE ~ "Not significant"
)


)

res$significance <- factor(
res$significance,
levels = c(
"Upregulated",
"Downregulated",
"Significant, smaller effect",
"Not significant"
)
)

cat("Gene classification complete.\n")
print(table(res$significance, useNA = "ifany"))
cat("\n")

# =========================================================

# 3. Extract significant genes for Table 2 and labels

# =========================================================

cat("STEP 3: Extracting significant genes...\n")

significant_genes <- res %>%
filter(
!is.na(padj),
!is.na(log2FoldChange),
padj < PADJ_CUTOFF
) %>%
select(
gene_symbol = gene_label,
ensembl_gene_id = ensembl_clean,
log2FoldChange,
padj,
classification = significance
) %>%
arrange(desc(abs(log2FoldChange)))

cat("Significant genes:", nrow(significant_genes), "\n")
print(significant_genes)
cat("\n")

write_csv(
significant_genes,
file.path(
OUTPUT_DIR,
"table_2_significant_genes_14dpi.csv"
)
)

# =========================================================

# 4. Generate Figure 6: Volcano plot

# =========================================================

cat("STEP 4: Generating Figure 6 volcano plot...\n")

label_genes <- res %>%
filter(
!is.na(padj),
!is.na(log2FoldChange),
padj < PADJ_CUTOFF
)

volcano_plot <- ggplot(
res %>% filter(!is.na(padj), !is.na(log2FoldChange)),
aes(
x = log2FoldChange,
y = neg_log10_padj,
color = significance
)
) +
geom_point(
alpha = 0.70,
size = 1.5
) +
geom_vline(
xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF),
linetype = "dashed",
linewidth = 0.4
) +
geom_hline(
yintercept = -log10(PADJ_CUTOFF),
linetype = "dashed",
linewidth = 0.4
) +
geom_label_repel(
data = label_genes,
aes(label = gene_label),
size = 3,
max.overlaps = Inf,
min.segment.length = 0,
box.padding = 0.4,
point.padding = 0.2,
show.legend = FALSE
) +
scale_color_manual(
values = c(
"Upregulated" = "red",
"Downregulated" = "blue",
"Significant, smaller effect" = "purple",
"Not significant" = "gray70"
),
name = "Gene classification"
) +
theme_minimal() +
labs(
title = "Volcano Plot: Male Spinal Cord at 14 Days Post-Injury",
subtitle = "Combined GSE306455 and GSE202166 analysis",
x = expression(log[2]~fold~change),
y = expression(-log[10]~adjusted~italic(p)*-value)
)

ggsave(
filename = file.path(
OUTPUT_DIR,
"figure_6_volcano_plot_14dpi.png"
),
plot = volcano_plot,
width = 10,
height = 7,
dpi = 300
)

cat("Saved Figure 6 volcano plot.\n\n")

# =========================================================

# 5. Generate Figure 7: Horizontal bar chart

# =========================================================

cat("STEP 5: Generating Figure 7 horizontal bar chart...\n")

bar_data <- significant_genes %>%
mutate(
gene_symbol = factor(
gene_symbol,
levels = gene_symbol[order(log2FoldChange)]
),


label_position = if_else(
  log2FoldChange >= 0,
  log2FoldChange + 0.12,
  log2FoldChange - 0.12
),

label_hjust = if_else(
  log2FoldChange >= 0,
  0,
  1
)


)

bar_plot <- ggplot(
bar_data,
aes(
x = gene_symbol,
y = log2FoldChange,
fill = classification
)
) +
geom_hline(
yintercept = 0,
linewidth = 0.4
) +
geom_col(
width = 0.7
) +
geom_text(
aes(
y = label_position,
label = sprintf("%.2f", log2FoldChange),
hjust = label_hjust
),
size = 3.5
) +
coord_flip() +
scale_fill_manual(
values = c(
"Upregulated" = "red",
"Downregulated" = "blue",
"Significant, smaller effect" = "purple"
),
name = "Classification"
) +
scale_y_continuous(
expand = expansion(mult = c(0.15, 0.20))
) +
theme_minimal() +
labs(
title = "Significant Genes at 14 Days Post-Injury",
subtitle = "Male spinal cord: combined GSE306455 and GSE202166 analysis",
x = "Gene",
y = expression(log[2]~fold~change)
)

ggsave(
filename = file.path(
OUTPUT_DIR,
"figure_7_significant_gene_log2FC_bar_chart_14dpi.png"
),
plot = bar_plot,
width = 9,
height = 6,
dpi = 300
)

cat("Saved Figure 7 horizontal bar chart.\n\n")

# =========================================================

# 6. Generate Table 2 as PNG

# =========================================================

cat("STEP 6: Generating Table 2...\n")

table_2_display <- significant_genes %>%
mutate(
log2FoldChange = round(log2FoldChange, 2),
padj = format(padj, scientific = TRUE, digits = 3)
) %>%
rename(
`Gene symbol` = gene_symbol,
`Ensembl ID` = ensembl_gene_id,
`log2 fold change` = log2FoldChange,
`Adjusted p-value` = padj,
Classification = classification
)

table_grob <- tableGrob(
table_2_display,
rows = NULL,
theme = ttheme_minimal(
base_size = 10,
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

table_png_file <- file.path(
OUTPUT_DIR,
"table_2_significant_genes_14dpi.png"
)

png(
filename = table_png_file,
width = 12,
height = 5,
units = "in",
res = 300
)

grid.newpage()
grid.draw(table_grob)

dev.off()

cat("Saved Table 2 as CSV and PNG.\n\n")

# =========================================================

# 7. Save summary statistics

# =========================================================

cat("STEP 7: Saving analysis summary...\n")

summary_table <- tibble(
metric = c(
"genes_tested",
"significant_genes_padj_0.05",
"upregulated_log2FC_gt_1",
"downregulated_log2FC_lt_minus_1",
"significant_smaller_effect_genes"
),
value = c(
nrow(res),
sum(res$padj < PADJ_CUTOFF, na.rm = TRUE),
sum(
res$padj < PADJ_CUTOFF &
res$log2FoldChange > LOG2FC_CUTOFF,
na.rm = TRUE
),
sum(
res$padj < PADJ_CUTOFF &
res$log2FoldChange < -LOG2FC_CUTOFF,
na.rm = TRUE
),
sum(
res$padj < PADJ_CUTOFF &
abs(res$log2FoldChange) <= LOG2FC_CUTOFF,
na.rm = TRUE
)
)
)

write_csv(
summary_table,
file.path(
OUTPUT_DIR,
"summary_14dpi_combined_analysis.csv"
)
)

print(summary_table)

cat("\n============================================\n")
cat("14 DPI VISUALIZATION PIPELINE COMPLETE\n")
cat("============================================\n")
