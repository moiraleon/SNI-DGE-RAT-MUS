setwd("~/Desktop/Professional/Development/Masters Data Analysis Scripts/SNI-DGE-RAT-MUS")

cat("\n============================================\n")
cat("Starting 7 DPI Shared-Gene Visualization Pipeline\n")
cat("============================================\n\n")

# =========================================================

# Package installation: run once if needed

# =========================================================
# options(repos = c(CRAN = "https://cloud.r-project.org"))


# install.packages(c(
#   "readr",
#   "dplyr",
#   "tidyr",
#   "tibble",
#   "ggplot2",
#   "VennDiagram",
#   "gridExtra"
# ))

# =========================================================

# Load libraries

# =========================================================

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(VennDiagram)
library(grid)
library(gridExtra)

cat("Libraries loaded successfully.\n\n")

# =========================================================

# User-defined settings

# =========================================================

# IMPORTANT:

# Use the COMPLETE DESeq2 results files here, not the files

# containing only significant genes.

FEMALE_RESULTS_FILE <- "final_matrix/PHASE_2/COMB_1/FEMALE_GSE241361/deseq2_results_SNI_SC_MUS_DPI_7_FEMALE_GSE241361_SNI_vs_CTRL.csv"

MALE_RESULTS_FILE <- "final_matrix/PHASE_2/COMB_1/MALE_GSE202166/deseq2_results_SNI_SC_MUS_DPI_7_MALE_GSE202166_SNI_vs_CTRL.csv"

OUTPUT_DIR <- "VIS_temporal_analysis/COMB_1_7DPI_Female_GSE241361_vs_Male_GSE202166"

PADJ_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 1

FEMALE_DATASET <- "GSE241361"
MALE_DATASET <- "GSE202166"

FEMALE_SEX <- "Female"
MALE_SEX <- "Male"

TISSUE <- "Spinal cord"
DPI <- 7

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

# Helper functions

# =========================================================

check_required_columns <- function(data, file_label) {

required_columns <- c(
"ensembl_gene_id",
"gene_symbol",
"log2FoldChange",
"padj"
)

missing_columns <- setdiff(required_columns, colnames(data))

if (length(missing_columns) > 0) {
stop(
"ERROR: ", file_label,
" is missing the following required columns:\n",
paste(missing_columns, collapse = ", ")
)
}
}

prepare_results <- function(file_path, dataset_label) {

data <- read_csv(file_path, show_col_types = FALSE)

check_required_columns(data, dataset_label)

data %>%
mutate(
# Removes Ensembl version suffixes, if present:
# ENSMUSG00000012345.1 becomes ENSMUSG00000012345
ensembl_clean = sub("[.].*$", "", ensembl_gene_id),


  # Uses Ensembl ID only if a gene symbol is missing
  gene_label = if_else(
    is.na(gene_symbol) | gene_symbol == "",
    ensembl_clean,
    gene_symbol
  ),
  
  dataset = dataset_label
)


}

count_upregulated <- function(data) {
sum(
data$padj < PADJ_CUTOFF &
data$log2FoldChange > LOG2FC_CUTOFF,
na.rm = TRUE
)
}

count_downregulated <- function(data) {
sum(
data$padj < PADJ_CUTOFF &
data$log2FoldChange < -LOG2FC_CUTOFF,
na.rm = TRUE
)
}

# =========================================================

# 1. Load DESeq2 results

# =========================================================

cat("STEP 1: Loading DESeq2 results...\n")

female_full <- prepare_results(
FEMALE_RESULTS_FILE,
FEMALE_DATASET
)

male_full <- prepare_results(
MALE_RESULTS_FILE,
MALE_DATASET
)

cat("Female results loaded:", nrow(female_full), "genes\n")
cat("Male results loaded:", nrow(male_full), "genes\n\n")

# =========================================================

# 2. Select statistically significant genes

# =========================================================

cat("STEP 2: Filtering statistically significant genes...\n")

female_sig <- female_full %>%
filter(
!is.na(padj),
!is.na(log2FoldChange),
padj < PADJ_CUTOFF
) %>%
arrange(padj) %>%
distinct(ensembl_clean, .keep_all = TRUE)

male_sig <- male_full %>%
filter(
!is.na(padj),
!is.na(log2FoldChange),
padj < PADJ_CUTOFF
) %>%
arrange(padj) %>%
distinct(ensembl_clean, .keep_all = TRUE)

cat("Female significant genes:", nrow(female_sig), "\n")
cat("Male significant genes:", nrow(male_sig), "\n\n")

# =========================================================

# 3. Identify shared significant genes

# =========================================================

cat("STEP 3: Identifying shared significant genes...\n")

female_ids <- unique(female_sig$ensembl_clean)
male_ids <- unique(male_sig$ensembl_clean)

shared_ids <- intersect(female_ids, male_ids)

female_shared <- female_sig %>%
filter(ensembl_clean %in% shared_ids) %>%
select(
ensembl_clean,
female_gene_symbol = gene_label,
female_log2FC = log2FoldChange,
female_padj = padj
)

male_shared <- male_sig %>%
filter(ensembl_clean %in% shared_ids) %>%
select(
ensembl_clean,
male_gene_symbol = gene_label,
male_log2FC = log2FoldChange,
male_padj = padj
)

shared_genes <- full_join(
female_shared,
male_shared,
by = "ensembl_clean"
) %>%
mutate(
gene_symbol = coalesce(female_gene_symbol, male_gene_symbol)
) %>%
select(
gene_symbol,
ensembl_gene_id = ensembl_clean,
female_log2FC,
female_padj,
male_log2FC,
male_padj
) %>%
arrange(gene_symbol)

cat("Shared significant genes:", nrow(shared_genes), "\n")
print(shared_genes)
cat("\n")

write_csv(
shared_genes,
file.path(OUTPUT_DIR, "shared_significant_genes_7dpi.csv")
)

# =========================================================

# 4. Generate Venn diagram

# =========================================================

cat("STEP 4: Generating Venn diagram...\n")

shared_gene_label <- if (nrow(shared_genes) == 0) {
"No shared genes"
} else {
paste(shared_genes$gene_symbol, collapse = "\n")
}

venn_output_file <- file.path(
OUTPUT_DIR,
"figure_4_venn_shared_significant_genes_7dpi.png"
)

png(
filename = venn_output_file,
width = 8,
height = 7,
units = "in",
res = 300
)

grid.newpage()

venn_grob <- draw.pairwise.venn(
area1 = length(female_ids),
area2 = length(male_ids),
cross.area = length(shared_ids),
category = c(
paste0(FEMALE_SEX, ": ", FEMALE_DATASET),
paste0(MALE_SEX, ": ", MALE_DATASET)
),
fill = c("#F8766D", "#00BFC4"),
alpha = c(0.50, 0.50),
lty = "blank",
cex = 1.2,
cat.cex = 1.0,
cat.pos = c(180, 0),
cat.dist = c(0.06, 0.06),
scaled = FALSE,
ind = FALSE
)

grid.draw(venn_grob)

# Adjust y-coordinate slightly if this text overlaps your Venn counts.

grid.text(
label = paste0(
"Shared genes (n = ", length(shared_ids), "):"
),
x = 0.5,
y = 0.39,
gp = gpar(fontsize = 10, fontface = "plain")
)

grid.text(
label = shared_gene_label,
x = 0.5,
y = 0.33,
gp = gpar(fontsize = 11, fontface = "italic")
)

grid.text(
label = "Significant genes defined as adjusted p < 0.05",
x = 0.5,
y = 0.06,
gp = gpar(fontsize = 9)
)

dev.off()

cat("Saved Venn diagram:\n")
cat(venn_output_file, "\n\n")

# =========================================================

# 5. Generate grouped bar chart for shared genes

# =========================================================

cat("STEP 5: Generating shared-gene log2FC bar chart...\n")

if (nrow(shared_genes) == 0) {

warning(
"No shared significant genes were found. ",
"The grouped bar chart will not be generated."
)

} else {

bar_data <- shared_genes %>%
select(
gene_symbol,
female_log2FC,
male_log2FC
) %>%
pivot_longer(
cols = c(female_log2FC, male_log2FC),
names_to = "dataset_group",
values_to = "log2FoldChange"
) %>%
mutate(
dataset_group = recode(
dataset_group,
female_log2FC = paste0(FEMALE_SEX, ": ", FEMALE_DATASET),
male_log2FC = paste0(MALE_SEX, ": ", MALE_DATASET)
),
label_y = log2FoldChange +
if_else(log2FoldChange >= 0, 0.18, -0.18)
)

preferred_gene_order <- c("Gpr151", "Atf3")

gene_order <- c(
intersect(preferred_gene_order, unique(bar_data$gene_symbol)),
setdiff(sort(unique(bar_data$gene_symbol)), preferred_gene_order)
)

bar_data <- bar_data %>%
mutate(
gene_symbol = factor(
gene_symbol,
levels = gene_order
)
)

#setting default colors
fill_colors <- setNames(
  c("#F8766D", "#00BFC4"),
  c(
    paste0(FEMALE_SEX, ": ", FEMALE_DATASET),
    paste0(MALE_SEX, ": ", MALE_DATASET)
  )
)

fold_change_plot <- ggplot(
bar_data,
aes(
x = gene_symbol,
y = log2FoldChange,
fill = dataset_group
)
) +
geom_hline(
yintercept = 0,
linewidth = 0.4
) +
geom_col(
position = position_dodge(width = 0.75),
width = 0.65
) +
geom_text(
aes(
y = label_y,
label = sprintf("%.2f", log2FoldChange)
),
position = position_dodge(width = 0.75),
size = 3.5
) +
scale_fill_manual(
values = fill_colors
) +
scale_y_continuous(
expand = expansion(mult = c(0.10, 0.20))
) +
theme_minimal() +
labs(
title = "Shared Significant Genes at 7 Days Post-Injury",
subtitle = "Log2 fold changes for genes significant in both datasets",
x = "Gene",
y = expression(log[2]~fold~change),
fill = "Dataset / sex"
)

ggsave(
filename = file.path(
OUTPUT_DIR,
"figure_5_shared_gene_log2FC_bar_chart_7dpi.png"
),
plot = fold_change_plot,
width = 7,
height = 5,
dpi = 300
)

cat("Saved grouped bar chart.\n\n")
}

# =========================================================

# 6. Generate Table 1 summary

# =========================================================

cat("STEP 6: Generating Table 1 summary...\n")

shared_gene_summary <- if (nrow(shared_genes) == 0) {
"None"
} else {
paste(shared_genes$gene_symbol, collapse = "; ")
}

table_1 <- tibble(
Dataset = c(FEMALE_DATASET, MALE_DATASET),
Sex = c(FEMALE_SEX, MALE_SEX),
Tissue = c(TISSUE, TISSUE),
DPI = c(DPI, DPI),
Genes_tested = c(
nrow(female_full),
nrow(male_full)
),
Significant_genes_padj_lt_0.05 = c(
nrow(female_sig),
nrow(male_sig)
),
Upregulated_padj_lt_0.05_log2FC_gt_1 = c(
count_upregulated(female_full),
count_upregulated(male_full)
),
Downregulated_padj_lt_0.05_log2FC_lt_minus_1 = c(
count_downregulated(female_full),
count_downregulated(male_full)
),
Shared_significant_genes = c(
shared_gene_summary,
shared_gene_summary
)
)

print(table_1)

write_csv(
table_1,
file.path(
OUTPUT_DIR,
"table_1_7dpi_single_dataset_DEG_summary.csv"
)
)

# Create a cleaner version specifically for the PNG figure
table_1_display <- table_1 %>%
  rename(
    `GEO dataset` = Dataset,
    Sex = Sex,
    Tissue = Tissue,
    `Days post-injury` = DPI,
    `Genes tested` = Genes_tested,
    `Significant genes\n(padj < 0.05)` = Significant_genes_padj_lt_0.05,
    `Upregulated genes\n(padj < 0.05; log2FC > 1)` =
      Upregulated_padj_lt_0.05_log2FC_gt_1,
    `Downregulated genes\n(padj < 0.05; log2FC < -1)` =
      Downregulated_padj_lt_0.05_log2FC_lt_minus_1,
    `Shared significant genes` = Shared_significant_genes
  )

table_grob <- tableGrob(
  table_1_display,
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
  "table_1_7dpi_single_dataset_DEG_summary.png"
)

png(
  filename = table_png_file,
  width = 18,
  height = 3.5,
  units = "in",
  res = 300
)

grid.newpage()
grid.draw(table_grob)

dev.off()

cat("Saved Table 1 summary as CSV and PNG.\n\n")

cat("============================================\n")
cat("7 DPI SHARED-GENE VISUALIZATION COMPLETE\n")
cat("============================================\n")
