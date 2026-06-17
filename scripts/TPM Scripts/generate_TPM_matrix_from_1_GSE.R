setwd("~/Desktop/Professional/Development/Masters Data Analysis Scripts/SNI-DGE-RAT-MUS")

# =========================================================
# User-defined settings
# =========================================================

ANALYSIS_NAME <- "SNI_SC_MUS_DPI_7_MALE_GSE202166"

METADATA_FILE <- paste0("metadata/", ANALYSIS_NAME, "_METADATA.CSV")

TPM_OUTPUT_FILE <- "final_matrix/gene_level_TPM_matrix.csv"
COUNTS_OUTPUT_FILE <- "final_matrix/gene_level_counts_matrix.csv"

cat("\n==============================\n")
cat("Starting Gene-Level Matrix Generation Pipeline\n")
cat("==============================\n\n")

library(readr)
library(dplyr)
library(tibble)
library(tximport)
library(biomaRt)

cat("Libraries loaded successfully.\n\n")

# =========================================================
# 1. Load sample metadata
# =========================================================

cat("STEP 1: Loading sample metadata...\n")

samples <- read_csv(METADATA_FILE)

cat("Metadata file:", METADATA_FILE, "\n")
cat("Metadata loaded.\n")
cat("Number of samples:", nrow(samples), "\n\n")

# =========================================================
# 2. Locate Salmon quant files
# =========================================================

cat("STEP 2: Locating Salmon quant.sf files...\n")

files <- file.path(
  "quant",
  samples$gse,
  samples$sample_id,
  "quant.sf"
)

names(files) <- samples$sample_id

cat("Checking file existence...\n")
print(file.exists(files))

missing_files <- files[!file.exists(files)]

if (length(missing_files) > 0) {
  stop(
    "\nERROR: Missing quant.sf files:\n",
    paste(missing_files, collapse = "\n")
  )
}

cat("All quant.sf files located successfully.\n\n")

# =========================================================
# 3. Build transcript-to-gene mapping
# =========================================================

cat("STEP 3: Connecting to Ensembl BioMart mirror...\n")

mart <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl",
  mirror = "asia"
)

cat("Connected to Ensembl.\n")
cat("Downloading transcript-to-gene mapping...\n")

tx2gene <- getBM(
  attributes = c(
    "ensembl_transcript_id_version",
    "ensembl_gene_id",
    "external_gene_name"
  ),
  mart = mart
)

cat("Transcript-to-gene mapping downloaded.\n")
cat("Rows in tx2gene:", nrow(tx2gene), "\n")

write_csv(
  tx2gene,
  "gene_mapping/ensembl_tx2gene_mapping.csv"
)

cat("Saved transcript-to-gene mapping.\n\n")

tx2gene_simple <- tx2gene %>%
  dplyr::select(
    ensembl_transcript_id_version,
    ensembl_gene_id
  ) %>%
  distinct()

# =========================================================
# 4. Import Salmon results
# =========================================================

cat("STEP 4: Importing Salmon quantification files...\n")

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene_simple,
  ignoreTxVersion = FALSE
)

cat("Salmon quantifications imported successfully.\n")
cat("Genes imported:", nrow(txi$abundance), "\n")
cat("Samples imported:", ncol(txi$abundance), "\n\n")

# =========================================================
# 5. Create TPM and counts matrices
# =========================================================

cat("STEP 5: Creating TPM and counts matrices...\n")

gene_tpm <- as.data.frame(txi$abundance) %>%
  rownames_to_column("ensembl_gene_id")

gene_counts <- as.data.frame(txi$counts) %>%
  rownames_to_column("ensembl_gene_id")

cat("TPM matrix dimensions:", dim(gene_tpm), "\n")
cat("Counts matrix dimensions:", dim(gene_counts), "\n\n")

# =========================================================
# 6. Add gene symbols
# =========================================================

cat("STEP 6: Adding gene symbols...\n")

gene_symbols <- tx2gene %>%
  dplyr::select(
    ensembl_gene_id,
    external_gene_name
  ) %>%
  distinct() %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    gene_symbol = first(na.omit(external_gene_name)),
    .groups = "drop"
  )

gene_tpm <- gene_tpm %>%
  left_join(gene_symbols, by = "ensembl_gene_id") %>%
  relocate(gene_symbol, .after = ensembl_gene_id)

gene_counts <- gene_counts %>%
  left_join(gene_symbols, by = "ensembl_gene_id") %>%
  relocate(gene_symbol, .after = ensembl_gene_id)

cat("Gene symbols added successfully.\n\n")

# =========================================================
# 7. Save TPM and counts matrices
# =========================================================

cat("STEP 7: Saving gene-level matrices...\n")

write_csv(
  gene_tpm,
  TPM_OUTPUT_FILE
)

write_csv(
  gene_counts,
  COUNTS_OUTPUT_FILE
)

cat("Saved TPM matrix to:", TPM_OUTPUT_FILE, "\n")
cat("Saved counts matrix to:", COUNTS_OUTPUT_FILE, "\n\n")

cat("==============================\n")
cat("PIPELINE COMPLETE\n")
cat("==============================\n")