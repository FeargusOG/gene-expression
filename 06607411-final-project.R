library(dplyr)
library(tidyr)

# Step 3. Read the RNA-seq file
data_rnaseq <- read.delim("brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt")

# Step 4. Read the patient data. Skip first 4 rows as they don't contain observations.
data_patient <- read.delim("brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt", skip = 4)

# Step 5. Read the Copy Number Aberrations Data
data_cna <- read.delim("brca_tcga_pan_can_atlas_2018/data_cna.txt")

# Step 6. Match the IDS
## First match the cols between RNASeq and CNA
common_cols <- intersect(colnames(data_rnaseq), colnames(data_cna))
data_rnaseq <- data_rnaseq[, common_cols]
data_cna <- data_cna[, common_cols]

## Filter rows where Hugo_Symbol equals "ERBB2"
erbb2_cna <- data_cna[data_cna$Hugo_Symbol == "ERBB2", ] |>
  select(-Hugo_Symbol, -Entrez_Gene_Id)
rownames(erbb2_cna) <- "ERBB2"

## Print the counts
##   ERBB2 <= 0 (Not Amplified)
##   ERBB2 > 0 (Amplified)
print(table(unlist(erbb2_cna)))

## Extract assay data, skipping the first two columns (gene identifiers)
assay <- round(as.matrix(data_rnaseq[, -c(1, 2)]))
rownames(assay) <- data_rnaseq[, 1]

# Step 7. Create Metadata
## Build metadata for analysis
metadata <- matrix(0, ncol = 1, nrow = ncol(assay))

## Loop over each column in the assay dataset.
for (i in seq_len(ncol(assay))) {
  ## Find the patient with the matching ID for this col in the assay
  patient_col <- erbb2_cna |> 
    select(colnames(assay)[i])
  ## Check the ERBB2 amplification level (greater than 0 is amplified)
  metadata[i, 1] <- ifelse(as.numeric(patient_col[1,1]) > 0, 1, 0)
}

# Replace NA values in metadata with 0
metadata[is.na(metadata)] <- 0
colnames(metadata) <- c("ERBB2_Amp")

# Prepare RNASeq data for DESeq2
assay[is.na(assay)] <- 0  # Replace NA values with zeros
assay[assay < 0] <- 0     # Set negative values to zero

# Filter out genes with low expression
smallest_group_size <- 3
keep <- rowSums(assay >= 10) >= smallest_group_size
assay <- assay[keep, ]

# Ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Check and install each package individually
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

library(DESeq2)

# differential expression analysis.
#   assay:    each row is a gene and each column a patient.
#   metadata: each row is a patient and the single column is a grouping factor.
# Therefore, the rows in metadata must match the columns in countData in both number and order!
#
# Basically, what this function does is say "here is a bunch of patients with their genome data. 
# Here is a grouping factor (could be age, could be ERB22+ expression, whatever). What patterns 
# can be see in their genome data from this grouping factor?"
dds <- DESeqDataSetFromMatrix(countData = assay,
                              colData = metadata,
                              design = ~ ERBB2_Amp)

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

# Display top 10 differentially expressed genes by adjusted p-value
top_genes <- res[order(res$padj)[1:10], ]
print(top_genes)


# Perform PCA
vsd <- vst(dds)
plotPCA(vsd, intgroup = c("ERBB2_Amp"))


# Generate heatmap for top differentially expressed genes
if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")
library(pheatmap)

top_DE <- order(res$padj)
vsd_DE <- assay(vsd)[top_DE[1:20],]

annotation_col <- data.frame(Early = as.matrix(metadata[, 1]))
rownames(annotation_col) <- colnames(vsd)

pheatmap(vsd_DE, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row",
         show_colnames = FALSE, show_rownames = TRUE, annotation_col = annotation_col)
