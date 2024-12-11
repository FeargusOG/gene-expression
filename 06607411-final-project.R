ensure_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

library(dplyr)
library(tidyr)

# 1. Download the dataset on: https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018
## Performed outside of this script....

# 2. Untar the folder and extract the files.
## Performed outside of this script....

# 3. Read the RNA-seq file: data_mrna_seq_v2_rsem.txt
data_rnaseq <- read.delim("brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt")

# 4. Read the Patient Data file: data_clinical_patient.txt
## Skip first 4 rows as they don't contain observations.
data_patient <- read.delim("brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt", skip = 4)

# 5. Read the Copy Number Aberrations Data: data_cna.txt
data_cna <- read.delim("brca_tcga_pan_can_atlas_2018/data_cna.txt")

# 6. Match the RNA-seq patient ids with the CNA ids and the Patient Data ids.

## In both data_rnaseq and data_cna the patient IDs are formatte with dots:
##   - TCGA.A2.A0D3.01
## Whereas, in data_patient, they are formatted with dashes:
##   - TCGA-A2-A0D3
## TBD!!! Not sure about the .01, or whether I even need to match them here...

# 7. Create metadata using the CNA level of ERBB2+
## The patients contained in RNASeq and CNA datasets don't match fully; because
## we want to examine the relationship between patients in these two datasets,
## lets drop (the very few) patients from the datasets that aren't in both.
common_cols <- intersect(colnames(data_rnaseq), colnames(data_cna))
data_rnaseq <- data_rnaseq[, common_cols]
data_cna <- data_cna[, common_cols]

## Extract the "ERBB2" row from the CNA dataset and filter it so that we just
## have the patient data.
erbb2_cna <- data_cna[data_cna$Hugo_Symbol == "ERBB2", ] |>
  select(-Hugo_Symbol, -Entrez_Gene_Id)
rownames(erbb2_cna) <- "ERBB2"

## Print the counts of each ERBB2 value.
##   - ERBB2 <= 0 (Not Amplified)
##   - ERBB2 > 0 (Amplified)
print(table(unlist(erbb2_cna)))

## Extract assay data, skipping the first two columns (gene identifiers).
assay <- round(as.matrix(data_rnaseq[, -c(1, 2)]))
rownames(assay) <- data_rnaseq[, 1]

## Create an empty matrix to store our metadata.
metadata <- matrix(0, ncol = 1, nrow = ncol(assay))

## Loop over each column in the assay dataset.
for (i in seq_len(ncol(assay))) {
  ## Find the patient with the matching ID in the CNA ERB22 dataset.
  patient_col <- erbb2_cna |> 
    select(colnames(assay)[i])
  ## Check the ERBB2 amplification level (greater than 0 is amplified)
  metadata[i, 1] <- ifelse(as.numeric(patient_col[1,1]) > 0, 1, 0)
}

## Replace NA values in metadata with 0
metadata[is.na(metadata)] <- 0
ERBB2_Amp_Col <- "ERBB2_Amp"
colnames(metadata) <- c(ERBB2_Amp_Col)

# 8 - Normalize data using DESeq2
## Prepare RNASeq data for DESeq2
assay[is.na(assay)] <- 0  # Replace NA values with zeros
assay[assay < 0] <- 0     # Set negative values to zero

## Filter out genes with low expression
## TODO Why these numbers?? Justify....
smallest_group_size <- 3
keep <- rowSums(assay >= 10) >= smallest_group_size
assay <- assay[keep, ]

## Ensure BiocManager and DESeq2 are installed
ensure_package("BiocManager")
ensure_package("DESeq2")
library(DESeq2)

## Create the differential expression analysis dataset.
##   assay:    each row is a gene and each column a patient.
##   metadata: each row is a patient and the single column is a grouping factor.
##
## Therefore, the rows in metadata must match the columns in countData in both 
## number and order!
##
## Basically, what this function does is say "here is a bunch of patients with 
## their genome data. Here is a grouping factor (could be age, could be ERB22+ 
## expression, whatever). What patterns can be see in their genome data from 
## this grouping factor?"
dds <- DESeqDataSetFromMatrix(countData = assay,
                              colData = metadata,
                              design = ~ ERBB2_Amp)

## Run DESeq2 differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

# 9. Obtain Differentially Expressed Genes
## Display top 10 differentially expressed genes by adjusted p-value
top_genes <- res[order(res$padj)[1:10], ]
print(top_genes)

# 10. Perform a Pathway Enrichment Analysis
## Gene Ontology (GO) enrichment analysis
## KEGG Pathway Enrichment Analysis

# 11. Get the variance stabilised transformed expression values.
vsd <- vst(dds)

# 12. With the vst values obtain a PCA plot and a heatmap.
## PCA Plot
plotPCA(vsd, intgroup = c(ERBB2_Amp_Col))

## Generate heatmap for top differentially expressed genes
ensure_package("pheatmap")
library(pheatmap)

top_DE <- order(res$padj)
vsd_DE <- assay(vsd)[top_DE[1:20],]

annotation_col <- data.frame(ERBB2_Amp = as.matrix(metadata[, 1]))
rownames(annotation_col) <- colnames(vsd)

pheatmap(vsd_DE, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row",
         show_colnames = FALSE, show_rownames = TRUE, annotation_col = annotation_col)
