print_banner("Preprocessing...")

library(dplyr)

process_cna_gene <- function(data_cna, gene_symbol) {
  gene_cna <- data_cna[data_cna$Hugo_Symbol == gene_symbol, ] |>
    dplyr::select(-Hugo_Symbol, -Entrez_Gene_Id)
  rownames(gene_cna) <- gene_symbol
  return(gene_cna)
}

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
## The patients contained in RNASeq and CNA datasets don't match fully; because
## we want to examine the relationship between patients in these two datasets,
## lets drop (the very few) patients from the datasets that aren't in both.
common_cols <- intersect(colnames(data_rnaseq), colnames(data_cna))
data_rnaseq <- data_rnaseq[, common_cols]
data_cna <- data_cna[, common_cols]

# 7. Create metadata using the CNA level of ERBB2+
cna_erbb2 <- process_cna_gene(data_cna, "ERBB2")

## Print the counts of each ERBB2 value.
##   - ERBB2 <= 0 (Not Amplified)
##   - ERBB2 > 0 (Amplified)
cat("ERBB2 CNA Values:")
print(table(unlist(cna_erbb2)))

## Extract assay data, skipping the first two columns (gene identifiers).
assay <- round(as.matrix(data_rnaseq[, -c(1, 2)]))
rownames(assay) <- data_rnaseq[, 1]

## Create an empty matrix to store our metadata.
metadata <- matrix(0, ncol = 1, nrow = ncol(assay))

## Loop over each column in the assay dataset.
for (i in seq_len(ncol(assay))) {
  ## Find the patient with the matching ID in the CNA ERB22 dataset.
  patient_col <- cna_erbb2 |> 
    dplyr::select(colnames(assay)[i])
  ## Check the ERBB2 amplification level (greater than 0 is amplified)
  metadata[i, 1] <- ifelse(as.numeric(patient_col[1,1]) > amplification_threshold, 1, 0)
}

## Replace NA values in metadata with 0
metadata[is.na(metadata)] <- 0
ERBB2_Amp_Col <- "ERBB2_Amp"
colnames(metadata) <- c(ERBB2_Amp_Col)
