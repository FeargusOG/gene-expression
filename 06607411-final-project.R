ensure_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

library(dplyr)
library(tidyr)
library(ggplot2)

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
  dplyr::select(-Hugo_Symbol, -Entrez_Gene_Id)
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
    dplyr::select(colnames(assay)[i])
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

## Filter out genes with low expression or expressed in few samples; this
## protects our analysis from spurious results due to noise etc.
smallest_group_size <- 3
keep <- rowSums(assay >= 10) >= smallest_group_size
assay <- assay[keep, ]

## Ensure BiocManager and DESeq2 are installed
ensure_package("BiocManager")
ensure_package("BiocParallel") # DESeq2 was slow! Let's use more cores....
ensure_package("DESeq2")
library(BiocParallel)
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
## Use one less core than the total, but ensure at least 1 worker
workers <- max(1, parallel::detectCores() - 2)
register(MulticoreParam(workers = workers))
dds <- DESeq(dds, parallel = TRUE)
res <- results(dds)

# 9. Obtain Differentially Expressed Genes
## Display top 10 differentially expressed genes by Fold Change.
## The genes listed will have either a positive or negative value for 
## `log2FoldChange` depending on whether they are up-regulated or 
## down-regulated in the presence of ERBB2+, respectively.
top_genes <- res[order(abs(res$log2FoldChange), decreasing = TRUE)[1:10],]
print(top_genes)

# 10. Perform a Pathway Enrichment Analysis
## So from our above DESeq analysis, we have found a bunch of genes that are
## either over or under expressed. Now with enrichment analysis we can determine
## what biological processes are invovled in these genes and how statistically
## significant the results are. For example, ERBB2 amplification might cause a
## bunch of other genes related to the cell cycle to become overexpressed. This 
## enrichment analysis would say "we're seeing a lot of cell cycle processes
## being impacted in patients with ERBB2 amplification; more than would be
## statistically likely by chance."
ensure_package("clusterProfiler")
library(clusterProfiler)
ensure_package("org.Hs.eg.db")
library(org.Hs.eg.db)
ensure_package("enrichplot")
library(enrichplot)

## Use the standard significance threshold to test the Benjamini-Hochberg
## adjusted p-values.
SIGNIFICANCE_LEVEL <- 0.05
res_sig = res[res$padj < SIGNIFICANCE_LEVEL,]
DE_over = rownames(res_sig[res_sig$log2FoldChange>0,])
DE_under = rownames(res_sig[res_sig$log2FoldChange<0,])

## Gene Ontology (GO) enrichment analysis
### Find genes that are over-expressed in ERBB2 amplified patients.
go_results_over = enrichGO(
  gene          = DE_over,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP", # Biological Process
  qvalueCutoff  = SIGNIFICANCE_LEVEL # Use a stricter threshold than 0.2 to 
                                     # provide more confidence in our results.
)

### Find genes that are under-expressed in ERBB2 amplified patients.
go_results_under = enrichGO(
  gene          = DE_under,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",  
  ont           = "BP",
  qvalueCutoff  = SIGNIFICANCE_LEVEL
)

print(head(go_results_over))
print(head(go_results_under))

### Plot results of GO enrichment analysis
dotplot(go_results_over) + 
  ggtitle("Gene Ontology Enrichment Over Expressed")

dotplot(go_results_under) + 
  ggtitle("Gene Ontology Enrichment Under Expressed")

## KEGG Pathway Enrichment Analysis
### Map gene symbol to gene entrez ID. 
### This is necessary for Reactome and Keggs.
gene_entrez_over <- bitr(
  DE_over,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

gene_entrez_under <- bitr(
  DE_under,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

### defaults to 'hsa' for organism i.e. home sapien
kegg_results_over =  enrichKEGG(
  gene          = gene_entrez_over[,2],
  qvalueCutoff  = SIGNIFICANCE_LEVEL
)

kegg_results_under =  enrichKEGG(
  gene          = gene_entrez_under[,2],
  qvalueCutoff  = SIGNIFICANCE_LEVEL
)

print(head(kegg_results_over))
print(head(kegg_results_under))

### Plot results of KEGG pathway enrichment analysis
dotplot(kegg_results_over) + 
  ggtitle("Kegg Pathway Enrichment Over Expressed")
dotplot(kegg_results_under) + 
  ggtitle("Kegg Pathway Enrichment Under Expressed")

## Reactome Pathway Enrichment Analysis
ensure_package("ReactomePA")
library(ReactomePA)
ensure_package("pathview")
library(pathview)

### defaults to 'human' for organism
reactome_results_over =  enrichPathway(
  gene          = gene_entrez_over[,2],
  qvalueCutoff  = SIGNIFICANCE_LEVEL,
)

reactome_results_under =  enrichPathway(
  gene          = gene_entrez_under[,2],
  qvalueCutoff  = SIGNIFICANCE_LEVEL,
)

print(head(reactome_results_over))
print(head(reactome_results_under))

### Plot results of Reactome pathway enrichment analysis
dotplot(reactome_results_over, showCategory=10) + 
  ggtitle("Reactome Pathway Enrichment Over Expressed")
dotplot(reactome_results_under, showCategory=10) + 
  ggtitle("Reactome Pathway Enrichment Under Expressed")

### TREE PLOTS
go_results_over_pw = pairwise_termsim(go_results_over)
treeplot(go_results_over_pw)+ ggtitle("GO Enrichment Over Expressed")

go_results_under_pw = pairwise_termsim(go_results_under)
treeplot(go_results_under_pw)+ ggtitle("GO Enrichment Under Expressed")

kegg_results_over_pw = pairwise_termsim(kegg_results_over)
treeplot(kegg_results_over_pw)+ ggtitle("KEGG Enrichment Over Expressed")

kegg_results_under_pw = pairwise_termsim(kegg_results_under)
treeplot(kegg_results_under_pw)+ ggtitle("KEGG Enrichment Under Expressed")

reactome_results_over_pw = pairwise_termsim(reactome_results_over)
treeplot(reactome_results_over_pw)+ ggtitle("Reactome Enrichment Over Expressed")

reactome_results_under_pw = pairwise_termsim(reactome_results_under)
treeplot(reactome_results_under_pw)+ ggtitle("Reactome Enrichment Under Expressed")

# # Extract similarity matrix as a data frame
# similarity_matrix <- reactome_results_under_pw@termsim
# similarity_df <- as.data.frame(as.table(similarity_matrix))
# # Filter for significant pathway relationships (e.g., similarity > 0.5)
# significant_relationships <- subset(similarity_df, Freq > 0.5)
# print(significant_relationships)

# 11. Get the variance stabilised transformed expression values.
vsd <- vst(dds)

# 12. With the vst values obtain a PCA plot and a heatmap.
## 12.1 PCA Plot
pcaData <- plotPCA(vsd, intgroup = c(ERBB2_Amp_Col), returnData = TRUE)
## Creating a factor here allows for a binary colour legend
pcaData$ERBB2_Amp <- factor(pcaData$ERBB2_Amp, 
                            levels = c(0, 1), 
                            labels = c("ERBB2-", "ERBB2+"))

ggplot(pcaData, aes(PC1, PC2, color = ERBB2_Amp)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("darkblue", "red")) +
  labs(color = "ERBB2+/-") +
  ggtitle("PCA of DEG for ERBB2 Amplification, Variance-Stabilised Data") +
  xlab("Principal Component 1") +
  ylab("Principal Component 2")

## 12.2 Generate heatmap for top differentially expressed genes
ensure_package("pheatmap")
library(pheatmap)

# # Order genes by absolute log2 fold change
# top_DE <- order(abs(res$log2FoldChange), decreasing = TRUE)
# vsd_DE <- assay(vsd)[top_DE[1:20],]

# Order genes by adjusted p-value (offers better clustering...)
top_DE <- order(res$padj)
vsd_DE <- assay(vsd)[top_DE[1:20],]

## Again, a factor allows the legend to be binary, not a gradient
annotation_col <- data.frame(ERBB2_Amp = factor(as.matrix(metadata[, 1]), 
                                                levels = c(0, 1), 
                                                labels = c("ERBB2-", "ERBB2+")))

# Define the annotation colours for ERBB2_Amp to match the PCA chart above.
annotation_colors <- list(
  ERBB2_Amp = c("ERBB2-" = "darkblue", "ERBB2+" = "red")
)

rownames(annotation_col) <- colnames(vsd)

pheatmap(vsd_DE, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row",
         show_colnames = FALSE, 
         show_rownames = TRUE, 
         annotation_col = annotation_col,
         annotation_colors = annotation_colors)

# 13. With the vst values of the DE genes generate an overall survival model using the glmnet package.
## TODO lets see if this aligns with this paper: "Delving into the Heterogeneity of Different Breast 
## Cancer Subtypes and the Prognostic Models Utilizing scRNA-Seq and Bulk RNA-Seq"

# SOURCES
## https://glmnet.stanford.edu/articles/Coxnet.html
## https://bookdown.org/staedler_n/highdimstats/survival-analysis.html#regularized-cox-regression

library(glmnet)
library(survival)

##
## Preprocessing for X
##
# Extract the rownames (gene IDs) of significant genes
sig_genes <- rownames(res_sig)

# Extract the assay data (matrix of vst values)
vsd_matrix <- assay(vsd)

# Subset the vst matrix for only significant genes
vsd_matrix_sig <- vsd_matrix[sig_genes, ]

# Transpose the matrix to have patients as rows and genes as columns
X_vsd_sig <- t(vsd_matrix_sig)

# Standardise IDs in the vsd dataset
X_vsd_sig_rownames <- rownames(X_vsd_sig) # Extract row names (patient IDs) from vsd
X_vsd_sig_rownames <- sub("\\.\\d+$", "", X_vsd_sig_rownames) # Remove trailing ".01"
X_vsd_sig_rownames <- gsub("\\.", "-", X_vsd_sig_rownames) # Replace '.' with '-'

# Update column names in the vsd object
rownames(X_vsd_sig) <- X_vsd_sig_rownames

##
## Preprocessing for Y
##
# Set the PATIENT_ID column as the rownames
Y_patient <- data.frame(
  patient_id = data_patient$PATIENT_ID,
  time = data_patient$OS_MONTHS,
  status = data_patient$OS_STATUS
)
rownames(Y_patient) <- Y_patient$patient_id

# Drop the patient_id column
Y_patient$patient_id <- NULL
# Ensure the OS_MONTHS column is numeric and > 0
Y_patient$time <- as.numeric(Y_patient$time)
Y_patient <- Y_patient[Y_patient$time > 0,]
# Conver the OS_STATUS to numeric 1 or 0
Y_patient$status <- ifelse(Y_patient$status == "1:DECEASED", 1, 
                           ifelse(Y_patient$status == "0:LIVING", 0, NA))

##
## Match X and Y
##
# The patients in each dataset don't match perfectly;
# Find the common rownames between X_vsd_sig and Y_patient
common_rownames <- intersect(rownames(X_vsd_sig), rownames(Y_patient))

# Subset X_vsd_sig and Y_patient to only include rows with common rownames
X <- X_vsd_sig[common_rownames, , drop = FALSE]
Y_patient <- Y_patient[common_rownames, , drop = FALSE]
Y <- Surv(time = Y_patient$time, event = Y_patient$status)


##
## Split training and test data
##
# Set seed for reproducibility
set.seed(1234)

# Get the indices for the training set (80% of the data)
train_ind <- sample(1:nrow(X), size = floor(0.8 * nrow(X)))

# Split X into training and test sets
X_train <- X[train_ind, , drop = FALSE]
X_test <- X[-train_ind, , drop = FALSE]

# Split Y into training and test sets
Y_train <- Y[train_ind]
Y_test <- Y[-train_ind]

# Confirm the dimensions
cat("Training set: ", nrow(X_train), "patients\n")
cat("Test set: ", nrow(X_test), "patients\n")


# Its slow! Let's paralelise it.
library(doParallel)
cl <- makeCluster(workers)
registerDoParallel(cl)

fit <- glmnet(X_train, Y_train, family = "cox", parallel = TRUE)
plot(fit)
cvfit <- cv.glmnet(X_train, Y_train, family = "cox", type.measure = "C", parallel = TRUE)
plot(cvfit)

# Stop the parelel stuff.
stopCluster(cl)

# Extract coefficients at the optimal lambda as a matrix
model_genes <- coef(cvfit, s = "lambda.min")
non_zero_genes <- as.matrix(model_genes)
non_zero_genes <- non_zero_genes[non_zero_genes[, 1] != 0, , drop = FALSE]
# Sort the non-zero genes by coefficient values
sorted_genes <- non_zero_genes[order(non_zero_genes[, 1], decreasing = TRUE), , drop = FALSE]
# Extract the top 10 genes for high risk (positive coefficients)
top_high_risk_genes <- head(sorted_genes, 10)
# Extract the top 10 genes for low risk (negative coefficients)
top_low_risk_genes <- tail(sorted_genes, 10)
print(top_high_risk_genes)
print(top_low_risk_genes)


# Split patients into groups (e.g., high/low risk) based on scores
library(survminer)
# Generate risk scores for the test set
risk_scores <- predict(cvfit, newx = X_test, s = "lambda.min")

dat_test <- data.frame(
  time = Y_test[, "time"], 
  status = Y_test[, "status"], 
  risk_group = ifelse(risk_scores > median(risk_scores), "High Risk", "Low Risk")
)

# Fit a survival curve using the test data
s_fit <- survfit(Surv(time, status) ~ risk_group, data = dat_test)

# Plot the survival curve
ggsurvplot(
  s_fit, 
  data = dat_test, 
  conf.int = TRUE, 
  # risk.table = TRUE, 
  title = "Survival Analysis by ERBB2+ Differential Gene Expression", 
  xlab = "Time (Months)"
)

##
### ERBB2+/- Survival Model
##
# Transpose the matrix to have patients as rows and genes as columns
transpose_erbb2_cna <- t(erbb2_cna)

# Standardise IDs in the transpose_erbb2_cna dataset
transpose_erbb2_cna_rownames <- rownames(transpose_erbb2_cna) # Extract row names (patient IDs)
transpose_erbb2_cna_rownames <- sub("\\.\\d+$", "", transpose_erbb2_cna_rownames) # Remove trailing ".01"
transpose_erbb2_cna_rownames <- gsub("\\.", "-", transpose_erbb2_cna_rownames) # Replace '.' with '-'

# Update rownames in the transpose_erbb2_cna object
rownames(transpose_erbb2_cna) <- transpose_erbb2_cna_rownames

# Merge Y_patient and transpose_erbb2_cna by patient ID (rownames)
Y_erbb2_amp <- merge(Y_patient, transpose_erbb2_cna, by = "row.names", all = FALSE)

# Set the merged rownames to match the patient IDs
rownames(Y_erbb2_amp) <- Y_erbb2_amp$Row.names
Y_erbb2_amp$Row.names <- NULL

# Ensure time and status columns are numeric
Y_erbb2_amp$time <- as.numeric(Y_erbb2_amp$time)  # Ensure time is numeric
Y_erbb2_amp$status <- as.numeric(Y_erbb2_amp$status)  # Ensure status is numeric (1 = event, 0 = censored)

# Create the erbb2_amp column as a binary categorical variable
dat_test_erbb2 <- data.frame(
  time = Y_erbb2_amp[, "time"], 
  status = Y_erbb2_amp[, "status"], 
  erbb2_amp = ifelse(Y_erbb2_amp[, "ERBB2"] > 0, "ERBB2+", "ERBB2-")  # Create binary column
)

# Fit a survival curve using the test data
s_fit <- survfit(Surv(time, status) ~ erbb2_amp, data = dat_test_erbb2)

# Plot the survival curve
ggsurvplot(
  s_fit, 
  data = dat_test_erbb2, 
  conf.int = TRUE, 
  # risk.table = TRUE, 
  title = "Kaplan-Meier Survival by ERBB2 Amplification",
  xlab = "Time (Months)",
  ylab = "Survival Probability",
  legend.labs = c("ERBB2-", "ERBB2+"), # Customise legend labels
  palette = c("darkblue", "red")       # Customise colours
)



