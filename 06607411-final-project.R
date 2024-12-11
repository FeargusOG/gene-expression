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

## Filter out genes with low expression
## TODO Why these numbers?? Justify....
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
## Display top 10 differentially expressed genes by adjusted p-value
top_genes <- res[order(res$padj)[1:10], ]
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

## TODO Why this number?? p-value?
res_sig = res[res$padj<0.05,]
DE_over = rownames(res_sig[res_sig$log2FoldChange>0,])
DE_under = rownames(res_sig[res_sig$log2FoldChange<0,])

## Gene Ontology (GO) enrichment analysis
### Find genes that are over-expressed in ERBB2 amplified patients.
go_results_over = enrichGO(
  gene          = DE_over,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP", # TODO why not use default MF
  pAdjustMethod = "BH", # TODO This is the default... can remove.
  pvalueCutoff  = 0.05, # TODO This is the default... can remove.
  qvalueCutoff  = 0.05 # TODO Why not use default 0.2?
)
print(head(go_results_over))

### Find genes that are under-expressed in ERBB2 amplified patients.
go_results_under = enrichGO(
  gene          = DE_under,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",  
  ont           = "BP", # TODO why not use default MF
  pAdjustMethod = "BH", # TODO This is the default... can remove.
  pvalueCutoff  = 0.05, # TODO This is the default... can remove.
  qvalueCutoff  = 0.05 # TODO Why not use default 0.2?
)

print(head(go_results_under))

### Plot results of GO enrichment analysis
dotplot(go_results_over, showCategory=10) + 
  ggtitle("Gene Ontology Enrichment Over Expressed")

dotplot(go_results_under, showCategory=10) + 
  ggtitle("Gene Ontology Enrichment Under Expressed")

## KEGG Pathway Enrichment Analysis
### Map symbol to entrez for Reactome and Keggs
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

kegg_results_over =  enrichKEGG(
  gene          = gene_entrez_over[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

kegg_results_under =  enrichKEGG(
  gene          = gene_entrez_under[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

print(head(kegg_results_over))
print(head(kegg_results_under))

### Plot results of KEGG pathway enrichment analysis
dotplot(kegg_results_over, showCategory=10) + 
  ggtitle("Kegg Pathway Enrichment Over Expressed")
dotplot(kegg_results_under, showCategory=10) + 
  ggtitle("Kegg Pathway Enrichment Under Expressed")

## Reactome Pathway Enrichment Analysis
ensure_package("ReactomePA")
library(ReactomePA)
ensure_package("pathview")
library(pathview)

reactome_results_over =  enrichPathway(
  gene          = gene_entrez_over[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
)

reactome_results_under =  enrichPathway(
  gene          = gene_entrez_under[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
)

print(head(reactome_results_over))
print(head(reactome_results_under))

### Plot results of Reactome pathway enrichment analysis
dotplot(reactome_results_over, showCategory=10) + 
  ggtitle("Reactome Pathway Enrichment Over Expressed")
dotplot(reactome_results_under, showCategory=10) + 
  ggtitle("Reactome Pathway Enrichment Under Expressed")

go_results_under_pw = pairwise_termsim(go_results_under)
treeplot(go_results_under_pw)+ ggtitle("GO Enrichment Under Expressed")

kegg_results_under_pw = pairwise_termsim(kegg_results_under)
treeplot(kegg_results_under_pw)+ ggtitle("KEGG Enrichment Under Expressed")

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
