print_banner("DEG Analysis...")

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
## expression, whatever). What patterns can be see in their genome exoression 
## data from this grouping factor?"
dds <- DESeqDataSetFromMatrix(countData = assay,
                              colData = metadata,
                              design = ~ ERBB2_Amp)

## Run DESeq2 differential expression analysis
## Use 2 less cores than the total, but ensure at least 1 worker
workers <- max(1, parallel::detectCores() - 2)
register(MulticoreParam(workers = workers))
dds <- DESeq(dds, parallel = TRUE)
dds_res <- results(dds)

# 9. Obtain Differentially Expressed Genes
## Display top 10 differentially expressed genes by Fold Change.
## The genes listed will have either a positive or negative value for 
## `log2FoldChange` depending on whether they are up-regulated or 
## down-regulated in the presence of ERBB2+, respectively.
degs_top_10 <- dds_res[order(abs(dds_res$log2FoldChange), decreasing = TRUE)[1:10],]
print(degs_top_10)

degs_significant = dds_res[dds_res$padj < SIGNIFICANCE_LEVEL,]
deg_overx = rownames(degs_significant[degs_significant$log2FoldChange>0,])
deg_underx = rownames(degs_significant[degs_significant$log2FoldChange<0,])

# Count total number of genes
total_genes <- nrow(dds_res)

# Count number of significant DEGs
num_significant_degs <- nrow(degs_significant)

# Count over- and under-expressed DEGs
num_over_expressed <- length(rownames(degs_significant[degs_significant$log2FoldChange > 0, ]))
num_under_expressed <- length(rownames(degs_significant[degs_significant$log2FoldChange < 0, ]))

# Print the result
cat("Number of significant DEGs found:", num_significant_degs, 
    "out of a total of", total_genes, "genes (",
    round(100 * num_significant_degs / total_genes, 2), "%)\n")
cat("  - Over-expressed genes:", num_over_expressed, "\n")
cat("  - Under-expressed genes:", num_under_expressed, "\n")



