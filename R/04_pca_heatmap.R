print_banner("PCA and HeatMap...")
# 11. Get the variance stabilised transformed expression values.
vsd <- vst(dds)

# 12. With the vst values obtain a PCA plot and a heatmap.
## 12.1 PCA Plot
pcaData <- plotPCA(vsd, intgroup = c(ERBB2_Amp_Col), returnData = TRUE)

## Creating a factor here allows for a binary colour legend
pcaData$ERBB2_Amp <- factor(pcaData$ERBB2_Amp, 
                            levels = ERBB2_Factor_Levels, 
                            labels = ERBB2_Factor_Labels)

print(ggplot(pcaData, aes(PC1, PC2, color = ERBB2_Amp)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(ERBB2_Negative_Colour, ERBB2_Positive_Colour)) +
  labs(color = "ERBB2+/-") +
  ggtitle("PCA of DEG for ERBB2 Amplification, Variance-Stabilised Data") +
  xlab("Principal Component 1") +
  ylab("Principal Component 2"))

## 12.2 Generate heatmap for top differentially expressed genes
ensure_package("pheatmap")
library(pheatmap)

# # Order genes by absolute log2 fold change
# degs_fold_ordered <- order(abs(dds_res$log2FoldChange), decreasing = TRUE)
# degs_vsd_fold_ordered <- assay(vsd)[degs_fold_ordered[1:20],]

# Order genes by adjusted p-value (offers better clustering...)
degs_p_ordered <- order(dds_res$padj)
degs_vsd_p_ordered <- assay(vsd)[degs_p_ordered[1:20],]

## Again, a factor allows the legend to be binary, not a gradient
annotation_col <- data.frame(ERBB2_Amp = factor(as.matrix(metadata[, 1]), 
                                                levels = ERBB2_Factor_Levels, 
                                                labels = ERBB2_Factor_Labels))
rownames(annotation_col) <- colnames(vsd)

# Define the annotation colours for ERBB2_Amp to match the PCA chart above.
annotation_colors <- list(
  ERBB2_Amp = c("ERBB2-" = ERBB2_Negative_Colour, 
                "ERBB2+" = ERBB2_Positive_Colour)
)

pheatmap(degs_vsd_p_ordered, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row",
         show_colnames = FALSE, 
         show_rownames = TRUE, 
         annotation_col = annotation_col,
         annotation_colors = list(ERBB2_Amp = c(
           "ERBB2-" = ERBB2_Negative_Colour, 
           "ERBB2+" = ERBB2_Positive_Colour))
         )
