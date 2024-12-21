print_banner("Pathway Enrichment...")

# 10. Perform a Pathway Enrichment Analysis
## So from our above DESeq analysis, we have found a bunch of genes that are
## either over or under expressed. Now with enrichment analysis we can determine
## what biological processes are invovled in these genes and how statistically
## significant the results are. For example, ERBB2 amplification might cause a
## bunch of other genes related to the cell cycle to become overexpressed. This 
## enrichment analysis would say "we're seeing a lot of cell cycle processes
## being impacted in patients with ERBB2 amplification; more than would be
## statistically likely by chance."
library(ggplot2)
ensure_package("clusterProfiler")
library(clusterProfiler)
ensure_package("org.Hs.eg.db")
library(org.Hs.eg.db)
ensure_package("enrichplot")
library(enrichplot)

print("Gene Ontology...")

## Gene Ontology (GO) enrichment analysis
### Find genes that are over-expressed in ERBB2 amplified patients.
go_results_over = enrichGO(
  gene          = deg_overx,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP", # Biological Process
  qvalueCutoff  = SIGNIFICANCE_LEVEL # Use a stricter threshold than 0.2 to 
  # provide more confidence in our results.
)

### Find genes that are under-expressed in ERBB2 amplified patients.
go_results_under = enrichGO(
  gene          = deg_underx,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",  
  ont           = "BP",
  qvalueCutoff  = SIGNIFICANCE_LEVEL
)


### Plot results of GO enrichment analysis
print(dotplot(go_results_over) + 
        ggtitle("Gene Ontology Enrichment Over Expressed"))

print(dotplot(go_results_under) + 
        ggtitle("Gene Ontology Enrichment Under Expressed"))

print(treeplot(pairwise_termsim(go_results_over)) + 
        ggtitle("GO Enrichment Over Expressed"))

print(treeplot(pairwise_termsim(go_results_under)) + 
        ggtitle("GO Enrichment Under Expressed"))

print("KEGG...")

## KEGG Pathway Enrichment Analysis
### Map gene symbol to gene entrez ID. 
### This is necessary for Reactome and Keggs.
deg_overx_entrez <- bitr(
  deg_overx,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

deg_underx_entrez <- bitr(
  deg_underx,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

### defaults to 'hsa' for organism i.e. home sapien
kegg_results_over =  enrichKEGG(
  gene          = deg_overx_entrez[,2],
  qvalueCutoff  = SIGNIFICANCE_LEVEL
)

kegg_results_under =  enrichKEGG(
  gene          = deg_underx_entrez[,2],
  qvalueCutoff  = SIGNIFICANCE_LEVEL
)

### Plot results of KEGG pathway enrichment analysis
print(dotplot(kegg_results_over) + 
        ggtitle("Kegg Pathway Enrichment Over Expressed"))

print(dotplot(kegg_results_under) + 
        ggtitle("Kegg Pathway Enrichment Under Expressed"))

print(treeplot(pairwise_termsim(kegg_results_over)) + 
        ggtitle("KEGG Enrichment Over Expressed"))

print(treeplot(pairwise_termsim(kegg_results_under)) + 
        ggtitle("KEGG Enrichment Under Expressed"))

print("Reactome...")

## Reactome Pathway Enrichment Analysis
ensure_package("ReactomePA")
library(ReactomePA)
ensure_package("pathview")
library(pathview)

### defaults to 'human' for organism
reactome_results_over =  enrichPathway(
  gene          = deg_overx_entrez[,2],
  qvalueCutoff  = SIGNIFICANCE_LEVEL,
)

reactome_results_under =  enrichPathway(
  gene          = deg_underx_entrez[,2],
  qvalueCutoff  = SIGNIFICANCE_LEVEL,
)

### Plot results of Reactome pathway enrichment analysis
print(dotplot(reactome_results_over, showCategory = 10) + 
        ggtitle("Reactome Pathway Enrichment Over Expressed"))

print(dotplot(reactome_results_under, showCategory = 10) + 
        ggtitle("Reactome Pathway Enrichment Under Expressed"))

print(treeplot(pairwise_termsim(reactome_results_over)) + 
        ggtitle("Reactome Enrichment Over Expressed"))

print(treeplot(pairwise_termsim(reactome_results_under)) + 
        ggtitle("Reactome Enrichment Under Expressed"))