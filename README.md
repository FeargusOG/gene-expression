# gene-expression

This repository contains R scripts for analysing breast cancer RNA-seq and CNA data, performing differential expression analysis, pathway enrichment, and survival modelling.

The code contains significant amounts of comments for greater detail, but here is an overview of the files!

## File Descriptions

### 00_utilities.R

This file contains utility functions used throughout the project:

- **ensure_package(pkg)**: Ensures the specified R package is installed and loaded.
- **print_banner(section_name)**: Prints a formatted banner to console, useful for workflow progress updates.

---

### 01_preprocessing.R

Handles preprocessing of input datasets:

- **read_data(patient, rnaseq, cna)**: Reads and processes the patient, RNA-seq, and CNA data files.
- **extract_cna_gene(data_cna, gene_symbol)**: Extract just a single gene from the CNA dataset.
- **create_amp_metadata(data_rnaseq, data_cna, gene, amp_thresh)**: Creates metadata for amplification status of a specified gene (e.g., ERBB2) based on a threshold. Returns assay data and amplification metadata.

---

### 02_DEG.R

Performs differential expression analysis:

- **run_deg_analysis(assay, metadata, group_col)**: Executes DESeq2-based differential expression analysis.
- **split_degs(dds_res)**: Categorises significant DEGs into over-expressed and under-expressed groups based on log2FoldChange.

---

### 03_pathway_enrichment.R

Conducts pathway enrichment analysis:

- **run_pathway_enrichment(deg_gene_symbols, gene_list_descr, p_thresh, q_thresh)**: Performs enrichment analysis for GO, KEGG, and Reactome pathways using specified gene lists and thresholds.
- **plot_pathway_enrichment_results(results, title)**: Visualises enrichment results via dotplots and treeplots, ensuring there are enough data points for clustering.

---

### 04_pca_heatmap.R

Generates PCA plots and heatmaps:

- Uses variance-stabilised transformed (VST) data to create:
  - **PCA plots**: Visualises variance across principal components, with colour coding based on amplification status.
  - **Heatmaps**: Highlights clustering of differentially expressed genes.

---

### 05_cox_regress.R

Implements survival modelling:

- **prepare_cox_reg_x_y(vsd, patient_data, gene_list)**: Prepares data for Cox regression analysis by aligning patient metadata and expression data.
- **split_training_data(X, Y, seed, split)**: Splits the dataset into training and test sets.
- **fit_cox_models(X_train, Y_train, plot_it)**: Fits Cox regression models using glmnet, with optional plotting.
- **make_cox_predictions(cox_cv_model, X_test, lambda)**: Generates risk scores for test data.
- **fit_survival_curves(Y_time, Y_status, group_by)**: Creates Kaplan-Meier survival curves for risk stratification.
- **extract_cox_model_genes(cox_cv_model, lambda)**: Extracts genes contributing to the Cox model's predictions.

---

### main.R

Coordinates the entire workflow
