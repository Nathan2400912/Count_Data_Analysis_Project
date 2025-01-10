# RNA-seq Results Analyzer
This repository contains the final assignment project for the **R Bioinformatics course (BF591)** at Boston University. The goal of this project is to create an interactive **RShiny application** for manipulating and visualizing RNA-seq experiment results.

## Background
This application is based on an RNA-seq study of Huntington's Disease (HD) detailed in the paper:
“RNA Sequence Analysis of Human Huntington Disease Brain Reveals an Extensive Increase in Inflammatory and Developmental Gene Expression.”
Huntington's Disease is a severe neurodegenerative disorder caused by an expanded CAG trinucleotide repeat in the Huntingtin (HTT) gene. While transcriptional dysregulation in the HD brain is well-documented, it remains incompletely understood. This study presents a genome-wide analysis of mRNA expression in the human prefrontal cortex from 20 HD patients and 49 neuropathologically normal controls, using next-generation sequencing.

For an overview of the application, you can watch the following introductory video:


https://github.com/user-attachments/assets/bf4b9d7c-bb35-4271-ab6a-e321bb090beb



If you need further clarification, continue reading this README for detailed guidance.

## Data Preparation
The application is structured into **four sections**, each requiring specific input files in CSV format. Detailed instructions on preparing these files are provided in **file_preparation.Rmd**. The input data for this app was obtained from **GSE64810**.

## Application features
The app has four main sections, each designed for a specific purpose:

### 1. Samples 
This section summarizes the metadata associated with the study. For this project, you can use the file study_metadata.csv. The Samples tab has three subsections:
- Summary: Displays the mean ± SE for quantitative metadata and lists distinct values for qualitative metadata.
- Table: Provides a sortable table of metadata columns. You can use selectors to customize which columns to display.
- Plots: Generates overlapping density plots for quantitative metadata, allowing users to select the variables to visualize.

### 2. Counts
This section visualizes the count data matrix. The input file should ideally be a normalized count matrix, as no additional normalization is performed. You can use the file normalized_counts.csv provided with the study. Two filter sliders are available:
- Variance Slider: Filters genes based on variance quartile thresholds.
- Non-Zero Slider: Filters genes with more than X zero-value samples.

Filters are applied by clicking the "Apply Filters" button. The Counts tab has four subsections:
- Summary: Shows the number of samples, genes, and the percentage of genes passing the current filters.
- Scatter Plot: Visualizes variance vs. median count or zero-value samples vs. median count.
- Heatmap: Displays a heatmap of the normalized count matrix using pheatmap.
- PCA: Performs PCA on the dataset and lets users select PCs to plot. Metadata and count data must be uploaded for PCA computation.

### 3. Differential Expression (DE)
This section visualizes differential expression results (e.g., from DESeq2). For this study, the DE results are saved in the file deseq2_results.csv. The DE tab includes:
- Table: A sortable table of DE results.
- Plots: Customizable plots using columns from the DE results file. Includes a p-value threshold slider to highlight points of interest. For instance, you can create a volcano plot with log2FoldChange on the x-axis and padj (log-transformed) on the y-axis.

### 4. Gene Set Enrichment Analysis (GSEA)
This section displays GSEA results but does not compute them. For this project, GSEA results were generated using fgsea with log2FoldChange values from DESeq2 results against the hallmark gene set from MSigDB. More details are available in file_preparation.Rmd. GSEA features:
- Top Pathways: Shows the top 10 positive and top 10 negative pathways based on the p-value threshold.
- Results: Displays a sortable table of GSEA results. Users can filter pathways (all, positive, or negative) and download the filtered results.
- Scatter Plot: Plots padj values against NES scores from the GSEA results.

## Summary
This README provides a complete walkthrough of the RNA-seq Results Analyzer app. For further clarification on data preparation or functionality, refer to **file_preparation.Rmd** or explore the app to better understand its capabilities.

Happy analyzing! 🎉



