# RNA-seq Results Analyzer
This repository contains a project for the final assignment of the R Bioinformatics course (BF591) at Boston University. The goal of this project is to create an RShiny application for the manipulation and visualization of results from an RNA-seq experiment.

## Background
This app is made based on a study on Huntingtons Disease that contains RNA-seq results. The paper used is as follows: RNA Sequence Analysis of Human Huntington Disease Brain Reveals an Extensive Increase in Inflammatory and Developmental Gene Expression. Huntington’s Disease (HD) is a devastating neurodegenerative disorder that is caused by an expanded CAG trinucleotide repeat in the Huntingtin (HTT) gene. Transcriptional dysregulation in the human HD brain has been documented but is incompletely understood. Here we present a genome-wide analysis of mRNA expression in human prefrontal cortex from 20 HD and 49 neuropathologically normal controls using next generation high-throughput sequencing. 

## Data Preparation
This app contains four different tabs with different functions each requiring different csv files. To see how the data was obtained and how the different CSV files were created, please refer to **file_preparation.Rmd**. In general, files were obtained from GSE64810.

## Tabs
This application has 4 sections serving different purposes. Below is a brief walkthrough of each section's use and functionality. 

### Samples 
This section is used to give a summary view of the metadata for the respective study. In this study specifically, the metadata file you can use is called **study_metadata.csv**. There are three specific tabs that you can use to display the metadata, here is a brief description of each tabs function.
  
  - Summary: Provides the mean +- SE of quantitative metadata. Also includes qualitative metadata but displays distinct values instead
  - Table: Sortable table of the columns of the metadata csv file. At the bottom are selectors where you can decide what to include in the table
  - Plots: Overlapping density plots of quantitative metadata, you can also choose what to plot

### Counts
The counts section is where you can visualize your count data by uploading a count data matrix. The uploaded is preferably a normalized matrix as no normalization steps are included in the background. You can try using the respective count matrix of this study in the file **normalized_counts.csv**. This section contains two sliders that you can use to control the output of the sliders. For one, there is a variance slider, which filters out genes based on the threshold variance quartile. Second, there is non-zero slider that filters out genes with greater than X amount of zero samples. To apply the filters, simply click apply filters under the sliders. Here is a breakdown of the function of the four tabs in this section:

  - Summary: Shows the number of samples, as well as the number and percent of genes passing the current filters applied
  - Scatter Plot: Scatterplots of either variance vs median count or number of zero samples vs median count. Basically visualization of the filters and the summary table in the previous tab
  - Heatmap: Heatmap of the normalized count matrix generated using pheatmap 
  - PCA: Allows you to plot whichever PCs you want to. PCA is computed in the background once the metadata and count data is uploaded. Once computed, you then can choose what PC to plot based on the dropdown of the x and y-axis options.

### DE (Differential Expression)
The DE section allows you to visualize any differential expression results that you may have. In this case, the differential expression results were generated using DESEq2 (preferred DE method for this program) and the DE results were downloaded into a csv file called **deseq2_results.csv**. There are two tabs in this DE section:

  - Table: A sortable table of the differential expression results
  - Plots: Allows you to plot anything that you want based on the columns of your differential expression reuslts file. There is also a p-value threshold slider that allows to color out points based on the p-value threshold of your interest. There are color selectors as well. For example, you can plot log2FoldChange in the x-axis and padj in the y-axis to generate a volcano plot. Note that the y-axis is log-transformed to aid in visualization.

## GSEA
The GSEA section is used to display GSEA results, it does not compute GSEA results for you. To generate the GSEA results for this study, I ran fgsea using log2FoldChange as the parameter from the DESEq2 results against hallmark gene set from MSigDB. For a more comprehensive view of how I obtained the results, please once again refer to **file_preparation.Rmd**. In this section, there is a p-value threshold slider for you to choose what results to exclude, the lower the slider the more it excludes. There are 3 tabs in this section, here is what they do:

  - Top Pathways: Displays the top positive 10 and negative 10 pathways. If there are less than 10 based on the exclusion criteria or if your results have less than 10, it will simply show the remaining ones.
  - Results: Sortable table of your GSEA results. You can select all the pathways or only positive/negative pathways using the radiobuttons under the slider. A download button is also included for you to donwload desired filtered results into a csv file.
  - Scatter Plot: Visualization of padj values against NES scores from the gsea results.

That is the complete rundown of the app.



