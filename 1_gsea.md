# Part 1: Running GSEA on in-vivo datasets to establish a baseline prediction.
Corresponding script: invivo_prediction.R

## GSEA
- GSEA, specifically the fgsea R package, requires passing in a list of ranked genes and pathways.
  - The datasets Nikos gave me were normalized count matrices, so they could not be passed into **DESeq2** (which expects raw count matrices) in order to calculate the logFC, stat, pvalue, etc statistics.
  - The **limma** package is an alternative for differential expression analysis which can handle log-transformed or normalized data. However, it was designed for data with continuous, normally distributed values (like log-transformed microarray intensities).
  - Then, the **voom** package converts the count data for RNA-seq into log-transformed values while accounting for the mean-variance relationship inherent in RNA-seq data.
    - Genes with low expression tend to have higher relative variance. Thus, voom applies precision weights to each observation, and then the transformed values and associated weights are fed into limma to allow it to perform linear modeling and differential expression analysis using the empirical Bayes method.