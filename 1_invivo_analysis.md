# Part 1: Running GSEA on in-vivo datasets to establish a baseline prediction.
Corresponding script: invivo_prediction.R

## GSEA
- GSEA, specifically the fgsea R package, requires passing in a list of ranked genes and pathways.
  - The datasets Nikos gave me were normalized count matrices, so they could not be passed into **DESeq2** (which expects raw count matrices) in order to calculate the logFC, stat, pvalue, etc statistics.
  - The **limma** package is an alternative for differential expression analysis which can handle log-transformed or normalized data. However, it was designed for data with continuous, normally distributed values (like log-transformed microarray intensities).
  - Then, the **voom** package converts the count data for RNA-seq into log-transformed values while accounting for the mean-variance relationship inherent in RNA-seq data.
    - Genes with low expression tend to have higher relative variance. Thus, voom applies precision weights to each observation, and then the transformed values and associated weights are fed into limma to allow it to perform linear modeling and differential expression analysis using the empirical Bayes method.
- [Perform GSEA starting with DESeq outputs](https://stephenturner.github.io/deseq-to-fgsea/)
- [Bioconductor's fgsea tutorial](https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html)
- [Bioconductor's DESeq2 tutorial](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [Biostatsquid's What is GSEA?](https://www.youtube.com/watch?v=egO7Lt92gDY)

- Is my data log-transformed?
  - Log transformations are defined only for positive values. Thus, since I have negative values, it was likely not log-transformed.

- limma differential expression using design (~ NAS + fibrosis).
- extract the t-statistic.
- fgsea using pathways (geneset_list containing gNAS and gFib), stat (t-statistic), nperm = 1000.

## Questions
- Which designs to choose from? NAS + fibrosis, NAS, fibrosis
- Is the data log-transformed normalized?
- Should I only use the t-statistic for fgsea analysis? Why or why not select any other stat?
- How many permutations should I run for gsea? What does nperm mean for gsea?
- What is the value under each coefficient when I run topTable(fit)?