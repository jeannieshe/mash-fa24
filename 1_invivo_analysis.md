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

## Changes to code 9/11
- The data are normalized (CPM, log2 transformed, and then centered)
- Nikos: you can use a code like above! where measurements is basically you X. And this will give you 2 matrices. One with rows the genesets and columns the sample that has Normalized Enrichment Scores (NES) and one with the same dimesnsions and p-values (for the later we dont care at this point so you can avoid using it if you want and for that reason you may use just 100 permutations)

- GSEA actually doesn't account for how the gene affects itself within the geneset. It doesn't account for the correlation of each gene with NAS and fibrosis, but it matters in combination which of them are upregulated and which are down. The way the markers are annotated, it matters if they are positively or negatively correlated with NAS and fibrosis

- Instead, let's run viper (an activity inference technique similar to GSEA). It requires a package called [decoupleR](https://saezlab.github.io/decoupleR/)

- If I want to use it in python, I will need to create a conda environment and install it there.
  - conda create -n myenv; conda activate myenv; conda install decoupler-py -c conda-forge
  
## Updates 10/1
- Instead of using GSEA, I will now be using VIPER. VIPER uses the regulon of the datasets and the given gene sets in order to predict a score for each patient.
- VIPER gives a direction of regulation in addition to the magnitude of regulation, and the score is in the range (-2, 2).

- Move datasets to Python to train on PLSR
- tuning:
- train on one dataset, Govaere, 10 fold cross validation (90% train, 10% test)
PLSR: try diff. num_lvs = [2, 4, 6, 8, 10] to find the latent variable hyperparameter that provides the highest performance
  - save the models and take the best one and make prediction on test set
  - save the actual data splits
  - or you can train with all of Govaere and test on the other two
  - as csv, after selecting the model:
  - train, validation, Hoang, Pantano, shuffled on x
  - y: pearson coefficient for fibrosis
  - shuffle: shuffle the genes but with the same names, 
  
- test sets on Hoang and Pantano
- then repeat tuning with Hoang and Pantano

- right now, focus on predicting clinical scores
later, after the model is established, see how well it will predict inferred scores

- for the data splits, try to match the distribution of the data points (NAS/ fibrosis scores)

- Let me perform an exercise on identifying the most optimal PLSR models.
  - First, let's perform a tuning/optimization algorithm.
  - Let's take our tuning dataset = Govaere.
    - Create a 10-fold cross validation dataset split where the train and test parts of the 


- Now that I have identified the most optimal PLSR models, I can work on PLSR, SVM, ElasticNet, KNN, random forest (how long)

## Updates 10/8
- Use sk.learn's StratifiedKFold function to split evenly on the NAS/Fibrosis scores as classes, with each iteration having mutually exclusive data points


## Updates 10/29
- created scatter plots to see how well viper scores predict the NAS and Fibrosis clinical scores
- a score of 0 means uncertainty, which isn't really well depicted here
- the Pearson correlation coefficients are not very high, but there is an obvious relationship there still

## Updates 11/5
- finished with the boxplots for PLSR training! number of LVs: 8 was best for all of the models
- moving onto KNN and ElasticNet today
- KNN:
  - choosing KNN Regression because we have a relationship between the Fib/NAS score, whereas classification would imply that the classes were independent
  - number of neighbors: 
    - Pantano: From 10-fold cross validation, 20 neighbors achieves an average Pearson coeff of 0.5833896804381291.
    - Govaere: From 10-fold cross validation, 5 neighbors achieves an average Pearson coeff of 0.5372950994900263.
    - Hoang: From 10-fold cross validation on Hoang data, 10 neighbors achieves an average Pearson coeff of 0.6693899907884641.
- ElasticNet:
  - choosing ElasticNet to try and combine Lasso and Ridge Regression tactics
  - ran into warnings when using GridSearchCV, specifically about ConvergenceWarnings
    - next time, try a max_iter=2000 and also tol=0.01
    - if it performs fine on the training data, it is okay:)
    - doesn't necessarily use gradient descnet to calculate, therefore we know that running it in a GPU is not necessarily any better

## Updates 11/12
- Finished ElasticNet training
- alpha and L1 ratio have different meanings
- cannot necessarily just simplify it to mean Lasso vs. Ridge regression
- Doing LinearSVR: it's really bad. 
  - Needed to use MultiOutputRegressor because of the two variables to predict
  - Mean Pearson correlation (Fibrosis): 0.3809897076428476
  - Mean Pearson correlation (NAS): 0.5101705543138401

## Updates 11/19
- Reconstruct the datatable and save as csv
- Don't trust my saved models because validation is higher than training?
- Only run RBF SVM because linear and poly take too long
- Figure out why the validation dataset is so different???

## Updates 11/26
- Finally reconstructed all of the Pearson correlation coefficients for the PLSR tuning
  - Pantano: LV 9
  - Govaere: LV 8
  - Hoang: LV 10
- Reconstructed coefficients for KNN and ElasticNet as well
- Working on RBF SVM
- Actually discovered what was happening with my saved models (incorrect data splitting
when using skf from sklearn)
- Satisfyingly debugged my own code! Found a few errors, namely labeling NAS scores as
Fib and vice versa
- Reconstructed the csv dataset for all of these coefficients!
- Choosing Govaere dataset from PLSR models because it has the most diverse dataset and
rather good performance

## Updates 12/3
- Created the ideal model using PLSR, Govaere, all training data
- Finished predicting NAS/Fib scores on in vitro data and also calculated VIPER scores
  - Plotted -> interesting positive correlation amongst NAS and Fib scores
