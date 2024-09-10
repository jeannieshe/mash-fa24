# fa24-urop part one
# using in vivo data, see if gsea will be a good predictor of clinical severity
# later using ML

setwd("~/Desktop/fa24-urop")

# load in libraries
library(tidyverse)
library(RColorBrewer)
library(fgsea)
library(limma)

# load in data (count normalized)
X_Govaere <- readRDS("datasets/X_Govaere.rds")
X_Hoang <- readRDS("datasets/X_Hoang.rds")
X_Pantano <- readRDS("datasets/X_Pantano.rds")

Y_Govaere <- readRDS("datasets/Y_Govaere.rds")
Y_Hoang <- readRDS("datasets/Y_Hoang.rds")
Y_Pantano <- readRDS("datasets/Y_Pantano.rds")

# initialize signatures
gNAS <- c("ARPC5", "YWHAH", "ARF4", "TNFRSF12A", "ADHFE1",
          "USP33", "CD52", "ACVR2B", "ING5", "ASB3", "IFI30",
          "ERVW-1", "YWHAZ", "ERBB3", "KPNA2", "COQ10B", "MAGI1",
          "MAPRE1", "ABCA6")

gFib <- c("TIMP1", "MYLI2B", "LUM", "ZNF395", "AKAP9",
          "ACTR2", "LGALS3", "MAPRE1", "FRK", "ANKRD28",
          "IGFBP7", "YWHAZ", "USP33", "CD59", "TAX1BP3",
          "FAM221A", "ADHFE1", "TNFRSF12A")
gene_sets <- list(
    "gNAS" = gNAS,
    "gFib" = gFib
)
# begin by performing limma for differential expression analysis on the datasets
X_data <- c("X_Govaere", "X_Hoang", "X_Pantano")
Y_data <- c("Y_Govaere", "Y_Hoang", "Y_Pantano")

for (ii in 1:length(X_data)) {
    X_name <- X_data[ii]
    Y_name <- Y_data[ii]
    X <- get(X_name)
    Y <- get(Y_name)

    X <- t(X)
    Y <- as.data.frame(Y)
    design <- model.matrix(~ NAS + fibrosis, data = Y)
    fit <- lmFit(X, design) # linear modelling
    fit <- eBayes(fit) # empirical Bayes moderation carried out to obtain precise estimates of gene-wise variability
    # for eBayes(fit), set robust=TRUE, trend=TRUE if log-transformed normalized data

    # summary(decideTests(fit))
    diffex_NAS <- topTable(fit, coef = "NAS", n = Inf)
    diffex_fib <- topTable(fit, coef = "fibrosis", n = Inf)

    # then run fgsea() with the pathways
    # should i just use the t-statistic? why or why not select any other stat?
    tstat_NAS <- diffex_NAS$t
    names(tstat_NAS) <- rownames(diffex_NAS)

    fgsea_NAS <- fgsea(
        pathways = gene_sets, 
        stats = tstat_NAS, 
        nperm = 1000) #error: there are ties in the preranked stats (1.96%)

    fgsea_Govaere <- fgsea_NAS[, c("pathway", "NES")]
    colnames(fgsea_Govaere) <- c("Pathway", "NES", "DiffEx")
    fgsea_Govaere$DiffEx <- c("NAS", "NAS")
}



