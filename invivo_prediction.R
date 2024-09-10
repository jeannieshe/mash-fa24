# fa24-urop part one
# using in vivo data, see if gsea will be a good predictor of clinical severity
# later using ML

# setwd("~/Desktop/fa24-urop")

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
datasets <- c("Govaere", "Hoang", "Pantano")
fgsea_NASt <- data.frame(pathway = character(), NES = numeric(), dataset = character(), stringsAsFactors = FALSE) #tstat from fgsea is NAS coefficient
fgsea_Fibt <- data.frame(pathway = character(), NES = numeric(), dataset = character(), stringsAsFactors = FALSE) #tstat from fgsea is Fibrosis coefficient

for (ii in 1:length(datasets)) {
    name <- datasets[ii]
    X <- get(sprintf("X_%s", name))
    Y <- get(sprintf("Y_%s", name))

    X <- t(X)
    Y <- as.data.frame(Y)
    
    # begin by performing limma for differential expression analysis on the datasets
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
    tstat_fib <- diffex_fib$t
    names(tstat_fib) <- rownames(diffex_fib)
    fgsea_fib <- fgsea(
      pathways = gene_sets, 
      stats = tstat_fib, 
      nperm = 1000) #error: there are ties in the preranked stats (1.96%)

    # store the NES values in a dataframe
    temp_NASdf <- fgsea_NAS[, c("pathway", "NES")]
    temp_NASdf$dataset <- name
    fgsea_NASt <- rbind(fgsea_NASt, temp_NASdf)
    
    temp_Fibdf <- fgsea_fib[, c("pathway", "NES")]
    temp_Fibdf$dataset <- name
    fgsea_Fibt <- rbind(fgsea_Fibt, temp_Fibdf)
}

saveRDS(fgsea_NASt, "datasets/fgsea_NASt.rds")
saveRDS(fgsea_Fibt, "datasets/fgsea_Fibt.rds")


