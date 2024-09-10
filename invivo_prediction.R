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
X_Govaere <- readRDS("X_Govaere.rds") # nolint
X_Hoang <- readRDS("X_Hoang.rds")
X_Pantano <- readRDS("X_Pantano.rds")

Y_Govaere <- readRDS("Y_Govaere.rds")
Y_Hoang <- readRDS("Y_Hoang.rds")
Y_Pantano <- readRDS("Y_Pantano.rds")

# initialize signatures
gNAS <- c("ARPC5", "YWHAH", "ARF4", "TNFRSF12A", "ADHFE1",
          "USP33", "CD52", "ACVR2B", "ING5", "ASB3", "IFI30",
          "ERVW-1", "YWHAZ", "ERBB3", "KPNA2", "COQ10B", "MAGI1",
          "MAPRE1", "ABCA6")

gFib <- c("TIMP1", "MYLI2B", "LUM", "ZNF395", "AKAP9",
          "ACTR2", "LGALS3", "MAPRE1", "FRK", "ANKRD28",
          "IGFBP7", "YWHAZ", "USP33", "CD59", "TAX1BP3",
          "FAM221A", "ADHFE1", "TNFRSF12A")

# begin by performing limma for differential expression analysis on the datasets
design <- model.matrix(~ NAS + fibrosis, data = Y_Govaere)
fit <- lmFit(X_Govaere, design)
fit <- eBayes(fit) # robust=TRUE, trend=TRUE if log-transformed normalized data
# then run fgsea() with the pathways

# need to set directory to github!