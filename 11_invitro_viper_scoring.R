setwd("~/Documents/GitHub/mash-fa24")

# load in libraries
library(tidyverse)
library(RColorBrewer)
library(fgsea)
library(limma)
library(decoupleR)
library(viper)

# Signatures
gNAS <- c("ARPC5", "YWHAH", "ARF4", "TNFRSF12A", "ADHFE1",
          "USP33", "CD52", "ACVR2B", "ING5", "ASB3", "IFI30",
          "ERVW-1", "YWHAZ", "ERBB3", "KPNA2", "COQ10B", "MAGI1",
          "MAPRE1", "ABCA6")

gFib <- c("TIMP1", "MYLI2B", "LUM", "ZNF395", "AKAP9",
          "ACTR2", "LGALS3", "MAPRE1", "FRK", "ANKRD28",
          "IGFBP7", "YWHAZ", "USP33", "CD59", "TAX1BP3",
          "FAM221A", "ADHFE1", "TNFRSF12A")

X_MPS <- data.table::fread("datasets/MPS/X_MPS.csv") %>% column_to_rownames('V1')
metadata <- read.csv("datasets/metadata_MPS.csv")

sig_regulon <- rbind(
  data.frame(source = "NAS", target = gNAS, mor = c(1,1,1,1,-1,-1,1,-1,-1,-1,1,-1,1,-1,1,1,-1,1,-1)),
  data.frame(source = "Fibrosis stage", target = gFib, mor = c(1,1,1,-1,-1,1,1,1,-1,-1,1,1,-1,1,1,-1,-1,-1)))

sig_viper <- decoupleR::run_viper(t(X_MPS), net = sig_regulon) # ran for each of 3 datasets
saveRDS(sig_viper, "datasets/MPS/Y_viper_MPS.rds") # change for different datasets
