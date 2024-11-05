# creating a scatter plot to determine how well viper correlates with the known 
# values, determined quantitatively by the pearson correlation coefficient

library(tidyverse)
library(RColorBrewer)
library(fgsea)
library(limma)
library(decoupleR)
library(viper)
library(ggplot2)
library(ggpubr)

setwd("~/Documents/GitHub/mash-fa24")

# load in data (count normalized)
X_Govaere <- readRDS("datasets/X_Govaere.rds")
X_Hoang <- readRDS("datasets/X_Hoang.rds")
X_Pantano <- readRDS("datasets/X_Pantano.rds")

Y_Govaere <- readRDS("datasets/Y_Govaere.rds")
Y_Hoang <- readRDS("datasets/Y_Hoang.rds")
Y_Pantano <- readRDS("datasets/Y_Pantano.rds")

viper_X_Govaere <- readRDS("datasets/viper_X_Govaere.rds")
viper_X_Hoang <- readRDS("datasets/viper_X_Hoang.rds")
viper_X_Pantano <- readRDS("datasets/viper_X_Pantano.rds")

# Govaere
viper_Govaere_Fib <- viper_X_Govaere %>%
  filter(source=="Fibrosis stage")
viper_Govaere_NAS <- viper_X_Govaere %>%
  filter(source=="NAS")

Govaere_clinical_viper_Fib <- data.frame(clinical=as.data.frame(Y_Govaere)$fibrosis, viper=viper_Govaere_Fib$score)
Govaere_clinical_viper_NAS <- data.frame(clinical=as.data.frame(Y_Govaere)$NAS, viper=viper_Govaere_NAS$score)

Govaere_Fib_scores <- ggscatter(data=Govaere_clinical_viper_Fib, 
                                x="clinical", y="viper", size=2, 
                                cor.coef=TRUE, cor.method="pearson", 
                                title="Govaere Fibrosis Score Correlation"
          )
Govaere_NAS_scores <- ggscatter(data=Govaere_clinical_viper_NAS, 
                            x="clinical", y="viper", size=2, 
                            cor.coef=TRUE, cor.method="pearson", 
                            title="Govaere NAS Score Correlation"
          )

ggsave("datasets/Govaere_NAS_Score_Correlation.png")

# Hoang
viper_Hoang_Fib <- viper_X_Hoang %>%
  filter(source=="Fibrosis stage")
viper_Hoang_NAS <- viper_X_Hoang %>%
  filter(source=="NAS")

Hoang_clinical_viper_Fib <- data.frame(clinical=as.data.frame(Y_Hoang)$fibrosis, viper=viper_Hoang_Fib$score)
Hoang_clinical_viper_NAS <- data.frame(clinical=as.data.frame(Y_Hoang)$NAS, viper=viper_Hoang_NAS$score)

Hoang_Fib_scores <- ggscatter(data=Hoang_clinical_viper_Fib, 
                                x="clinical", y="viper", size=2, 
                                cor.coef=TRUE, cor.method="pearson", 
                                title="Hoang Fibrosis Score Correlation"
)
Hoang_NAS_scores <- ggscatter(data=Hoang_clinical_viper_NAS, 
                                x="clinical", y="viper", size=2, 
                                cor.coef=TRUE, cor.method="pearson", 
                                title="Hoang NAS Score Correlation"
)

ggsave("datasets/Hoang_NAS_Score_Correlation.png")

# Pantano
viper_Pantano_Fib <- viper_X_Pantano %>%
  filter(source=="Fibrosis stage")
viper_Pantano_NAS <- viper_X_Pantano %>%
  filter(source=="NAS")

Pantano_clinical_viper_Fib <- data.frame(clinical=as.data.frame(Y_Pantano)$fibrosis, viper=viper_Pantano_Fib$score)
Pantano_clinical_viper_NAS <- data.frame(clinical=as.data.frame(Y_Pantano)$NAS, viper=viper_Pantano_NAS$score)

Pantano_Fib_scores <- ggscatter(data=Pantano_clinical_viper_Fib, 
                              x="clinical", y="viper", size=2, 
                              cor.coef=TRUE, cor.method="pearson", 
                              title="Pantano Fibrosis Score Correlation"
)
Pantano_NAS_scores <- ggscatter(data=Pantano_clinical_viper_NAS, 
                              x="clinical", y="viper", size=2, 
                              cor.coef=TRUE, cor.method="pearson", 
                              title="Pantano NAS Score Correlation"
)

ggsave("datasets/Pantano_NAS_Score_Correlation.png")
