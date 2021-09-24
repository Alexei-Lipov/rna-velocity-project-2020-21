if(!require(dplyr)) install.packages("dplyr",repos = "http://cran.us.r-project.org")
if(!require(Seurat)) install.packages("Seurat",repos = "http://cran.us.r-project.org")
if(!require(patchwork)) install.packages("Seurat",repos = "http://cran.us.r-project.org")
if(!require(tidyverse)) install.packages("tidyverse",repos = "http://cran.us.r-project.org")
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
rm(list=ls())
gc()
pdf(file="/data/phar-ta-heart/lina3770/e125qc.pdf")  
############## E125 ##############



print("Loading e1251s")
e1251s_raw <- read.csv(file='/data/phar-ta-heart/lina3770/data/Data/Spliced_Unspliced/E12_5_1/E12_5_1_Spliced.csv')
print("Loading e1252s")
e1252s_raw <- read.csv(file='/data/phar-ta-heart/lina3770/data/Data/Spliced_Unspliced/E12_5_2/E12_5_2_Spliced.csv')


#get rid of the 'rownames' from the second batch, and then 
#removes duplicated gne names, then checks this has worked.
e1252s_raw <- e1252s_raw[,-1]
e125s_raw <- cbind(e1251s_raw,e1252s_raw)
e125s_raw <- e125s_raw[!duplicated(e125s_raw[,1]),]

rownames(e125s_raw) <- e125s_raw[,1]
e125s_raw <- e125s_raw[,-1]




e125s <- CreateSeuratObject(counts = e125s_raw, project = "E125S")
print("Finding mitochondrial counts")
e125s[["percent.mt"]] <- PercentageFeatureSet(e125s, pattern = "^mt-")

#change any NaNs to 0
e125s[["percent.mt"]][is.na(e125s[["percent.mt"]])] <- 0


VlnPlot(e125s, features = c("nFeature_RNA"), ncol = 1, pt.size = 0.3) 

VlnPlot(e125s, features = c("nFeature_RNA"), ncol = 1, pt.size = 0.3) + scale_y_continuous(limits = c(0,2000))

VlnPlot(e125s, features = c("nCount_RNA"), ncol = 1, pt.size = 0.3)

VlnPlot(e125s, features = c("nCount_RNA"), ncol = 1, pt.size = 0.3) + scale_y_continuous(limits = c(0,12500))

VlnPlot(e125s, features = c("percent.mt"), ncol = 1, pt.size = 0.3) 

VlnPlot(e125s, features = c("percent.mt"), ncol = 1, pt.size = 0.3) + scale_y_continuous(limits = c(0,10))

FeatureScatter(e125s, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.3)
FeatureScatter(e125s, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.3)

dev.off()