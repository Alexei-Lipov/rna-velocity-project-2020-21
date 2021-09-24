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
pdf(file="/data/phar-ta-heart/lina3770/e105qc.pdf")  
############## E105 ##############

print("Loading e1051s")
e1051s_raw <- read.csv(file='/data/phar-ta-heart/lina3770/data/Data/Spliced_Unspliced/E10_5_1/E10_5_1_Spliced.csv')
print("Loading e1052s")
e1052s_raw <- read.csv(file='/data/phar-ta-heart/lina3770/data/Data/Spliced_Unspliced/E10_5_2/E10_5_2_Spliced.csv')
print("Debug 1")

#get rid of the 'rownames' from the second batch, and then 
#removes duplicated gne names, then checks this has worked.
e1052s_raw <- e1052s_raw[,-1]
print("Debug 2")
e105s_raw <- cbind(e1051s_raw,e1052s_raw)
print("Debug 3")
ls()
rm(e1051s_raw, e1052s_raw)
e105s_raw <- e105s_raw[!duplicated(e105s_raw[,1]),]
print("Debug 4")

rownames(e105s_raw) <- e105s_raw[,1]
print("Debug 5")
e105s_raw <- e105s_raw[,-1]
print("Debug 6")
ls()
e105s_raw <- as.matrix(e105s_raw)
print("Debug 6.5")
e105s_raw <- as.sparse(e105s_raw)

ls()
e105s <- CreateSeuratObject(counts = e105s_raw, project = "E105S")
print("Debug 7")
rm(e105s_raw)
print("Finding mitochondrial counts")
e105s[["percent.mt"]] <- PercentageFeatureSet(e105s, pattern = "^mt-")
print("Debug 8")

#change any NaNs to 0
e105s[["percent.mt"]][is.na(e105s[["percent.mt"]])] <- 0


VlnPlot(e105s, features = c("nFeature_RNA"), ncol = 1) 

VlnPlot(e105s, features = c("nFeature_RNA"), ncol = 1) + scale_y_continuous(limits = c(0,2000))

VlnPlot(e105s, features = c("nCount_RNA"), ncol = 1)

VlnPlot(e105s, features = c("nCount_RNA"), ncol = 1) + scale_y_continuous(limits = c(0,12500))

VlnPlot(e105s, features = c("percent.mt"), ncol = 1) 

VlnPlot(e105s, features = c("percent.mt"), ncol = 1) + scale_y_continuous(limits = c(0,10))

FeatureScatter(e105s, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(e105s, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")



dev.off()