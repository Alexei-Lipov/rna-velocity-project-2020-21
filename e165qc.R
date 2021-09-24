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
pdf(file="/data/phar-ta-heart/lina3770/e165qc.pdf")
############## E165 ##############



print("Loading e1651s")
e1651s_raw <- read.csv(file='/data/phar-ta-heart/lina3770/data/Data/Spliced_Unspliced/E16_5_1/E16_5_1_Spliced.csv')
print("Loading e1652s")
e1652s_raw <- read.csv(file='/data/phar-ta-heart/lina3770/data/Data/Spliced_Unspliced/E16_5_2/E16_5_2_Spliced.csv')


#get rid of the 'rownames' from the second batch, and then 
#removes duplicated gne names, then checks this has worked.
e1652s_raw <- e1652s_raw[,-1]
e165s_raw <- cbind(e1651s_raw,e1652s_raw)
rm(e1651s_raw, e1652s_raw)
e165s_raw <- e165s_raw[!duplicated(e165s_raw[,1]),]

rownames(e165s_raw) <- e165s_raw[,1]
e165s_raw <- e165s_raw[,-1]

e165s <- CreateSeuratObject(counts = e165s_raw, project = "E165S")
print("Finding mitochondrial counts")
e165s[["percent.mt"]] <- PercentageFeatureSet(e165s, pattern = "^mt-")

#change any NaNs to 0
e165s[["percent.mt"]][is.na(e165s[["percent.mt"]])] <- 0


VlnPlot(e165s, features = c("nFeature_RNA"), ncol = 1) 

VlnPlot(e165s, features = c("nFeature_RNA"), ncol = 1) + scale_y_continuous(limits = c(0,2000))

VlnPlot(e165s, features = c("nCount_RNA"), ncol = 1)

VlnPlot(e165s, features = c("nCount_RNA"), ncol = 1) + scale_y_continuous(limits = c(0,12500))

VlnPlot(e165s, features = c("percent.mt"), ncol = 1) 

VlnPlot(e165s, features = c("percent.mt"), ncol = 1) + scale_y_continuous(limits = c(0,10))

FeatureScatter(e165s, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(e165s, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

dev.off()