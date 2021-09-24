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
pdf(file="/data/phar-ta-heart/lina3770/e145qc.pdf")  
############## E145 ##############




print("Loading e1451s")
e1451s_raw <- read.csv(file='/data/phar-ta-heart/lina3770/data/Data/Spliced_Unspliced/E14_5_1/E14_5_1_Spliced.csv')
print("Loading e1452s")
e1452s_raw <- read.csv(file='/data/phar-ta-heart/lina3770/data/Data/Spliced_Unspliced/E14_5_2/E14_5_2_Spliced.csv')


#get rid of the 'rownames' from the second batch, and then 
#removes duplicated gne names, then checks this has worked.
e1452s_raw <- e1452s_raw[,-1]
e145s_raw <- cbind(e1451s_raw,e1452s_raw)
e145s_raw <- e145s_raw[!duplicated(e145s_raw[,1]),]

rownames(e145s_raw) <- e145s_raw[,1]
e145s_raw <- e145s_raw[,-1]




e145s <- CreateSeuratObject(counts = e145s_raw, project = "E145S")
print("Finding mitochondrial counts")
e145s[["percent.mt"]] <- PercentageFeatureSet(e145s, pattern = "^mt-")

#change any NaNs to 0
e145s[["percent.mt"]][is.na(e145s[["percent.mt"]])] <- 0


VlnPlot(e145s, features = c("nFeature_RNA"), ncol = 1) 

VlnPlot(e145s, features = c("nFeature_RNA"), ncol = 1) + scale_y_continuous(limits = c(0,2000))

VlnPlot(e145s, features = c("nCount_RNA"), ncol = 1)

VlnPlot(e145s, features = c("nCount_RNA"), ncol = 1) + scale_y_continuous(limits = c(0,12500))

VlnPlot(e145s, features = c("percent.mt"), ncol = 1) 

VlnPlot(e145s, features = c("percent.mt"), ncol = 1) + scale_y_continuous(limits = c(0,10))

FeatureScatter(e145s, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(e145s, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

dev.off()