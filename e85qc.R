if(!require(dplyr)) install.packages("dplyr",repos = "http://cran.us.r-project.org")
if(!require(Seurat)) install.packages("Seurat",repos = "http://cran.us.r-project.org")
if(!require(patchwork)) install.packages("Seurat",repos = "http://cran.us.r-project.org")
if(!require(tidyverse)) install.packages("tidyverse",repos = "http://cran.us.r-project.org")
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)


rm(list=ls())        # this and line below clears the memory & "collects the garbage"
gc()
pdf(file="/data/phar-ta-heart/lina3770/e85qc.pdf")  
############## E85 ##############

print("Loading e851s")
e851s_raw <- read.csv(file='/data/phar-ta-heart/lina3770/data/Data/Spliced_Unspliced/E8_5_1/E8_5_1_Spliced.csv')
print("Loading e852s")
e852s_raw <- read.csv(file='/data/phar-ta-heart/lina3770/data/Data/Spliced_Unspliced/E8_5_2/E8_5_2_Spliced.csv')


#get rid of the 'rownames' from the second batch, and then 
#removes duplicated gne names, then checks this has worked.
e852s_raw <- e852s_raw[,-1]
e85s_raw <- cbind(e851s_raw,e852s_raw)
e85s_raw <- e85s_raw[!duplicated(e85s_raw[,1]),]

rownames(e85s_raw) <- e85s_raw[,1]
e85s_raw <- e85s_raw[,-1]




e85s <- CreateSeuratObject(counts = e85s_raw, project = "E85S")
print("Finding mitochondrial counts")
e85s[["percent.mt"]] <- PercentageFeatureSet(e85s, pattern = "^mt-")

#change any NaNs to 0
e85s[["percent.mt"]][is.na(e85s[["percent.mt"]])] <- 0


VlnPlot(e85s, features = c("nFeature_RNA"), ncol = 1) 

VlnPlot(e85s, features = c("nFeature_RNA"), ncol = 1) + scale_y_continuous(limits = c(0,2000))

VlnPlot(e85s, features = c("nCount_RNA"), ncol = 1)

VlnPlot(e85s, features = c("nCount_RNA"), ncol = 1) + scale_y_continuous(limits = c(0,12500))

VlnPlot(e85s, features = c("percent.mt"), ncol = 1) 

VlnPlot(e85s, features = c("percent.mt"), ncol = 1) + scale_y_continuous(limits = c(0,10))

FeatureScatter(e85s, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(e85s, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

dev.off()