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
rm(list=ls())
gc()
pdf(file="/data/phar-ta-heart/lina3770/e105qc.pdf")  
############## E105 ##############

print("Loading e1051s")
e1051s_raw <- read.csv(file='/data/phar-ta-heart/lina3770/data/Data/Spliced_Unspliced/E10_5_1/E10_5_1_Spliced.csv')
print("Loading e1052s")
e1052s_raw <- read.csv(file='/data/phar-ta-heart/lina3770/data/Data/Spliced_Unspliced/E10_5_2/E10_5_2_Spliced.csv')


#get rid of the 'rownames' from the second batch, and then 
#removes duplicated gne names, then checks this has worked.
e1052s_raw <- e1052s_raw[,-1]
e105s_raw <- cbind(e1051s_raw,e1052s_raw)
e105s_raw <- e105s_raw[!duplicated(e105s_raw[,1]),]

rownames(e105s_raw) <- e105s_raw[,1]
e105s_raw <- e105s_raw[,-1]




e105s <- CreateSeuratObject(counts = e105s_raw, project = "E105S")
print("Finding mitochondrial counts")
e105s[["percent.mt"]] <- PercentageFeatureSet(e105s, pattern = "^mt-")

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


VlnPlot(e125s, features = c("nFeature_RNA"), ncol = 1) 

VlnPlot(e125s, features = c("nFeature_RNA"), ncol = 1) + scale_y_continuous(limits = c(0,2000))

VlnPlot(e125s, features = c("nCount_RNA"), ncol = 1)

VlnPlot(e125s, features = c("nCount_RNA"), ncol = 1) + scale_y_continuous(limits = c(0,12500))

VlnPlot(e125s, features = c("percent.mt"), ncol = 1) 

VlnPlot(e125s, features = c("percent.mt"), ncol = 1) + scale_y_continuous(limits = c(0,10))

FeatureScatter(e125s, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(e125s, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

dev.off()
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
rm(list=ls())
gc()
pdf(file="/data/phar-ta-heart/lina3770/NRCMqc.pdf")
############## NRCM ##############



print("Loading NRCM1s")
NRCM1s_raw <- read.csv(file='/data/phar-ta-heart/lina3770/data/Data/NRCM.Spliced/NRCM/NRCM_Spliced.csv')
print("Loading NRCM2s")
NRCM2s_raw <- read.csv(file='/data/phar-ta-heart/lina3770/data/Data/NRCM.Spliced/NRCM2/NRCM2_Spliced.csv')


#get rid of the 'rownames' from the second batch, and then 
#removes duplicated gne names, then checks this has worked.
NRCM2s_raw <- NRCM2s_raw[,-1]
NRCMs_raw <- cbind(NRCM1s_raw,NRCM2s_raw)
NRCMs_raw <- NRCMs_raw[!duplicated(NRCMs_raw[,1]),]

rownames(NRCMs_raw) <- NRCMs_raw[,1]
NRCMs_raw <- NRCMs_raw[,-1]




NRCMs <- CreateSeuratObject(counts = NRCMs_raw, project = "NRCMS")
print("Finding mitochondrial counts")
NRCMs[["percent.mt"]] <- PercentageFeatureSet(NRCMs, pattern = "^mt-")

#change any NaNs to 0
NRCMs[["percent.mt"]][is.na(NRCMs[["percent.mt"]])] <- 0


VlnPlot(NRCMs, features = c("nFeature_RNA"), ncol = 1) 

VlnPlot(NRCMs, features = c("nFeature_RNA"), ncol = 1) + scale_y_continuous(limits = c(0,2000))

VlnPlot(NRCMs, features = c("nCount_RNA"), ncol = 1)

VlnPlot(NRCMs, features = c("nCount_RNA"), ncol = 1) + scale_y_continuous(limits = c(0,12500))

VlnPlot(NRCMs, features = c("percent.mt"), ncol = 1) 

VlnPlot(NRCMs, features = c("percent.mt"), ncol = 1) + scale_y_continuous(limits = c(0,10))

FeatureScatter(NRCMs, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(NRCMs, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")



dev.off()
rm(list=ls())
gc()

