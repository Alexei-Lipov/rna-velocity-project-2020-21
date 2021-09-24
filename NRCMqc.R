if(!require(dplyr)) install.packages("dplyr",repos = "http://cran.us.r-project.org")
if(!require(Seurat)) install.packages("Seurat",repos = "http://cran.us.r-project.org")
if(!require(patchwork)) install.packages("Seurat",repos = "http://cran.us.r-project.org")
if(!require(tidyverse)) install.packages("tidyverse",repos = "http://cran.us.r-project.org")
if(!require(Matrix)) install.packages("Matrix",repos = "http://cran.us.r-project.org")
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(Matrix)

rm(list=ls())
gc()
pdf(file="/data/phar-ta-heart/lina3770/NRCMqc.pdf")
############## NRCM ##############




Echonk <- readMM("/data/phar-ta-heart/corp2591/CMs/stages/Neo/matrix.mtx")
Echonk_cell_ids <- read_tsv("/data/phar-ta-heart/corp2591/CMs/stages/Neo/barcodes.tsv", col_names=FALSE)$X1
Echonk_cell_ids <- paste0("Neo.",Echonk_cell_ids)
Echonk_genes <- read_tsv("/data/phar-ta-heart/corp2591/CMs/stages/Neo/features.tsv", col_names=FALSE)
Echonk_gene_ids <- Echonk_genes$X1
rownames(Echonk) <- Echonk_gene_ids
colnames(Echonk) <- Echonk_cell_ids

annots_clust <- readRDS("/data/phar-ta-heart/corp2591/CMs/Echonko_rownames_annots.csv")
print(annots_clust[1:10,2])
print(annots_clust[1:10,1])

CMs <- Echonk
CMs <- CMs[(intersect(rownames(CMs),annots_clust[,2])),]
print("Debug 1")
for(i in 1:length(rownames(CMs))){
  
  for(j in 1:length(annots_clust[,2])){
    if(rownames(CMs)[i] == annots_clust[j,2]){
       rownames(CMs)[i] <- annots_clust[j,1]
    }
  }
}






#print("Loading NRCM1s")
#NRCM1s_raw <- read.csv(file='/data/phar-ta-heart/lina3770/data/Data/NRCM.Spliced/NRCM/NRCM_Spliced.csv')
#print("Loading NRCM2s")
#NRCM2s_raw <- read.csv(file='/data/phar-ta-heart/lina3770/data/Data/NRCM.Spliced/NRCM2/NRCM2_Spliced.csv')


#get rid of the 'rownames' from the second batch, and then 
#removes duplicated gne names, then checks this has worked.
#NRCM2s_raw <- NRCM2s_raw[,-1]
#NRCMs_raw <- cbind(NRCM1s_raw,NRCM2s_raw)
#NRCMs_raw <- NRCMs_raw[!duplicated(NRCMs_raw[,1]),]

#rownames(NRCMs_raw) <- NRCMs_raw[,1]
#NRCMs_raw <- NRCMs_raw[,-1]



print("Debug 2")
NRCMs <- CreateSeuratObject(counts = CMs, project = "NRCMS")
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