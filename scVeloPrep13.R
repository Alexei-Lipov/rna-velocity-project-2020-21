## THIS VERSION IS FOR UNTRUNCATED LARGER DATASET WITH NO SCTRANSFORM (NO BATCH REGRESSION) - THIS IS VERSION IS MEANT TO BE PUT THROUGH SCANPY PROCESSING

#ntasks-per-node=16 is an SBATCH that I removed
#ulimit -s 2560000000 is part of slurm_submission_r that I removed

library(BiocManager)
#BiocManager::install("BiocGenerics", lib="/data/phar-ta-heart/lina3770/RLib")
#BiocManager::install("glmGamPoi", lib="/data/phar-ta-heart/lina3770/RLib") 
library(glmGamPoi)
# Install devtools from CRAN
#install.packages("devtools", lib="/data/phar-ta-heart/lina3770/RLib", repos='http://cran.us.r-project.org')
# Use devtools to install hdf5r and loomR from GitHub
#devtools::install_github(repo = "hhoeflin/hdf5r", lib="/data/phar-ta-heart/lina3770/RLib")
#devtools::install_github(repo = "mojaveazure/loomR", ref = "develop", lib="/data/phar-ta-heart/lina3770/RLib")
# Load loomR
library(loomR)
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(data.table)
library(sctransform)
library(BiocManager)
library(sctransform)
#library(scater)

# to get past memory limit error in sctransform when it comes to E10.5
options(future.globals.maxSize = 256000 * 1024^2)

setwd("/data/phar-ta-heart/lina3770/data/Data")
print(paste("working directory:", getwd()))

#file_list <- list.files(pattern = "*Unspliced.csv", recursive=TRUE)

stage_list <- c("E8","E10","E12","E14","E16","NRCM")

# I have hard coded the number of rows - this is bad, should be changed at some point somehow
corrected_unspliced <- matrix(, nrow = 27933, ncol = 0)

#print(list.files(pattern = "E14.*Unspliced", recursive=TRUE))


for (stage in stage_list){
  # finds files with stage, e.g. E14, and Unspliced in filename, i.e. gives the two batches of the stage
  file_list <- list.files(pattern = paste0(stage,".*Unspliced"), recursive=TRUE)

  print(paste("reading in", file_list[1]))
  batch_1 <- fread(file = file_list[1])
  print(paste("reading in", file_list[2]))
  batch_2 <- fread(file = file_list[2])
  
  print("converting from dataframe to matrix")
  batch_1 <- as.matrix(batch_1)
  batch_2 <- as.matrix(batch_2)
  
  rownames(batch_1) <- batch_1[,1]
  batch_1 <- batch_1[,-1]
  rownames(batch_2) <- batch_2[,1]
  batch_2 <- batch_2[,-1]
  
  batch_1 <- batch_1[unique(rownames(batch_1)),]
  batch_2 <- batch_2[unique(rownames(batch_2)),]
  
  class(batch_1) <- "numeric"
  class(batch_2) <- "numeric"
  
  print("dimension of batch 1:")
  print(dim(batch_1))
  print("dimension of batch 2:")
  print(dim(batch_2))
  
  stage <- cbind(batch_1,batch_2)
  
  print("concatenated to give matrix of dimension:")
  print(dim(stage))

  batch_labels <- c(rep(1, dim(batch_1)[2]), rep(2, dim(batch_2)[2]))
  
  # creates cell attributes dataframe (needed for sctransform:vst)
  cell_attr <- data.frame(n_umi = colSums(stage), n_gene = colSums(stage > 0))
  rownames(cell_attr) <- colnames(stage)
  cell_attr$batch_var <- batch_labels
  
  # Go through each row and determine if a value is zero (we are finding cells with zero n_umi & n_gene)
  row_sub = apply(cell_attr, 1, function(row) all(row !=0 ))
  
  # Filter out the cells with zero n_umi & zero n_gene
  cell_attr <- cell_attr[row_sub,]
  stage <- stage[,row_sub]
  
  # not sure if I need this part (and latent_var =, in sct) or not
  #cell_attr$umi_per_gene <- cell_attr$n_umi / cell_attr$n_gene
  #cell_attr$log_umi_per_gene <- log10(cell_attr$umi_per_gene)      # used to be log10, had NaN errors so changed it to log1p
  
  #possible needed argument: , latent_var = 'log_umi_per_gene'
  #vst_out <- sctransform::vst(stage, cell_attr = cell_attr, batch_var = "batch_var", return_corrected_umi = TRUE, method="glmGamPoi")
  
  #corrected_umi <- vst_out$umi_corrected
  corrected_umi <- stage
  print("dimension of the corrected umi count matrix after batch regressing using sctransform:")
  print(dim(corrected_umi))
  
  # if number of genes is not the same then we pad with zeros in order to be able to concatenate after this
  if (dim(corrected_unspliced)[1] > dim(corrected_umi)[1]){
  
  difference = dim(corrected_unspliced)[1] - dim(corrected_umi)[1]
  
  addition <- matrix(0L, nrow=difference, ncol=dim(corrected_umi)[2])
  
  corrected_umi <- rbind(corrected_umi, addition)
  
  } else if (dim(corrected_unspliced)[1] < dim(corrected_umi)[1]){
  
  difference = dim(corrected_umi)[1] - dim(corrected_unspliced)[1]
  
  addition <- matrix(0L, nrow=difference, ncol=dim(corrected_unspliced)[2])
  
  corrected_unspliced <- rbind(corrected_unspliced, addition)
  
  }
  
  corrected_unspliced <- cbind(corrected_unspliced, corrected_umi)
  print("dimension of the corrected_unspliced matrix (where all batches & stages are appended to):")
  print(dim(corrected_unspliced))
  print("debugging rownames error")
  #print(rownames(corrected_unspliced))
  
}

####################

print("read in ootlist_Spliced.csv (new larger spliced dataset)")
ootlist_Spliced <- readRDS("/data/phar-ta-heart/lina3770/ootlist_Spliced.csv")

new_spliced <- matrix(, nrow = 27477, ncol = 0)

for (i in c(1,2,3,4,5,6)){
  print(paste("Stage",str(i)))
  ootlist_temp <- data.matrix(ootlist_Spliced[[i]])
  print(dim(ootlist_temp))
  
  rownames(ootlist_temp) <- row.names(ootlist_Spliced[[i]])
  colnames(ootlist_temp) <- colnames(ootlist_Spliced[[i]])
  
  new_spliced <- cbind(new_spliced, ootlist_temp)
  print(dim(new_spliced))
  }

###################

#saveRDS(corrected_unspliced, 'corrected_unspliced.RDS')

#print("read in spliced")
#data_spliced <-  readRDS('../Lipov.csv')

# make the cell labels consistent
colnames(new_spliced) <- sub("^(.*)_(.*)-(.*)$", "\\1_\\3:\\2", colnames(new_spliced))
colnames(corrected_unspliced) <- sub("_","",colnames(corrected_unspliced))
colnames(corrected_unspliced) <- sub("NRCM","Neo_",colnames(corrected_unspliced))
colnames(corrected_unspliced) <- sub("Neo_:","Neo_1:",colnames(corrected_unspliced))

print("column names of new_spliced after correction")
print(colnames(new_spliced[1:10,1:10]))
print(colnames(new_spliced[1:10,160000:160010]))

print("column names of corrected_unspliced after correction")
print(colnames(corrected_unspliced[1:10,1:10]))

print("dim of corrected_unspliced (check to see if number of cells is less than for new_spliced)")
print(dim(corrected_unspliced))

print("dim of new_spliced")
print(dim(new_spliced))

print("dimension of corrected_unspliced before & after filtering according to cells of the spliced, respectively:")
print(dim(corrected_unspliced))
# filter the unspliced according to the cells of the spliced
corrected_unspliced <- corrected_unspliced[,intersect(colnames(new_spliced),colnames(corrected_unspliced))]
print(dim(corrected_unspliced))

#change any NaNs to 0
corrected_unspliced[is.na(corrected_unspliced)] <- 0



#print(rownames(corrected_unspliced))
#corrected_unspliced <- corrected_unspliced[1:625,]
#print(rownames(corrected_unspliced))

print(dim(corrected_unspliced))


# creates cell attributes dataframe (needed for sctransform:vst)
cell_attr <- data.frame(n_umi = colSums(as.array(corrected_unspliced)), n_gene = colSums(as.array(corrected_unspliced) > 0))
rownames(cell_attr) <- colnames(corrected_unspliced)

# Go through each row and determine if a value is zero (we are finding cells with zero n_umi & n_gene)
row_sub = apply(cell_attr, 1, function(row) all(row !=0 ))

# Filter out the cells with zero n_umi & zero n_gene
cell_attr <- cell_attr[row_sub,]
corrected_unspliced <- corrected_unspliced[,row_sub]
print("after filtering out cells with zero n_umi & zero n_gene")
print(dim(corrected_unspliced))

write.csv(corrected_unspliced, "spliced_unprocessed_untruncated_150221.csv")
write.csv(new_spliced, "unspliced_unprocessed_untruncated_150221.csv")


#unspliced <- CreateSeuratObject(counts = corrected_unspliced, project = "unspliced", meta.data=NULL)
#rm(corrected_unspliced)
#spliced <- CreateSeuratObject(counts = data_spliced, project = "spliced", meta.data=NULL)
#rm(data_spliced)

#spliced[["percent.mt"]] <- PercentageFeatureSet(spliced, pattern = "^mt-")
#unspliced[["percent.mt"]] <- PercentageFeatureSet(unspliced, pattern = "^mt-")


#spliced[["nCount_RNA"]] <- spliced@meta.data$nCount_RNA
#unspliced[["nCount_RNA"]] <- unspliced@meta.data$nCount_RNA

#, method = 'glmGamPoi'
#print("Performing SCTransform on spliced dataset of dimension:")
#print(dim(spliced@assays$RNA@counts))
#spliced_test <- SCTransform(spliced, vars.to.regress = c("percent.mt","nCount_RNA"), verbose = TRUE, min_cells=1, return.only.var.genes=FALSE)
#rm(spliced)
#print("Performing SCTransform on unspliced dataset of dimension:")
#print(dim(unspliced@assays$RNA@counts))
#print(unspliced@assays$RNA@counts[1:3,1:3])
#print(max(unspliced@meta.data$nCount_RNA))
#print(min(unspliced@meta.data$nCount_RNA))
#unspliced_test <- SCTransform(unspliced, vars.to.regress = c("percent.mt","nCount_RNA"), verbose = TRUE, min _cells=1, return.only.var.genes=FALSE)
#rm(unspliced)
#bay_out <- bayNorm(spliced@assays[["RNA"]]@counts,mean_version = TRUE,verbose=TRUE)

#write.csv(spliced_test[["SCT"]]@data, "spliced_processed_untruncated_1.csv")
#write.csv(unspliced_test[["SCT"]]@data, "unspliced_processed_untruncated_1.csv")
#saveRDS(spliced_test[["SCT"]]@data, file = "spliced_processed_truncated.RDS")
#saveRDS(unspliced_test[["SCT"]]@data, file = "unspliced_processed_truncated.RDS"








###### possibly useful code snippets ####

#print(class(ootlist_Spliced))
#print(class(ootlist_Spliced[[1]]))


#print(ootlist_Spliced[[1]][1:10,1:10])

#print("read in old spliced")
#data_spliced <-  readRDS('../Lipov.csv')

#print("colnames new_spliced")
#print(colnames(new_spliced[1:10,1:10]))
#print(colnames(new_spliced[1:10,160000:160010]))
#print("colnames old spliced")
#print(colnames(data_spliced[1:10,1:10]))
#print(colnames(new_spliced[1:10,23000:23010]))


#print(esf
#print(ootlist_Spliced[[1]])
#print("converting to matrix, first time")
#ootlist_test <- data.matrix(ootlist_Spliced[[1]])


#print(row.names(ootlist_Spliced[[1]][1:10,1:10]))
#print(colnames(ootlist_Spliced[[1]][1:10,1:10]))

#rownames(ootlist_test) <- row.names(ootlist_Spliced[[1]])
#colnames(ootlist_test) <- colnames(ootlist_Spliced[[1]])

#print(rownames(ootlist_test[1:10,1:10]))
#print(colnames(ootlist_test[1:10,1:10]))



#print(class(ootlist_test))
#print(dim(ootlist_test))
#print(ootlist_test[1:10,1:10])
#print(rownames(ootlist_test[1:10]))
#print(colnames(ootlist_test[1:10]))

#print(ootlist_test[1,1:10])
#print(ootlist_test[1:10,1])



#print(lengths(ootlist_Spliced))
#print(lengths(ootlist_Spliced[1]))

#print("debug")
#new_spliced <- c(as.matrix(ootlist_Spliced[[1]]), as.matrix(ootlist_Spliced[[2]]))
#print(dim(new_spliced[1]))
#print(dim(new_spliced[2]))
#new_spliced <- rep(NA, 6)
#print(new_spliced)
#print(class(new_spliced))

#print(new_spliced[1])
#new_spliced[1] <- as.matrix(ootlist_Spliced[[1]])
#print(new_spliced)
#print(class(new_spliced))
#print(dim(new_spliced))


#for (i in c(1,2,3,4,5,6)){  
#  print(ootlist)
#  print(ootlist[[i]])
#  print(dim(as.matrix(ootlist_Spliced[[i]])))
#  print(class(as.matrix(ootlist_Spliced[[i]])))
#  ootlist[[i]] <- as.matrix(ootlist_Spliced[[i]])
#  print(class(ootlist[[i]]))
#  print(dim(ootlist[[i]]))
#  rownames(ootlist[[i]]) <- ootlist[[i]][,1]
#  ootlist[[i]] <- ootlist[[i]][,-1]
#  print(rownames(ootlist[[i]][1:10]))
  
#  colnames(ootlist[[i]]) <- ootlist[[i]][1,]
#  ootlist[[i]] <- ootlist[[i]][-1,]
#  print(colnames(ootlist[[i]][1:10]))
  
  
#  }

#print("corrected_unspliced")
#print(colnames(corrected_unspliced))
#print(rownames(corrected_unspliced))
#print("ootlist")
#print(colnames(ootlist_Spliced[1]))
#print(rownames(ootlist_Spliced[1]))
#print(len(ootlist_Spliced))

#print("read in old spliced")
#data_spliced <-  readRDS('../Lipov.csv')