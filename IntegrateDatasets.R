# BRCA 1 skull metastasis (own data 10X) + 4 primary + 4 paired lymph node metastases(GSE225600)

setwd('/BCBM/0.rawdata')
# library packages ====
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)

# Integration 2 datasets ====
## load GSE225600 data
gse225600 <- readRDS('GSE225600_rm_Doublets.rds')

## load own data
ownData <- readRDS('ownData_rm_Doublets.rds')

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list(gse225600,ownData),
                                      nfeatures = 5000)

inte.anchors <- FindIntegrationAnchors(object.list = list(gse225600,ownData), 
                                       anchor.features = features,
                                       dims = 1:50) 

# this command creates an 'integrated' data assay
bcbm <- IntegrateData(anchorset = inte.anchors,
                      normalization.method = "LogNormalize") # "SCT"
barcodes <- bcbm@assays[["RNA"]]@data@Dimnames[[2]] %>% strsplit('-') %>% 
  unlist() %>% 
  matrix(ncol = 2,byrow = T) %>% as.data.frame()

barcodes$V3 <- gsub('L','LN',barcodes$V2)
barcodes$V3 <- gsub('T','PT',barcodes$V3)
barcodes$V3 <- gsub('1','BoM',barcodes$V3)
table(barcodes$V3)

# add barcode annotation
bcbm <- AddMetaData(bcbm,barcodes$V3,col.name = 'Patient')
#37333 features across 34626 samples

tumorType <- bcbm$Patient %>% as.character()
tumorType[grep('LN',tumorType)] <- 'LN'
tumorType[grep('PT',tumorType)] <- 'PT'
tumorType <- factor(tumorType)
bcbm<- AddMetaData(bcbm,tumorType,col.name = 'tumorType')

saveRDS(bcbm,'BCBM_Raw_SeuratObj_seuratIntegration.rds')

bcbm_sub.7 <- subset(bcbm, Patient %in% c('BoM','LN2','LN3','LN7','PT2','PT6','PT7'))
#40333 features across 34375 cells; 7 samples

table(bcbm_sub.7$Patient)
#BoM  LN2  LN3  LN7  PT2  PT6  PT7 
#6208 4181 4315 7244 2904 3390 6133

saveRDS(bcbm_sub.7,'BCBM_sub7samples_SeuratObj.rds')

rm(list = ls());gc()

# 7 samples scale data ====
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
bcbm <- readRDS('BCBM_sub7samples_SeuratObj.rds')
#DefaultAssay(bcbm) <- "integrated"

# Run the standard workflow for visualization and clustering
bcbm <- ScaleData(bcbm, verbose = T)

# PCA 降维 ====
bcbm <- RunPCA(bcbm, npcs = 30, 
               verbose = T,
               features = VariableFeatures(object = bcbm))

## 确定数据集维度
bcbm <- JackStraw(bcbm, num.replicate = 100)
bcbm <- ScoreJackStraw(bcbm, dims = 1:20)
#JackStrawPlot(bcbm, dims = 1:20)
pdf('inteData_PCA_elbow_7samples.pdf',width = 6,height = 4,onefile = F)
ElbowPlot(bcbm) # 20
dev.off()

bcbm <- RunUMAP(bcbm, reduction = "pca", dims = 1:20)

# Cluster ====
bcbm <- FindNeighbors(bcbm,#reduction = "harmony",
                      dims = 1:20)

## set resolution
bcbm <- FindClusters(bcbm, resolution = .3)

# UMAP & tSNE 降维 ====
bcbm <- RunUMAP(bcbm, dims = 1:20)
bcbm <- RunTSNE(bcbm, dims = 1:20)

## plot
pdf('inteData_UMAP.pdf',width = 20,height = 5.5,onefile = F)
DimPlot(bcbm, reduction = "umap", label = TRUE, 
        group.by = c("seurat_clusters", "Patient", "orig.ident"))
dev.off()
pdf('inteData_tSNE.pdf',width = 20,height = 5.5,onefile = F)
DimPlot(bcbm, reduction = "tsne", label = TRUE, 
        group.by = c("seurat_clusters", "Patient", "orig.ident"))
dev.off()

saveRDS(bcbm,'BCBM_inteData7samples_dimRedu.rds')
#40333 features across 34375 cells; 7 samples
