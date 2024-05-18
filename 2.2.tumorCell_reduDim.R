setwd('BCBM/2.2.tumor_reduceDim')
# library packages ====
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)

# load data ====
bcbm <- readRDS('../1.anno/BCBM_7samples_commonCellTypes.rds')
epi <- subset(bcbm, commonCellTypes == 'Epithelial')
#40333 features across 3768 samples

## normalize ====
epi <- NormalizeData(object = epi, assay = 'RNA',
                     normalization.method = "LogNormalize",
                     scale.factor = 10000,
                     margin = 1, #If performing CLR normalization, normalize across features (1) or cells (2)
                     verbose = TRUE) 

# Run the standard workflow for visualization and clustering
epi <- ScaleData(epi, verbose = T)

# PCA 降维 ====
epi <- RunPCA(epi, npcs = 30, 
              verbose = T,
              #features = VariableFeatures(object = epi)
              )

## 确定数据集维度
epi <- JackStraw(epi, num.replicate = 100)
epi <- ScoreJackStraw(epi, dims = 1:20)
#JackStrawPlot(epi, dims = 1:20)
pdf('epi_PCA_elbow.pdf',width = 6,height = 4,onefile = F)
ElbowPlot(epi) # 20
dev.off()

# Cluster ====
epi <- FindNeighbors(epi,reduction = "pca", 
                     dims = 1:20)

epi <- FindClusters(epi, resolution = .1) %>% 
  FindClusters( resolution = .2) %>% FindClusters( resolution = .3) %>%
  FindClusters( resolution = .4) %>% FindClusters( resolution = .5) %>% 
  FindClusters( resolution = .6) %>% FindClusters( resolution = .7) %>% 
  FindClusters( resolution = .8) %>% FindClusters( resolution = .9) %>% 
  FindClusters( resolution = 1) %>% FindClusters( resolution = 1.1) %>% 
  FindClusters( resolution = 1.2) %>% FindClusters( resolution = 1.3) %>% 
  FindClusters( resolution = 1.4)

## plot
library(clustree)
pdf('epi_clustree.pdf',width = 15,height = 11,onefile = F)
clustree(epi)
dev.off()

## set resolution
epi <- FindClusters(epi, resolution = .8)

# UMAP & tSNE 降维 ====
epi <- RunUMAP(epi, dims = 1:20)
epi <- RunTSNE(epi, dims = 1:20)

## plot
pdf('epi_UMAP.pdf',width = 13,height = 4,onefile = F)
DimPlot(epi, reduction = "umap", label = TRUE, 
        group.by = c("seurat_clusters", "tumorType", "orig.ident"))
dev.off()
pdf('epi_tSNE.pdf',width = 13,height = 4,onefile = F)
DimPlot(epi, reduction = "tsne", label = TRUE, 
        group.by = c("seurat_clusters", "tumorType", "orig.ident"))
dev.off()

table(epi$tumorType,epi$seurat_clusters)

## plot cell ratio ====
library(plyr)
library(ggplot2)

plotData <- table(epi$tumorType,epi$seurat_clusters) %>% data.frame()# %>% 
#reshape2::dcast(Var1~Var2,value.var = 'Freq')
colnames(plotData) <- c('Patient','Cluster','Freq')
plotData$Patient <- factor(plotData$Patient, levels = c('PT','LN','BoM'))

plotData <- ddply(plotData, "Cluster", transform,
                  percent = Freq / sum(Freq) * 100)

pdf('tumorCell_Patient_Percent.pdf',width = 7,height = 5,onefile = F)
ggplot(plotData, aes(x=Cluster, y=percent, fill=Patient)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('#ff7e67','#ffc93c','#0dceda','#6ef3d6',
                             '#a8e6cf','#fdffab','#ffd3b6','#ffaaa5',
                             '#11999e','#3f72af','#ffc7c7','#9896f1',
                             '#fbafaf','#f2c6b4','#f3e8cb','#99e1e5',
                             '#ff9a3c','#ff6f3c','#155263','#a6d0e4',
                             '#c86b85','#e6a4b4','#f3d7ca','#bfcfff',
                             '#a5dee5','#e0f9b5','#ffcfdf','#e8ffe8',
                             '#a6fff2','#74f9ff','#00e0ff','#a56cc1',
                             '#15b7b9','#10ddc2','#f57170','#ff9de2',
                             '#878ecd','#b9bbdf','#dde7f2','#dff4f3'))+
  
  theme_bw()+
  labs(x='Cluster',y='Percent',fill='Patient')
dev.off() 

plotData <- ddply(plotData, "Patient", transform,
                  percent = Freq / sum(Freq) * 100)

pdf('tumorCell_Cluster_Percent.pdf',width = 7,height = 5,onefile = F)
ggplot(plotData, aes(x=Patient, y=percent, fill=Cluster)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('#a8e6cf','#fdffab','#ffd3b6','#ffaaa5',
                             '#11999e','#3f72af','#ffc7c7','#9896f1',
                             '#ff7e67','#0dceda','#6ef3d6','#a6d0e4',
                             '#fbafaf','#f2c6b4','#f3e8cb','#99e1e5',
                             '#ffc93c','#ff9a3c','#ff6f3c','#155263',
                             '#c86b85','#e6a4b4','#f3d7ca','#bfcfff',
                             '#a5dee5','#e0f9b5','#ffcfdf','#e8ffe8',
                             '#a6fff2','#74f9ff','#00e0ff','#a56cc1',
                             '#15b7b9','#10ddc2','#f57170','#ff9de2',
                             '#878ecd','#b9bbdf','#dde7f2','#dff4f3'))+
  
  theme_bw()+
  labs(x='Patient',y='Percent',fill='Cluster')
dev.off() 

# find markers ====
allmarkers <- FindAllMarkers(
  object = epi, assay = NULL,
  only.pos= TRUE,min.pct = .25, logfc.threshold = 0.25) %>% 
  group_by(cluster)

table(allmarkers$cluster)
#0   1   2   3   4   5   6   7   8   9  10  11  12  13  14 
#17 379 282  50 165 273 502 196 191 234 196 148 544 486 352 

## add Patien cell types metaData ====
epiSubtypes <- ifelse(epi$seurat_clusters %in% c(0,3,5,6,12,13),'BoM',
                      ifelse(epi$seurat_clusters %in% c(1,2,4,7:10),'PT',
                             ifelse(epi$seurat_clusters %in% c(11),'PT+LN',
                                    'PT+BoM')
                      )) %>% as.factor()

epi <- AddMetaData(epi,epiSubtypes,col.name = 'epiSubtypes')
table(epi$epiSubtypes)
table(epi$epiSubtypes,epi$Patient)
# BoM   BoM+PT  PT   PT+LN 
#2336     40   2411    155 

### plot ====
plotData <- table(epi$tumorType,epi$epiSubtypes) %>% data.frame()# %>% 
#reshape2::dcast(Var1~Var2,value.var = 'Freq')
colnames(plotData) <- c('Patient','Cluster','Freq')
plotData$Patient <- factor(plotData$Patient, levels = c('PT','LN','BoM'))
plotData$Cluster <- factor(plotData$Cluster, levels = c('PT','PT+LN','PT+BoM','BoM'))

plotData <- ddply(plotData, "Cluster", transform,
                  percent = Freq / sum(Freq) * 100)

pdf('tumorCell_subtype_Patient_Percent.pdf',width = 7,height = 5,onefile = F)
ggplot(plotData, aes(x=Cluster, y=percent, fill=Patient)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('#11999e','#3f72af','#f57170','#ffc7c7',
                             '#ff7e67','#ffc93c','#0dceda','#6ef3d6',
                             '#a8e6cf','#fdffab','#ffd3b6','#ffaaa5',
                             '#9896f1',
                             '#a6d0e4',
                             '#fbafaf','#f2c6b4','#f3e8cb','#99e1e5',
                             '#ffc93c','#ff9a3c','#ff6f3c','#155263',
                             '#c86b85','#e6a4b4','#f3d7ca','#bfcfff',
                             '#a5dee5','#e0f9b5','#ffcfdf','#e8ffe8',
                             '#a6fff2','#74f9ff','#00e0ff','#a56cc1',
                             '#15b7b9','#10ddc2','#f57170','#ff9de2',
                             '#878ecd','#b9bbdf','#dde7f2','#dff4f3'))+
  
  theme_bw()+
  labs(x='Cluster',y='Percent',fill='Patient')
dev.off() 

plotData <- ddply(plotData, "Patient", transform,
                  percent = Freq / sum(Freq) * 100)

pdf('tumorCell_subtype_Cluster_Percent.pdf',width = 7,height = 5,onefile = F)
ggplot(plotData, aes(x=Patient, y=percent, fill=Cluster)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('#a8e6cf','#fdffab','#ffd3b6','#ffaaa5',
                             '#11999e','#3f72af','#ffc7c7','#9896f1',
                             '#ff7e67','#0dceda','#6ef3d6','#a6d0e4',
                             '#fbafaf','#f2c6b4','#f3e8cb','#99e1e5',
                             '#ffc93c','#ff9a3c','#ff6f3c','#155263',
                             '#c86b85','#e6a4b4','#f3d7ca','#bfcfff',
                             '#a5dee5','#e0f9b5','#ffcfdf','#e8ffe8',
                             '#a6fff2','#74f9ff','#00e0ff','#a56cc1',
                             '#15b7b9','#10ddc2','#f57170','#ff9de2',
                             '#878ecd','#b9bbdf','#dde7f2','#dff4f3'))+
  
  theme_bw()+
  labs(x='Patient',y='Percent',fill='Cluster')
dev.off() 

pdf('epi_full_tSNE.pdf',width = 11,height = 8,onefile = F)
DimPlot(epi, reduction = "tsne", label = TRUE, 
        group.by = c("epiSubtypes","seurat_clusters", "tumorType", "orig.ident"))
dev.off()

pdf('epi_full_UMAP.pdf',width = 11,height = 8,onefile = F)
DimPlot(epi, reduction = "umap", label = TRUE, 
        group.by = c("epiSubtypes","seurat_clusters", "tumorType", "orig.ident"))
dev.off()

pdf('epi_subType_UMAP.pdf',width = 8,height = 6.5,onefile = F)
DimPlot(epi, reduction = "umap", label = TRUE, 
        group.by = c("epiSubtypes"))
dev.off()

saveRDS(epi,'BCBM_epi7samples_dimRedu.rds')
#40333 features across 4942 cells; 7 samples
write.table(allmarkers,'allmarkers_reso0.8.txt',sep = '\t',quote = F,row.names = F)
