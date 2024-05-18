setwd('/Users/icu/Desktop/BCBM/1.anno')
# library packages ====
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)

# load data ====
immune <- readRDS('BCBM_immuneCells.rds')
#DefaultAssay(immune) <- "integrated"

# normalize ====
immune <- NormalizeData(object = immune, assay = 'RNA',
                        normalization.method = "LogNormalize",
                        scale.factor = 10000,
                        margin = 1, #If performing CLR normalization, normalize across features (1) or cells (2)
                        verbose = TRUE) # 
## find variable features ====
immune <- FindVariableFeatures(
  object = immune,
  assay = 'RNA',
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = "auto",
  num.bin = 20,
  binning.method = "equal_width",
  nfeatures = 5000,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  verbose = TRUE
) 

## 线性变化 ====
immune <- ScaleData(immune, features = rownames(immune))

## PCA降维 ====
immune <- RunPCA(immune, features = VariableFeatures(object = immune))

## 确定数据集维度
immune <- JackStraw(immune, num.replicate = 100)
immune <- ScoreJackStraw(immune, dims = 1:20)
ElbowPlot(immune) # 20

immune <- FindClusters(immune, resolution = 2.2)
# Look at cluster IDs of the first 5 cells
head(Idents(immune), 5)
table(immune$seurat_clusters)

#pdf('immuneCell_UMAP.pdf',width = 10,height = 5,onefile = F)
#DimPlot(immune, reduction = "umap", label = TRUE, 
#        group.by = c("Patient",'seurat_clusters'))
#dev.off()
#pdf('immuneCell_tSNE.pdf',width = 10,height = 5,onefile = F)
#DimPlot(immune, reduction = "tsne", label = TRUE, 
#        group.by = c("Patient",'seurat_clusters'))
#dev.off()

saveRDS(immune, 'BCBM_immuneCells_reso22.rds')

# Immune cells annotation ===
allMarkers.list <- list(
  'Cell Type' = c('PTPRC','CD3D','CD3E', #T cell
                  'CD8A','CD8B', #CD8+ T
                  'CD4', #CD4+ T
                  'GNLY','KLRD1','KLRK1', #NK
                  'CD79A','CD79B','MS4A1', #B cell
                  'CD14','CD16','CD11b','CD11B', #Monocytes
                  'CD68','CD163','CD206', #Macrophages
                  'CD15','CD66b','CD66B', #Neutrophils
                  'CD123','CD203c','CD203C', #Basophils
                  'CD11c','CD11C','CD141','BDCA3', #Dendritic cells
                  'CD117', #Mast Cells
                  'LYZ','CD68','ITGAX') %>% #Myeloid
    unique() %>% .[which(. %in% immune@assays$RNA@data@Dimnames[[1]])],
  'Treg' = c('IKZF2','IL2RA','FOXP3') %>% 
    unique() %>% .[which(. %in% immune@assays$RNA@data@Dimnames[[1]])],
  'Naive' = c('CCR7','LEF1','SELL','TCF7') %>% 
    unique() %>% .[which(. %in% immune@assays$RNA@data@Dimnames[[1]])],
  'Cytotoxic' = c('IFNG','GZMA','GZMB','GZMK','PRF1','NKG7') %>% 
    unique() %>% .[which(. %in% immune@assays$RNA@data@Dimnames[[1]])],
  'Co-stimulatory' = c('TNFRSF9','TNFRSF14','ICOS','CD28') %>% 
    unique() %>% .[which(. %in% immune@assays$RNA@data@Dimnames[[1]])]
)
## plot
library(ggplot2)
library(RColorBrewer)

## dot plot
#dir.create('../3.1.immune_anno')
#setwd('../3.1.immune_anno')
pdf('immuneCellMarkers.pdf',width = 16,height = 9,onefile = F)
DotPlot(immune, assay = 'RNA',
        features = allMarkers.list,group.by = "seurat_clusters",
        cols = c("#f5af19", "#f12711"))+
  RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  theme(legend.position="bottom")
dev.off()

### violin plot
names(allMarkers.list)
lapply(names(allMarkers.list), function(cellType){
  pdf(paste0('immune_',cellType,'_vlnPlot.pdf'),width = 15,height = 15,onefile = F)
  print(
    VlnPlot(immune, assay = 'RNA',
            features = allMarkers.list[[cellType]], pt.size = 0)
  )
  dev.off()
})

# find markers ====
allmarkers <- FindAllMarkers(
  object = immune, assay = 'integrated',
  only.pos= TRUE,min.pct = .25, logfc.threshold = 0.25) %>% 
  group_by(cluster)

## add immune cell types metaData ====
immuneTypes <- ifelse(immune$seurat_clusters %in% c(0,6,10,19,27,25,7,37),'Cytotoix NK-T',
                      ifelse(immune$seurat_clusters %in% c(1:4,8,9,13,16:18,20:22,29,30,33),
                             'Naive T',
                             ifelse(immune$seurat_clusters %in% c(5,11,15,23,24,31,32,35,36),
                                    'B',
                                    ifelse(immune$seurat_clusters %in% c(12,28),
                                           'Regulatory T',
                                           ifelse(immune$seurat_clusters %in% c(14),
                                                  'Myeloid',
                                                         ifelse(immune$seurat_clusters %in% c(26),
                                                                'TNFAIP2-MAFB','UBE2C-DLGAP5'))
                                    )
                             ))
) %>% factor(levels = c('UBE2C-DLGAP5','TNFAIP2-MAFB',
                        'Myeloid','B','Cytotoix NK-T','Regulatory T','Naive T'))

immune <- AddMetaData(immune,immuneTypes,col.name = 'immuneTypes')
table(immune$immuneTypes)
table(immune$immuneTypes,immune$tumorType)

## plot

pdf('immuneCellMarkers_mergeClusters.pdf',width = 12,height = 3.5,onefile = F)
DotPlot(immune, assay = 'RNA',
        features = allMarkers.list,group.by = "immuneTypes",
        cols = c("#f5af19", "#f12711"))+
  #coord_flip()+
  RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  theme(legend.position="bottom")
dev.off()

pdf('immuneType_UMAP.pdf',width = 17,height = 5,onefile = F)
DimPlot(immune, reduction = "umap", label = T, 
        group.by = c("tumorType", "immuneTypes",'seurat_clusters'))
dev.off()
pdf('immuneType_tSNE.pdf',width = 17,height = 5,onefile = F)
DimPlot(immune, reduction = "tsne", label = T, 
        group.by = c("tumorType", "immuneTypes",'seurat_clusters'))
dev.off()

saveRDS(immune, 'BCBM_7samples_immuneTypes_reso22.rds')

# plot cell ratio ====
library(plyr)
library(ggplot2)

plotData <- table(immune$tumorType,immune$immuneTypes) %>% data.frame()# %>% 
#reshape2::dcast(Var1~Var2,value.var = 'Freq')
colnames(plotData) <- c('Patient','Cluster','Freq')
plotData$Patient <- factor(plotData$Patient, levels = c('PT','LN','BoM'))
plotData$Cluster <- factor(plotData$Cluster, levels = c('Naive T','Regulatory T',
                                                        'Cytotoix NK-T','B','Myeloid',
                                                        'UBE2C-DLGAP5','TNFAIP2-MAFB'))

plotData <- ddply(plotData, "Cluster", transform,
                  percent = Freq / sum(Freq) * 100)

pdf('tumorCell_Patient_Percent.pdf',
    width = 9,height = 5,onefile = F)
print(ggplot(plotData, aes(x=Cluster, y=percent, fill=Patient)) +
        geom_bar(stat="identity") +
        scale_fill_manual(values=c('#a8e6cf','#fdffab','#ffd3b6','#ffaaa5',
                                   '#11999e','#3f72af','#ffc7c7','#9896f1',
                                   '#ff7e67','#0dceda','#6ef3d6','#a6d0e4',
                                   '#fbafaf','#f2c6b4','#f3e8cb','#99e1e5',
                                   '#ffc93c','#ff9a3c','#ff6f3c','#155263',
                                   '#c86b85','#e6a4b4','#f3d7ca','#bfcfff',
                                   '#a5dee5','#e0f9b5','#ffcfdf','#ffc7c7'))+
        
        theme_bw()+
        labs(x='Cluster',y='Percent',fill='Patient'))
dev.off() 

plotData <- ddply(plotData, "Patient", transform,
                  percent = Freq / sum(Freq) * 100)

pdf('tumorCell_Cluster_Percent.pdf'
    ,width = 7,height = 5,onefile = F)
print(ggplot(plotData, aes(x=Patient, y=percent, fill=Cluster)) +
        geom_bar(stat="identity") +
        scale_fill_manual(values=c('#a8e6cf','#fdffab','#ffd3b6','#ffaaa5',
                                   '#11999e','#3f72af','#ffc7c7','#9896f1',
                                   '#ff7e67','#0dceda','#6ef3d6','#a6d0e4',
                                   '#fbafaf','#f2c6b4','#f3e8cb','#99e1e5',
                                   '#ffc93c','#ff9a3c','#ff6f3c','#155263',
                                   '#c86b85','#e6a4b4','#f3d7ca','#bfcfff',
                                   '#a5dee5','#e0f9b5','#ffcfdf','#ffc7c7'))+
        
        theme_bw()+
        labs(x='Patient',y='Percent',fill='Cluster'))
dev.off() 
