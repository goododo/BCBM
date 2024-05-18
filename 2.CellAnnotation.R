setwd('BCBM/1.anno')
# library packages ====
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)

# load data ====
bcbm <- readRDS('../0.rawdata/BCBM_inteData7samples_dimRedu.rds')
#DefaultAssay(bcbm) <- "integrated"

allMarkers.list <- list(
  ImmuneCells = c('PTPRC','CD3E','CD3D','CD2','TRAC','TRBC2', #T
                  'MS4A1','CD79A','CD79B','CD19','CD74', #B
                  'CD68','FCGR3A', #NK
                  'CD14','LYZ') %>% unique() %>% #myeloid
    .[which(. %in% bcbm@assays$RNA@data@Dimnames[[1]])],
  EpithelialCells = c('KRT18','EPCAM','KRT19','CD24','CLDN4','MUC1',
                      'KRT8','KRT7',
                      'CALD1','AGRN','CTNNB1') %>% unique() %>% 
    .[which(. %in% bcbm@assays$RNA@data@Dimnames[[1]])],
  EndothelialCells = c('PECAM1','VWF','KDR','EFNB2',
                       'ENG','PLVAP','CDH5','RAMP2') %>% unique() %>% 
    .[which(. %in% bcbm@assays$RNA@data@Dimnames[[1]])],
  Fibroblasts = c('PDGFRB','ACTA2','THY1',
                  'S100A4',
                  'COL1A1','FN1','CD90','POSTN','FAP','DCN','PDGFRA','LUM','RUNX2','COL3A1') %>%
    unique() %>% .[which(. %in% bcbm@assays$RNA@data@Dimnames[[1]])]
)
## plot
library(ggplot2)
library(RColorBrewer)

### dot plot
pdf('commonCellMarkers.pdf',width = 13,height = 5,onefile = F)
DotPlot(bcbm, assay = 'RNA',
        features = allMarkers.list,group.by = "seurat_clusters",
        cols = c("#f8fcfb", "#2193b0"),
        col.min = 0)+
  RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
dev.off()

### violin plot
names(allMarkers.list)
lapply(names(allMarkers.list), function(cellType){
  pdf(paste0(cellType,'_vlnPlot.pdf'),width = 15,height = 15,onefile = F)
  print(
    VlnPlot(bcbm, assay = 'RNA',
            features = allMarkers.list[[cellType]], pt.size = 0)
  )
  dev.off()
})

## add common cell types metaData ====
commonTypes <- ifelse(bcbm$seurat_clusters %in% c(0:3,7,9:11),'ImmuneCell',
                      ifelse(bcbm$seurat_clusters %in% c(4,8,12),'Epithelial',
                             ifelse(bcbm$seurat_clusters %in% c(5),'Endothelial',
                                    ifelse(bcbm$seurat_clusters %in% c(6,13),
                                           'Fibroblast','Unknown'))
                      )) %>% factor(levels = c('ImmuneCell','Epithelial','Endothelial','Fibroblast'))

bcbm <- AddMetaData(bcbm,commonTypes,col.name = 'commonCellTypes')
table(bcbm$commonCellTypes)
table(bcbm$commonCellTypes,bcbm$Patient)

## plot
pdf('commonCellType_UMAP.pdf',width = 20,height = 5.5,onefile = F)
DimPlot(bcbm, reduction = "umap", label = TRUE, 
        group.by = c("Patient", "commonCellTypes",'seurat_clusters'))
dev.off()
pdf('commonCellType_tSNE.pdf',width = 20,height = 5.5,onefile = F)
DimPlot(bcbm, reduction = "tsne", label = TRUE, 
        group.by = c("Patient", "commonCellTypes",'seurat_clusters'))
dev.off()

## markers
allMarkers <- c('PTPRC','CD3E','CD3D','CD2','TRAC','TRBC2', #T
                'MS4A1','CD79A','CD79B','CD19','CD74', #B
                'CD68','FCGR3A', #NK
                'CD14','LYZ',
                'KRT18','EPCAM','KRT19','CD24','CLDN4','MUC1',
                'KRT8','KRT7',
                'CALD1','AGRN','CTNNB1',
                'PECAM1','VWF','KDR','EFNB2',
                'ENG','PLVAP','CDH5','RAMP2',
                'PDGFRB','ACTA2','THY1',
                'S100A4',
                'COL1A1','FN1','CD90','POSTN','FAP','DCN','PDGFRA',
                'LUM','RUNX2','COL3A1') %>% unique() %>% 
    .[which(. %in% bcbm@assays$RNA@data@Dimnames[[1]])]

pdf('commonCellMarkers_mergeClusters.pdf',width = 5,height = 10,onefile = F)
DotPlot(bcbm, assay = 'RNA',
        features = allMarkers,group.by = "commonCellTypes",
        cols = c("#f5af19", "#f12711"),
        #col.min = -1
        )+
  RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  coord_flip()+
  theme(legend.position="right")
dev.off()

saveRDS(bcbm, 'BCBM_7samples_commonCellTypes.rds')

## plot cell type percent of each patient ====
library(plyr)
library(ggplot2)

plotData <- table(bcbm$tumorType,bcbm$commonCellTypes) %>% data.frame()# %>% 
#reshape2::dcast(Var1~Var2,value.var = 'Freq')
colnames(plotData) <- c('Patient','Cluster','Freq')
plotData$Patient <- factor(plotData$Patient, levels = c('PT','LN','BoM'))

plotData <- ddply(plotData, "Patient", transform,
                  percent = Freq / sum(Freq) * 100)

pdf('commonCell_Percent_Patient.pdf',width = 7,height = 5,onefile = F)
ggplot(plotData, aes(x=Patient, y=percent, fill=Cluster)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('#a8e6cf','#fdffab','#ffd3b6','#ffaaa5'))+
  
  theme_bw()+
  labs(x='Tumor Type',y='Percent',fill='Cell Type')
dev.off() 

# subclusters ====
epi <- subset(bcbm, commonCellTypes == 'Epithelial')
#40333 features across 4942 samples
saveRDS(epi, 'BCBM_epiCells.rds')

immune <- subset(bcbm, commonCellTypes == 'ImmuneCell')
#40333 features across 24795 samples
saveRDS(immune, 'BCBM_immuneCells.rds')

endo <- subset(bcbm, commonCellTypes == 'Endothelial')
#40333 features across 2288 samples
saveRDS(endo, 'BCBM_endoCells.rds')

fibro <- subset(bcbm, commonCellTypes == 'Fibroblast')
#40333 features across 2350 samples
saveRDS(fibro, 'BCBM_fibroCells.rds')
