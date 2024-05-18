setwd('/Users/icu/Desktop/BCBM/4.State1_epi_immune_cellChat')
# library packages ====
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)

# load data ====
epi <- readRDS('../2.4.monocle/epi_withPseudotime.RDS')
epi <- subset(epi, pseudotimeState == 1)
immune <- readRDS('../3.1.immune_anno/BCBM_7samples_immuneTypes_reso22.rds')
immune <- subset(immune, immuneTypes %in% 
                   c('Naive T','B','Regulatory T','Cytotoix NK-T','Myeloid'))

# cell chat ====
#devtools::install_github("sqjin/CellChat")
library(CellChat)
library(patchwork)
library(ggplot2)
library(ggalluvial)
library(svglite)
options(stringsAsFactors = FALSE)

## Import ligand-receptor interaction database ====
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

lapply(unique(epi$tumorType),function(p){
  epiData <- subset(epi, tumorType == p)
  immuneData <- subset(immune, tumorType == p)
  
  # Part I: Data input & processing and initialization of CellChat object ====
  ## prepare inputs ====
  data.input1 <- epiData@assays$RNA@data
  data.input2 <- immuneData@assays$RNA@data
  data.input <- cbind(data.input1,data.input2)
  
  identity1 <- data.frame(group = rep('Epithelial',ncol(epiData)),
                          type = rep('Carcinoma',ncol(epiData)),
                          row.names = names(epiData$tumorType))
  identity2 <- data.frame(group = 'Immune',
                          type = immuneData$immuneTypes,
                          row.names = names(immuneData$immuneTypes))
  identity <- rbind(identity1,identity2) # create a dataframe consisting of the cell labels
  table(identity$group,identity$type) # check the cell labels
  
  rm(data.input1,data.input2,identity1,identity2)
  ## create cell chat obj ====
  cellchat <- createCellChat(object = data.input, meta = identity, group.by = "type")
  
  ### Add cell information into *meta* slot of the object  (Optional)
  #cellchat <- addMeta(cellchat, meta = identity, meta.name = 'labels')
  cellchat <- setIdent(cellchat, ident.use = "type") # set "labels" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
  
  # set the used database in the object
  cellchat@DB <- CellChatDB.use
  
  ## Preprocessing the expression data ====
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
  # cellchat <- projectData(cellchat, PPI.human)
  
  # Part II: Inference of cell-cell communication network ====
  ## Compute the communication probability and infer cellular communication network ====
  cellchat <- computeCommunProb(cellchat)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 3)
  
  ## Infer the cell-cell communication at a signaling pathway level ====
  cellchat <- computeCommunProbPathway(cellchat)
  
  ## Calculate the aggregated cell-cell communication network ====
  cellchat <- aggregateNet(cellchat)
  
  saveRDS(cellchat,paste0('tumorType_',p,'_cellChat.RDS'))
  
})

# combind all cell chat res ====
grep('_cellChat.RDS',dir(),value = T)
object.list <- list(
  PT = readRDS(grep('_cellChat.RDS',dir(),value = T)[2]),
  #LN = readRDS(grep('_cellChat.RDS',dir(),value = T)[2]),
  BoM = readRDS(grep('_cellChat.RDS',dir(),value = T)[1])
)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# Visualization ====
## Compare the total number of interactions and interaction strength ====
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1:2)#,group.levels = c('PT','LN','BoM')
)
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1:2),#group.levels = c('PT','LN','BoM'),
                           measure = 'weight'
)
gg1 + gg2

## Compare the overall information flow of each signaling pathway ====
gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2),
               color.use = c('PT'='#dcedc2','LN'='#ffd3b5','BoM'='#ffaaa6'),
               stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2),
               color.use = c('PT'='#dcedc2','LN'='#ffd3b5','BoM'='#ffaaa6'),
               stacked = F, do.stat = TRUE)
pdf('informationFlow.pdf',width = 7,height = 6,onefile = F)
gg1 + gg2
dev.off()

## net up & down ====
## Chord diagram (cell) ====
pos.dataset = 'BoM'
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets",
                                       pos.dataset = pos.dataset,
                                       features.name = features.name,
                                       only.pos = FALSE,
                                       thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = .05)

#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = pos.dataset,
                              ligand.logFC = 0.2, receptor.logFC = NULL)
pdf(paste0("Up-regulated signaling pathways - ",
           names(object.list)[2],'.pdf'),width = 6,height = 6)
netVisual_chord_cell(object.list[[2]],signaling = NULL,
                     net = net.up, 
                     slot.name = 'netP',
                     color.use = c('Carcinoma'='#ff7e67','Cytotoix NK-T'='#fdffab',
                                   'B'='#ffd3b6','Naive T'='#ffaaa5',
                                   'Regulatory T'='#11999e','Myeloid'='#e0f9b5',
                                   'TNFAIP2-MAFB'='#ffc7c7','UBE2C-DLGAP5'='#a6d0e4'),
                     title.name = paste0("Up-regulated signaling pathways - ",
                                         names(object.list)[2]), legend.pos.x = 5)
dev.off()

pdf(paste0("Up-regulated signaling pathways from carcinoma - ",
           names(object.list)[2],'.pdf'),width = 10,height = 10)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP',
                     net = net.up,
                     sources.use = 'Carcinoma',
                     #targets.use = 'Carcinoma',
                     #remove.isolate = T,
                     show.legend = T,
                     lab.cex = 0.5,
                     #small.gap=0, big.gap=10,
                     color.use = c('Carcinoma'='#ff7e67','Cytotoix NK-T'='#fdffab',
                                   'B'='#ffd3b6','Naive T'='#ffaaa5',
                                   'Regulatory T'='#11999e','Myeloid'='#e0f9b5'),
                     title.name = paste0("Up-regulated signaling pathways from carcinoma - ",
                                         names(object.list)[2]), legend.pos.x = 5)
dev.off()

pdf(paste0("Up-regulated signaling pathways target carcinoma - ",
           names(object.list)[2],'.pdf'),width = 10,height = 10)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP',
                     net = net.up,
                     #sources.use = 'Carcinoma',
                     targets.use = 'Carcinoma',
                     #remove.isolate = T,
                     show.legend = T,
                     lab.cex = 0.5,
                     #small.gap=0, big.gap=10,
                     color.use = c('Carcinoma'='#ff7e67','Cytotoix NK-T'='#fdffab',
                                   'B'='#ffd3b6','Naive T'='#ffaaa5',
                                   'Regulatory T'='#11999e','Myeloid'='#e0f9b5'),
                     title.name = paste0("Up-regulated signaling pathways target carcinoma - ",
                                         names(object.list)[2]), legend.pos.x = 5)
dev.off()

intersect(object.list[[1]]@netP$pathways,object.list[[2]]@netP$pathways) %>% 
  intersect(object.list[[2]]@netP$pathways)

union(object.list[[1]]@netP$pathways,object.list[[2]]@netP$pathways) %>% 
  union(object.list[[2]]@netP$pathways)

pdf('netVisual_bubble_MHC-I_CD99.pdf',width = 6,height = 4,onefile = F)
netVisual_bubble(cellchat,
                 #sources.use = c(1:6),
                 #targets.use = c(2:4),
                 signaling = c("MHC-I","CD99"),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('PT'='#66c6ba', 'BoM'='#e67a7a'),
                 grid.on = F,
                 comparison = c(1:2),angle.x = 90, remove.isolate = T, return.data = F)
#> Comparing communications on a merged object
dev.off()

## specific in BoM (bubble plot)
pdf('netVisual_bubble_BoM_fromCarcinoma.pdf',width = 3.6,height = 4,onefile = F)
netVisual_bubble(object.list[[2]],
                 sources.use = 'Carcinoma',
                 #targets.use = 'Carcinoma',
                 signaling = c('MHC-II','FN1','GALECTIN','MK'),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('PT'='#66c6ba','BoM'='#e67a7a'),
                 grid.on = F,
                 angle.x = 90, remove.isolate = T, return.data = F)
dev.off()

pdf('netVisual_bubble_BoM_targetCarcinoma.pdf',width = 3.6,height = 4,onefile = F)
netVisual_bubble(object.list[[2]],
                 #sources.use = 'Carcinoma',
                 targets.use = 'Carcinoma',
                 signaling = c('SPP1','FN1','GALECTIN','GRN'),
                 sort.by.target = T,
                 #color.heatmap = "Spectral",
                 color.text = c('PT'='#66c6ba','BoM'='#e67a7a'),
                 grid.on = F,
                 angle.x = 90, remove.isolate = T, return.data = F)
dev.off()

## specific in BoM (chord)
interestPW <- c('MHC-II','FN1','GALECTIN','MK')
lapply(interestPW, function(pw){
  pdf(paste0('fromCarcinoma_',pw,'.pdf'),width = 4,height = 4,onefile = F)
  print(
    netVisual_chord_gene(object.list[[2]], slot.name = 'net',
                         net = net.up,
                         signaling = pw,
                         sources.use = 'Carcinoma',
                         #targets.use = 'Carcinoma',
                         #remove.isolate = T,
                         show.legend = F,
                         lab.cex = 0.5,
                         #small.gap=0, big.gap=10,
                         color.use = c('Carcinoma'='#ff7e67','Cytotoix NK-T'='#fdffab',
                                       'B'='#ffd3b6','Naive T'='#ffaaa5',
                                       'Regulatory T'='#11999e','Myeloid'='#e0f9b5'),
                         title.name = paste0(names(object.list)[2],' - ',pw),
                         legend.pos.x = 5)
  )
  dev.off()
})

interestPW <- c('SPP1','FN1','GALECTIN','GRN','CEACAM')
lapply(interestPW, function(pw){
  pdf(paste0('targetCarcinoma_',pw,'.pdf'),width = 4,height = 4,onefile = F)
  print(
    netVisual_chord_gene(object.list[[2]], slot.name = 'net',
                         net = net.up,
                         signaling = pw,
                         #sources.use = 'Carcinoma',
                         targets.use = 'Carcinoma',
                         #remove.isolate = T,
                         show.legend = F,
                         lab.cex = 0.5,
                         #small.gap=0, big.gap=10,
                         color.use = c('Carcinoma'='#ff7e67','Cytotoix NK-T'='#fdffab',
                                       'B'='#ffd3b6','Naive T'='#ffaaa5',
                                       'Regulatory T'='#11999e','Myeloid'='#e0f9b5'),
                         title.name = paste0(names(object.list)[2],' - ',pw),
                         legend.pos.x = 5)
  )
  dev.off()
})
