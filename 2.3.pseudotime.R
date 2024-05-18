setwd('BCBM/2.4.monocle')
# library packages ====
library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(data.table)
#BiocManager::install('monocle')
library(monocle)
library(data.table)

epi.all <- readRDS('../2.2.tumor_reduceDim/BCBM_epi7samples_dimRedu.rds')

epi <- subset(epi.all, tumorType == 'LN', invert = T)
saveRDS(epi, 'epi.RDS')

#(1) count表达矩阵 ====
expr_matrix = GetAssayData(epi, layer = "data")
expr_matrix[1:3,1:3]

#(2) cell meta注释信息 ====
p_data <- epi@meta.data 
head(p_data)
pd <- new('AnnotatedDataFrame', data = p_data) 

#(3) gene meta注释信息 ====
f_data <- data.frame(gene_short_name = row.names(epi),
                     row.names = row.names(epi))
head(f_data)
fd <- new('AnnotatedDataFrame', data = f_data)

#构建cds对象 ====
cds_pre <- newCellDataSet(expr_matrix,
                          phenoData = pd,
                          featureData = fd,
                          expressionFamily = negbinomial.size())
##预处理
# Add Size_Factor文库因子 ====
cds_pre <- estimateSizeFactors(cds_pre)
cds_pre$Size_Factor %>% head()
# [1] 0.5018784 0.6026831 0.4911545 0.6798951 0.7013429 0.9244000

# 计算基因表达量的离散度 ====
cds_pre <- estimateDispersions(cds_pre)

## variable(high dispersion) gene by Seurat ====
gene_sle <- VariableFeatures(epi)

### 标记所选择的基因 ====
cds <- setOrderingFilter(cds_pre, gene_sle)

#降维(关键步骤)
cds <- reduceDimension(cds, method = 'DDRTree')

#排序,得到轨迹分化相关的若干State
cds <- orderCells(cds)

### plot ====
pdf('monocle_tumorType_variable.pdf',height = 4,width = 10,onefile = F)
plot_cell_trajectory(cds, color_by = "Pseudotime") +
  facet_wrap("~tumorType", nrow = 1)
dev.off()

pData(cds) %>% colnames()

# plot ====

pdf('monocle_State_variable.pdf',height = 4,width = 4,onefile = F)
plot_cell_trajectory(cds, color_by = "State")
dev.off()

pdf('monocle_epiSubtype_variable.pdf',height = 4,width = 4.2,onefile = F)
plot_cell_trajectory(cds, color_by = "tumorType")
dev.off()

pdf('monocle_epiSubtype_pseudotime_variable.pdf',height = 4,width = 6,onefile = F)
plot_cell_trajectory(cds, color_by = "Pseudotime") +
  facet_wrap("~tumorType", nrow = 1)
dev.off()

#pdf('monocle_tumorType.pdf',height = 5,width = 6.5,onefile = F)
plot_cell_trajectory(cds, color_by = "tumorType")
#dev.off()

pdf('monocle_clusters.pdf',#height = 5,width = 5.5,
    onefile = F)
plot_cell_trajectory(cds, color_by = "seurat_clusters") +
  facet_wrap("~seurat_clusters", nrow = 3)
dev.off()

pdf('monocle_clusters_BoM.pdf',#height = 5,width = 5.5,
    onefile = F)
plot_cell_trajectory(cds[,which(cds$epiSubtypes == 'BoM')],
                     color_by = "seurat_clusters") +
  facet_wrap("~seurat_clusters", nrow = 3)
dev.off()

#pdf('monocle_Pseudotime.pdf',height = 5,width = 6.5,onefile = F)
plot_cell_trajectory(cds, color_by = "Pseudotime")
#dev.off()

# 鉴定轨迹分化相关基因 =====
diff_pseudo <- differentialGeneTest(cds[gene_sle,], cores = 1, 
                                    fullModelFormulaStr = "~sm.ns(Pseudotime)")
head(diff_pseudo)
table(diff_pseudo$qval<0.05)
#FALSE  TRUE 
#3481  1519

## 选取其中最显著的进行可视化 ====
diff_pseudo_gene <- diff_pseudo %>% 
  dplyr::arrange(qval) %>% 
  rownames() %>% head(50)

p <- plot_pseudotime_heatmap(cds[diff_pseudo_gene,], 
                             num_clusters = 7,  # default 6
                             return_heatmap=T)

pdf('pseudotime_heatmap.pdf',width = 3,height = 4,onefile = F)
p
dev.off()

## 获得具体每个cluster的组成基因 ====
pseudotime_clusters <- cutree(p$tree_row, k = 7) %>% 
  data.frame(gene = names(.), cluster = . )
head(pseudotime_clusters)
table(pseudotime_clusters$cluster)
# 1  2  3  4  5  6  7 
#15  6  6  5  2 15  1

## 分支点基因变化情况 ====
epi <- AddMetaData(epi,cds$State,col.name = 'pseudotimeState')
Idents(epi) <- epi$pseudotimeState

allMarkers <- FindAllMarkers(object = epi, assay = NULL,
                             only.pos= TRUE,min.pct = .25, logfc.threshold = 0.25) %>% 
  group_by(cluster)

table(allMarkers$cluster)
#  1    2    3 
# 769 1323  779 

allMarkers.sig <- subset(allMarkers, avg_log2FC > 1.1 & p_val_adj < 0.05)

saveRDS(allMarkers.sig,'allMarkersSig.RDS')
saveRDS(allMarkers,'allMarkers.RDS')
saveRDS(epi,'epi_withPseudotime.RDS')
