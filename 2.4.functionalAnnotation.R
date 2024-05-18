setwd('/Users/icu/Desktop/BCBM/2.5.func_pseudotime')
# library packages ====
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(ggplot2)

epi <- readRDS('../2.3.monocle/epi.RDS')
allmarkers <- readRDS('../2.3.monocle/allMarkersSig.RDS')

# functional enrichment ====
library(clusterProfiler)
library(org.Hs.eg.db)
ids <- bitr(allmarkers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
allmarkers <- merge(allmarkers,ids,by.x='gene',by.y='SYMBOL')

genes <- allmarkers %>% subset(avg_log2FC > 1)
genes <- split(genes$ENTREZID,genes$cluster)

# KEGG ====
compare_kegg <- compareCluster(genes,
                                  fun = "enrichKEGG",
                                  organism = "hsa", pvalueCutoff = 0.05)

#pdf('cluster_KEGG.pdf',width = 6,height = 10,onefile = F)
dotplot(compare_kegg,
        showCategory = 6)+ 
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 0.5, hjust = 0.5
  ))
#dev.off()

compare_kegg.df <- compare_kegg@compareClusterResult %>%
  subset(p.adjust < 0.05)

library(RColorBrewer)
library(wesanderson)

pdf('cluster_KEGG.pdf',width = 6,height = 6,onefile = F)
ggplot(compare_kegg.df,aes(x=Cluster,y=Description,color=p.adjust,
                              size=Count))+
  geom_point()+
  #facet_grid(rows = NULL,cols = vars(Cluster),scales = "free",switch = 'x') +
  scale_color_continuous(low='#f12711',high='#f5af19')+
  
  theme_classic()+
  ylab(NULL)
dev.off()

# GO ====
compare_go <- compareCluster(genes,
                             fun = "enrichGO",
                             OrgDb = "org.Hs.eg.db",
                             ont = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.01,
                             qvalueCutoff = 0.05)

compare_go.df <- compare_go@compareClusterResult %>%
  subset(p.adjust < 0.05)

library(RColorBrewer)
library(wesanderson)

pdf('cluster_GO.pdf',width = 6,height = 6,onefile = F)
ggplot(compare_go.df,aes(x=Cluster,y=Description,color=p.adjust,
                            size=Count))+
  geom_point()+
  #facet_grid(rows = NULL,cols = vars(Cluster),scales = "free",switch = 'x') +
  scale_color_continuous(low='#654ea3',high='#eaafc8')+
  
  theme_classic()+
  ylab(NULL)
dev.off()

# hallmarks AUC ====
library(UCell)
library(irGSEA)

epi <- irGSEA.score(object = epi, assay = "RNA",
                    slot = "data", seeds = 123, ncores = 1,
                    min.cells = 3, min.feature = 0,
                    custom = F, geneset = NULL, msigdb = T,
                    species = "Homo sapiens", category = "H",  
                    subcategory = NULL, geneid = "symbol",
                    method = c("AUCell", "UCell", "singscore",
                               "ssgsea"),
                    aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                    kcdf = 'Gaussian')
Seurat::Assays(epi)

saveRDS(epi,'epi_irGSEAscore.RDS')

## 整合差异基因集
result.dge <- irGSEA.integrate(object = epi,
                               group.by = "pseudotimeState",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea"))

#pdf('irGSEA_heatmap.pdf',width = 10,height = 7,onefile = F)
irGSEA.heatmap(object = result.dge,
               method = "RRA",
               top = 50,
               show.geneset = NULL)
#dev.off()

#library(pheatmap)
library(reshape2)
library(reshape)
library(hrbrthemes)
plotData <- subset(result.dge[["RRA"]],pvalue<0.05)# %>% .[,-5] %>% cast(Name~cluster)
colnames(plotData) <- c('Name','pValue','Cluster','Direction','Method')
plotData$Name <- substr(plotData$Name,10,nchar(plotData$Name))
plotData$Name <- gsub('-',' ',plotData$Name)

plotData$Cluster <- factor(plotData$Cluster, levels = c('PT','LN','BoM'))

library(Hmisc)
plotData$Direction <- capitalize(plotData$Direction)

pdf('hallmarks_heatmap.pdf',width = 4,height = 4,onefile = F)
ggplot(plotData,aes(x=Cluster,y=Name,fill=Direction))+
  geom_tile() +
  scale_fill_manual(values=c('Up'="#f6416c", 'Down'="#00b8a9")) +
  
  theme_classic()+
  xlab(NULL)+ylab(NULL)+
  #coord_flip()+ 
  theme(legend.position = 'top'#,axis.text.x = element_text(angle = 70,vjust = .5)
  )
dev.off() 
