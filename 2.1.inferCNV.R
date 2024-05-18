setwd('BCBM/2.1.epi_inferCNV')
# library packages ====
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)

# load data ====
epi <- readRDS('../1.anno/BCBM_epiCells.rds')
fibro <- readRDS('../1.anno/BCBM_fibroCells.rds')

# inferCNV ====
#install.packages("rjags")
#BiocManager::install("infercnv")
library(infercnv)
library(rtracklayer)

## gene reference file
geneRef <- rtracklayer::import('gencode.v44.annotation.gtf') %>% as.data.frame() %>% 
  subset(type == 'gene' & #seqnames != 'chrX' & 
           seqnames != 'chrY' & seqnames != 'chrM' &
           gene_type %in% c('lncRNA','miRNA','protein_coding'))
table(duplicated(geneRef$gene_name))

geneRef2 <- geneRef[,c('gene_name','seqnames','start','end','width','gene_type')] %>% 
  subset(gene_name %in% rownames(epi))
table(duplicated(geneRef2$gene_name))

### delete duplicated genes
dup_genes <- geneRef2$gene_name[duplicated(geneRef2$gene_name)] %>% unique()

library(data.table)
addgenes <- lapply(dup_genes, function(x){
  dup_rows <- geneRef2[which(geneRef2$gene_name %in% x),] %>% arrange(desc(width))
  each <- dup_rows[1,]
  return(each)
}) %>% rbindlist()

geneRef3 <- geneRef2[-which(geneRef2$gene_name %in% dup_genes),]
geneRef3 <- rbind(geneRef3, addgenes) %>% .[,1:4] %>% arrange(seqnames,start)
table(duplicated(geneRef3$gene_name))

rownames(geneRef3) <- geneRef3$gene_name

write.table(geneRef3,'geneOrderingFile.txt',quote = F,sep = '\t',row.names = T)

## sample annotation file
sampleAnno <- data.frame(row.names = c(colnames(epi),colnames(fibro)),
                         patient = c(epi$Patient,fibro$Patient),
                         group = c(rep('carcinoma',ncol(epi)),rep('fibro',ncol(fibro)))
                         #group = c(epi$seurat_clusters,rep('fibro',ncol(fibro)))
                         )

write.table(sampleAnno,'cellAnnotations.txt',quote = F,sep = '\t',row.names = T)

sampleAnno <- data.frame(
  row.names = rownames(sampleAnno),
  group = sampleAnno$group
)

## raw count data
rawCount1 <- epi@assays$RNA@counts %>% .[geneRef3$gene_name,] %>% as.matrix()
rawCount2 <- fibro@assays$RNA@counts %>% .[geneRef3$gene_name,] %>% as.matrix()
rawCount <- cbind(rawCount1,rawCount2)
dim(rawCount)# 4871 7292
dim(rawCount1)# 4871 4942
dim(rawCount2)# 4871 2350
rm(rawCount1,rawCount2)

## infer CNV
geneRef3 <- read.table('geneOrderingFile.txt',sep = '\t',row.names = 1)
geneRef3 <- geneRef3[,-1]
geneRef3$start <- as.numeric(geneRef3$start)
geneRef3$end <- as.numeric(geneRef3$end)

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=rawCount, 
                                    annotations_file=sampleAnno,
                                    #delim="\t",
                                    gene_order_file=geneRef3,
                                    ref_group_names="fibro")

# perform infercnv operations to reveal cnv signal
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,  
                              out_dir="fibro_output",  
                              cluster_by_groups=T,   
                              denoise=T, 
                              HMM=T, 
                              HMM_type = "i3")
