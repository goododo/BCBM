setwd('/Users/icu/Desktop/BCBM/2.5.clinical')
# library packages ====
library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(data.table)

# load data ====
## pseudotime state DEGs
epi <- readRDS('../2.4.monocle/epi_withPseudotime.RDS')

allMarkers <- readRDS('../2.4.monocle/allMarkers.RDS')
#  1    2    3 
# 769 1323  779 

allMarkers.sig <- subset(allMarkers, avg_log2FC > 1.5 & p_val_adj < 0.05)
table(allMarkers.sig$cluster)
#1   2   3 
#278 737 171

# decision tree ====
library(rpart)
#library(randomForest)
epiExp <- epi %>% GetAssayData() %>% as.matrix() %>% t() %>% as.data.frame()
epiExp <- cbind(
  class=epi$pseudotimeState[rownames(epiExp)],
  epiExp)
dt <- rpart(class~.,data = epiExp,method = 'class')

dt.list <- lapply(1:3, function(state){
  epiExp.dt <- epiExp
  epiExp.dt$class <- ifelse(epiExp$class == state,1,0)
  dt <- rpart(class~.,data = epiExp.dt,method = 'class')
  dt.prune <- prune(dt, cp=0.01)
  
  return(dt.prune)
})

## plot decision tree ====
plot(dt.list[[1]],compress=T,margin=0.2)
text(dt.list[[1]],cex=1.5)

library(rpart.plot)

layout(matrix(c(1, 2,  # First, second
                3), # and third plot
              nrow = 1,
              ncol = 3,
              byrow = TRUE))
prp(dt.list[[1]],type=4,extra=2,digits=3)
prp(dt.list[[2]],type=4,extra=2,digits=3)
prp(dt.list[[3]],type=4,extra=2,digits=3)

layout(matrix(c(1), nrow = 1, ncol = 1,byrow = TRUE))

## markers
dt.markers <- data.frame(
  gene=c(dt.list[[1]]$frame$var,dt.list[[2]]$frame$var,dt.list[[3]]$frame$var),
  state=c(rep(1,nrow(dt.list[[1]]$frame)),rep(2,nrow(dt.list[[2]]$frame)),
          rep(3,nrow(dt.list[[3]]$frame))) %>% as.factor(),
  complexity=c(dt.list[[1]]$frame$complexity,dt.list[[2]]$frame$complexity,
               dt.list[[3]]$frame$complexity)
) %>% subset(gene!='<leaf>')

saveRDS(dt.markers,'decisionTree_markers.RDS')
saveRDS(dt.list,'decisionTree_res.RDS')

# cBioPortal ====
clinical <- read.table('cBioPortal/combined_study_clinical_data.tsv',header = T,sep = '\t')
exp1 <- read.table(paste0('cBioPortal/',subset(dt.markers,state==1)$gene[1],'.txt'),
                   header = T,sep = '\t') %>% 
  subset(SAMPLE_ID %in% clinical$Sample.ID) %>% .[,-1] %>% 
  cbind(.,gene=subset(dt.markers,state==1)$gene[1])
exp2 <- read.table(paste0('cBioPortal/',subset(dt.markers,state==1)$gene[2],'.txt'),
                   header = T,sep = '\t') %>% 
  subset(SAMPLE_ID %in% clinical$Sample.ID) %>% .[,-1] %>% 
  cbind(.,gene=subset(dt.markers,state==1)$gene[2])
exp3 <- read.table(paste0('cBioPortal/',subset(dt.markers,state==1)$gene[3],'.txt'),
                   header = T,sep = '\t') %>% 
  subset(SAMPLE_ID %in% clinical$Sample.ID) %>% .[,-1] %>% 
  cbind(.,gene=subset(dt.markers,state==1)$gene[3])
colnames(exp1)[2] <- colnames(exp2)[2] <- colnames(exp3)[2] <- 'Exp'

library(reshape2)
expMarkers <- rbind(exp1,exp2,exp3)
expMarkers$Exp <- as.numeric(expMarkers$Exp)
expMarkers <-  dcast(expMarkers,gene~SAMPLE_ID,value.var = 'Exp')

rm(exp1,exp2,exp3,expMarkers,clinical)

# survival data ====
## TCGA BRCA (UCSC Xena)
probe <- read.table('gencode.v22.annotation.gene.probeMap',header = T)
colnames(probe)[1] <- 'Ensembl_ID'

counts <- read.table('TCGA-BRCA.htseq_counts.tsv',header = T,check.names = F)
counts <- merge(probe[,1:2],counts,by='Ensembl_ID')
exp <- counts[match(dt.markers$gene,counts$gene),-1] %>%
  group_by(gene) %>%
  summarise_all(mean) %>% as.data.frame() %>% na.omit()

rownames(exp) <- exp$gene
exp <- exp[,-1]

### survival data
surv <- read.table('TCGA-BRCA.survival.tsv',header = T) %>%
  subset(., OS.time >= 30) %>% .[grep('-11',.$sample,invert = T),]

### exp
exp <- exp[,match(surv$sample,colnames(exp)) %>% na.omit()]

### surv data
surv <- subset(surv, sample %in% colnames(exp))
surv_km <- data.frame(
  row.names = surv$sample,
  sample = surv$sample,
  event = surv$OS,
  time = surv$OS.time)
surv_10years <- surv_km
surv_10years$event[which(surv_10years$time > 365*10)] <- 0
surv_10years$time[which(surv_10years$time > 365*10)] <- 365*10

# set class ====
sampleType <- apply(exp, 2, function(x){
  scoreEach <- x[dt.markers[1:3,1]] * as.numeric(dt.markers[1:3,3])
  score <- sum(scoreEach)
})
sampleType <- ifelse(sampleType>=median(sampleType),"State1", "nonState1")

surv_km$Class <- sampleType
table(surv_km$Class)
surv_10years$Class <- sampleType
table(surv_10years$Class)

# survival K-M plot ====
library(survival)
library(survminer)

pairwise_survdiff(Surv(time, event) ~ Class,
                  data = surv_10years, p.adjust.method = 'none')

# 10 years
fit1 <- survfit(Surv(time, event) ~ Class, data= surv_10years)

pdf('10years_KMcurv.pdf',width = 6,height = 5,onefile = F)
ggsurvplot(fit1, data = surv_10years,
           pval=TRUE,
           pval.coord = c(0, 0.25),
           pval.size = 6,
           pval.method=TRUE,
           pval.method.coord = c(0, 0.35),
           conf.int=T, 
           
           risk.table=TRUE,
           tables.height = .3, 
           tables.col='strata',
           risk.table.pos='out',
           
           palette=c('#f6416c','#00b8a9'), 
           
           xscale='d_y',
           #title="Kaplan-Meier Curve for OS ", 
           xlab='Time(Year)',
           legend.labs=c("State1", "nonState1"), 
           legend.title="State",  
           surv.median.line = "hv" 
)
dev.off()

pairwise_survdiff(Surv(time, event) ~ Class,
                  data = surv_km, p.adjust.method = 'none')

# 总体生存曲线
fit1 <- survfit(Surv(time, event) ~ Class, data= surv_km)

pdf('all_KMcurv.pdf',width = 6,height = 5,onefile = F)
ggsurvplot(fit1, data = surv_km,
           pval=TRUE, 
           pval.coord = c(0, 0.25),
           pval.size = 6,
           pval.method=TRUE,
           pval.method.coord = c(0, 0.35),
           conf.int=T, 
           
           risk.table=TRUE,
           tables.height = .3, 
           tables.col='strata',
           risk.table.pos='out',
           
           palette=c('#f6416c','#00b8a9'), 
           
           xscale='d_y',
           #title="Kaplan-Meier Curve for OS ", 
           xlab='Time(Year)',
           legend.labs=c("State1", "nonState1"), 
           legend.title="State",  
           surv.median.line = "hv" 
)
dev.off()

# expression level ====
library(ggplot2)
library(ggsignif)
library(patchwork)

## plot data
plotData <- exp[dt.markers$gene,] %>% cbind(gene=rownames(.),.) %>% melt()
plotData$Class <- sampleType[plotData$variable] %>%
  factor(levels = c('nonState1','State1'))
colnames(plotData) <- c('Gene','Sample','Value','State')

plotData$geneType <- ifelse(plotData$Gene %in% dt.markers[1:3,]$gene,1,
                            ifelse(plotData$Gene %in% dt.markers[4:7,]$gene,2,3)) %>%
  factor(levels = c('Normal','nonState1','State1'))

## without Normal samples
pdf('expLevel_withoutNormal.pdf',width = 6,height = 3,onefile = F)
ggplot(subset(plotData,geneType==1),
       aes(x=State,y=Value,color=State))+
  geom_boxplot(linetype=1)+
  #geom_jitter()+
  scale_color_manual(values = c('Normal'='#ffd460', 'nonState1' = "#00b8a9",
                                'State1' = "#f6416c"))+
  
  geom_signif(comparisons = list(c('nonState1','State1')),
              test = 'wilcox.test', # wilcox.test, t.test
              step_increase = c(0,.05,.07),
              map_signif_level = F)+
  
  
  facet_wrap(~Gene)+
  
  theme_classic()#+
#theme(axis.text.y = element_blank(),axis.title.x = element_blank());p1
dev.off()

## with Normal samples
exp <- counts[match(dt.markers$gene,counts$gene),-1] %>%
  group_by(gene) %>%
  summarise_all(mean) %>% as.data.frame() %>% na.omit()

rownames(exp) <- exp$gene
exp <- exp[,-1]

exp <- cbind(exp[,grep('-11',colnames(exp))],
             exp[match(surv$sample,colnames(exp)) %>% na.omit()])

sampleType <- apply(exp, 2, function(x){
  scoreEach <- x[dt.markers[1:3,1]] * as.numeric(dt.markers[1:3,3])
  score <- sum(scoreEach)
})
sampleType <- ifelse(sampleType>=median(sampleType[grep('-11',colnames(exp),invert = T)]),
                     'State1','nonState1')
sampleType[grep('-11',colnames(exp))] <- 'Normal'
table(sampleType)

names(sampleType) <- colnames(exp)

plotData <- exp[dt.markers$gene,] %>% cbind(gene=rownames(.),.) %>% melt()
plotData$Class <- sampleType[plotData$variable] %>%
  factor(levels = c('Normal','nonState1','State1'))
colnames(plotData) <- c('Gene','Sample','Value','State')

plotData$geneType <- ifelse(plotData$Gene %in% dt.markers[1:3,]$gene,1,
                            ifelse(plotData$Gene %in% dt.markers[4:7,]$gene,2,3)) %>%
  factor(levels = c('Normal','nonState1','State1'))

pdf('expLevel_withNormal.pdf',width = 8,height = 4,onefile = F)
ggplot(subset(plotData,geneType==1),
       aes(x=State,y=Value,color=State))+
  geom_boxplot(linetype=1)+
  #geom_jitter()+
  scale_color_manual(values = c('Normal'='#ffd460', 'nonState1' = "#00b8a9",
                                'State1' = "#f6416c"))+
  
  geom_signif(comparisons = list(c('Normal','nonState1'),
                                 c('State1','Normal'),
                                 c('nonState1','State1')),
              test = 'wilcox.test', # wilcox.test, t.test
              step_increase = c(0,.05,.07),
              map_signif_level = T)+
  
  
  facet_wrap(~Gene)+
  
  theme_classic()#+
#theme(axis.text.y = element_blank(),axis.title.x = element_blank());p1
dev.off()

pheno <- read.csv('TCGA-BRCA.GDC_phenotype.csv',header = T) %>% 
  .[grep('-11',.$submitter_id.samples,invert = T),]
# recurrence ====
plotData <- data.frame(
  Sample=pheno$submitter_id.samples,
  Status=pheno[,66]
) %>% subset(Status != '' & Sample %in% colnames(exp))
table(plotData$Status)

plotData$Class <- factor(sampleType[plotData$Sample],levels = c('Normal','nonState1','State1'))
table(plotData$Class)

plotData <- table(plotData$Status,plotData$Class) %>% data.frame()# %>% 
#reshape2::dcast(Var1~Var2,value.var = 'Freq')
colnames(plotData) <- c('Status','Class','Freq')
plotData$Status <- factor(plotData$Status, levels = c('NO','YES'))
plotData$Class <- factor(plotData$Class, levels = c('nonState1','State1'))

library(dplyr)
library(plyr)
plotData <- ddply(plotData, "Status", transform,
                  percent = Freq / sum(Freq) * 100)

plotData1 <- plotData
plotData1$Type <- 'Recurrence'

pdf('recurrence_barplot.pdf',width = 3.3,height = 4,onefile = F)
ggplot(plotData, aes(x=Status, y=percent, fill=Class)) +
  geom_bar(stat = "identity") +
  #coord_flip()+
  scale_fill_manual(values = c('Normal'="#e8d3ff", 'nonState1'="#c6cfff", 'State1'='#deecff'))+
  
  theme_classic()
dev.off()

# stage ====
library(Hmisc)
plotData <- data.frame(
  Sample=pheno$submitter_id.samples,
  Status=pheno[,119]
) %>% subset(Status %nin% c('','not reported','stage x') & Sample %in% colnames(exp))
table(plotData$Status)

plotData$Status[grep('stage iv',plotData$Status)] <- 'Stage IV'
plotData$Status[grep('stage iii',plotData$Status)] <- 'Stage III'
plotData$Status[grep('stage ii',plotData$Status)] <- 'Stage II'
plotData$Status[grep('stage i',plotData$Status)] <- 'Stage I'

plotData$Class <- factor(sampleType[plotData$Sample],levels = c('nonState1','State1'))
table(plotData$Class)

plotData <- table(plotData$Status,plotData$Class) %>% data.frame()# %>% 
#reshape2::dcast(Var1~Var2,value.var = 'Freq')
colnames(plotData) <- c('Status','Class','Freq')
plotData$Status <- factor(plotData$Status, levels = c('Stage I','Stage II','Stage III','Stage IV'))
plotData$Class <- factor(plotData$Class, levels = c('nonState1','State1'))

plotData <- ddply(plotData, "Status", transform,
                  percent = Freq / sum(Freq) * 100)

plotData2 <- plotData
plotData2$Type <- 'Stage'

pdf('stage_barplot.pdf',width = 5,height = 4,onefile = F)
ggplot(plotData, aes(x=Status, y=percent, fill=Class)) +
  geom_bar(stat = "identity") +
  #coord_flip()+
  scale_fill_manual(values = c('Normal'="#e8d3ff", 'nonState1'="#c6cfff", 'State1'='#deecff'))+
  
  theme_classic()
dev.off()

# tumor size ====
library(Hmisc)
plotData <- data.frame(
  Sample=pheno$submitter_id.samples,
  Status=pheno[,72]
) %>% subset(Status %nin% c('','TX') & Sample %in% colnames(exp))
table(plotData$Status)

plotData$Status[grep('T1',plotData$Status)] <- 'T1'
plotData$Status[grep('T2',plotData$Status)] <- 'T2'
plotData$Status[grep('T3',plotData$Status)] <- 'T3'

plotData$Status[grep('T4',plotData$Status)] <- 'T4'

plotData$Class <- factor(sampleType[plotData$Sample],levels = c('nonState1','State1'))
table(plotData$Class)

plotData <- table(plotData$Status,plotData$Class) %>% data.frame()# %>% 
#reshape2::dcast(Var1~Var2,value.var = 'Freq')
colnames(plotData) <- c('Status','Class','Freq')
plotData$Status <- factor(plotData$Status, levels = c('T1','T2','T3','T4'))
plotData$Class <- factor(plotData$Class, levels = c('nonState1','State1'))

plotData <- ddply(plotData, "Status", transform,
                  percent = Freq / sum(Freq) * 100)

plotData3 <- plotData
plotData3$Type <- 'TumorSize'

pdf('tumorSize_barplot.pdf',width = 5,height = 4,onefile = F)
ggplot(plotData, aes(x=Status, y=percent, fill=Class)) +
  geom_bar(stat = "identity") +
  #coord_flip()+
  scale_fill_manual(values = c('Normal'="#e8d3ff", 'nonState1'="#c6cfff", 'State1'='#deecff'))+
  
  theme_classic()
dev.off()

# lymph node ====
library(Hmisc)
plotData <- data.frame(
  Sample=pheno$submitter_id.samples,
  Status=pheno[,71]
) %>% subset(Status %nin% c('','NX') & Sample %in% colnames(exp))
table(plotData$Status)

plotData$Status[grep('N0',plotData$Status)] <- 'N0'
plotData$Status[grep('N1|N2|N3',plotData$Status)] <- 'N1'

plotData$Class <- factor(sampleType[plotData$Sample],levels = c('nonState1','State1'))
table(plotData$Class)

plotData <- table(plotData$Status,plotData$Class) %>% data.frame()# %>% 
#reshape2::dcast(Var1~Var2,value.var = 'Freq')
colnames(plotData) <- c('Status','Class','Freq')
plotData$Status <- factor(plotData$Status, levels = c('N0','N1','N2','N3'))
plotData$Class <- factor(plotData$Class, levels = c('nonState1','State1'))

plotData <- ddply(plotData, "Status", transform,
                  percent = Freq / sum(Freq) * 100)

plotData4 <- plotData
plotData4$Type <- 'LymphNode'

pdf('lymphNode_barplot.pdf',width = 5,height = 4,onefile = F)
ggplot(plotData, aes(x=Status, y=percent, fill=Class)) +
  geom_bar(stat = "identity") +
  #coord_flip()+
  scale_fill_manual(values = c('Normal'="#e8d3ff", 'nonState1'="#c6cfff", 'State1'='#deecff'))+
  
  theme_classic()
dev.off()

# metastasis ====
library(Hmisc)
plotData <- data.frame(
  Sample=pheno$submitter_id.samples,
  Status=pheno[,70]
) %>% subset(Status %nin% c('','MX') & Sample %in% colnames(exp))
table(plotData$Status)

plotData$Status[grep('M0',plotData$Status)] <- 'M0'
plotData$Status[grep('M1',plotData$Status)] <- 'M1'

plotData$Class <- factor(sampleType[plotData$Sample],levels = c('nonState1','State1'))
table(plotData$Class)

plotData <- table(plotData$Status,plotData$Class) %>% data.frame()# %>% 
#reshape2::dcast(Var1~Var2,value.var = 'Freq')
colnames(plotData) <- c('Status','Class','Freq')
plotData$Status <- factor(plotData$Status, levels = c('M0','M1'))
plotData$Class <- factor(plotData$Class, levels = c('nonState1','State1'))

plotData <- ddply(plotData, "Status", transform,
                  percent = Freq / sum(Freq) * 100)

plotData5 <- plotData
plotData5$Type <- 'Metastasis'

pdf('lymphNode_barplot.pdf',width = 5,height = 4,onefile = F)
ggplot(plotData, aes(x=Status, y=percent, fill=Class)) +
  geom_bar(stat = "identity") +
  #coord_flip()+
  scale_fill_manual(values = c('Normal'="#e8d3ff", 'nonState1'="#c6cfff", 'State1'='#deecff'))+
  
  theme_classic()
dev.off()

# bind clinical features ====
plotData.bind <- rbind(plotData1,plotData2,plotData3,plotData4,plotData5)
plotData.bind$Type <- factor(plotData.bind$Type,
                             levels = c('Recurrence',"Stage","LymphNode","TumorSize","Metastasis"))

pdf('clinicalBind_barplot.pdf',width = 9,height = 4,onefile = F)
ggplot(plotData.bind, aes(x=Status, y=percent, fill=Class)) +
  geom_bar(stat = "identity") +
  #coord_flip()+
  scale_fill_manual(values = c('Normal'="#e8d3ff", 'nonState1'="#c6cfff", 'State1'='#deecff'))+
  
  theme_classic()
dev.off()
