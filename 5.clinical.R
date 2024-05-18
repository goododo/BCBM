setwd('/Users/icu/Desktop/BCBM/5.clinical')
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

## TCGA BRCA (UCSC Xena)
probe <- read.table('gencode.v22.annotation.gene.probeMap',header = T)
colnames(probe)[1] <- 'Ensembl_ID'

counts <- read.table('TCGA-BRCA.htseq_counts.tsv',header = T,check.names = F)
counts <- merge(probe[,1:2],counts,by='Ensembl_ID')
exp <- counts[match(allMarkers.sig$gene,counts$gene),-1] %>%
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

# surv ====
## library packages ====
library(dplyr)
library(survival)
library(plyr)

# univariate Cox ====
uni_cox <- function(surv.data, gene){
  gene.exp <- exp[gene,] %>% as.numeric()
  gene.group <- ifelse(gene.exp > median(gene.exp),2,1) %>% as.factor(.)
  names(gene.group) <- colnames(exp)
  surv.data$group <- gene.group[match(surv.data$sample,names(gene.group))]
  surv <- Surv(time = surv.data$time, event = surv.data$event)
  surv.data$surv<-with(surv.data,surv) 
  geneCox<-coxph(surv~group,data = surv.data)
  geneCox<-summary(geneCox) 
  CI<-paste0(round(geneCox$conf.int[,3:4],2),collapse = "-") 
  Pvalue<-round(geneCox$coefficients[,5],3) 
  HR<-round(geneCox$coefficients[,2],2)
  
  ### 查看一下结果
  Unicox<-data.frame("Characteristics"=gene,"Hazard Ratio"=HR,
                     "CI95"=CI,"P value"=Pvalue)
  return(Unicox)
}

library(data.table)
## All year
uniCox.all <- lapply(rownames(exp), function(marker){
  uniCox <- uni_cox(surv_km,marker)
}) %>% rbindlist(.)
uniCox.all <- uniCox.all[order(uniCox.all$P.value),]

uniCox.markers <- subset(allMarkers.sig, gene %in%
                           subset(uniCox.all, P.value < 0.05)$Characteristics)
table(uniCox.markers$cluster)

## 10 years
uniCox.10year <- lapply(rownames(exp), function(marker){
  uniCox <- uni_cox(surv_10years,marker)
}) %>% rbindlist(.)
uniCox.10year <- uniCox.10year[order(uniCox.10year$P.value),]

uniCox10year.markers <- subset(allMarkers.sig, gene %in%
                           subset(uniCox.10year, P.value < 0.05)$Characteristics)
table(uniCox10year.markers$cluster)

# decision tree ====
library(rpart)
#library(randomForest)
epiExp <- epi[rownames(exp),] %>% GetAssayData() %>% as.matrix() %>% t() %>% as.data.frame()
epiExp <- cbind(
  class=epi$pseudotimeState[rownames(epiExp)],
  epiExp)
dt <- rpart(class~.,data = epiExp,method = 'class')
#rf <- randomForest(epiExp[,2:ncol(epiExp)],epiExp[,1]) # error rate 太高

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

# ssGSEA ====
library(GSVA)

dtMarkers.list<- split(dt.markers[,1], dt.markers[,2])

ssgseaRes <- gsva(as.matrix(exp), dtMarkers.list,
                  method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
sampleType <- apply(ssgseaRes, 2, function(x){
  ifelse(x[1]==max(x),1,ifelse(x[2]==max(x),2,3))
})

surv_km$Class <- sampleType
surv_10years$Class <- sampleType

# survival K-M plot ====
library(survival)
library(survminer)

pairwise_survdiff(Surv(time, event) ~ Class,
                  data = surv_10years, p.adjust.method = 'none')

# 10 years生存曲线
fit1 <- survfit(Surv(time, event) ~ Class, data= surv_10years)

pdf('10years_KMcurv.pdf',width = 6,height = 6,onefile = F)
ggsurvplot(fit1, data = surv_10years,
           pval=TRUE, 
           #pval.coord = c(0, 0.2),
           pval.size =5,
           pval.method=TRUE,
           conf.int=T, 
           risk.table=TRUE,
           risk.table.height = 0.25,
           palette=c('#00b8a9','#f6416c','#ffd460'), 
           #title="Kaplan-Meier Curve for OS ",
           legend.labs=c("State1", "State2","State3"), 
           legend.title="State", 
           surv.median.line = "hv" 
)
dev.off()

pairwise_survdiff(Surv(time, event) ~ Class,
                  data = surv_km, p.adjust.method = 'none')

# 总体生存曲线
fit1 <- survfit(Surv(time, event) ~ Class, data= surv_km)

pdf('all_KMcurv.pdf',width = 6,height = 6,onefile = F)
ggsurvplot(fit1, data = surv_km,
           pval=TRUE, 
           #pval.coord = c(0, 0.2),
           pval.size =5,
           pval.method=TRUE,
           conf.int=T, 
           risk.table=TRUE,
           risk.table.height = 0.25,
           palette=c('#00b8a9','#f6416c','#ffd460'), 
           #title="Kaplan-Meier Curve for OS ",
           legend.labs=c("State1", "State2","State3"), 
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
plotData$Class <- sampleType[plotData$variable] %>% as.factor()
colnames(plotData) <- c('Gene','Sample','Value','State')

plotData$geneType <- ifelse(plotData$Gene %in% c("ZNF831","CTLA4","GIMAP7"),1,
                            ifelse(plotData$Gene %in% c("ARHGAP24","EMP3","DNAJC18"),2,3)) %>%
  factor(levels = c(1:3))

## with p value
p1 <- ggplot(subset(plotData,geneType==1),
             aes(x=State,y=Value,color=State))+
  geom_boxplot(linetype=1)+
  #geom_jitter()+
  scale_color_manual(values = c('1' = "#00b8a9", '2' = "#f6416c", '3' = '#ffd460'))+
  
  geom_signif(comparisons = list(c(1,2),
                                 c(3,1),
                                 c(2,3)),
              test = 'wilcox.test', # wilcox.test, t.test
              step_increase = c(.05,.1,.11),
              map_signif_level = T)+
  
  
  facet_wrap(~Gene)+
  
  theme_classic()#+
  #theme(axis.text.y = element_blank(),axis.title.x = element_blank());p1

p2 <- ggplot(subset(plotData,geneType==2),
             aes(x=State,y=Value,color=State))+
  geom_boxplot(linetype=1)+
  #geom_jitter()+
  scale_color_manual(values = c('1' = "#00b8a9", '2' = "#f6416c", '3' = '#ffd460'))+
  
  geom_signif(comparisons = list(c(1,2),
                                 c(3,1),
                                 c(2,3)),
              test = 'wilcox.test', # wilcox.test, t.test
              step_increase = c(.05,.1,.11),
              map_signif_level = T)+
  
  
  facet_wrap(~Gene)+
  
  theme_classic()#+
  #theme(axis.text.y = element_blank(),axis.title.x = element_blank());p2

p3 <- ggplot(subset(plotData,geneType==3),
             aes(x=State,y=Value,color=State))+
  geom_boxplot(linetype=1)+
  #geom_jitter()+
  scale_color_manual(values = c('1' = "#00b8a9", '2' = "#f6416c", '3' = '#ffd460'))+
  
  geom_signif(comparisons = list(c(1,2),
                                 c(3,1),
                                 c(2,3)),
              test = 'wilcox.test', # wilcox.test, t.test
              step_increase = c(.05,.1,.11),
              map_signif_level = T)+
  
  facet_wrap(~Gene,nrow = 1)+
  
  theme_classic()#+
  #theme(axis.text.y = element_blank(),axis.title.x = element_blank());p3

library(lemon)
pdf('expLevel_withPvalue.pdf',width = 13.2,height = 4.5,onefile = F)
grid_arrange_shared_legend(p1, p2, p3, ncol = 3, nrow = 1, position='top')
dev.off()

## without p value
p1 <- ggplot(subset(plotData,geneType==1),
             aes(x=State,y=Value,color=State))+
  geom_boxplot(linetype=1)+
  #geom_jitter()+
  scale_color_manual(values = c('1' = "#00b8a9", '2' = "#f6416c", '3' = '#ffd460'))+
  
  facet_wrap(~Gene)+
  
  theme_classic()#+
  #theme(axis.text.y = element_blank(),axis.title.x = element_blank());p1

p2 <- ggplot(subset(plotData,geneType==2),
             aes(x=State,y=Value,color=State))+
  geom_boxplot(linetype=1)+
  #geom_jitter()+
  scale_color_manual(values = c('1' = "#00b8a9", '2' = "#f6416c", '3' = '#ffd460'))+
  
  facet_wrap(~Gene)+
  
  theme_classic()#+
  #theme(axis.text.y = element_blank(),axis.title.x = element_blank());p2

p3 <- ggplot(subset(plotData,geneType==3),
             aes(x=State,y=Value,color=State))+
  geom_boxplot(linetype=1)+
  #geom_jitter()+
  scale_color_manual(values = c('1' = "#00b8a9", '2' = "#f6416c", '3' = '#ffd460'))+
  
  facet_wrap(~Gene,nrow = 1)+
  
  theme_classic()#+
  #theme(axis.text.y = element_blank(),axis.title.x = element_blank());p3

pdf('expLevel_withoutPvalue.pdf',width = 11.2,height = 3,onefile = F)
grid_arrange_shared_legend(p1, p2, p3, ncol = 3, nrow = 1, position='top')
dev.off()

## ssgsea score
plotData <- ssgseaRes %>% melt()
plotData$Class <- sampleType[plotData$Var2] %>% factor(levels = 1:3)
colnames(plotData) <- c('State','Sample','Value','Class')
plotData$State <- factor(plotData$State)

pdf('ssgseaScore_withPvalue.pdf',width = 8,height = 3.5,onefile = F)
ggplot(plotData,
       aes(x=State,y=Value,color=State))+
  geom_boxplot(linetype=1)+
  #geom_jitter()+
  scale_color_manual(values = c('1' = "#00b8a9", '2' = "#f6416c", '3' = '#ffd460'))+
  
  geom_signif(comparisons = list(c(1,2),
                                 c(3,1),
                                 c(2,3)),
              test = 'wilcox.test', # wilcox.test, t.test
              step_increase = c(.05,.1,.11),
              map_signif_level = T)+
  
  facet_wrap(~Class)+
  
  theme_classic()#+
#theme(axis.text.y = element_blank(),axis.title.x = element_blank())
dev.off()

pdf('ssgseaScore_withoutPvalue.pdf',width = 6,height = 2,onefile = F)
ggplot(plotData,
       aes(x=Class,y=Value,color=State))+
  geom_boxplot(linetype=1)+
  #geom_jitter()+
  scale_color_manual(values = c('1' = "#00b8a9", '2' = "#f6416c", '3' = '#ffd460'))+
  
  theme_classic()#+
#theme(axis.text.y = element_blank(),axis.title.x = element_blank())
dev.off()

pheno <- read.csv('TCGA-BRCA.GDC_phenotype.csv',header = T)
# recurrence ====
plotData <- data.frame(
  Sample=pheno$submitter_id.samples,
  Status=pheno[,66]
) %>% subset(Status != '' & Sample %in% colnames(exp))
table(plotData$Status)

plotData$Class <- factor(sampleType[plotData$Sample],levels = 1:3)
table(plotData$Class)

plotData <- table(plotData$Status,plotData$Class) %>% data.frame()# %>% 
#reshape2::dcast(Var1~Var2,value.var = 'Freq')
colnames(plotData) <- c('Status','Class','Freq')
plotData$Status <- factor(plotData$Status, levels = c('YES','NO'))
plotData$Class <- factor(plotData$Class, levels = c(1:3))

plotData <- ddply(plotData, "Status", transform,
                  percent = Freq / sum(Freq) * 100)

plotData1 <- plotData
plotData1$Type <- 'Recurrence'

pdf('recurrence_barplot.pdf',width = 3.3,height = 4,onefile = F)
ggplot(plotData, aes(x=Status, y=percent, fill=Class)) +
  geom_bar(stat = "identity") +
  #coord_flip()+
  scale_fill_manual(values = c('1' = "#e8d3ff", '2' = "#c6cfff", '3' = '#deecff'))+
  
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

plotData$Class <- factor(sampleType[plotData$Sample],levels = 1:3)
table(plotData$Class)

plotData <- table(plotData$Status,plotData$Class) %>% data.frame()# %>% 
#reshape2::dcast(Var1~Var2,value.var = 'Freq')
colnames(plotData) <- c('Status','Class','Freq')
plotData$Status <- factor(plotData$Status, levels = c('Stage I','Stage II','Stage III','Stage IV'))
plotData$Class <- factor(plotData$Class, levels = c(1:3))

plotData <- ddply(plotData, "Status", transform,
                  percent = Freq / sum(Freq) * 100)

plotData2 <- plotData
plotData2$Type <- 'Stage'

pdf('stage_barplot.pdf',width = 5,height = 4,onefile = F)
ggplot(plotData, aes(x=Status, y=percent, fill=Class)) +
  geom_bar(stat = "identity") +
  #coord_flip()+
  scale_fill_manual(values = c('1' = "#e8d3ff", '2' = "#c6cfff", '3' = '#deecff'))+
  
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
plotData$Status[grep('N1',plotData$Status)] <- 'N1'
plotData$Status[grep('N2',plotData$Status)] <- 'N2'
plotData$Status[grep('N3',plotData$Status)] <- 'N3'

plotData$Class <- factor(sampleType[plotData$Sample],levels = 1:3)
table(plotData$Class)

plotData <- table(plotData$Status,plotData$Class) %>% data.frame()# %>% 
#reshape2::dcast(Var1~Var2,value.var = 'Freq')
colnames(plotData) <- c('Status','Class','Freq')
plotData$Status <- factor(plotData$Status, levels = c('N0','N1','N2','N3'))
plotData$Class <- factor(plotData$Class, levels = c(1:3))

plotData <- ddply(plotData, "Status", transform,
                  percent = Freq / sum(Freq) * 100)

plotData3 <- plotData
plotData3$Type <- 'LymphNode'

pdf('lymphNode_barplot.pdf',width = 5,height = 4,onefile = F)
ggplot(plotData, aes(x=Status, y=percent, fill=Class)) +
  geom_bar(stat = "identity") +
  #coord_flip()+
  scale_fill_manual(values = c('1' = "#e8d3ff", '2' = "#c6cfff", '3' = '#deecff'))+
  
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

plotData$Class <- factor(sampleType[plotData$Sample],levels = 1:3)
table(plotData$Class)

plotData <- table(plotData$Status,plotData$Class) %>% data.frame()# %>% 
#reshape2::dcast(Var1~Var2,value.var = 'Freq')
colnames(plotData) <- c('Status','Class','Freq')
plotData$Status <- factor(plotData$Status, levels = c('T1','T2','T3','T4'))
plotData$Class <- factor(plotData$Class, levels = c(1:3))

plotData <- ddply(plotData, "Status", transform,
                  percent = Freq / sum(Freq) * 100)

plotData4 <- plotData
plotData4$Type <- 'TumorSize'

pdf('tumorSize_barplot.pdf',width = 5,height = 4,onefile = F)
ggplot(plotData, aes(x=Status, y=percent, fill=Class)) +
  geom_bar(stat = "identity") +
  #coord_flip()+
  scale_fill_manual(values = c('1' = "#e8d3ff", '2' = "#c6cfff", '3' = '#deecff'))+
  
  theme_classic()
dev.off()

# bind clinical features ====
plotData.bind <- rbind(plotData1,plotData2,plotData3,plotData4)
plotData.bind$Type <- factor(plotData.bind$Type,
                             levels = c('Recurrence',"Stage","LymphNode","TumorSize"))

pdf('clinicalBind_barplot.pdf',width = 9,height = 4,onefile = F)
ggplot(plotData.bind, aes(x=Status, y=percent, fill=Class)) +
  geom_bar(stat = "identity") +
  #coord_flip()+
  scale_fill_manual(values = c('1' = "#e8d3ff", '2' = "#c6cfff", '3' = '#deecff'))+
  
  theme_classic()
dev.off()
