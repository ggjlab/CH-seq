## load required libraries
source('./slice.R')
library(Seurat)
library(reshape2)
library(ggplot2)
library(ggsignif)
library(RColorBrewer)

## load the kappa similarity matrix of human genes
load("./hs_km.Rda")

## Data prepraration
load('../WT/Split/Stomach/Stomach.rdata')
table(Idents(pbmcsub))
head(pbmcsub[[]])
annouse<-FetchData(pbmcsub,vars = c('ident','organ','batch'))

## Normalization
dgeuse<-as.data.frame(GetAssayData(pbmcsub,slot = 'counts'))

cpm <- apply(dgeuse,2, function(x) (x/sum(x))*10000) 
cpm<-as.data.frame(cpm)

## axo2human gene translation
geneuse<-colsplit(rownames(cpm), pattern = '-', names = c('c1','c2','c3'))
geneuse<-cbind(rownames(cpm), geneuse)
geneuse<-geneuse[!is.na(geneuse$c2),]
geneuse<-cbind(geneuse, new=toupper(geneuse$c3))
geneuse<-geneuse[!duplicated(geneuse$new),]

cpm<-cpm[geneuse$`rownames(cpm)`,]
rownames(cpm)<-geneuse$new

###### ------ Start SLICE ------ ######
sc <- construct(exprmatrix=cpm, 
                cellidentity=annouse$ident,
                projname="Meta")

## bootstrap calculation of scEntropy
sc <- getEntropy(sc, km=km,                             # use the pre-computed kappa similarity matrix of human genes
                 calculation="bootstrap",               # choose the bootstrap calculation
                 B.num=10,                              # default:100 iterations, can be reduced for efficiency when data is large
                 exp.cutoff=1,                          # the threshold for expressed genes
                 B.size=1000,                           # the size of bootstrap sample
                 clustering.k=floor(sqrt(1000/2)),      # the number of functional clusters  
                 random.seed=123)                    # set the random seed to reproduce the results in the paper

## extract out entropy data
annouse<-FetchData(pbmcsub,vars = c('organ','batch'))
entropy<-sc@entropies
annouse<-cbind(annouse, Entropy = entropy$scEntropy.bootstrap)
colnames(annouse)[1:3]<-c('Organ','Batch','Entropy')
# saving
save(annouse, file = './Wt_Stomach.RData')

###### -------- Downstream analysis -------- ######

## load entropy data 
# META
setwd('./META/')
filename<-list.files(pattern = '.RData')
datause_meta<-data.frame()
for (i in 1:length(filename)) {
  load(filename[i])
  rownames(annouse)<-paste(annouse$Organ, rownames(annouse), sep = '_')
  datause_meta<-rbind(datause_meta, annouse)
}
# WT
setwd('../WT/')
filename<-list.files(pattern = '.RData')
datause_wt<-data.frame()
for (i in 1:length(filename)) {
  load(filename[i])
  rownames(annouse)<-paste(annouse$Organ, rownames(annouse), sep = '_')
  datause_wt<-rbind(datause_wt, annouse)
}

## data preprocessing
datause<-rbind(datause_wt, datause_meta)
table(datause$Organ)
table(datause$Batch)

datause$Batch[grep(datause$Batch, pattern = 'Meta')]<-'META'
datause$Batch[grep(datause$Batch, pattern = 'WT')]<-'WT'
table(datause$Batch)

## split-organs for comparison
organuse<-setdiff(unique(datause$Organ), c('Gill','Pancreas','Prostate'))
datause_split<-datause[datause$Organ %in% organuse,]

compare<-list(c('WT','META'))
ggplot(datause_split, aes(x=Batch, y=Entropy,fill=Batch)) +
  geom_boxplot(color="black", notch = FALSE)+
  scale_fill_manual(values = brewer.pal(2,'Set1'))+
  geom_signif(comparisons = compare,
              step_increase = 0.1, vjust=0.6,
              textsize = 3,
              map_signif_level = F,
              test = 'wilcox.test')+
  theme(axis.text.x = element_text(colour = 'black',angle = 45,hjust = 1),
        axis.text.y = element_text(colour = 'black'),
        panel.border = element_rect(colour = 'black',fill = NA),
        panel.background = element_blank(),
        legend.position = 0)+
  facet_wrap(vars(Organ),nrow = 2)

## META compare to WT
ggplot(datause, aes(x=Batch, y=Entropy,fill=Batch)) +
  geom_boxplot(color="black", notch = FALSE)+
  scale_fill_manual(values = brewer.pal(2,'Set1'))+
  geom_signif(comparisons = compare,
              step_increase = 0.1, vjust=0.6,
              textsize = 3,
              map_signif_level = F,
              test = 'wilcox.test')+
  theme(axis.text.x = element_text(colour = 'black',angle = 45,hjust = 1),
        axis.text.y = element_text(colour = 'black'),
        panel.border = element_rect(colour = 'black',fill = NA),
        panel.background = element_blank(),
        legend.position = 0)

