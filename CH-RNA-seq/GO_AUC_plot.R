library(AUCell)
library(GSEABase)
library(data.table)
library(doMC)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(RColorBrewer)

setwd("E:\\AXO\\GO_AUC")

test<-read.table("go.txt")
test<-unique(test$V1)
test<-test[grep("GO",test)]
gmtFile='mmusculus.GO_BP.name.gmt'
geneSets=getGmt(gmtFile)
genset<-geneSets[test]
as_matrix <- function(mat){
  
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}



pbmc<-subset(pbmc,features=rownames(aause))
meta<-pbmc@meta.data
aa<-pbmc@assays$RNA@counts
rownames(aa)<-aause$`5`
pbmc<-CreateSeuratObject(aa,meta.data = meta)
Idents(pbmc)<-pbmc$tissue

anno2<-unique(pbmc$tissue)
table(pbmc$tissue)

i=21
sub<-subset(pbmc,idents=anno2[i])
exprMat<-sub@assays$RNA@counts
aucellRankings=AUCell_buildRankings(exprMat, nCores=1, plotStats=FALSE)
# abline(v=500, col="skyblue3", lwd=5, lty=3) # aucellRankings@nGenesDetected["5%"]
regulonAUC=AUCell_calcAUC(genset, aucellRankings, aucMaxRank=127, nCores=1) # aucellRankings@nGenesDetected["5%"]
regulonMatix=getAUC(regulonAUC)
aucdata=as.data.frame(regulonMatix)

for (i in 22:30) {
  sub<-subset(pbmc,idents=anno2[i])
  exprMat<-sub@assays$RNA@counts
  aucellRankings=AUCell_buildRankings(exprMat, nCores=1, plotStats=FALSE)
  regulonAUC=AUCell_calcAUC(genset, aucellRankings, aucMaxRank=127, nCores=1) # aucellRankings@nGenesDetected["5%"]
  regulonMatix=getAUC(regulonAUC)
  regulonMatix=as.data.frame(regulonMatix)
  aucdata=merge(aucdata,regulonMatix,by='row.names',all=T)
  aucdata[is.na(aucdata)]=0
  row.names(aucdata)<-aucdata$Row.names
  aucdata<-aucdata[,-1]
}
save(aucdata,file = "auc3.rdata")


######plot#######


col_flg1<-c("#279E68","#1F77B4","#FF7F0E")
col_flg2<-c("#279E68","#1F77B4")
gmtFile='mmusculus.GO_BP.name.gmt'
geneSets=getGmt(gmtFile)

genesets<-NULL
for(i in 1:length(geneSets)){
  temp<-data.frame(description=geneSets[[i]]@shortDescription,
                   GOid=geneSets[[i]]@setName)
  genesets<-rbind(genesets,temp)
}

rownames(genesets)<-genesets$GOid

epigo<-c("GO:0030216","GO:0006614","GO:0030162","GO:0032787","GO:0002366","GO:0010817","GO:0043277","GO:0071621")

endogo<-c("GO:0006413","GO:0030198","GO:0043588","GO:0048732","GO:0050673","GO:0001568","GO:0009611","GO:2000147","GO:0030036","GO:0043542")

strgo<-c("GO:0030198","GO:0045055","GO:0048514","GO:0001503","GO:0007423","GO:0061448","GO:0007517","GO:0006099")

secgo<-c("GO:0031960","GO:0070633","GO:0032787","GO:0006979","GO:0007586","GO:0009235","GO:0016126","GO:0002446")
  
j=1
  for(j in 1:length(epigo)){
    go<-epigo[j]
    temp.geneset<-genesets[go,]
    title<-paste0(temp.geneset$GOid," ",temp.geneset$description)
    result<-as.data.frame(Epi@assays$RNA@counts)
    result<-as.data.frame(t(result[go,]))
    colnames(result) <- "value"
    result$stage <- Epi$type
    result$GO<-title
    result$stage<-factor(result$stage,levels=c("wt","meta"))
    aa<-compare_means(value ~ stage, data = result)
    my_comparisons <- list( c("wt", "meta"))
    p1<-ggplot(result,aes(stage,value,fill=stage))+geom_boxplot(width=0.5)  +
      stat_compare_means(label ="p.signif", method = "t.test",comparisons = my_comparisons)+labs(x=title, y= 'Activity')
    p1<-p1+theme_bw() +theme(
      axis.title.x = element_text(size = 13),
      axis.title.y = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      title = element_text(size = 6),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"))+scale_fill_manual(values = col_flg)
    p1
    ggsave(paste0("E:\\AXO\\GO_AUC\\plot\\Epi\\",j,"_final.pdf"),w=4,h=4)
    save(p1,file = paste0("E:\\AXO\\GO_AUC\\plot\\Epi\\",j,".rdata"))
  }

j=1
for(j in 1:length(endogo)){
  go<-endogo[j]
  temp.geneset<-genesets[go,]
  title<-paste0(temp.geneset$GOid," ",temp.geneset$description)
  result<-as.data.frame(Endo@assays$RNA@counts)
  result<-as.data.frame(t(result[go,]))
  colnames(result) <- "value"
  result$stage <- Endo$type
  result$GO<-title
  result$stage<-factor(result$stage,levels=c("wt","meta"))
  aa<-compare_means(value ~ stage, data = result)
  my_comparisons <- list( c("wt", "meta"))
  p1<-ggplot(result,aes(stage,value,fill=stage))+geom_boxplot(width=0.5)  +
    stat_compare_means(label ="p.signif", method = "t.test",comparisons = my_comparisons)+labs(x=title, y= 'Activity')
  p1<-p1+theme_bw() +theme(
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    title = element_text(size = 6),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"))+scale_fill_manual(values = col_flg)
  p1
  ggsave(paste0("E:\\AXO\\GO_AUC\\plot\\Endo\\",j,"_final.pdf"),w=4,h=4)
  save(p1,file = paste0("E:\\AXO\\GO_AUC\\plot\\Endo\\",j,".rdata"))
}

j=1
for(j in 1:length(secgo)){
  go<-secgo[j]
  temp.geneset<-genesets[go,]
  title<-paste0(temp.geneset$GOid," ",temp.geneset$description)
  result<-as.data.frame(Sec@assays$RNA@counts)
  result<-as.data.frame(t(result[go,]))
  colnames(result) <- "value"
  result$stage <- Sec$type
  result$GO<-title
  result$stage<-factor(result$stage,levels=c("wt","meta"))
  aa<-compare_means(value ~ stage, data = result)
  my_comparisons <- list( c("wt", "meta"))
  p1<-ggplot(result,aes(stage,value,fill=stage))+geom_boxplot(width=0.5)  +
    stat_compare_means(label ="p.signif", method = "t.test",comparisons = my_comparisons)+labs(x=title, y= 'Activity')
  p1<-p1+theme_bw() +theme(
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    title = element_text(size = 6),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"))+scale_fill_manual(values = col_flg)
  p1
  ggsave(paste0("E:\\AXO\\GO_AUC\\plot\\Sec\\",j,"_final.pdf"),w=4,h=4)
  save(p1,file = paste0("E:\\AXO\\GO_AUC\\plot\\Sec\\",j,".rdata"))
}

j=1
for(j in 1:length(strgo)){
  go<-strgo[j]
  temp.geneset<-genesets[go,]
  title<-paste0(temp.geneset$GOid," ",temp.geneset$description)
  result<-as.data.frame(Str@assays$RNA@counts)
  result<-as.data.frame(t(result[go,]))
  colnames(result) <- "value"
  result$stage <- Str$type
  result$GO<-title
  result$stage<-factor(result$stage,levels=c("wt","meta"))
  aa<-compare_means(value ~ stage, data = result)
  my_comparisons <- list( c("wt", "meta"))
  p1<-ggplot(result,aes(stage,value,fill=stage))+geom_boxplot(width=0.5) + 
    stat_compare_means(label ="p.signif", method = "t.test",comparisons = my_comparisons)+labs(x=title, y= 'Activity')
  p1<-p1+theme_bw() +theme(
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    title = element_text(size = 6),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"))+scale_fill_manual(values = col_flg)
  p1
  ggsave(paste0("E:\\AXO\\GO_AUC\\plot\\Str\\",j,"_final.pdf"),w=4,h=4)
  save(p1,file = paste0("E:\\AXO\\GO_AUC\\plot\\Str\\",j,".rdata"))
}



