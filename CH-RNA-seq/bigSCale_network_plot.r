setwd("F:\\")
library(Seurat)
library(bigSCale)
library(Matrix)
library(MASS)
library(tidyverse)
library(reshape2)

data<-readRDS("dge.RDS")
anno<-read.csv("ann.csv",row.names = 1)
META<-anno[anno$Type=="META",]
WT<-anno[anno$Type=="WT",]
METAdata<-data[,META$pseudo.id]
WTdata<-data[,WT$pseudo.id]

gene.names<-rownames(METAdata)
results.META=compute.network(expr.data = METAdata,gene.names = gene.names)

gene.names<-rownames(WTdata)
results.WT=compute.network(expr.data = WTdata,gene.names = gene.names)

aa<-list()

aa[[1]]<-results.META
aa[[2]]<-results.WT
output<-homogenize.networks(input.networks = aa)


setwd("E:\\AXO\\Network/")
write.csv(results.META$centrality,file = "META_centrality.csv")
write.csv(results.WT$centrality,file = "wt_centrality.csv")

library(igraph)
g<-layout_with_fr(results.WT$graph,niter =10000,grid = "nogrid")
METAg<-g
wtg<-g

aa<-cluster_fast_greedy(results.WT$graph)
aa<-membership(aa) 
aa<-as.matrix(aa)
aa<-as.data.frame(aa)
write.csv(aa,file = "wt_cluster.csv")

library(RColorBrewer)
col_flg<-colorRampPalette(brewer.pal(8,"Set1"))(length(unique(aa[,1])))
for (i in 1:length(aa$V1)) {
  aa$V1[i]<-col_flg[as.numeric(aa$V1[i])]
}

pdf("META_top50.pdf",w=10,h=10)
plot(results.META$graph,layout=METAg,vertex.size=log(degree(results.META$graph))*1.5,vertex.label = ifelse(degree(results.META$graph) > 10000,names(V(results.META$graph)),NA),edge.width=0.1,vertex.color=ifelse(results.META$centrality$PAGErank > 0.002289029,"Red","White"))
dev.off()
pdf("META_group.pdf",w=10,h=10)
plot(results.META$graph,layout=METAg,vertex.size=log(degree(results.META$graph))*1.5,vertex.label = ifelse(degree(results.META$graph) > 10000,names(V(results.META$graph)),NA),edge.width=0.1,vertex.color=aa$V1)
dev.off()
pdf("wt_top50.pdf",w=10,h=10)
plot(results.WT$graph,layout=wtg,vertex.size=log(degree(results.WT$graph))*1.5,vertex.label = ifelse(degree(results.WT$graph) > 10000,names(V(results.WT$graph)),NA),edge.width=0.1,vertex.color=ifelse(results.WT$centrality$PAGErank > 0.00228775,"Red","White"))
dev.off()
pdf("wt_group.pdf",w=10,h=10)
plot(results.WT$graph,layout=wtg,vertex.size=log(degree(results.WT$graph))*1.5,vertex.label = ifelse(degree(results.WT$graph) > 10000,names(V(results.WT$graph)),NA),edge.width=0.1,vertex.color=aa$V1)
dev.off()
gene<-c("Isx","Cdx1","Cdx2","Trim31","Morc1","Nr1i2","Atoh1","Rbfox2","Ndn","Lrrc19","Osr1","Brca1","Prrx1","Ddr2","Cacng7","Cenpk","Rgma","Nsun2","Gli3","Larp6","Mixl1","Pygo1","Tacc3","Gas1","Pdgfra","Wnt5a","Esrp1","Aurka","Traip","Vsx2","Fermt1","Il1f9","Zfp334","Exosc2","Ssc5d","Ccnb1","Gli2","Iapp","Il17d","Lyar","Klf17","Sfrp1","Ckap2","Eif2s1","Plk1","Cdc45","Lum","Megf8","Rpl26","Uhrf1")
gene<-toupper(gene)

pdf("wt_gene.pdf",w=10,h=10)
plot(results.WT$graph,layout=wtg,vertex.size=log(degree(results.WT$graph))*1.5,vertex.label = ifelse(names(V(results.WT$graph)) %in%gene,names(V(results.WT$graph)),NA),edge.width=0.1,vertex.color=ifelse(results.WT$centrality$PAGErank>0.00228774,"Red","White"))
dev.off()
pdf("META_gene.pdf",w=10,h=10)
plot(results.META$graph,layout=METAg,vertex.size=log(degree(results.META$graph))*1.5,vertex.label = ifelse(names(V(results.META$graph)) %in%gene,names(V(results.META$graph)),NA),edge.width=0.1,vertex.color=ifelse(results.META$centrality$PAGErank>0.002289028,"Red","White"))
dev.off()




save(METAg,results.META,aa,file = "results_META.rdata")
save(wtg,results.WT,aa,file = "results_wt.rdata")

########order100gene
ab<-results.META$centrality$PAGErank
ab<-ab[order(-as.numeric(ab))]
ab[50]
ab<-results.WT$centrality$PAGErank
ab<-ab[order(-as.numeric(ab))]
ab[50]



write.csv(results_wt$centrality,file = "wt_centrality.csv")
write.csv(aa,file = "wt_cluster.csv")
plot(1:21,rep(1,21),col= col_flg,pch=16,cex=2)



aa<-cluster_fast_greedy(results_meta$graph)
aa<-membership(aa) 
aa<-as.array(aa2)
aa<-as.matrix(aa)
aa<-as.data.frame(aa)
library(RColorBrewer)
col_flg<-colorRampPalette(brewer.pal(8,"Set1"))(length(unique(aa[,1])))
for (i in 1:length(aa$V1)) {
  aa$V1[i]<-col_flg[as.numeric(aa$V1[i])]
}
write.csv(results_meta$centrality,file = "meta_centrality.csv")
results_wt$graph
edge_density(results_meta$graph, loops = FALSE)




############split cluster
setwd("F:\\META\\Split\\")
{
load("results_META.rdata")
load("results_wt.rdata")
aa<-cluster_fast_greedy(results.WT$graph)
aa<-membership(aa) 
aa<-as.matrix(aa)
aa<-as.data.frame(aa)

dir.create("wt_split")
setwd("./wt_split")

for (i in 1:15) {
  
  col_flg<-colorRampPalette(brewer.pal(8,"Set1"))(length(unique(aa[,1])))
  col_flg[1:length(unique(aa[,1]))]<-"#FFFFFF"
  col_flg2<-colorRampPalette(brewer.pal(8,"Set1"))(length(unique(aa[,1])))
  col_flg[i]<-col_flg2[i]
  for (j in 1:length(aa$V1)) {
    aa$V2[j]<-col_flg[as.numeric(aa$V1[j])]
  }
  pdf(paste0("wt_group_",i,".pdf"),w=10,h=10)
  plot(results.WT$graph,layout=wtg,vertex.size=log(degree(results.WT$graph))*1.5,vertex.label = ifelse(degree(results.WT$graph) > 10000,names(V(results.WT$graph)),NA),edge.width=0.1,vertex.color=aa$V2)
  dev.off()
}


setwd("../")

aa<-cluster_fast_greedy(results.META$graph)
aa<-membership(aa) 
aa<-as.matrix(aa)
aa<-as.data.frame(aa)

dir.create("META_split")
setwd("./META_split")

for (i in 1:15) {
  
  col_flg<-colorRampPalette(brewer.pal(8,"Set1"))(length(unique(aa[,1])))
  col_flg[1:length(unique(aa[,1]))]<-"#FFFFFF"
  col_flg2<-colorRampPalette(brewer.pal(8,"Set1"))(length(unique(aa[,1])))
  col_flg[i]<-col_flg2[i]
  for (j in 1:length(aa$V1)) {
    aa$V2[j]<-col_flg[as.numeric(aa$V1[j])]
  }
  pdf(paste0("META_group_",i,".pdf"),w=10,h=10)
  plot(results.META$graph,layout=METAg,vertex.size=log(degree(results.META$graph))*1.5,vertex.label = ifelse(degree(results.META$graph) > 10000,names(V(results.META$graph)),NA),edge.width=0.1,vertex.color=aa$V2)
  dev.off()
}
}
