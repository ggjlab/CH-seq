install.packages("readxl")
library(readxl)
setwd("G:\\iMT\\itermedia")
epi<-read_excel("epistart.xlsx")
endo<-read_excel("endostart.xlsx")
gm<-read_excel("gmstart.xlsx")
xenneu<-read_excel("xenneustart.xlsx")
epigene<-epi$MyList
endogene<-endo$MyList
gmgene<-gm$MyList
xenneugene<-xenneu$MyList

epigene<-data.frame(epigene,epigene)
colnames(epigene)<-c("Gene","Epi")
epigene$Epi<-"Epi"
endogene<-data.frame(endogene,endogene)
colnames(endogene)<-c("Gene","Endo")
endogene$Endo<-"Endo"
gmgene<-data.frame(gmgene,gmgene)
colnames(gmgene)<-c("Gene","Gm")
gmgene$Gm<-"Gm"
xenneugene<-data.frame(xenneugene,xenneugene)
colnames(xenneugene)<-c("Gene","Xenneu")
xenneugene$Xenneu<-"Xenneu"

gene<-merge(epigene,endogene,by = "Gene",all = T)
gene<-merge(gene,gmgene,by = "Gene",all = T)
gene<-merge(gene,xenneugene,by = "Gene",all = T)
gene[is.na(gene)]<-0
gene$freq <- rep("0")
?setdiff
for (i in 1:498) {
  gene$freq[i] <- length(setdiff(as.character(gene[i,]),"0"))-1
}
write.csv(gene,"4.csv")
