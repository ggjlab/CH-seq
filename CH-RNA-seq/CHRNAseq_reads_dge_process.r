# define the report folder from sci-RNA-seq pipeline
report_folder ="I:\\923BGIAXO\\axo"

# define the output folder for output the df_cell, df_gene and gene_count matrix
output_folder = "I:\\923BGIAXO\\axo"

suppressMessages(library(Matrix))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

combine_exon_intron <- function (df_gene, gene_count) 
{
  gene_count_exon = gene_count[df_gene$exon_intron == "exon", 
  ]
  gene_count_intron = gene_count[df_gene$exon_intron == "intron", 
  ]
  if (nrow(gene_count_exon) == nrow(gene_count_intron)) {
    gene_count_combine = gene_count_exon + gene_count_intron
  }
  else {
    gene_count_combine = gene_count_exon[-nrow(gene_count_exon), 
    ] + gene_count_intron
    gene_count_combine = rbind(gene_count_combine, gene_count_exon[nrow(gene_count_exon), 
    ])
  }
  return(gene_count_combine)
}

sciRNAseq_gene_count_summary <- function (gene_count_folder) {
  gene_matrix = paste0(gene_count_folder,'/',gene_count_folder, "/count.MM")
  df_gene = paste0(gene_count_folder,'/', gene_count_folder,"/gene_name_annotate.txt")
  df_cell = paste0(gene_count_folder,'/',gene_count_folder, "/cell_annotate.txt")
  df_report = paste0(gene_count_folder,'/',gene_count_folder, "/report.MM")
  report_annotate = paste0(gene_count_folder,'/',gene_count_folder, "/report_annotate.txt")
  df_gene = read.csv(df_gene, header = F)
  df_cell = read.csv(df_cell, header = F)
  gene_matrix = read.csv(gene_matrix, header = F)
  colnames(df_gene) = c("gene_id", "gene_type", "exon_intron", 
                        "gene_name", "index")
  colnames(df_cell) = c("sample", "index")
  rownames(df_gene) = df_gene$gene_id
  rownames(df_cell) = df_cell$cell_name
  gene_count = sparseMatrix(i = gene_matrix$V1, j = gene_matrix$V2, 
                            x = gene_matrix$V3)
  df_gene = df_gene[1:nrow(gene_count), ]
  rownames(gene_count) = df_gene$gene_id
  colnames(gene_count) = df_cell$cell_name
  gene_count = combine_exon_intron(df_gene, gene_count)
  df_gene = df_gene %>% filter(exon_intron == "exon")
  reportMM = read.csv(df_report, header = F)
  df_report = sparseMatrix(i = reportMM$V1, j = reportMM$V2, 
                           x = reportMM$V3)
  df_report = as.matrix(t(df_report))
  df_report_annotate = read.csv(report_annotate, header = F)
  colnames(df_report) = df_report_annotate$V2
  df_report = data.frame(df_report)
  df_report["index"] = as.numeric(rownames(df_report))
  df_cell_combine = inner_join(df_cell, df_report, by = "index")
  df_cell_combine["all_exon"] = df_cell_combine$X.Perfect.intersect.exon.match + 
    df_cell_combine$X.Nearest.intersect.exon.match + df_cell_combine$X.Perfect.combine.exon.match + 
    df_cell_combine$X.Nearest.combine.exon.match
  df_cell_combine["all_intron"] = df_cell_combine$X.Perfect.intersect.gene.match + 
    df_cell_combine$X.Nearest.intersect.gene.match + df_cell_combine$X.Perfect.combine.gene.match + 
    df_cell_combine$X.Nearest.combine.gene.match
  df_cell_combine["all_reads"] = df_cell_combine$all_exon + 
    df_cell_combine$all_intron + df_cell_combine$X.No.match
  df_cell_combine["unmatched_rate"] = df_cell_combine$X.No.match/df_cell_combine$all_reads
  df_cell = df_cell_combine %>% select(sample, unmatched_rate)
  df_cell$UMI_count = df_cell_combine$all_exon + df_cell_combine$all_intron
  df_gene = df_gene %>% select(gene_id, gene_type, gene_name)
  return(list(df_cell, df_gene, gene_count))
}

setwd("I:\\923BGIAXO\\axo")

report_folder<-list.files(pattern = "*COL")
i=1

length(report_folder)
for(i in c(1:length(report_folder))){
  result = sciRNAseq_gene_count_summary(report_folder[i])
  df_cell = result[[1]]
  df_gene = result[[2]]
  df_cell$sample<-as.character(df_cell$sample)
  df_cell$sample<-paste(report_folder[i],df_cell$sample,sep = "_")
  df_gene$gene_id<-as.character(df_gene$gene_id)
  df_gene$gene_name<-as.character(df_gene$gene_name)
  rownames(df_gene)<-df_gene$gene_id
  gene_count = result[[3]]
  cellname<-df_cell$sample
  geneid<-as.character(df_gene$gene_id)
  dataset<-data.frame(as.matrix(gene_count))
  colnames(dataset)<-cellname[1:dim(dataset)[2]]
  geneid_u<-df_gene$gene_id[!duplicated(df_gene$gene_name)]
  df_gene1<-df_gene[geneid_u,]
  dataset<-data.frame(dataset[rownames(df_gene1),],row.names = df_gene1$gene_name)
  dataset1<-na.omit(dataset)
  dataset1<-dataset1[,colSums(dataset1)>10]
  dataset1<-dataset1[rowSums(dataset1)>=3,]
  write.table(dataset1,file = paste(report_folder[i],"_ht_dge.txt",sep = "/"))
  message(report_folder[i])
}



tongji<-data.frame(matrix(ncol = 8,nrow = length(report_folder)))
colnames(tongji)<-c("COL","CELLNUMBER","UMI","GENE","READS","Mean_Gene","Mean_UMI","Mean_Read")

for(i in c(1:length(report_folder))){
  col<-report_folder[i]
  filename<-paste(col,"_ht_dge.txt",sep = "/")
  aa<-read.table(filename,header = T,row.names = 1)
  reads<-read.csv(paste0(col,"/",col,"/",col,"_reads.txt"),header = F,stringsAsFactors = F)
  reads$V1<-paste(col,reads$V1,sep = "_")
  su<-data.frame(UMI=colSums(aa),Gene=colSums(aa>0))
  su<-su[order(su$UMI,decreasing = T),]
  su_50UMI<-su[su$UMI>=500,]
  reads_500UMI<-reads[reads$V1%in%rownames(su_50UMI),]
  tongji$COL[i]<-col
  tongji$CELLNUMBER[i]<-dim(su_50UMI)[1]
  tongji$UMI[i]<-as.character(list(summary(su_50UMI$UMI)))
  tongji$GENE[i]<-as.character(list(summary(su_50UMI$Gene)))
  tongji$READS[i]<-as.character(list(summary(reads_500UMI$V2)))
  tongji$Mean_Gene[i]<-mean(su_50UMI$Gene)
  tongji$Mean_UMI[i]<-mean(su_50UMI$UMI)
  tongji$Mean_Read[i]<-mean(reads_500UMI$V2)
  message(col)
}



tongji<-na.omit(tongji)
write.table(tongji,file = "tongji_reads_500UMI.csv",quote = F,sep = "\t",row.names = F)

tongji<-data.frame(matrix(ncol = 8,nrow = length(report_folder)))
colnames(tongji)<-c("COL","CELLNUMBER","UMI","GENE","READS","Mean_Gene","Mean_UMI","Mean_Read")

for(i in c(1:length(report_folder))){
  col<-report_folder[i]
  filename<-paste(col,"_ht_dge.txt",sep = "/")
  aa<-read.table(filename,header = T,row.names = 1)
  reads<-read.csv(paste0(col,"/",col,"/",col,"_reads.txt"),header = F,stringsAsFactors = F)
  reads$V1<-paste(col,reads$V1,sep = "_")
  su<-data.frame(UMI=colSums(aa),Gene=colSums(aa>0))
  su<-su[order(su$UMI,decreasing = T),]
  su_50UMI<-su[su$UMI>=200,]
  reads_500UMI<-reads[reads$V1%in%rownames(su_50UMI),]
  tongji$COL[i]<-col
  tongji$CELLNUMBER[i]<-dim(su_50UMI)[1]
  tongji$UMI[i]<-as.character(list(summary(su_50UMI$UMI)))
  tongji$GENE[i]<-as.character(list(summary(su_50UMI$Gene)))
  tongji$READS[i]<-as.character(list(summary(reads_500UMI$V2)))
  tongji$Mean_Gene[i]<-mean(su_50UMI$Gene)
  tongji$Mean_UMI[i]<-mean(su_50UMI$UMI)
  tongji$Mean_Read[i]<-mean(reads_500UMI$V2)
  message(col)
}



tongji<-na.omit(tongji)
write.table(tongji,file = "tongji_reads_200UMI.csv",quote = F,sep = "\t",row.names = F)

