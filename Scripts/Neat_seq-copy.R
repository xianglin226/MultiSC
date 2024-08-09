library(rhdf5)
library(Seurat)
library(openxlsx)
data <- H5Fopen("/Users/siqijiang/Desktop/src/realdata/GSE178707_neatseq_lane2.h5")
setwd("/Users/siqijiang/Desktop/画图/10.23")

ADT <- as.data.frame(data$X2)
mRNA <- as.data.frame(data$X1)
ATAC <- as.data.frame(data$X4)
Peaks <- as.data.frame(data$X3)

##Norm后再转置的，目的是为了加行和列的名，因为gene name 重复的，不转置加不上去
mRNA <- as.data.frame(t(NormalizeData(mRNA, normalization.method = "LogNormalize", scale.factor = 10000)))
ATAC <- as.data.frame(t(NormalizeData(ATAC, normalization.method = "LogNormalize", scale.factor = 10000)))
ADT <- as.data.frame(t(NormalizeData(ADT, normalization.method = "LogNormalize", scale.factor = 10000)))
Peaks <- as.data.frame(t(NormalizeData(Peaks, normalization.method = "LogNormalize", scale.factor = 10000)))


cell_name <- data$Barcodes
ADT_name <- data$ADT
gene_name <- data$Genes
ATAC_name <- data$GeneFromPeaks

rownames(ADT) <- cell_name
colnames(ADT) <- ADT_name

rownames(mRNA) <- cell_name
colnames(mRNA) <- gene_name

rownames(ATAC) <- cell_name
colnames(ATAC) <- ATAC_name


mRNA <- mRNA[-1648,] ## delete the 1647th cell indexed from 0
label <- read.csv('1_pred.csv', header =F) ## read label file
names(label) <- c('label')
m_l <- cbind(mRNA,label) ## concat mRNA info with label info

##According to the label information, filter out the cells we need. Here we need class 0 and class 5
clu_low <- rownames(m_l[which(m_l$label==0),])
clu_high <- rownames(m_l[which(m_l$label==5),])

atac_l <- ATAC[which(rownames(ATAC)%in%clu_low),]
atac_h <- ATAC[which(rownames(ATAC)%in%clu_high),]
mrna_l <- mRNA[which(rownames(mRNA)%in%clu_low),]
mrna_h <- mRNA[which(rownames(mRNA)%in%clu_high),]
adt_l <- ADT[which(rownames(ADT)%in%clu_low),]
adt_h <- ADT[which(rownames(ADT)%in%clu_high),]

##Below, we start to calculate the P_value and logFC value
#Calculate DE for mRNA

List_P_value <- list()
List_logFC <- list()

for(i in 1:36601) ##totally have 36601 genes
{
  res <- wilcox.test(mrna_h[,i],mrna_l[,i],alternative = 'two.sided')
  List_P_value[[i]] <- res$p.value
  res2 <- (mean(mrna_h[,i])-mean(mrna_l[,i]))
  List_logFC[[i]] <- res2
}

p_adjust <- p.adjust(List_P_value, method = "BH")

df_Pvalue <- data.frame(matrix(unlist(List_P_value), nrow=36601, byrow=T))
df_logFC <- data.frame(matrix(unlist(List_logFC), nrow=36601, byrow=T))
df_final <- cbind(df_Pvalue,df_logFC)

ADT_mRNA <- cbind(gene_name,df_final)
ADT_mRNA <- cbind(ADT_mRNA[,1:2],p_adjust,ADT_mRNA[,3])
names(ADT_mRNA) <- c("Gene_name","P_Value","Adjust_P","LogFC")

#Calculate DE for ATAC

List_P_value2 <- list()
List_logFC2 <- list()

for(i in 1:20010) ##totally have 20010 cells
{
  result <- wilcox.test(atac_h[,i],atac_l[,i],alternative = 'two.sided')
  List_P_value2[[i]] <- result$p.value
  result2 <- (mean(atac_h[,i])-mean(atac_l[,i]))
  List_logFC2[[i]] <- result2
}

p_adjust2 <- p.adjust(List_P_value2, method = "BH")

df_Pvalue2 <- data.frame(matrix(unlist(List_P_value2), nrow=20010, byrow=T))
df_logFC2 <- data.frame(matrix(unlist(List_logFC2), nrow=20010, byrow=T))
df_final2 <- cbind(df_Pvalue2,df_logFC2)

ADT_ATAC <- cbind(ATAC_name,df_final2)
ADT_ATAC <- cbind(ADT_ATAC[,1:2],p_adjust2,ADT_ATAC[,3])
names(ADT_ATAC) <- c("GeneFromPeaks_name","P_Value","Adjust_P","LogFC")

#install.packages("openxlsx") # 如果没有openxlsx包，运行该命令
#library(openxlsx)
sheets = list("FOXP3_mRNA" = ADT_mRNA,"FOXP3_ATAC" = ADT_ATAC)
write.xlsx(sheets,"FOXP3_two.sided.xlsx")


