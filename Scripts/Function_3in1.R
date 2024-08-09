library(rhdf5)
library(openxlsx)
library(tidyverse)
library(mediation)

data <- H5Fopen("/Users/siqijiang/Desktop/src/realdata/GSE178707_neatseq_lane1.h5")  ## Input file

##Totally have 3 functions here:
##DE(data)
##ADT_ATAC(data)
##ADT1_ADT2(data)


##Read dataset
DE <- function(x){
  
  ADT <- as.data.frame(t(data$X2))
  mRNA <- as.data.frame(t(data$X1))
  ATAC <- as.data.frame(t(data$X4))
  Peaks <- as.data.frame(t(data$X3))
  
  cell_name <- data$Barcodes
  ADT_name <- data$ADT
  gene_name <- data$Genes
  ATAC_name <- data$GeneFromPeaks
  Peaks_name <- data$Peaks
  
  rownames(ADT) <- cell_name
  colnames(ADT) <- ADT_name
  rownames(mRNA) <- cell_name
  colnames(mRNA) <- gene_name
  rownames(ATAC) <- cell_name
  colnames(ATAC) <- ATAC_name
  rownames(Peaks) <- cell_name
  colnames(Peaks) <- Peaks_name
  
  #Normalizing the data
  ADT <- as.data.frame(NormalizeData(ADT, normalization.method = "LogNormalize", scale.factor = 10000))
  mRNA <- as.data.frame(NormalizeData(mRNA, normalization.method = "LogNormalize", scale.factor = 10000))
  ATAC <- as.data.frame(NormalizeData(ATAC, normalization.method = "LogNormalize", scale.factor = 10000))
  Peaks <- as.data.frame(NormalizeData(Peaks, normalization.method = "LogNormalize", scale.factor = 10000))
  
  mRNA <- mRNA[-5000,] ## delete the 4999th cell indexed from 0
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
  
  for(i in 1:36601) ##totally have 36601 cells
  {
    res <- wilcox.test(mrna_h[,i],mrna_l[,i],alternative = 'greater')
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
    result <- wilcox.test(atac_h[,i],atac_l[,i],alternative = 'greater')
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
  return(ADT_ATAC)
}

ADT_ATAC <- function(x){
  
  ADT <- as.data.frame(t(data$X2))
  mRNA <- as.data.frame(t(data$X1))
  ATAC <- as.data.frame(t(data$X4))
  Peaks <- as.data.frame(t(data$X3))
  
  #Normalizing the data
  ADT <- as.data.frame(NormalizeData(ADT, normalization.method = "LogNormalize", scale.factor = 10000))
  mRNA <- as.data.frame(NormalizeData(mRNA, normalization.method = "LogNormalize", scale.factor = 10000))
  ATAC <- as.data.frame(NormalizeData(ATAC, normalization.method = "LogNormalize", scale.factor = 10000))
  Peaks <- as.data.frame(NormalizeData(Peaks, normalization.method = "LogNormalize", scale.factor = 10000))
  
  cell_name <- data$Barcodes
  ADT_name <- data$ADT
  gene_name <- data$Genes
  ATAC_name <- data$GeneFromPeaks
  Peaks_name <- data$Peaks
  
  rownames(ADT) <- cell_name
  colnames(ADT) <- ADT_name
  rownames(mRNA) <- cell_name
  colnames(mRNA) <- gene_name
  rownames(ATAC) <- cell_name
  colnames(ATAC) <- ATAC_name
  rownames(Peaks) <- cell_name
  colnames(Peaks) <- Peaks_name
  
  ##MLR
  #Find the names of those genes where mRNA and ATAC intersect
  name <- intersect(x=gene_name,y=ATAC_name)
  
  List_t1 <- list()
  List_t2 <- list()
  List_t3 <- list()
  List_p1 <- list()
  List_p2 <- list()
  List_p3 <- list()
  
  for(i in name)
  {
    x1<- ADT[,c("GATA3")]   
    x2 <- ATAC[[i]]
    y <- mRNA[[i]]
    res <- lm(y~x1+x2+x1*x2)
    results <- summary(res)
    table <- results$coefficients
    j <- nrow(table)
    if(j==2){
      List_t1[[i]] <- results$coefficients[2,3]
      List_p1[[i]] <- results$coefficients[2,4]
      List_t2[[i]] <- 0
      List_p2[[i]] <- 0
      List_t3[[i]] <- 0
      List_p3[[i]] <- 0
      
    } else if(j==3) {
      List_t1[[i]] <- results$coefficients[2,3]
      List_p1[[i]] <- results$coefficients[2,4]
      List_t2[[i]] <- results$coefficients[3,3]
      List_p2[[i]] <- results$coefficients[3,4]
      List_t3[[i]] <- 0
      List_p3[[i]] <- 0
    } else if(j==4){
      List_t1[[i]] <- results$coefficients[2,3]
      List_p1[[i]] <- results$coefficients[2,4]
      List_t2[[i]] <- results$coefficients[3,3]
      List_p2[[i]] <- results$coefficients[3,4]
      List_t3[[i]] <- results$coefficients[4,1]
      List_p3[[i]] <- results$coefficients[4,4]
    }
  }
  
  t1 <- data.frame(matrix(unlist(List_t1), nrow=18026, byrow=T))
  t2 <- data.frame(matrix(unlist(List_t2), nrow=18026, byrow=T))
  t3 <- data.frame(matrix(unlist(List_t3), nrow=18026, byrow=T))
  p1 <- data.frame(matrix(unlist(List_p1), nrow=18026, byrow=T))
  p2 <- data.frame(matrix(unlist(List_p2), nrow=18026, byrow=T))
  p3 <- data.frame(matrix(unlist(List_p3), nrow=18026, byrow=T))
  
  final <- cbind(t1,t2,t3,p1,p2,p3)
  names(final) <- c("t1","t2","t3","p1","p2","p3")
  rownames(final) <- name
  final <- cbind(name,final)
  return(final)
}


ADT1_ADT2 <- function(x){
  
ADT <- as.data.frame(t(data$X2))
mRNA <- as.data.frame(t(data$X1))
ATAC <- as.data.frame(t(data$X4))
Peaks <- as.data.frame(t(data$X3))

#Normalizing the data
ADT <- as.data.frame(NormalizeData(ADT, normalization.method = "LogNormalize", scale.factor = 10000))
mRNA <- as.data.frame(NormalizeData(mRNA, normalization.method = "LogNormalize", scale.factor = 10000))
ATAC <- as.data.frame(NormalizeData(ATAC, normalization.method = "LogNormalize", scale.factor = 10000))
Peaks <- as.data.frame(NormalizeData(Peaks, normalization.method = "LogNormalize", scale.factor = 10000))

cell_name <- data$Barcodes
ADT_name <- data$ADT
gene_name <- data$Genes
ATAC_name <- data$GeneFromPeaks
Peaks_name <- data$Peaks

rownames(ADT) <- cell_name
colnames(ADT) <- ADT_name

rownames(mRNA) <- cell_name
colnames(mRNA) <- gene_name

rownames(ATAC) <- cell_name
colnames(ATAC) <- ATAC_name

rownames(Peaks) <- cell_name
colnames(Peaks) <- Peaks_name

##MLR

name <- unique(gene_name)

List_t1 <- list()
List_t2 <- list()
List_t3 <- list()
List_p1 <- list()
List_p2 <- list()
List_p3 <- list()

for(i in name)
{
  x1 <- ADT[,c("GATA3")]   
  x2 <- ADT[,c("FOXP3")]
  y <- mRNA[[i]]
  res <- lm(y~x1+x2+x1*x2)
  results <- summary(res)
  table <- results$coefficients
  j <- nrow(table)
  if(j==2){
    List_t1[[i]] <- results$coefficients[2,3]
    List_p1[[i]] <- results$coefficients[2,4]
    List_t2[[i]] <- 0
    List_p2[[i]] <- 0
    List_t3[[i]] <- 0
    List_p3[[i]] <- 0
    
  } else if(j==3) {
    List_t1[[i]] <- results$coefficients[2,3]
    List_p1[[i]] <- results$coefficients[2,4]
    List_t2[[i]] <- results$coefficients[3,3]
    List_p2[[i]] <- results$coefficients[3,4]
    List_t3[[i]] <- 0
    List_p3[[i]] <- 0
  } else if(j==4){
    List_t1[[i]] <- results$coefficients[2,3]
    List_p1[[i]] <- results$coefficients[2,4]
    List_t2[[i]] <- results$coefficients[3,3]
    List_p2[[i]] <- results$coefficients[3,4]
    List_t3[[i]] <- results$coefficients[4,1]
    List_p3[[i]] <- results$coefficients[4,4]
  }
}

t1 <- data.frame(matrix(unlist(List_t1), nrow=36591, byrow=T))
t2 <- data.frame(matrix(unlist(List_t2), nrow=36591, byrow=T))
t3 <- data.frame(matrix(unlist(List_t3), nrow=36591, byrow=T))
p1 <- data.frame(matrix(unlist(List_p1), nrow=36591, byrow=T))
p2 <- data.frame(matrix(unlist(List_p2), nrow=36591, byrow=T))
p3 <- data.frame(matrix(unlist(List_p3), nrow=36591, byrow=T))

final <- cbind(t1,t2,t3,p1,p2,p3)
names(final) <- c("t1","t2","t3","p1","p2","p3")
rownames(final) <- name
final <- cbind(name,final)
return(final)
}

