library(rhdf5)
library(Seurat)
library(openxlsx)
library(fgsea)
library(DOSE)
data <- H5Fopen("/Users/siqijiang/Desktop/src/realdata/GSE178707_neatseq_lane1.h5")
setwd("/Users/siqijiang/Desktop/画图/")


ADT <- as.data.frame(data$X2)
mRNA <- as.data.frame(data$X1)
ATAC <- as.data.frame(data$X4)

##Norm后再转置的，目的是为了加行和列的名，因为gene name 重复的，不转置加不上去
mRNA <- as.data.frame(t(NormalizeData(mRNA, normalization.method = "LogNormalize", scale.factor = 10000)))
ATAC <- as.data.frame(t(NormalizeData(ATAC, normalization.method = "LogNormalize", scale.factor = 10000)))

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

#install.packages("openxlsx") # 如果没有openxlsx包，运行该命令
#library(openxlsx)
sheets = list("GATA3_mRNA" = ADT_mRNA,"GATA3_ATAC" = ADT_ATAC)
write.xlsx(sheets,"GATA3.xlsx")


################################################################################
##Enrich分析

react <- read.gmt("c2.cp.reactome.v7.5.1.symbols.gmt")
kegg <- read.gmt("c2.cp.kegg.v7.5.1.symbols.gmt")

out.list <- list()
adts <- c()
for (i in c(1)) {
  dat1 <- read.xlsx("GATA3.xlsx", sheet = 1)
  dat1  <- dat1[!is.na(dat1$P_Value),]
  dat1 <- dat1[order(dat1$P_Value),]
  dat1 <- dat1[!duplicated(dat1$Gene_name),]
  dat1.f <- dat1[dat1$P_Value<0.05,] #mRNA
  dat2 <- read.xlsx("GATA3.xlsx", sheet = 1+i)
  dat2  <- dat2[!is.na(dat2$P_Value),]
  dat2 <- dat2[order(dat2$P_Value),]
  dat2 <- dat2[!duplicated(dat2$GeneFromPeaks_name),]
  dat2.f <- dat2[dat2$P_Value<0.05,] #ATAC
  hit <- intersect(dat2.f$GeneFromPeaks_name, dat1.f$Gene_name) # overlapping for mRNA and ATAC
  
  direct1 <-dat1.f[dat1.f$Gene_name%in%hit,4]
  direct2 <-dat2.f[dat2.f$GeneFromPeaks_name%in%hit,4]
  adt <- ("GATA3")
  out <- data.frame(Gene=hit, ADT=adt, DEgene=direct1*-1, DEatac=direct2*-1)
  out.list[[1]] <- out
  adts <- c(adts,adt)
}

out.list <- out.list[c(1)]
names(out.list) <- adts
write.xlsx(out.list, "Neat_hint.xlsx")

en.df.list<-list()
en.list <- list()
#for (i in 1) {
genelist <- out.list[[1]][,1]
if(length(genelist)>100){
  en <- enricher(genelist, TERM2GENE=react,pvalueCutoff = 0.05,qvalueCutoff=0.2)
  en.df <- as.data.frame(en)
  en.df.list[[1]] <- en.df
  names(en.df.list)[1] <- adts[1]
  en.list[[1]] <- en
  names(en.list)[1] <- adts[1]
}

write.xlsx(en.df.list, "Neat-seq_lane1_enrich_react.xlsx")

library(enrichplot)
library(ggupset)
library(cowplot)
#for (i in 1) {
  adt <- names(en.list[1])
  en <- en.list[[1]]
  p1 <- dotplot(en,font.size=15) 
  p2 <- emapplot(pairwise_termsim(en))
  tiff(filename = paste0("Neatseq_react_lane1_",adt,".tiff"),res = 300, height = 2000, width = 5000, compression = "lzw")
  print(plot_grid(p1,p2, ncol=2))
  dev.off()
#}

mode(en)

library(ReactomePA)

out2 <- enrichPathway(genelist, pvalueCutoff=0.05)



enrich.test<-function(gseaFile="c2.cp.reactome.v7.5.1.symbols.gmt", DEs, EEs, method=c("Hypergeometric", "Binomial", "Fisher", "All"), alternative=c("two.sided", "greater", "less"), sorting=TRUE, n.total.down=10,n.total.up=10000) {
  #gseaFile: pathway denfition file in gmt format(tab separated), for each line: ID\tPathway Description\tGene1\tGene2\t...GeneK 
  #DEs: differentially expressed genes
  #EEs: equally expressed genes
  #method: choose one of the three methods or All; 
  #alternative: for Fisher test only; by default two-sided; note Fisher "greater" == Hypergeometric,
  #sorting: output the results by significance in decreasing order
  #n.total.down, up: only pathways with   n.total.down<= #genes <= n.total.up will be considered
  #             some times, the rlts for a pathway with small size may not be stable, and too large may be uninformative
  c3.tft=readLines(gseaFile)
  c3.tft.rlts<-matrix(unlist(lapply(c3.tft, getPWp, DEs=DEs, method=method, alternative=alternative, EEs=EEs)),byrow=T,nrow=length(c3.tft))
  method <- match.arg(method)
  if(method!="All") {
    tmp<-data.frame(Name=as.character(c3.tft.rlts[,1]),
                    pVal=as.numeric(c3.tft.rlts[,2]),
                    nPos=as.numeric(c3.tft.rlts[,3]),
                    nTotal=as.numeric(c3.tft.rlts[,4]),
                    nSize=as.numeric(c3.tft.rlts[,5]),
                    url=as.character(c3.tft.rlts[,6]),
                    hits=as.character(c3.tft.rlts[,7]),
                    stringsAsFactors=FALSE
    )
    
    c3.tft.rlts<-tmp[tmp[,"nTotal"]>=n.total.down & tmp[,"nTotal"]<=n.total.up,]
    
    qVal=p.adjust(as.numeric(c3.tft.rlts[,"pVal"]),method="BH");
    pBonf=p.adjust(as.numeric(c3.tft.rlts[,"pVal"]),method="bonferroni");
    rtns=data.frame(c3.tft.rlts[,1:2],qVal, pBonf, c3.tft.rlts[,-(1:2)],stringsAsFactors=FALSE)
    if (sorting) {
      rtns<-rtns[order(as.numeric(rtns[,"pVal"])),]
    }
  } else {
    
    tmp<-data.frame(Name  =as.character(c3.tft.rlts[,1]),
                    H.pval  =as.numeric(c3.tft.rlts[,2]),
                    B.pval  =as.numeric(c3.tft.rlts[,3]),
                    F.pval  =as.numeric(c3.tft.rlts[,4]),                    
                    nPos  =as.numeric(c3.tft.rlts[,5]),
                    nTotal=as.numeric(c3.tft.rlts[,6]),
                    nSize =as.numeric(c3.tft.rlts[,7]),
                    url   =as.character(c3.tft.rlts[,8]),
                    hits  =as.character(c3.tft.rlts[,9]),
                    stringsAsFactors=FALSE
    )
    
    c3.tft.rlts<-tmp[tmp[,"nTotal"]>=n.total.down & tmp[,"nTotal"]<=n.total.up,]
    
    H.qval=p.adjust(as.numeric(c3.tft.rlts[,"H.pval"]),method="BH");
    H.pBonf=p.adjust(as.numeric(c3.tft.rlts[,"H.pval"]),method="bonferroni");
    
    B.qval=p.adjust(as.numeric(c3.tft.rlts[,"B.pval"]),method="BH");
    B.pBonf=p.adjust(as.numeric(c3.tft.rlts[,"B.pval"]),method="bonferroni");
    
    F.qval=p.adjust(as.numeric(c3.tft.rlts[,"F.pval"]),method="BH");
    F.pBonf=p.adjust(as.numeric(c3.tft.rlts[,"F.pval"]),method="bonferroni");
    
    rtns=data.frame(c3.tft.rlts[,1:2], H.qval, H.pBonf, c3.tft.rlts[3], B.qval, B.pBonf, c3.tft.rlts[4], F.qval, F.pBonf, c3.tft.rlts[,-(1:4)],stringsAsFactors=FALSE)
    if (sorting) {
      rtns<-rtns[order(as.numeric(rtns[,"H.pval"])),]
    }
  }
  data.frame(index=rownames(rtns), rtns,stringsAsFactors=FALSE)
  #Output
  #name: ID
  #p.val: enrichment P value
  #qVal: fdr level
  #pBonf: Bonferroni adjusted P value
  #nPos:  #genes in the pathway are DE
  #nTotal: #DE or EE genes that can be found in the pathway
  #nSize:  #genes in the pathway, nTotal<=nSize
  #url: Pathway despcription
  #hits: the DE genes in the pathway, comma-separated  
}

allgenes <- dat1$Gene_name
ees <- allgenes[!allgenes%in%hit]
des <- hit
results <- enrich.test(gseaFile = "c2.cp.reactome.v7.5.1.symbols.gmt",
                  DEs=des,EEs=ees,method = "Hypergeometric")

write.xlsx(
  results,
  "Neat-seq_lane1_enrich.xlsx",
  colNames = TRUE,
  rowNames = TRUE,

)


library(dplyr)
library(ggplot2)
library(ggrepel)


#按照PValue从低到高排序[升序]
KEGG_dataset <- arrange(KEGG_dataset,KEGG_dataset[,4])
#Pathway列最好转化成因子型，否则作图时ggplot2会将所有Pathway按字母顺序重排序
#将Pathway列转化为因子型
KEGG_dataset$Term <- factor(KEGG_dataset$Term,levels = rev(KEGG_dataset$Term))

#图片背景设定
mytheme <- theme(axis.title=element_text(face="bold", size=14,colour = 'black'), #坐标轴标题
                 axis.text=element_text(face="bold", size=14,colour = 'black'), #坐标轴标签
                 axis.line = element_line(size=0.5, colour = 'black'), #轴线
                 panel.background = element_rect(color='black'), #绘图区边框
                 legend.key = element_blank() #关闭图例边框
)

#绘制KEGG气泡图
p <- ggplot(KEGG_dataset,aes(x=Gene.ratio,y=Term,colour=-1*log10(PValue),size=Count))+
  geom_point()+
  scale_size(range=c(2, 8))+
  scale_colour_gradient(low = "blue",high = "red")+
  theme_bw()+
  ylab("KEGG Pathway Terms")+
  xlab("Gene Ratio")+
  labs(color=expression(-log[10](PValue)))+
  theme(legend.title=element_text(size=14), legend.text = element_text(size=14))+
  theme(axis.title.y = element_text(margin = margin(r = 50)),axis.title.x = element_text(margin = margin(t = 20)))+
  theme(axis.text.x = element_text(face ="bold",color="black",angle=0,vjust=1))
plot <- p+mytheme
plot
#保存图片
ggsave(plot,filename = "KEGG.pdf",width = 10,height = 6,dpi=300)
ggsave(plot,filename = "KEGG.png",width = 10,height = 6,dpi=300)