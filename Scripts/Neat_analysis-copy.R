library(openxlsx)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
setwd("G:/My Drive/Papers/Neat-seq/ADT_grouping")
react <- read.gmt("c2.cp.reactome.v7.5.1.symbols.gmt")
kegg <- read.gmt("c2.cp.kegg.v7.5.1.symbols.gmt")

out.list <- list()
adts <- c()
for (i in c(1,3,5,7,9,11,13)) {
dat1 <- read.xlsx("ADT+Means.xlsx", sheet = i)
dat1  <- dat1[!is.na(dat1$P_Value),]
dat1 <- dat1[order(dat1$P_Value),]
dat1 <- dat1[!duplicated(dat1$Gene_name),]
dat1.f <- dat1[dat1$P_Value<0.01,]
dat2 <- read.xlsx("ADT+Means.xlsx", sheet = i+1)
dat2  <- dat2[!is.na(dat2$P_Value),]
dat2 <- dat2[order(dat2$P_Value),]
dat2 <- dat2[!duplicated(dat2$GeneFromPeaks_name),]
dat2.f <- dat2[dat2$P_Value<0.01,]
hit <- intersect(dat2.f$GeneFromPeaks_name, dat1.f$Gene_name)
direct1 <-dat1.f[dat1.f$Gene_name%in%hit,5]
direct2 <-dat2.f[dat2.f$GeneFromPeaks_name%in%hit,5]
adt <- dat1$ADT[1]
out <- data.frame(Gene=hit, ADT=adt, DEgene=direct1*-1, DEatac=direct2*-1)
out.list[[i]] <- out
adts <- c(adts,adt)
}

out.list <- out.list[c(1,3,5,7,9,11,13)]
names(out.list) <- adts
write.xlsx(out.list, "Neat-seq_lane1_Hitgenes.xlsx")

en.df.list<-list()
en.list <- list()
for (i in 1:7) {
  genelist <- out.list[[i]][,1]
  if(length(genelist)>100){
  en <- enricher(genelist, TERM2GENE=kegg)
  en.df <- as.data.frame(en)
  en.df.list[[i]] <- en.df
  names(en.df.list)[i] <- adts[i]
  en.list[[i]] <- en
  names(en.list)[i] <- adts[i]
  }
}

#en.list <- en.list[-c(1,4)]
#en.df.list <- en.df.list[-c(1,4)]
write.xlsx(en.df.list, "Neat-seq_lane1_enrich_kegg.xlsx")

library(enrichplot)
library(ggupset)
library(cowplot)
for (i in 1:5) {
adt <- names(en.list[i])
en <- en.list[[i]]
p1 <- dotplot(en) + ggtitle(adt)
p2 <- barplot(en) + ggtitle(adt)
p3 <- upsetplot(en)
p4 <- emapplot(pairwise_termsim(en))
tiff(filename = paste0("Neatseq_kegg_lane1_",adt,".tiff"),res = 300, height = 5000, width = 7000, compression = "lzw")
print(plot_grid(p1,p3,p2,p4, ncol=2))
dev.off()
}

#############################################################################################################################
out.list <- list()
adts <- c()
for (i in c(1,3,5,7,9,11,13)) {
  dat1 <- read.xlsx("ADT+Means_line2.xlsx", sheet = i)
  dat1  <- dat1[!is.na(dat1$P_Value),]
  dat1 <- dat1[order(dat1$P_Value),]
  dat1 <- dat1[!duplicated(dat1$Gene_name),]
  dat1.f <- dat1[dat1$P_Value<0.01,]
  dat2 <- read.xlsx("ADT+Means_line2.xlsx", sheet = i+1)
  dat2  <- dat2[!is.na(dat2$P_Value),]
  dat2 <- dat2[order(dat2$P_Value),]
  dat2 <- dat2[!duplicated(dat2$GeneFromPeaks_name),]
  dat2.f <- dat2[dat2$P_Value<0.01,]
  hit <- intersect(dat2.f$GeneFromPeaks_name, dat1.f$Gene_name)
  direct1 <-dat1.f[dat1.f$Gene_name%in%hit,5]
  direct2 <-dat2.f[dat2.f$GeneFromPeaks_name%in%hit,5]
  adt <- dat1$ADT[1]
  out <- data.frame(Gene=hit, ADT=adt, DEgene=direct1*-1, DEatac=direct2*-1)
  out.list[[i]] <- out
  adts <- c(adts,adt)
}

out.list <- out.list[c(1,3,5,7,9,11,13)]
names(out.list) <- adts
write.xlsx(out.list, "Neat-seq_lane2_Hitgenes.xlsx")

en.df.list<-list()
en.list <- list()
for (i in 1:7) {
  genelist <- out.list[[i]][,1]
  if(length(genelist)>100){
    en <- enricher(genelist, TERM2GENE=kegg)
    en.df <- as.data.frame(en)
    en.df.list[[i]] <- en.df
    names(en.df.list)[i] <- adts[i]
    en.list[[i]] <- en
    names(en.list)[i] <- adts[i]
  }
}

en.list <- en.list[-c(1,4)]
en.df.list <- en.df.list[-c(1,4)]
write.xlsx(en.df.list, "Neat-seq_lane2_enrich_kegg.xlsx")

library(enrichplot)
library(ggupset)
library(cowplot)
for (i in 1:5) {
  adt <- names(en.list[i])
  en <- en.list[[i]]
  p1 <- dotplot(en) + ggtitle(adt)
  p2 <- barplot(en) + ggtitle(adt)
  p3 <- upsetplot(en)
  p4 <- emapplot(pairwise_termsim(en))
  tiff(filename = paste0("Neatseq_kegg_lane2_",adt,".tiff"),res = 300, height = 5000, width = 7000, compression = "lzw")
  print(plot_grid(p1,p3,p2,p4, ncol=2))
  dev.off()
}

#check overlap
for (i in 1:7) {
  dat1 <- read.xlsx("Neat-seq_lane1_Hitgenes.xlsx", sheet = i)
  dat2 <- read.xlsx("Neat-seq_lane2_Hitgenes.xlsx", sheet = i)
  print(adts[i])
  print(length(intersect(dat1[,1],dat2[,1])))
}

for (i in 1:5) {
  dat1 <- read.xlsx("Neat-seq_lane1_enrich_kegg.xlsx", sheet = i)
  dat2 <- read.xlsx("Neat-seq_lane2_enrich_kegg.xlsx", sheet = i)
  print(names(en.list)[i])
  print(length(intersect(dat1[,1],dat2[,1])))
}