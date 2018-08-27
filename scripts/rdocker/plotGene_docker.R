args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(grid)
library(reshape2)
library(ggplot2) 

gene<-args[1]
dataDir<-args[2]
outputDir<-args[3]

load(paste0(dataDir,"/tcga.Kallisto.counts.Robj"))
load(paste0(dataDir,"/tcga.Kallisto.TPM.Robj"))
load(paste0(dataDir,"/metadata.Robj"))
load(paste0(dataDir,"/geneNamecrossref.Robj"))

geneName<- as.character(geneNamecrossref[match(gene,geneNamecrossref$geneID),]$geneName)

textsize<-12
colourScale<-scale_color_manual(values = c("normal"="#93C572","primary(blood)"="#E62020","primary"="#00B7EB","metastatic"="#7851A9","recurrent"="#FCC200"))

gene.TPM<-as.data.frame(t(tcga.Kallisto.TPM[gene,]))
colnames(gene.TPM)<-"TPM"
gene.TPM$ID<-rownames(gene.TPM)
gene.counts<-as.data.frame(t(tcga.Kallisto.counts[gene,]))
colnames(gene.counts)<-"count"
gene.counts$ID<-rownames(gene.counts)
metadata<-metadata[match(gene.TPM$ID,metadata$ID),]
gene.info<-Reduce(function(...) merge(..., all = TRUE, by = c("ID")), list(metadata,gene.TPM,gene.counts))
gene.info$logTPM<-log2(gene.info$TPM+0.01)
gene.info = gene.info %>% filter(type!="CELLC")

write.table(gene.info,file=paste0(outputDir,"/",geneName,".",gene,".txt"),quote = F,row.names = F,sep="\t")

pdf(paste0(outputDir,"/",geneName,".",gene,".pdf"),width = 12,height=4)

p<-ggplot(gene.info, aes(y=logTPM, x=factor(tumor),color=factor(type,levels = c("NT","TP","TM","TR","TB"),labels = c("normal","primary","metastatic","recurrent","primary(blood)")))) +
  geom_boxplot() +
  geom_point(pch = 21, size=1, position = position_jitterdodge()) +
  theme_bw() + ggtitle(paste(geneName,",",gene)) +
  labs(x="",y="log2(TPM)")+
  theme(text = element_text(size=textsize),panel.grid = element_blank(),panel.border = element_blank(),legend.title=element_blank(),axis.line = element_line(color = 'grey'),axis.ticks.x = element_blank(),axis.ticks.y=element_line(color = 'grey'),axis.text.x = element_text(angle=90,hjust = 0.5,vjust = 0.5)) +
  colourScale
print(p)

dev.off()
