#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05        
qvalueFilter=0.05         

setwd("/Users/lcxexcellent/Desktop/R_files/TAA/Enrichment")         
rt=read.table("id.txt",sep="\t",header=T,check.names=F)    
rt=rt[is.na(rt[,"entrezID"])==F,]                          
colnames(rt)[1]="Gene"
gene=rt$entrezID

colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$Gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG,file="KEGG.txt",sep="\t",quote=F,row.names = F)

showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

pdf(file="KEGG_barplot.pdf",width = 9,height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)
dev.off()

pdf(file="KEGG_bubble.pdf",width = 9,height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)
dev.off()




