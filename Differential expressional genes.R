#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")
#install.packages("ggplot2")

library(ggplot2)
library(limma)
library(pheatmap)
logFCfilter=1                   
adjPfilter=0.05                   
expFile="geneMatrix.txt"          
conFile="sample1.txt"             
treatFile="sample2.txt"           
setwd("/Users/lcxexcellent/Desktop/R_files/TAA/dff")      

rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#data=log2(data+1)      
data=normalizeBetweenArrays(data)


sample1=read.table(conFile,sep="\t",header=F,check.names=F)
sample2=read.table(treatFile,sep="\t",header=F,check.names=F)
conData=data[,as.vector(sample1[,1])]
treatData=data[,as.vector(sample2[,1])]
rt=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)


Type=c(rep("con",conNum),rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)


allDiff=topTable(fit2,adjust='fdr',number=100000)
write.table(allDiff,file="GEO_all.xls",sep="\t",quote=F)


diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & P.Value < adjPfilter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file="GEO_diff.xls",sep="\t",quote=F,col.names=F)
write.table(diffSigOut,file="GEO_diff.txt",sep="\t",quote=F,col.names=F)


geneNum=50
diffSig=diffSig[order(as.numeric(as.vector(diffSig$logFC))),]
diffGeneName=as.vector(rownames(diffSig))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(geneNum*2)){
    hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
    hmGene=diffGeneName
}
hmExp=rt[hmGene,]
Type=c(rep("intima-media",conNum),rep("Adventitia",treatNum))
names(Type)=colnames(rt)
Type=as.data.frame(Type)
pdf(file="GEO_heatmap.pdf",height=6,width=6)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("royalblue4", "grey92", "tomato4"))(100),
         cluster_cols =T,
         show_colnames = F,
         border=F,
         scale="row",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=8)
dev.off()


Significant=ifelse((allDiff$P.Value<adjPfilter & abs(allDiff$logFC)>logFCfilter), ifelse(allDiff$logFC>logFCfilter,"Up","Down"), "Not")

p = ggplot(allDiff, aes(logFC, -log10(P.Value)))+
    geom_point(aes(col=Significant))+
    scale_color_manual(values=c("blue", "grey", "red"))+
    geom_hline(yintercept=-log10(0.05),linetype=4)+
    geom_vline(xintercept=c(-1,1),linetype=4)+
    labs(title = " ")+
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p=p+theme_bw()

pdf("GEO_vol.pdf",width=5.5,height=5)
print(p)
dev.off()

