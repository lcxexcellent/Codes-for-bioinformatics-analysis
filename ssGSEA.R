
#install.packages("BiocManager")
#BiocManager::install("GSVA")
#BiocManager::install("GSEABase")


library(GSVA)
library(limma)
library(GSEABase)

inputFile="geneMatrix_AAdv.txt"    
gmtFile="immune.gmt"    
setwd("/Users/lcxexcellent/Desktop/R_files/TAA/Immune/30.ssGSEA")        


rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]
geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())

ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}

ssgseaOut=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
write.table(ssgseaOut,file="ssgseaOut.tcga.txt",sep="\t",quote=F,col.names=F)
