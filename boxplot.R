
options(stringsAsFactors=F)
library(limma)
library(ggpubr)
library(reshape2)

riskFile="geoRisk.txt"       
scoreFile="ssgseaOut.AAdv.txt"        
setwd("/Users/lcxexcellent/Desktop/R_files/TAA/Immune/31.scoreCor")               
data=read.table(scoreFile,sep="\t",header=T,check.names=F,row.names=1)     
#group=sapply(strsplit(colnames(data),"\\-"),"[",4)
#data=data[,group==0]
#colnames(data)=gsub("(.*?)\\-(.*?)\\-.*","\\2",colnames(data))
data=avereps(t(data))


risk=read.table(riskFile,header=T,sep="\t",row.names=1,check.names=F)

sameSample=intersect(row.names(data),row.names(risk))
data=data[sameSample,]
risk=risk[sameSample,]
rt=cbind(data,risk[,c("riskScore","risk")])
rt=rt[,-(ncol(rt)-1)]

immCell=c("aDCs","B_cells","CD8+_T_cells","DCs","iDCs","Macrophages",
          "Mast_cells","Neutrophils","NK_cells","pDCs","T_helper_cells",
          "Tfh","Th1_cells","Th2_cells","TIL","Treg")
rt1=rt[,c(immCell,"risk")]
data=melt(rt1,id.vars=c("risk"))
colnames(data)=c("Risk","Type","Score")
data$Risk=factor(data$Risk, levels=c("low","high"))
p=ggboxplot(data, x="Type", y="Score", color = "Risk",
     ylab="Score",add = "none",xlab="",palette = c("blue","red") )
p=p+rotate_x_text(50)
pdf(file="immCell.boxplot.cgga.pdf",width=7,height=6)         
p+stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()


immFunction=c("APC_co_inhibition","APC_co_stimulation","CCR",
          "Check-point","Cytolytic_activity","HLA","Inflammation-promoting",
          "MHC_class_I","Parainflammation","T_cell_co-inhibition",
          "T_cell_co-stimulation","Type_I_IFN_Reponse","Type_II_IFN_Reponse")
rt1=rt[,c(immFunction,"risk")]
data=melt(rt1,id.vars=c("risk"))
colnames(data)=c("Risk","Type","Score")
data$Risk=factor(data$Risk, levels=c("low","high"))
p=ggboxplot(data, x="Type", y="Score", color = "Risk",
     ylab="Score",add = "none",xlab="",palette = c("blue","red") )
p=p+rotate_x_text(50)
pdf(file="immFunction.boxplot.cgga.pdf",width=7,height=6)          
p+stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()


