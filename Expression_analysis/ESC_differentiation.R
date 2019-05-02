####Inputs####
##Set working directory
setwd('~/todd/Expression_analysis/')

##load R packages
library(ggplot2)

##Load expression data set from Hon et al
hon=read.delim("F:/Todd_et_al_git_hub/Data/RNAseq/ESC_diff_Hon_GSE48519.txt")[,c(1,6:7,18:19)]

##Load HiC ROI-Gene interaction lookup from Pairing_Gothic_Data_to_Expression.R
es.got=read.delim("F:/Todd_et_al_git_hub/Data/Output/ES_gothic_gene_lookup.txt")[,c(1,17:20)]



####Script####
hon[,2:5]=log2(hon[,2:5]+0.01)
hon$WT=rowMeans(hon[,2:3])
hon$d6_WT=rowMeans(hon[,4:5])


es.got=es.got[es.got$Gene_ID%in%hon$gene,]
hon=hon[,c(1,6:7)]
colnames(hon)[1]="Gene_ID"
mer=merge(es.got,hon,by="Gene_ID")

mer$min=mer$WT>(-1)|mer$d6_WT>(-1)
mer$l_FC=mer$d6_WT-mer$WT

rede=data.frame(mer[mer$ES_REDE_only,6:9],class=rep("REDE"),type=rep(T))
nede=data.frame(mer[mer$ES_NEDE_only,6:9],class=rep("NEDE"),type=rep(F))
nonenh=data.frame(mer[mer$ES_nonEnhTE_only,6:9],class=rep("NonEnh"),type=rep(F))
noenh=data.frame(mer[mer$No_ES_Enh_only,6:9],class=rep("NoEnh"),type=rep(F))
all=data.frame(mer[,6:9],class=rep("All"),type=rep(F))
#tog=rbind(rede,nede,nonenh,noenh,all)
tog=rbind(all,nonenh,rede,nede)

## Plot data as in Figure 3 figure supplement 2 A
gg=ggplot(tog[tog$min,],aes(x=class,y=l_FC))+geom_boxplot(outlier.shape = NA)
gg+coord_cartesian(ylim = c(-4,4))+theme_bw()
