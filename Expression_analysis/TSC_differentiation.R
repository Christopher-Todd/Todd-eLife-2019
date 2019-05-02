####Inputs####
##set working directory
setwd('~/todd/Expression_analysis/')

##load R packages
library(ggplot2)

##Load expression data set from Latos et al
rna.seq=read.delim('F:/Todd_et_al_git_hub/Data/RNAseq/Latos_RNA_logFPKM.txt')

##Load HiC ROI-Gene interaction lookup from Pairing_Gothic_Data_to_Expression.R
ts.got=read.delim("F:/Todd_et_al_git_hub/Data/Output/TS_gothic_gene_lookup.txt")[,c(1,17:20)]

####Script####
rna.seq$TSC_WT_mean=rowMeans(rna.seq[,6:10])
rna.seq$diffTSC_d3_mean=rowMeans(rna.seq[,14:16])
rna.seq=rna.seq[,c(1,17:18)]

hist(rna.seq$TSC_WT_mean)
min=rna.seq$TSC_WT_mean>(-1)|rna.seq$diffTSC_d3_mean>(-1)
rna.seq=rna.seq[min,]
colnames(rna.seq)[1]="Gene_ID"


mer=merge(ts.got,rna.seq,by="Gene_ID")
mer$l_FC=mer$diffTSC_d3_mean-mer$TSC_WT_mean

rede=data.frame(mer[mer$TS_REDE_only,6:8],class=rep("REDE"),type=rep(T))
nede=data.frame(mer[mer$TS_NEDE_only,6:8],class=rep("NEDE"),type=rep(F))
nonenh=data.frame(mer[mer$TS_nonEnhTE_only,6:8],class=rep("NonEnh"),type=rep(F))
noenh=data.frame(mer[mer$No_TS_Enh_only,6:8],class=rep("NoEnh"),type=rep(F))
all=data.frame(mer[,6:8],class=rep("All"),type=rep(F))
tog=rbind(rede,nede,nonenh,noenh,all)
tog=rbind(all,nonenh,rede,nede)

## Plot data as in Figure 3 figure supplement 2 B
gg=ggplot(tog,aes(x=class,y=l_FC))+geom_boxplot(outlier.shape = NA)
gg+coord_cartesian(ylim = c(-3,3))+theme_bw()
plot(gg)


aov.df=aov(l_FC~class,tog)
summary(aov.df)
print(summary(glht(aov.df, linfct = mcp(class = "Tukey"))))
