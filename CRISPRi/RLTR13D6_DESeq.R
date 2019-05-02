####Inputs####
##Set working directory
setwd('~/todd/CRISPRi/')

#load R packages
library('DESeq2')
library(ggplot2)


#Quantified RNA-seq data collected from CRISPRi experiment
reads = read.table('RLTR13D6_RNA_counts.txt', header = T, row.names = 1)


####Script####
##Prepare the data
colData = DataFrame(sample=factor(rep(c('empty','sg5.8','sg1.4'),each=3)),
 replicate=factor(rep(c('r4d8','r4d5','r3d5'),3)))
colData1 = DataFrame(sample=factor(rep(c('empty','sgrna'),each=3)),
 replicate=factor(rep(c('r4d8','r4d5','r3d5'),2)))


##Expression values (VST)

dds = DESeqDataSetFromMatrix(reads,colData,design=~sample+replicate)
vsd = varianceStabilizingTransformation(dds)
expr = assay(vsd)
##Output normalised expression values
write.table(expr,"RLTR13D6_RNA_vsd.txt",col.names = T,row.names = T,quote = F)



##Clustering

d = dist(t(expr))
plot(hclust(d))


##PCA

pca = prcomp(expr)
plot(pca$rotation[,1],pca$rotation[,2],
 col=rep(c('blue','grey','tomato'),each=3),pch=19)
text(pca$rotation[,1],pca$rotation[,2],rownames(pca$rotation),pos=2,cex=0.5)


##Get DE genes

dds1 = DESeqDataSetFromMatrix(reads[,1:6],colData1,design=~sample+replicate)
dds1 = DESeq(dds1)
res1 = results(dds1,contrast=c('sample','sgrna','empty'))
sg5.8 = as.data.frame(subset(res1,padj<0.05))
sg5.8 = sg5.8[order(sg5.8$log2FoldChange,decreasing=T),]

dds2 = DESeqDataSetFromMatrix(reads[,c(1:3,7:9)],colData1,design=~sample+replicate)
dds2 = DESeq(dds2)
res2 = results(dds2,contrast=c('sample','sgrna','empty'))
sg1.4 = as.data.frame(subset(res2,padj<0.05))
sg1.4 = sg1.4[order(sg1.4$log2FoldChange,decreasing=T),]


##Check expression

plot.gene = function(gene,log=FALSE) {
	id = rownames(expr)==gene
	sample = colData$sample
	
	if (log) {
		data = expr[id,]
	} else {
		data = 2^expr[id,]
	}
	
	av = tapply(data,sample,mean)
	h=barplot(av,main=gene,ylab='Relative expression',col='lightblue',ylim=c(0,max(data)*1.05),las=3)
	for (i in 1:nlevels(sample)) {
		samp = sample==levels(sample)[i]
		points(rep(h[i],sum(samp)),data[samp],pch=19,cex=0.5)
	}
	return(data)
}

plot.gene('Tdrd12',log = F)
plot.gene('Hook3',log=T)
plot.gene('Spp1',log=T)


expr$WT.exp=rowMeans(expr[,2:4])
expr$Set1.exp=rowMeans(expr[,8:10])
expr$Set2.exp=rowMeans(expr[,5:7])

expr$Set1.exp.lfc=expr$Set1.exp-expr$WT.exp
expr$Set2.exp.lfc=expr$Set2.exp-expr$WT.exp

set1.sig.genes=c("Tdrd12","Spp1")
set2.sig.genes=c("Tdrd12","Hook3")

##Plot expression with CRISPRi experiment as seen in Figure 5 C-D

ggplot(expr,aes(x=WT.exp,y=Set1.exp.lfc))+geom_point()+
  geom_point(data = expr[expr$X%in%set1.sig.genes,],aes(x=WT.exp,y=Set1.exp.lfc),colour="red")+
  geom_text(data = expr[expr$X%in%set1.sig.genes,],aes(x=WT.exp,y=Set1.exp.lfc,label=X),colour="red")+
  ylim(c(-1,1))

ggplot(expr,aes(x=WT.exp,y=Set2.exp.lfc))+geom_point()+
  geom_point(data = expr[expr$X%in%set2.sig.genes,],aes(x=WT.exp,y=Set2.exp.lfc),colour="red")+
  geom_text(data = expr[expr$X%in%set2.sig.genes,],aes(x=WT.exp,y=Set2.exp.lfc,label=X),colour="red")+
  ylim(c(-1,1))

