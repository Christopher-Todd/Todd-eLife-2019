####Inputs####
##Set working directory
setwd("~/todd")

##Output files from Pairing_Gothic_Data_to_Expression.R
ES.gene.look=read.delim("ES_gothic_gene_lookup.txt")
TS.gene.look=read.delim("TS_gothic_gene_lookup.txt")

#load R packages
library(ggplot2)



####Script####
##Plotting ESC and TSC gene expression
ES.rede=ES.gene.look[ES.gene.look$ES_REDE_only,]
ES.rede$class=rep("ES_REDE_only")
ES.nede=ES.gene.look[ES.gene.look$ES_NEDE_only,]
ES.nede$class=rep("ES_NEDE_only")
ES.nonEnhTE=ES.gene.look[ES.gene.look$ES_nonEnhTE_only,]
ES.nonEnhTE$class=rep("ES_nonEnhTE")
ES.noEnh=ES.gene.look[ES.gene.look$No_ES_Enh_only,]
ES.noEnh$class=rep("ES_noEnh")
ES.all=ES.gene.look[ES.gene.look$All,]
ES.all$class=rep("All")
ES.tog=rbind(ES.rede,ES.nede,ES.nonEnhTE,ES.noEnh,ES.all)


es.exp=ggplot(ES.tog,aes(x=class,y=ES_mean))+geom_boxplot()+geom_violin()
plot(es.exp)
rel.exp=ggplot(ES.tog,aes(x=class,y=Rel_exp))+geom_boxplot()
plot(rel.exp)



TS.rede=TS.gene.look[TS.gene.look$TS_REDE_only,]
TS.rede$class=rep("TS_REDE_only")
TS.nede=TS.gene.look[TS.gene.look$TS_NEDE_only,]
TS.nede$class=rep("TS_NEDE_only")
TS.nonEnhTE=TS.gene.look[TS.gene.look$TS_nonEnhTE_only,]
TS.nonEnhTE$class=rep("TS_nonEnhTE")
TS.noEnh=TS.gene.look[TS.gene.look$No_TS_Enh_only,]
TS.noEnh$class=rep("TS_noEnh")
TS.all=TS.gene.look[TS.gene.look$All,]
TS.all$class=rep("All")
TS.tog=rbind(TS.rede,TS.nede,TS.nonEnhTE,TS.noEnh,TS.all)

TS.exp=ggplot(TS.tog,aes(x=class,y=TS_mean))+geom_boxplot()+geom_violin()
plot(TS.exp)
rel.exp=ggplot(TS.tog,aes(x=class,y=Rel_exp))+geom_boxplot()
plot(rel.exp)

