#### linking HiC data with Gene expression
####Inputs####
##Set work directory
setwd('~/todd/PCHiC/')

##Set R object containing HiC data
hic = '~/todd/Data/PCHiC/gothic_interactions_all_xy'

##generated lookup file from Bait_lookup_files.R for which HiC bait fragments represent which gene
bait.file="~/todd/Data/PCHiC/bait_RNAseq_look.txt"
##generated lookup file from Bait_lookup_files.R for which HiC target fragments represent which ROI
agg.look = "~/todd/Data/Output/ROI_classes_within_Target_agg_lookup.txt"

##ESC/TSC Expression data from Cambuli et al
gene.file<-"~/todd/Data/RNA-seq/Cambuli_RNA_logFPKM.txt"


#load R packages
library(ggplot2)


####Script####
##Load processed HiC data
hic.data<- local({
  load(hic)
  stopifnot(length(ls())==1)
  environment()[[ls()]]
})

#Format ESC PCHiC interactions
es<-hic.data$ESC[hic.data$ES.pp==F,]
es$name.bait<-unlist(lapply(strsplit(es$name.bait,split = "-|,"),function(x){x[[1]]}))
es$target.ID = paste(es$chr.target,es$start.target,es$end.target,sep="-")
es$unique<-hic.data$ES.specific[hic.data$ES.pp==F]
es$bait.ID= paste(es$chr.bait,es$start.bait,es$end.bait,sep="-")

#Format TSC PCHiC interactions
ts<-hic.data$TSC[hic.data$TS.pp==F,]
ts$name.bait<-unlist(lapply(strsplit(ts$name.bait,split = "-|,"),function(x){x[[1]]}))
ts$target.ID = paste(ts$chr.target,ts$start.target,ts$end.target,sep="-")
ts$unique<-hic.data$TS.specific[hic.data$TS.pp==F]
ts$bait.ID= paste(ts$chr.bait,ts$start.bait,ts$end.bait,sep="-")


agg = read.delim(agg.look)

for(n in 1:2){

  if(n==1){type=es}else{type=ts}
type.agg<-merge(type,agg,by="target.ID",all.x=T)
head(type.agg)
type.agg<-type.agg[,c(13:21)]
type.agg[is.na(type.agg)]<-0
get.agg=function(df,ID.name){
for(i in 1:(ncol(df)-1)){
  class.agg<-aggregate(df[,(1+i)],list(df[,1]),sum)
  colnames(class.agg)<-c(ID.name,colnames(df)[(1+i)])
  if(i==1){id.agg=class.agg}else{id.agg<-merge(id.agg,class.agg,by=ID.name)}
}
return(id.agg)}
bait.agg=get.agg(type.agg,"Bait_ID")


if(n==1){write.table(bait.agg,"ES_gothic_bait_agg_look.txt",sep="\t",quote = F,col.names = T,row.names = F)}
if(n==2){write.table(bait.agg,"TS_gothic_bait_agg_look.txt",sep="\t",quote = F,col.names = T,row.names = F)}

gene.id<-read.delim(gene.file)[,c(1,6:9)]
colnames(gene.id)[1]<-"Gene_ID"
bait.look=read.delim(bait.file,as.is = T,stringsAsFactors = F,h=F)
colnames(bait.look)=c("Bait_ID","Gene_ID")
gene.mer=merge(gene.id,bait.look,by="Gene_ID")


bait.mer<-merge(gene.mer,bait.agg,by="Bait_ID",all.x=T)
bait.mer[is.na(bait.mer)]=0
gene.info=unique(bait.mer[,c(2:6)])
int.info=bait.mer[,c(2,7:ncol(bait.mer))]
gene.agg=get.agg(int.info,"Gene_ID")
gene.mer=merge(gene.agg,gene.info,by="Gene_ID")

##now I've linked baits to genes perform second aggregate function 

gene.mer$ES_mean<-(gene.mer$ES_E14+gene.mer$ES_J1)/2
gene.mer$TS_mean<-(gene.mer$TS_EGFP+gene.mer$TS_Rs26)/2
gene.mer$Rel_exp<-(gene.mer$ES_mean-gene.mer$TS_mean)

if(n==1){
  gene.mer$ES_REDE_only=gene.mer$ESC_REDE>0&gene.mer$ESC_NEDE==0
  gene.mer$ES_NEDE_only=gene.mer$ESC_REDE==0&gene.mer$ESC_NEDE>0
  gene.mer$ES_nonEnhTE_only=gene.mer$ESC_REDE==0&gene.mer$ESC_NEDE==0&gene.mer$ESC_NonEnh>0
  gene.mer$No_ES_Enh_only=gene.mer$ESC_REDE==0&gene.mer$ESC_NEDE==0
  gene.mer$All=rep(T)
}
if(n==2){
  gene.mer$TS_REDE_only=gene.mer$TSC_REDE>0&gene.mer$TSC_NEDE==0
  gene.mer$TS_NEDE_only=gene.mer$TSC_REDE==0&gene.mer$TSC_NEDE>0
  gene.mer$TS_nonEnhTE_only=gene.mer$TSC_REDE==0&gene.mer$TSC_NEDE==0&gene.mer$TSC_NonEnh>0
  gene.mer$No_TS_Enh_only=gene.mer$TSC_REDE==0&gene.mer$TSC_NEDE==0
  gene.mer$All=rep(T)
}
if(n==1){write.table(gene.mer,"ES_gothic_gene_lookup.txt",sep="\t",quote = F,col.names = T,row.names = F)}
if(n==2){write.table(gene.mer,"TS_gothic_gene_lookup.txt",sep="\t",quote = F,col.names = T,row.names = F)}
}


