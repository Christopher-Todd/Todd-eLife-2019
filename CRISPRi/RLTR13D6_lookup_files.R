#Script takes CasOFFinder results for predicted Cas9 gRNA hits and intersects with target TE coordinates
#CasOFFinder analysis was performed using webtool at http://www.rgenome.net/cas-offinder/

####Inputs####
##Set working directory
setwd("~/todd/CRISPRi/")

#load functions
source('~/todd/Functions/intersectBed.R')
##path to bedtools
bedtools = '~/bedtools2/bin/'

#path to HiC data
hic='~/todd/Data/PCHiC/gothic_interactions_all_xy'

#Coordinates for ROIs
ROI=read.delim("~/todd/Data/Output/ROI_coord.txt")

#Annotations of ESC TEs
es.te=read.delim("~/todd/Data/Annotations/ESC_TEs_mm10.txt")



####Script####
cas.off.results=list.files(path="./",pattern = "Cas_OFFinder")
group1=c("guide1","guide2","guide3","guide4")
group2=c("guide5","guide6","guide7","guide8")



for(i in 1:length(cas.off.results)){
guide.hits=read.delim(cas.off.results[i])
which.guide=strsplit(cas.off.results[i],split = "_Cas")[[1]][1]
perfect.hits=guide.hits[guide.hits$Mismatches==0&guide.hits$Bulge.Size==0,]
start=c()
forward=perfect.hits$Direction=="+"
reverse=perfect.hits$Direction=="-"  
start[forward]=perfect.hits$Position[forward]+1
start[reverse]=perfect.hits$Position[reverse]-20
end=c()
end[forward]=perfect.hits$Position[forward]+20
end[reverse]=perfect.hits$Position[reverse]-1

bed=data.frame(Chr=perfect.hits$Chromosome,start,end,guide=rep(which.guide))
if(i==1){tog.bed=bed}else{tog.bed=rbind(tog.bed,bed)}
}

group1.bed=tog.bed[tog.bed$guide%in%group1,]
group2.bed=tog.bed[tog.bed$guide%in%group2,]


ROI.ES.REDE=ROI[ROI$class=="ESC_REDE",]


RLTR13D6=es.te[es.te$repName=="RLTR13D6",]

##To ensure correct matching of names
REDE.ids=paste(unlist(lapply(strsplit(as.character(ROI.ES.REDE$ID_class),split = "_"),FUN = function(x){x[1]})),"_",sep="")
ids=paste(RLTR13D6$ID,"_",sep="")

RLTR13D6$Enh=ids%in%REDE.ids

RLTR13D6.bed=RLTR13D6[,c(1:3,6)]

RLTR13D6$group1.hits=intersectBed(RLTR13D6.bed,group1.bed,opt.string=" -c ",path.to.bedtools=bedtools)[,5]
RLTR13D6$group2.hits=intersectBed(RLTR13D6.bed,group2.bed,opt.string=" -c ",path.to.bedtools=bedtools)[,5]

write.table(RLTR13D6,"RLTR13D6_CRISPRi_hits_lookup.txt",sep="\t",col.names = T,row.names = F,quote = F)

hic.data<- local({
  load(hic)
  stopifnot(length(ls())==1)
  environment()[[ls()]]
}) 
es<-hic.data$ESC[hic.data$ES.pp==F,]
es$target.ID<-paste(es$chr.target,es$start.target,es$end.target,sep="-")
es.target.bed=es[,c("chr.target","start.target","end.target","target.ID")]

RLTR.target=unique(intersectBed(RLTR13D6.bed,es.target.bed,opt.string=" -wa -wb ",path.to.bedtools=bedtools)[,c(4,8)])
colnames(RLTR.target)=c("TE_ID","target.ID")

##Output with RLTR13D6 elements, the genes they interact with, and number of predicted sgRNA hits
write.table(RLTR.target,"RLTR13D6_target_look.txt",sep="\t",col.names = T,row.names = F,quote = F)
