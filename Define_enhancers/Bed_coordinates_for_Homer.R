#####Inputs####
##set working directory
setwd('~/todd')

##set distance for flanking regions
flank<-1000

##ES and TS Enhancer files from enhancerTEs.R
es.enh = read.delim('Data/Output/ESC_enhancer_TEs.txt',as.is=T,h=F)
ts.enh = read.delim('Data/Output/TSC_enhancer_TEs.txt',as.is=T,h=F)

##ES and TS TE annotations
es.te = read.delim('Data/Annotations/ESC_TEs_mm10.txt',as.is=T)
ts.te = read.delim('Data/Annotations/TSC_TEs_mm10.txt',as.is=T)

####Script####
###Generate lookup tables with all TE elements, denote those which are enhancers, and add flanking regions

colnames(es.enh)=c("chr","start","end","strand","repName","ID","mappability","ID_class")
colnames(ts.enh)=c("chr","start","end","strand","repName","ID","mappability","ID_class")

##get non enh TEs
es.nonenh=es.te[!(es.te$ID %in% es.enh$ID),]
es.nonenh$ID_class=paste(es.nonenh$ID,"noclass",sep="_")
ts.nonenh=ts.te[!(ts.te$ID %in% ts.enh$ID),]
ts.nonenh$ID_class=paste(ts.nonenh$ID,"noclass",sep="_")


##put enh and non enh into one df

tog.df=rbind(es.enh,es.nonenh,ts.enh,ts.nonenh)
tog.df$Plot.start=tog.df$start-flank
tog.df$Plot.end=tog.df$end+flank

##Output
write.table(tog.df,"Data/Output/TE_coordinates_lookup.txt",sep='\t',quote=F,row.names=F,col.names=T)
