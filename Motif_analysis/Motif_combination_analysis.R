####Inputs####
##set working directory
setwd('~/todd/Motif_analysis')

##Set files for Enhancer and Non-Enhancer TE coordinates
es.enh=read.delim("~/todd/Data/Output/ESC_enhancer_TEs.txt",h=F)
ts.enh=read.delim("~/todd/Data/Output/TSC_enhancer_TEs.txt",h=F)
non.enh=read.delim("~/todd/Data/Output/NonEnhTE_candidates.txt")
#Set files for ES and TS TE coordinates
es.te=read.delim("~/todd/Data/Annotations/ESC_TEs_mm10.txt")
ts.te=read.delim("~/todd/Data/Annotations/TSC_TEs_mm10.txt")

#Load R packages
library(RColorBrewer)
library(scales)
library(ggplot2)
library(ggrepel)

#List the TE classes of interest
class.list=c("RLTR9E","RLTR13D6","RLTR13D5","RLTR13B1")

#Set the minimum size cut off (as proportion of longest Enh copy length)
#0.6 = 60%
cutoff=0.6

#FIMO analysis result directories from Clustal_and_Motif_analysis.R
hocomoco.dir="~/todd/Motif_analysis/hocomoco/"
uniprobe.dir="~/todd/Motif_analysis/uniprobe/"


####Script#####

te.tog=rbind(es.te,ts.te)
#get size for each element
te.tog$size=te.tog$end-te.tog$start
te.tog=te.tog[,c("repName","ID","size")]

for(i in 1:length(class.list)){
  class=class.list[i]
  
#load FIMO results
enh.h=read.delim(paste(hocomoco.dir,class,"_enh_fimo.txt",sep=""))
non.h=read.delim(paste(hocomoco.dir,class,"_non_fimo.txt",sep=""))
enh.u=read.delim(paste(uniprobe.dir,class,"_enh_fimo.txt",sep=""))
non.u=read.delim(paste(uniprobe.dir,class,"_non_fimo.txt",sep=""))

#set IDs for desired motifs within HOCOMOCO and Uniprobe datasets
id.oct=c("PO5F1_MOUSE.H11MO.0.A","PO5F1_MOUSE.H11MO.1.A")
id.nan=c("NANOG_MOUSE.H11MO.0.A","NANOG_MOUSE.H11MO.1.A")
id.sox=c("SOX2_MOUSE.H11MO.0.A","SOX2_MOUSE.H11MO.1.A")
id.eom=c("UP00068_1","UP00068_2")
id.elf=c("UP00409_1","ELF5_MOUSE.H11MO.0.A","ELF5_MOUSE.H11MO.1.A")
id.cdx=c("CDX2_MOUSE.H11MO.0.A","CDX2_MOUSE.H11MO.0.A")

mof.tog=rbind(enh.h[,c("X..motif_id","sequence_name")],enh.u[,c("X..motif_id","sequence_name")],non.h[,c("X..motif_id","sequence_name")],non.u[,c("X..motif_id","sequence_name")])

#Generate vectors of element names which were found to contain selected motifs
mot.oct=mof.tog$sequence_name[mof.tog$X..motif_id%in%id.oct]
mot.nan=mof.tog$sequence_name[mof.tog$X..motif_id%in%id.nan]
mot.sox=mof.tog$sequence_name[mof.tog$X..motif_id%in%id.sox]
mot.eom=mof.tog$sequence_name[mof.tog$X..motif_id%in%id.eom]
mot.elf=mof.tog$sequence_name[mof.tog$X..motif_id%in%id.elf]
mot.cdx=mof.tog$sequence_name[mof.tog$X..motif_id%in%id.cdx]

enh.tog=rbind(es.enh,ts.enh)
enh.tog$size=enh.tog[,3]-enh.tog[,2]


class.enh=enh.tog[enh.tog$V5==class,c(5,8)]

class.non=non.enh[non.enh$repName==class,c(5,8)]
colnames(class.enh)=c("repName","ID_class")
class.enh$type=rep("Enh")
class.non$type=rep("NonEnh")

class.tog=rbind(class.enh,class.non)

#check to see if element contains motif
class.tog$oct=class.tog$ID_class%in%mot.oct
class.tog$nan=class.tog$ID_class%in%mot.nan
class.tog$sox=class.tog$ID_class%in%mot.sox
class.tog$eom=class.tog$ID_class%in%mot.eom
class.tog$elf=class.tog$ID_class%in%mot.elf
class.tog$cdx=class.tog$ID_class%in%mot.cdx

class.tog$oct.nan=class.tog$oct&class.tog$nan
class.tog$oct.sox=class.tog$oct&class.tog$sox
class.tog$nan.sox=class.tog$nan&class.tog$sox
class.tog$oct.sox.nan=class.tog$oct&class.tog$sox&class.tog$nan

class.tog$eom.elf=class.tog$eom&class.tog$elf
class.tog$eom.cdx=class.tog$eom&class.tog$cdx
class.tog$elf.cdx=class.tog$elf&class.tog$cdx
class.tog$eom.elf.cdx=class.tog$eom&class.tog$elf&class.tog$cdx



class.te.tog=te.tog[grep(class,te.tog$repName),]

#Perform size based cutoff
class.size=max(class.te.tog$size)
class.cutoff=class.size*cutoff
cutoff.ids=paste(as.character(class.te.tog$ID[class.te.tog$size>class.cutoff]),"_",sep="")

class.tog=class.tog[grep(paste(cutoff.ids,collapse = "|"),class.tog$ID_class),]


##Calculate the percentage of elements which contain motif
df=data.frame()
colnames(class.tog)
tf.cols=colnames(class.tog)[4:ncol(class.tog)]
for(i in 1:length(tf.cols)){
col=tf.cols[i]
enh.num=sum(class.tog$type=="Enh")
non.num=sum(class.tog$type=="NonEnh")
enh.count=sum(class.tog[class.tog$type=="Enh",col])
non.count=sum(class.tog[class.tog$type=="NonEnh",col])
enh.per=(enh.count/enh.num)*100
non.per=(non.count/non.num)*100
row=data.frame(motif=col,enh.per,non.per)
df=rbind(df,row)
}


##Plot Motif percentages
lim.upper=round((max(df[c(1:3,7:10),2:3])*1.2)+5,-1)

es.gg=ggplot(df[c(1:3,7:10),],aes(x=enh.per,y=non.per))+geom_point()+
  geom_text(aes(label=motif),size=3,nudge_x = 1,angle=30)+
  geom_text_repel(aes(label=motif))+
  ylab("% NonEnh")+xlab("% Enh")+
  ylim(c(0,lim.upper))+xlim(c(0,lim.upper))+geom_abline(intercept = 0,slope = 1)+
  ggtitle(class)+theme_bw()
pdf(paste(class,"ES_Motif_comb.pdf"),width = 6,height = 6)
plot(es.gg)
dev.off()

lim.upper=round((max(df[c(3:6,11:14),2:3])*1.2)+5,-1)

ts.gg=ggplot(df[c(3:6,11:14),],aes(x=enh.per,y=non.per))+geom_point()+
  geom_text(aes(label=motif),size=3,nudge_x = 1,angle=30)+
  geom_text_repel(aes(label=motif))+
  ylab("% NonEnh")+xlab("% Enh")+
  ylim(c(0,lim.upper))+xlim(c(0,lim.upper))+geom_abline(intercept = 0,slope = 1)+ggtitle(class)
pdf(paste(class,"TS_Motif_comb.pdf"),width = 6,height = 6)
plot(ts.gg)
dev.off()
}
