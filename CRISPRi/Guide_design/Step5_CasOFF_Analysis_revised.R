####  User Inputs  ####

#Set working directory
setwd("~/Guide_design")

##import functions
source('~/todd/Functions/intersectBed.R')
#path to bedtools
bedtools = '~/bedtools2/bin/'

#path to repeatmasker dataset (BED format)
repmask=read.delim("/data/Blizard-BrancoLab/Chris/Todd_et_al/Data/Annotations/RepeatMasker.txt",h=F)
#path to exon dataset (BED format)
exon.coord=read.delim("NCBI_RefSeq_mm10_exon.txt",h=F)

#set name of Target repeat family as Identified in Repeatmasker dataset
target.repeat.name="RLTR13D5"
#list of names for similar repeat families
similar.repeat.names=c("RLTR13")

#path to CasOFFinder output
casoff.output.name = "cas_off_output.txt"


#set output file name
output.name = "Repeat_hits.txt"

####  Script  ####

off=read.delim(casoff.output.name,h=F)

#get off-target coordinates using start and orientation
off.start=c()
off.end=c()
for(i in 1:nrow(off)){
  if(off[i,5]=="+"){
    off.start[i]=off[i,3]
    off.end[i]=off[i,3]+20
  }else{
    off.start[i]=off[i,3]-20
    off.end[i]=off[i,3]
  }
}

guide.seq=substr(as.character(off[,1]),1,20)

#formated bed coords for all CasOFF target hits
off.df=data.frame(chr=off[,2],start=off.start,end=off.end,name=guide.seq,mismatches=off[,6],strand=off[,5])
off.df$ID=paste("Hit_",seq(1:nrow(off.df)),sep="")

target.reps=grepl(target.repeat.name,repmask$V4)
similar.reps=grepl(paste(similar.repeat.names,collapse = "|"),repmask$V4)
target.rep.coord=repmask[target.reps,]
similar.rep.coord=repmask[(similar.reps==T&target.reps==F),]
non.target.rep=repmask[(similar.reps==F&target.reps==F),]

target.check=intersectBed(off.df[,c(1:3,7)],target.rep.coord,opt.string="-u",path.to.bedtools=bedtools)[,4]
similar.check=intersectBed(off.df[,c(1:3,7)],similar.rep.coord,opt.string="-u",path.to.bedtools=bedtools)[,4]
non.target.rep.check=intersectBed(off.df[,c(1:3,7)],non.target.rep,opt.string="-u",path.to.bedtools=bedtools)[,4]
non.rep.check=intersectBed(off.df[,c(1:3,7)],repmask,opt.string="-v",path.to.bedtools=bedtools)[,4]
exon.check=intersectBed(off.df[,c(1:3,7)],exon.coord,opt.string="-u",path.to.bedtools=bedtools)[,4]

off.df$target.hit=off.df$ID%in%target.check
off.df$similar.hit=off.df$ID%in%similar.check
off.df$non.target.rep.hit=off.df$ID%in%non.target.rep.check
off.df$non.rep.hit=off.df$ID%in%non.rep.check
off.df$exon.hit=off.df$ID%in%exon.check

unique.guides=unique(off.df$name)

for(i in 1:length(unique.guides)){
  guide=unique.guides[i]
  df.sub=off.df[off.df$name==guide&off.df$mismatches==0,]
  
  total.hits=nrow(df.sub)
  correct.hits=sum(df.sub$target.hit)
  related.hits=sum(df.sub$similar.hit)
  non.target.rep.hits=sum(df.sub$non.target.rep.hit)
  non.rep.hits=sum(df.sub$non.rep.hit)
  exon.hits=sum(df.sub$exon.hit)
  
  df.row=data.frame(guide,total.hits,correct.hits,related.hits,non.target.rep.hits,non.rep.hits,exon.hits)
  if(i==1){hit.df=df.row}else{hit.df=rbind(hit.df,df.row)}
  print(paste(i,"/",length(unique.guides)))
}

write.table(hit.df,output.name,sep="\t",col.names = T,row.names = F,quote = F)