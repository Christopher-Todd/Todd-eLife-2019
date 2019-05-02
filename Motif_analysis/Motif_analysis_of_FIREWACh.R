####Inputs####
#set working directory
setwd('~/todd/Motif_analysis')

#load functions and set path to bedtools
source("~/todd/Functions/intersectBed.R")
bedtools="~/bedtools2/bin/"

##Files containing FIREWACh +ve results and full analysed library which overlap TE elements
FIRE=read.delim("~/todd/Data/Annotations/FIREWACh_mm10.txt",h=F)
FIRE.lib=read.delim("~/todd/Data/Annotations/library_NFRs_mm10.txt",h=F)

##File containing coordinates of ROIs
ROI=read.delim("~/todd/Data/Output/ROI_coord.txt")

##ESC TF Motif database
motif.database = "~/motif_databases/MOUSE/chen2008.meme"

##Path to required MEME and Homer tools
path.to.ame = "~/meme/bin/ame"
path.to.fimo = "~/meme/bin/fimo"
path.to.homer = "~/Homer/bin/homerTools"
path.to.homerformat.genome = "~/Homer/data/genomes/mm10/"


####Script####
ROI=ROI[!grepl(pattern = "NonEnh",x = ROI$class),]
ROI.bed=ROI[,c(1,7:8,4)]


##Check which ROIs were positive for FIREWACh enhancer activity
FIRE.TE=intersectBed(ROI.bed,FIRE,opt.string=" -u",path.to.bedtools=bedtools)[,4]
FIRE.lib.TE=intersectBed(ROI.bed,FIRE.lib,opt.string=" -u",path.to.bedtools=bedtools)[,4]
FIRE.coords=unique(ROI[ROI$ID_class%in%FIRE.TE,c(1,7:8,4:6)])
FIRElib.coords=unique(ROI[ROI$ID_class%in%FIRE.lib.TE,c(1,7:8,4:6)])





##get sequences
temp.enh.file<-tempfile()
write.table(FIRE.coords,file=temp.enh.file,quote=F,sep='\t',col.names=F,row.names=F)
command = paste(path.to.homer,' extract ',temp.enh.file,' ',path.to.homerformat.genome,' -fa > ',"ESC_TE_FIRE.fa",sep='')
try(system(command))
unlink(temp.enh.file)
temp.non.file=tempfile()
write.table(FIRElib.coords,file=temp.non.file,quote=F,sep='\t',col.names=F,row.names=F)
command = paste(path.to.homer,' extract ',temp.non.file,' ',path.to.homerformat.genome,' -fa > ',"ESC_TE_FIRElib.fa",sep='')
try(system(command))
unlink(temp.non.file)

##FIMO
for(x in c("ESC_TE_FIRE","ESC_TE_FIRElib")){
    command = paste(path.to.fimo,' --o ',x,'_fimo ',motif.database,' ',x,'.fa',sep='')
    try(system(command))
    try(system(paste('mv ',x,'_fimo/fimo.txt ./',x,'_fimo.txt',sep='')))
    try(system(paste('rm -r ',x,'_fimo/',sep='')))
}
##AME
files=c("ESC_TE_FIRE","ESC_TE_nonFIRE")
  command = paste(path.to.ame,' --o ESC_TE_ame --control ',files[2],'.fa  ' ,files[1],'.fa',' ',motif.database,sep='')
  try(system(command))
  ##move and rename file
  try(system(paste('mv ',"ESC_TE",'_ame/ame.txt ./',"ESC_TE",'_ame.txt',sep='')))
  try(system(paste('rm -r ',"ESC_TE",'_ame/',sep='')))
  command = paste(path.to.ame,' --o ESC_TE_ame_nocontrol ' ,files[1],'.fa',' ',motif.database,sep='')
  try(system(command))
  ##move and rename file
  try(system(paste('mv ',"ESC_TE",'_ame_nocontrol/ame.txt ./',"ESC_TE",'_ame_nocontrol.txt',sep='')))
  try(system(paste('rm -r ',"ESC_TE",'_ame_nocontrol/',sep='')))




