####Inputs####
##Set working directory
setwd('~/todd/DHS/')

#location of ATAC-seq files from 
#Smith et al PMID:28959968
#Wu et al PMID:27309802
Smith.DNase.dir = '~/todd/Data/DNase_peaks/Smith_et_al_mm10/'
Wu.DNase.dir = '~/todd/Data/DNase_peaks/Wu_et_al_mm10/'

#Bed file of ROIs
bed = read.delim('~/todd/Data/Output/ROI_coord.txt')

#Classes of interest
class.list<-c("ESC_REDE","ESC_NonEnh","ESC_NEDE","TSC_REDE","TSC_NonEnh","TSC_NEDE")

#Load functions and path to bedtools
source('~/todd/Functions/intersectBed.R')
bedtools='~/bedtools2/bin/'


####Script####

Smith.bedfiles = list.files(path=Smith.DNase.files,pattern='.bed')
Wu.bedfiles = list.files(path=Wu.DNase.files,pattern='.bed')

Smith.bedfiles=paste(Smith.DNase.dir,Smith.bedfiles,sep="")
Wu.bedfiles=paste(Wu.DNase.dir,Wu.bedfiles,sep="")
DNAse.bedfiles=c(Smith.bedfiles,Wu.bedfiles)



##Getting logical value for each DNase dataset for if peaks overlap with each element in the class
DNase.log=list()
for(n in 1:length(DNAse.bedfiles)){
  intersect.bed = intersectBed(bed[,c(1,7:8,4)],read.delim(DNAse.bedfiles[n],h=F),opt.string='-c -F 0.1',path.to.bedtools=bedtools)
  DNase.log[[n]]=intersect.bed[,ncol(intersect.bed)]>0}

##putting this logical value into a matrix and setting col and row names
num = lapply(DNase.log, as.numeric)
mat = matrix(unlist(num),ncol=length(num))
rownames(mat)<-bed[,4]
colnames(mat)<-DNAse.bedfiles

##Getting numbers of each class which overlap DNAse-sensitive sites
class.numbers=list()
for(i in 1:length(class.list)){
  class=class.list[i]
  IDs=bed[bed$class==class,"ID_class"]
  mat.sub=mat[IDs,]
  col.sums=c(colSums(mat.sub))
  class.numbers[[i]]=col.sums
}
output.df=as.data.frame(do.call(rbind,class.numbers))
rownames(output.df)=class.list

##Output
write.table(output.df,"DNAse_preimplantation_output.txt",sep="\t",col.names = T,row.names = T,quote = F)