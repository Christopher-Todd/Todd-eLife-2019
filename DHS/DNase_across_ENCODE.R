####Inputs####
#set working directory
setwd('~/todd/DHS/')

#location of DNase files as downloaded by download_ENCODE_DHS.R
DNase.files = '~/todd/Data/DNase_peaks/ENCODE_DNase/'

#Bed file
bed = read.delim('~/todd/Data/Output/ROI_coord.txt')


#Classes of interest
class.list<-c("ESC_REDE","ESC_NonEnh","ESC_NEDE","TSC_REDE","TSC_NonEnh","TSC_NEDE")


#Load functions and path to bedtools
source('~/todd/Functions/intersectBed.R')
bedtools='~/bedtools2/bin/'
#Load R packages
library(ComplexHeatmap,lib= '~/Rpackages/')
library(circlize,lib= '~/Rpackages/')

##Set file which contains the order of tissues to be plotted
dnase.order = read.delim('~/todd/Data/DNase_peaks/DNase_order.txt',h=F)



####Script####
DNase.bedfiles = list.files(path=DNase.files,pattern='.bed')


##Getting logical value for each DNase dataset for if peaks overlap with each element in the class
DNase.log=list()
for(n in 1:length(DNase.bedfiles)){
##Interesect ROIs with DNAse peaks with at least 10% of ROI overlapping the DNAse peak
intersect.bed = intersectBed(bed[,c(1,7:8,4)],read.delim(paste(DNase.files,DNase.bedfiles[n],sep=""),h=F),opt.string='-c -f 0.1',path.to.bedtools=bedtools)
DNase.log[[n]]=intersect.bed[,ncol(intersect.bed)]>0}

##putting this logical value into a matrix and setting col and row names
num = lapply(DNase.log, as.numeric)
mat = matrix(unlist(num),ncol=length(num))
rownames(mat)<-bed[,4]
colnames(mat)<-DNase.bedfiles



###Pooling tissue replicates into single columns###

col.split<-strsplit(colnames(mat),split = "_ENCFF")
name<-unlist(lapply(X = col.split,FUN = function(x){x[[1]]}))
colnames(mat)<-name

#taking mean of logical values to get measure of in how many replicates for each tissue is this element associated with a DNAse peak 
agg.list<-list()
for(i in 1:nrow(mat)){
  mat.row<-mat[i,]
  agg<-aggregate(x = unlist(mat.row),list(colnames(mat)),mean)
  agg.list[[i]]<-t(agg$x)}
tog<-do.call(rbind,agg.list)

#replacing the lost col and row names again
agg<-aggregate(unlist(mat[1,]),list(colnames(mat)),mean)
colnames(tog)<-agg$Group.1
rownames(tog)<-rownames(mat)


##Function for drawing the heatmap
draw.heatmap = function(input.mat,name,width=4,height=10,col) {
  ##saturation point for the heatmap
  sat=0.5
  ##make heatmap object
  colramp = colorRamp2(c(0,sat),c('white',col))
  hm.plot= Heatmap(input.mat,cluster_rows=T,cluster_columns=F,
                   show_row_names=F,show_row_dend = FALSE,show_column_names=F,show_heatmap_legend=F,col=colramp)
  jpeg(paste(name,'heatmap.jpg',sep='_'),w=width,h=height,units='in',res=600)
  draw(hm.plot)
  dev.off()
  
}

for(i in 1:length(class.list)){
  class.bed = bed[bed[,9] %in% class.list[i],1:4]
  class.mat = tog[row.names(tog) %in% class.bed[,4],]
  class.mat2 = class.mat[,match(dnase.order[,1],colnames(class.mat))]
  class.mat3 = cbind(rownames(class.mat2),class.mat2)
  colnames(class.mat3)[1]<-"ID"
  ##output the dataframe of values
  write.table(class.mat3,paste(class.list[i],"_DNAse.txt",sep=""),quote=F,col.names = T,row.names = F,sep="\t")
  ##output the heatmap as a jpeg
  draw.heatmap(input.mat = class.mat2,name = paste(class.list[i],"_DNase",sep=""),col='green')
  
  #get the order of elements as plotted within the heatmap
  x=Heatmap(class.mat2,cluster_rows=T,cluster_columns=F,
            show_row_names=F,show_row_dend = FALSE,show_column_names=F,show_heatmap_legend=F)
  new.mat=class.mat3[rev(unlist(row_order(x))),]
  write.table(new.mat,paste(class.list[i],"_new_mat_order.txt",sep=""),sep="\t",quote=F,col.names = T,row.names = F)
  
}
