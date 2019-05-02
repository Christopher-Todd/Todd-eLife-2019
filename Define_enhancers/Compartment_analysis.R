####Inputs####
##set working directory
setwd("~/todd/")

#load functions and path to bedtools
source("~/todd/Functions/intersectBed.R")
bedtools="~/bedtools2/bin/"

#load R packages
library(ggplot2)

#File containing coordinates of A and B compartments within ESCs
compartments=read.delim("~/todd/Data/Annotations/A_B_compartments.txt",h=F)

#Coordinates of ROIs
ROI=read.delim("~/todd/Data/Output/ROI_coord.txt")

####Script####
ROI.bed=ROI[,c(1,7:8,4)]

ROI.com=intersectBed(ROI.bed,compartments,opt.string=" -wa -wb",path.to.bedtools=bedtools)[,c(4,8)]
colnames(ROI.com)=c("ID_class","compartment")

ROI.com.A=ROI.com[ROI.com$compartment=="A",]
ROI.com.B=ROI.com[ROI.com$compartment=="B",]

class.list=unique(ROI$class)

##Get percentage of elements within each type of compartment for each class
df=list()
for(i in 1:length(class.list)){
  class=class.list[i]
  ROI.sub=ROI[ROI$class==class,]
  n.elements=nrow(ROI.sub)
  n.A=sum(ROI.sub$ID_class%in%ROI.com.A$ID_class)
  n.B=sum(ROI.sub$ID_class%in%ROI.com.B$ID_class)
  percent.A=(n.A/n.elements)*100
  percent.B=(n.B/n.elements)*100
  df[[i]]=data.frame(class,percent.A,percent.B)
}
df.tog=do.call(rbind,df)

class.to.plot=c("ESC_REDE","ESC_NonEnh","ESC_NEDE")
plot.df=df.tog[df.tog$class%in%class.to.plot,]


##Plot percentages as seen in Figure 3 figure supplement 1 B
ggplot(plot.df,aes(x=class,y=percent.A,colour=class))+geom_bar()
ggplot(plot.df,aes(x=class,y=percent.B,colour=class))+geom_bar()