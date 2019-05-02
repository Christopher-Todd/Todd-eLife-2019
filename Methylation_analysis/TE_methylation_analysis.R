####Inputs####
#Set working directory
setwd("~/todd/Methylation/")

#load R packages
library(ggplot2)

##Methylation data processed by Seqmonk
#min CpGs = 3
#min read count = 10
#liftOver to mm10 assembly
meth=read.delim("./Hon_TE_BSpipeline_mm10.txt")


####Script####

##Plots as seen in Figure 2 D
gg=ggplot(meth,aes(y=WT_MethylC,x=class))+geom_boxplot()
plot(gg)

gg=ggplot(meth,aes(y=WT_TAB,x=class))+geom_boxplot()
plot(gg)

