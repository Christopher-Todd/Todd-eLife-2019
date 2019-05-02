####Inputs####
##set working directory
setwd('~/todd/')

##path to bedtools
bedtools = '~/bedtools2/bin/'

##import functions
source('~/todd/Functions/intersectBed.R')

##RepeatMasker database
RepMask = read.delim('~/todd/Data/Annotations/RepeatMasker.txt',h=F)

#RefSeq TSS lookup
RefSeq = read.delim('~/todd/Data/Annotations/RefSeq_mRNA.txt')

##ChIP-seq/DNase-seq data
es.atac = read.delim('Data/DNase_peaks/ES_ATAC_peaks.npf',as.is=T,header=F)
es.k4me3 = read.delim('Data/ChIP_peaks/ES_H3K4me3_peaks.broadPeak',as.is=T,header=F)
es.k27 = read.delim('Data/ChIP_peaks/ES_H3K27ac_peaks.broadPeak',as.is=T,header=F)
es.nanog = read.delim('Data/ChIP_peaks/ES_Nanog_peaks.narrowPeak',as.is=T,header=F)
es.oct4 = read.delim('Data/ChIP_peaks/ES_Oct4_peaks.narrowPeak',as.is=T,header=F)
es.sox2 = read.delim('Data/ChIP_peaks/ES_Sox2_peaks.narrowPeak',as.is=T,header=F)

ts.k4me3 = read.delim('Data/ChIP_peaks/TS_H3K4me3_peaks.broadPeak',as.is=T,header=F)
ts.k27 = read.delim('Data/ChIP_peaks/TS_H3K27ac_peaks.broadPeak',as.is=T,header=F)
ts.cdx2 = read.delim('Data/ChIP_peaks/TS_CDX2_peaks.narrowPeak',as.is=T,header=F)
ts.elf5 = read.delim('Data/ChIP_peaks/TS_ELF5_peaks.narrowPeak',as.is=T,header=F)
ts.eomes = read.delim('Data/ChIP_peaks/TS_EOMES_peaks.narrowPeak',as.is=T,header=F)
ts.atac = read.delim('Data/DNase_peaks/TS_ATAC_peaks.npf',as.is=T,header=F)


##Determine promoter flanking regions
prom.flank<-500
prom.cord<-unique(data.frame(RefSeq$chrom, RefSeq$TSS-prom.flank, RefSeq$TSS+prom.flank))

##Determine size of enhancer elements 
size=2500

####Script####
##Begin with ATAC peaks which are not either in promoter or repeat regions
print(paste("number es atac:",nrow(es.atac)))
es.ex = intersectBed(es.atac,prom.cord,opt.string='-v',path.to.bedtools=bedtools)
print(paste("number minus prom:",nrow(es.ex)))
es.ex = intersectBed(es.ex,RepMask,opt.string='-v',path.to.bedtools=bedtools)
print(paste("number minus repmask:",nrow(es.ex)))
es.ex = intersectBed(es.ex,es.k4me3,opt.string='-v',path.to.bedtools=bedtools)
print(paste("number minus h3k4me3:",nrow(es.ex)))

print(paste("number ts atac:",nrow(ts.atac)))
ts.ex = intersectBed(ts.atac,prom.cord,opt.string='-v',path.to.bedtools=bedtools)
print(paste("number minus prom:",nrow(ts.ex)))
ts.ex = intersectBed(ts.ex,RepMask,opt.string='-v',path.to.bedtools=bedtools)
print(paste("number minus repmask:",nrow(ts.ex)))
ts.ex = intersectBed(ts.ex,ts.k4me3,opt.string='-v',path.to.bedtools=bedtools)
print(paste("number minus h3k4me3:",nrow(ts.ex)))

##intersect peaks with H3K27ac and a TF

es.enh = intersectBed(es.ex,es.k27,opt.string='-u',path.to.bedtools=bedtools)
print(paste("number es ex plus K27:",nrow(es.enh)))
es.enh = intersectBed(es.enh,rbind(es.sox2,es.oct4,es.nanog),opt.string='-u',path.to.bedtools=bedtools)
print(paste("number plus TF:",nrow(es.enh)))
ts.enh = intersectBed(ts.ex,ts.k27,opt.string='-u',path.to.bedtools=bedtools)
print(paste("number ts ex plus K27:",nrow(ts.enh)))
ts.enh = intersectBed(ts.enh,rbind(ts.cdx2,ts.elf5,ts.eomes),opt.string='-u',path.to.bedtools=bedtools)
print(paste("number plus TF:",nrow(ts.enh)))

es.enh[,4]<-paste(paste("ESCnonTE",seq(1:nrow(es.enh)),sep=""),'NEDE',sep="_")
ts.enh[,4]<-paste(paste("TSCnonTE",seq(1:nrow(ts.enh)),sep=""),'NEDE',sep="_")
es.enh<-es.enh[,1:4]
ts.enh<-ts.enh[,1:4]
es.enh[,5]<-rep(0)
ts.enh[,5]<-rep(0)

bed<-rbind(es.enh,ts.enh)
mid=(bed[,2]+bed[,3])/2
bed[,6]<-round(mid-(size/2))
bed[,7]<-round(mid+(size/2))
colnames(bed)<-c("Chr","Peak.start","Peak.end","ID","score","Plot.start","Plot.end")

#removing incomplete/unmapped chromosome regions
bed2<-bed[grepl('\\_',bed[,1])==F,]

##Output
write.table(es.enh,'Data/Output/ESC_enhancer_nonTEs.txt',sep='\t',quote=F,row.names=F,col.names=F)
write.table(ts.enh,'Data/Output/TSC_enhancer_nonTEs.txt',sep='\t',quote=F,row.names=F,col.names=F)
write.table(bed2[,1:7],'Data/Output/NonTE_Enhancer_coord.txt',sep="\t",quote=F,row.names = F,col.names = T)
