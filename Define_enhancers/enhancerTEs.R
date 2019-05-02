#####Inputs####
##set working directory
setwd('~/todd')

##path to bedtools
bedtools = '~/bedtools2/bin/'

##import functions
source('~/todd/Functions/intersectBed.R')

#RefSeq TSS lookup
RefSeq = read.delim('~/todd/Data/Annotations/RefSeq_mRNA.txt')

##TE data
es.te = read.delim('~/todd/Data/Annotations/ESC_TEs_mm10.txt',as.is=T)
ts.te = read.delim('~/todd/Data/Annotations/TSC_TEs_mm10.txt',as.is=T)

##ChIP-seq/DNase-seq data
es.k4me1 = read.delim('~/todd/Data/ChIP_peaks/ES_H3K4me1_peaks.broadPeak',as.is=T,header=F)
es.k4me3 = read.delim('~/todd/Data/ChIP_peaks/ES_H3K4me3_peaks.broadPeak',as.is=T,header=F)
es.k27 = read.delim('~/todd/Data/ChIP_peaks/ES_H3K27ac_peaks.broadPeak',as.is=T,header=F)
es.nanog = read.delim('~/todd/Data/ChIP_peaks/ES_Nanog_peaks.narrowPeak',as.is=T,header=F)
es.oct4 = read.delim('~/todd/Data/ChIP_peaks/ES_Oct4_peaks.narrowPeak',as.is=T,header=F)
es.sox2 = read.delim('~/todd/Data/ChIP_peaks/ES_Sox2_peaks.narrowPeak',as.is=T,header=F)
es.dhs = read.delim('~/todd/Data/DNase_peaks/ES_ATAC_peaks.npf',as.is=T,header=F)

ts.k4me1 = read.delim('~/todd/Data/ChIP_peaks/TS_H3K4me1_peaks.broadPeak',as.is=T,header=F)
ts.k4me3 = read.delim('~/todd/Data/ChIP_peaks/TS_H3K4me3_peaks.broadPeak',as.is=T,header=F)
ts.k27 = read.delim('~/todd/Data/ChIP_peaks/TS_H3K27ac_peaks.broadPeak',as.is=T,header=F)
ts.cdx2 = read.delim('~/todd/Data/ChIP_peaks/TS_CDX2_peaks.narrowPeak',as.is=T,header=F)
ts.elf5 = read.delim('~/todd/Data/ChIP_peaks/TS_ELF5_peaks.narrowPeak',as.is=T,header=F)
ts.eomes = read.delim('~/todd/Data/ChIP_peaks/TS_EOMES_peaks.narrowPeak',as.is=T,header=F)
ts.dhs = read.delim('~/todd/Data/DNase_peaks/TS_ATAC_peaks.npf',as.is=T,header=F)

####Script####
##Determine promoter flanking regions
prom.flank<-500
prom.cord<-unique(data.frame(RefSeq$chrom, RefSeq$TSS-prom.flank, RefSeq$TSS+prom.flank))

##Exclude elements near TSS and overlaping promoter associated H3K4me3

es.ex = intersectBed(es.te,es.k4me3,opt.string='-v',path.to.bedtools=bedtools)
es.ex = intersectBed(es.ex,prom.cord,opt.string='-v',path.to.bedtools=bedtools)
ts.ex = intersectBed(ts.te,ts.k4me3,opt.string='-v',path.to.bedtools=bedtools)
ts.ex = intersectBed(ts.ex,prom.cord,opt.string='-v',path.to.bedtools=bedtools)


##intersect TEs with H3K27ac and ATAC


es.enh = intersectBed(es.ex,es.k27,opt.string='-u',path.to.bedtools=bedtools)
es.enh = intersectBed(es.enh,es.dhs,opt.string='-u',path.to.bedtools=bedtools)

ts.enh = intersectBed(ts.ex,ts.k27,opt.string='-u',path.to.bedtools=bedtools)
ts.enh = intersectBed(ts.enh,ts.dhs,opt.string='-u',path.to.bedtools=bedtools)


##intersect with at least one TF

es.enh = intersectBed(es.enh,rbind(es.nanog,es.oct4,es.sox2),opt.string='-u',path.to.bedtools=bedtools)
ts.enh = intersectBed(ts.enh,rbind(ts.cdx2,ts.elf5,ts.eomes),opt.string='-u',path.to.bedtools=bedtools)

es.enh$ID = paste(es.enh[,6],"REDE",sep="_")
ts.enh$ID = paste(ts.enh[,6],"REDE",sep="_")

##write enhancer files
##Output
write.table(es.enh,'Data/Output/ESC_enhancer_TEs.txt',sep='\t',quote=F,row.names=F,col.names=F)
write.table(ts.enh,'Data/Output/TSC_enhancer_TEs.txt',sep='\t',quote=F,row.names=F,col.names=F)


