####Inputs####
##set working directory
setwd('~/todd')

##path to bedtools
bedtools = '~/bedtools2/bin/'

##import functions
source('Functions/intersectBed.R')

#File containing coordinates for HiC baits
bait.file = '~/todd/Data/PCHiC/mm10.baitmap'
#Files containing the coordinates for various gene annotations
rna.seq.file = '~/todd/Data/RNA-seq/Cambuli_RNA_logFPKM.txt'
ensmbl.genes = '~/todd/Data/Annotations/ENSMBL_genes_lookup.txt'
ensmbl.ncrna = '~/todd/Data/Annotations/ENSMBL_ncrna_lookup.txt'

####Script####

baits=read.delim(bait.file,h=F)
rna=read.delim(rna.seq.file)[,1:5]
e.genes=read.delim(ensmbl.genes)
e.ncrna=read.delim(ensmbl.ncrna)

#format the baits file
colnames(baits)=c("chr","start","end","score","names")
baits$chr=paste("chr",baits$chr,sep="")
baits$bait.id=paste(baits$chr,baits$start,baits$end,sep="-")
baits=baits[,c(1:3,6)]

#write.table(baits,"mm10_baits.txt",sep="\t",col.names = F,row.names = F,quote = F)
##will use intersectBed function to see what regions this overlaps with


##generate a bed list of promoter regions for Cambuli_RNA 
prom.start=c()
prom.end=c()
rna$Chromosome=paste("chr",rna$Chromosome,sep="")
for(i in 1:nrow(rna)){
  if(rna$Strand[i]=="-"){
    prom.start[i]=rna$End[i]-500
    prom.end[i]=rna$End[i]+500
  }
  if(rna$Strand[i]=="+"){
    prom.start[i]=rna$Start[i]-500
    prom.end[i]=rna$Start[i]+500
  }
}
rna.df=data.frame(chr=rna$Chromosome,p.start=prom.start,p.end=prom.end,name=rna$Probe)

##generate a bed list of promoter regions for ENSMBL_genes 
prom.start=c()
prom.end=c()
#remove unmapped regions and MT 
e.genes=e.genes[!grepl(e.genes$chr,pattern = "M"),]
for(i in 1:nrow(e.genes)){
  if(e.genes$strand[i]=="-"){
    prom.start[i]=e.genes$end[i]-500
    prom.end[i]=e.genes$end[i]+500
  }
  if(e.genes$strand[i]=="+"){
    prom.start[i]=e.genes$start[i]-500
    prom.end[i]=e.genes$start[i]+500
  }
}
e.genes.df=unique(data.frame(chr=e.genes$chr,p.start=prom.start,p.end=prom.end,name=e.genes$id))

##generate a bed list of promoter regions for ENSMBL_ncrna 
prom.start=c()
prom.end=c()
e.ncrna=e.ncrna[!grepl(e.ncrna$chr,pattern = "M"),]
for(i in 1:nrow(e.ncrna)){
 if(e.ncrna$strand[i]=="-"){
    prom.start[i]=e.ncrna$end[i]-500
    prom.end[i]=e.ncrna$end[i]+500
 }
  if(e.ncrna$strand[i]=="+"){
    prom.start[i]=e.ncrna$start[i]-500
    prom.end[i]=e.ncrna$start[i]+500
  }
}
e.ncrna.df=unique(data.frame(chr=e.ncrna$chr,p.start=prom.start,p.end=prom.end,name=e.ncrna$id))



#intersect promoter regions with bait map
rna.look=intersectBed(baits,rna.df,opt.string='-wa -wb',path.to.bedtools=bedtools)
e.genes.look=intersectBed(baits,e.genes.df,opt.string='-wa -wb',path.to.bedtools=bedtools)
e.ncrna.look=intersectBed(baits,e.ncrna.df,opt.string='-wa -wb',path.to.bedtools=bedtools)

#write results to file
write.table(rna.look[,c(4,8)],"./Data/PCHiC/bait_RNAseq_look.txt",sep="\t",col.names = F,row.names = F,quote = F)
write.table(unique(e.genes.look[,c(4,8)]),"./Data/PCHiC/bait_ENSMBL_genes_look.txt",sep="\t",col.names = F,row.names = F,quote = F)
write.table(unique(e.ncrna.look[,c(4,8)]),"./Data/PCHiC/bait_ENSMBL_ncrna_look.txt",sep="\t",col.names = F,row.names = F,quote = F)


