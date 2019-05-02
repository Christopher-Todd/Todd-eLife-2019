####Input####
##Set working directory
setwd("~/todd/Motif_analysis/")

##Set paths to Homertools and Genome
path.to.homertools="~/Homer/bin/homerTools"
path.to.homer.format.genome="~/Homer/data/genomes/mm10/"
##Set paths to Clustal Omega
path.to.clustalo="~/clustalo"
##Set paths to MEME tools
path.to.fimo="~/meme/bin/fimo"
path.to.ame="~/meme/bin/ame"

##Set files for Enhancer and Non-Enhancer TE coordinates
es.enh=read.delim("~/todd/Data/Output/ESC_enhancer_TEs.txt",h=F)
ts.enh=read.delim("~/todd/Data/Output/TSC_enhancer_TEs.txt",h=F)
non.enh=read.delim("~/todd/Data/Output/NonEnhTE_candidates.txt")

##list of TE classes to be analysed
rltr.list=c("RLTR9E","RLTR9D","RLTR9A3","RLTR13B1","RLTR13B2","RLTR13B3","RLTR13B4","RLTR13D6","RLTR13D5")

##set directory and database for HOCOMOCO motifs
directory.for.hocomoco='./hocomoco'
hocomoco.database="~/motif_databases/MOUSE/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme"
##set directory and database for Uniprobe motifs
directory.for.uniprobe='./uniprobe'
uniprobe.database="~/motif_databases/MOUSE/uniprobe_mouse.meme"

####Script####
enh=rbind(es.enh,ts.enh)

##Function for Clustal Omega sequence clustering 
clustalo = function(ltr) {
  
  ##Get fasta for enh and nonenh groups
  enh.sub=enh[grep(enh[,5],pattern = ltr),c(1:3,8,7,4)]
  non.sub=non.enh[grep(non.enh$repName,pattern = ltr),c(1:3,8,7,4)]
  
  temp.enh.file<-tempfile()
  write.table(enh.sub,file=temp.enh.file,quote=F,sep='\t',col.names=F,row.names=F)
  command = paste(path.to.homertools,' extract ',temp.enh.file,' ',path.to.homer.format.genome,' -fa > ',ltr,"_enh",".fa",sep='')
  try(system(command))
  unlink(temp.enh.file)
  temp.non.file=tempfile()
  write.table(non.sub,file=temp.non.file,quote=F,sep='\t',col.names=F,row.names=F)
  command = paste(path.to.homertools,' extract ',temp.non.file,' ',path.to.homer.format.genome,' -fa > ',ltr,"_non",".fa",sep='')
  try(system(command))
  unlink(temp.non.file)
  
  ##cutoff for getting only mostly full copies
  fa=scan(paste(ltr,'_enh.fa',sep=''), sep='\n', character()) 
  id=seq(from = 1,to = length(fa),by = 2)
  seq=seq(from = 2,to = length(fa),by = 2)
  max=max(unlist(lapply(fa[seq], nchar)))
  cutoff=max*0.6
  get.fa.file=function(seq){
    log=nchar(fa[seq])>cutoff
    file.id=list()
    for(i in 1:length(log)){
      if(log[i]==T){
        file.id[[i]]=c(fa[id[i]],fa[seq[i]])
      }else{next}
    }
    new.file=unlist(file.id)
    return(new.file)}
  write.table(get.fa.file(seq), paste(ltr,"_enh.fa",sep=''),sep = "\n",quote = F,col.names = F,row.names = F)
  fa=scan(paste(ltr,'_non.fa',sep=''), sep='\n', character()) 
  id=seq(from = 1,to = length(fa),by = 2)
  seq=seq(from = 2,to = length(fa),by = 2)
  write.table(get.fa.file(seq),paste(ltr,"_non.fa",sep=""),sep="\n",quote = F,col.names = F,row.names = F)
  
  ltr.subs=c(paste(ltr,"_enh",sep=""),paste(ltr,"_non",sep=""))
  for(l.sub in ltr.subs){
  ##clustalo round 1
  
  command = paste(path.to.clustalo,' -i ',l.sub,'.fa -t DNA -o ',l.sub,'_omega.fa --threads=2 --output-order=tree-order --force',sep='')
  try(system(command))
  
  ##read for sequence filtering
  fasta = scan(paste(l.sub,'omega.fa',sep='_'), sep='\n', character())
  
  ##Merge sequences
  newseq = grep('>',fasta)
  
  seq = character(length(newseq))
  for (i in 1:length(newseq)) {
    first = newseq[i]+1
    if (i==length(newseq)) {
      last = length(fasta)
    } else {
      last = newseq[i+1]-1
    }
    seq[i] = toupper(paste(fasta[first:last],collapse=''))
  }
  names(seq) = unlist(lapply(strsplit(fasta[newseq],'[> ]'),function(x) x[2]))
  
  bp = strsplit(seq,split='')
  bp.mat = matrix(unlist(bp),nrow=length(bp),byrow=T)
  
  ##Count gap size in sliding window
  n = 5
  gap.size = matrix(nrow=nrow(bp.mat),ncol=ncol(bp.mat)-n+1)
  for (i in 1:ncol(gap.size)) {
    window = bp.mat[,i:(i+n-1)]
    gap.size[,i] = rowSums(window=='-')
  }
  
  ##Get frequent gaps
  gap.freq = colSums(gap.size==n)
  frequent = gap.freq>=0.97*length(bp)
  
  ##Remove sequences making gaps
  not.empty = gap.size[,frequent]==0
  spurious = rowSums(not.empty)>0
  filtered = seq[!spurious]
  
  
  #Write
  file = paste(l.sub,'filt.fa',sep='_')
  unlink(file)
  for (i in 1:length(filtered)) {
    write(paste('>',names(filtered)[i],sep=''),file,append=T)
    write(filtered[i],file,append=T)
  }
  
  ##clustalo round 2
  command = paste(path.to.clustalo,' -i ',l.sub,'_filt.fa -t DNA -o ',l.sub,'_omega2.fa --threads=2 --output-order=tree-order --force',sep='')
  try(system(command))
}

}
##Function for FIMO motif identification
fimo=function(ltr){
  ltr.subs=c(paste(ltr,"_enh",sep=""),paste(ltr,"_non",sep=""))
  for(l.sub in ltr.subs){
    command = paste(path.to.fimo,' --o ',l.sub,'_fimo ',motif.db," ../",l.sub,'.fa',sep='')
    try(system(command))
    try(system(paste('mv ',l.sub,'_fimo/fimo.txt ./',l.sub,'_fimo.txt',sep='')))
    try(system(paste('rm -r ',l.sub,'_fimo/',sep='')))
  }}
##Function for AME motif enrichment analysis
ame=function(ltr){
  ltr.subs=c(paste(ltr,"_enh",sep=""),paste(ltr,"_non",sep=""))
  
  command = paste(path.to.ame,' --o ',ltr,'_ame --control ../',ltr.subs[2],'.fa  ../' ,ltr.subs[1],'.fa ',motif.db,sep='')
  try(system(command))
  ##move and rename file
  try(system(paste('mv ',ltr,'_ame/ame.txt ./',ltr,'_ame.txt',sep='')))
  try(system(paste('rm -r ',ltr,'_ame/',sep='')))
  

}



##Iterate Clustalo function for each class
for (ltr in rltr.list) clustalo(ltr)

##Perform analysis with hocomoco database
setwd(directory.for.hocomoco)
motif.db=hocomoco.database

##Iterate FIMO function for each class
for (ltr in rltr.list) fimo(ltr)
##Iterate AME function for each class
for (ltr in rltr.list) ame(ltr)

##Perform analysis with uniprobe database
setwd(directory.for.uniprobe)
motif.db=uniprobe.database

##Iterate FIMO function for each class
for (ltr in rltr.list) fimo(ltr)
##Iterate AME function for each class
for (ltr in rltr.list) ame(ltr)
