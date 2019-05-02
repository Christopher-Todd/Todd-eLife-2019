####  User Inputs  ####

#Set working directory
setwd("~/Guide_design")

#path to input fasta file
repeat.fa.file="~/Guide_design/RLTR13D5.fa"

#set output file name
output.name = "all_guide_hit_mat.txt"


####  Script  ####

##getting and formating guides from fasta file
fasta = scan(repeat.fa.file, sep='\n', character()) 
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

##reverse compliment sequence function
rev.comp<-function(dna){
  #reverse compliment function
  seq<-strsplit(dna,split="")[[1]]
  seq.rev<-rev(seq)
  seq.rev<-paste(seq.rev,collapse = "")
  sub<-gsub("C","g",seq.rev)
  sub<-gsub("G","c",sub)
  sub<-gsub("A","t",sub)
  sub<-gsub("T","a",sub)
  revcom<-toupper(sub)
  return(revcom)}


#getting all the possible guide options

#loop through each fasta seq
seq.guides=list()
for(n in 1:length(seq)){
rep.seq=seq[n]
#get forward seqs
for.rep.seq=unlist(gregexpr("GG",rep.seq))
n.guides=length(for.rep.seq)
for.seq.guides=c()
for(i in 1:n.guides){
pam.pos=for.rep.seq[i]  
for.seq.guides[i]=substr(rep.seq,pam.pos-21,pam.pos-2)
}
#cut to only full length guides
for.seq.guides=for.seq.guides[nchar(for.seq.guides)==20]

#get reverse seqs
rev.rep=rev.comp(rep.seq)
rev.rep.seq=unlist(gregexpr("GG",rev.rep))
n.guides=length(rev.rep.seq)
rev.seq.guides=c()
for(i in 1:n.guides){
  pam.pos=rev.rep.seq[i]  
  rev.seq.guides[i]=substr(rev.rep,pam.pos-21,pam.pos-2)
}
#cut to only full length guides
rev.seq.guides=rev.seq.guides[nchar(rev.seq.guides)==20]

seq.guides[[n]]=unique(append(for.seq.guides,rev.seq.guides))
}
rm(for.rep.seq,for.seq.guides,i,last,n,first,n.guides,newseq,pam.pos,rep.seq,rev.rep,rev.seq.guides,rev.rep.seq)
all.guide.options=unique(unlist(seq.guides))
names(seq.guides)=names(seq)

guide.present=list()
for(i in 1:length(all.guide.options)){
  guide.present[[i]]=grepl(all.guide.options[i],seq.guides)
  print(paste("    ",i,"/",length(all.guide.options),"      ",i/length(all.guide.options)*100,"%","    ",sep=" "))
}

mat=do.call(rbind,guide.present)
colnames(mat)=names(seq)
rownames(mat)=all.guide.options

write.table(mat,output.name,sep="\t",col.names = T,row.names = T,quote = F)



