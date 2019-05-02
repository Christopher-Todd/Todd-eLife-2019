####  User Inputs  ####

#Set working directory
setwd("~/Guide_design")

#set input file path
guide.mat.input = "all_guide_hit_mat.txt"

#set path to genome directory
#Genome downloaded from UCSC by "make_mm10.sh"
genome = "/data/home/hmx579/scratch/Genomes/Mouse/mm10/"

#set output name (CasOFF input file name)
output.name = "cas_off_input.txt"

####  Script  ####

bp.mismatch=2
min.hits=5

mat=read.delim(guide.mat.input)
mat.sums=rowSums(mat)

guides.for.consideration = row.names(mat[mat.sums>min.hits,])

motif="NNNNNNNNNNNNNNNNNNNNNGG"
guides=paste(row.names(mat),"NNN ",bp.mismatch,sep="")
head(guides)

cas.off.input=c(genome,motif,guides)
head(cas.off.input)

write.table(cas.off.input,output.name,sep="\t",col.names = F,row.names = F,quote = F)

