####  User Inputs  ####

#Set working directory
setwd("~/Guide_design")

#set path to homertools directory
path.to.homertools="~/Homer/bin/homerTools"

#set path to genome file
path.to.genome="~/Homer/data/genomes/mm10/"

#Set input BED coord file
bed.file="RLTR13D5.txt"

#Set output FASTA file name
output.name="RLTR13D5.fa"


####  Script  ####
command = paste(path.to.homertools,' extract ',bed.file, path.to.genome, ' -fa > ',output.name,sep=' ')
try(system(command))