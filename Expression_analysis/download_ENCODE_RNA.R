#set working directory
setwd('~/todd/')


##read file metadata
meta = read.delim('Data/RNA-seq/ENCODE_Gingeras_RNA_files.tsv',as.is=T)


##select mm10 gene quantification files
is.mm10 = meta$Assembly=='mm10'
is.gene = meta$Output.type=='gene quantifications'


##download files
write(meta$File.download.URL[is.mm10 & is.gene],'temp_files.txt')
try(system('xargs -n 1 curl -O -L < temp_files.txt'))
unlink('temp_files.txt')


##rename files with sample name
gene.meta = meta[is.mm10 & is.gene,]
fname = paste(gene.meta$File.accession,'tsv',sep='.')
new.name = paste(gsub(' ','_',gene.meta$Biosample.term.name),fname,sep='_')

file.rename(fname,new.name)
