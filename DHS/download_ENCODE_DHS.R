
setwd('~/todd/DHS')


##read file metadata

meta = read.delim('~/todd/Data/DNase_peaks/ENCODE_DNase_files.txt',as.is=T)


##select mm10 bed peak files

bed = meta[grep('bed.gz',meta$File.download.URL),]
is.mm10 = bed$Assembly=='mm10'
is.peaks = bed$Output.type=='peaks'


##download files

write(bed$File.download.URL[is.mm10 & is.peaks],'temp_files.txt')
try(system('xargs -n 1 curl -O -L < temp_files.txt'))
unlink('temp_files.txt')


##rename files with sample name

bed.meta = bed[is.mm10 & is.peaks,]
fname = paste(bed.meta$File.accession,'bed.gz',sep='.')
new.name = paste(gsub(' ','_',bed.meta$Biosample.term.name),fname,sep='_')

file.rename(fname,new.name)
