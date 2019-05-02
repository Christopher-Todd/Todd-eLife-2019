##Generic R implementation of Homer annotatePeak function


annotatetrend <- function(bed.file,tag,ref.genome=' mm10',path.to.homer='') {
 out = tempfile()
  
  command = paste('perl ',path.to.homer,'annotatePeaks.pl ',bed.file,' ',ref.genome, ' -size given -hist 100 -d ',tag,' > ',out,sep='')
  cat(command,'\n')
  try(system(command))
  res = read.delim(out,header=T,as.is=T)
  
  unlink(out)
  return(res)
}
