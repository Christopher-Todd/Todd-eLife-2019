
##Generic R implementation of Homer annotatePeak function


annotatepeaks <- function(bed.file,tag,ref.genome=' mm10',path.to.homer='') {
 out = tempfile()
  
  command = paste('perl ',path.to.homer,'annotatePeaks.pl ',bed.file,' ',ref.genome, '-size given -hist 100 -ghist -d ',tag,' > ',out,sep='')
  cat(command,'\n')
  try(system(command))
  res = read.table(out,header=T,as.is=T)
  
  unlink(out)
  return(res)
}