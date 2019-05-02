
##Genetic R implementation of intersectBed function


intersectBed <- function(a,b,opt.string='-wa',path.to.bedtools='') {
	a.file = tempfile()
	b.file = tempfile()
	write.table(a,file=a.file,quote=F,sep='\t',col.names=F,row.names=F)
	write.table(b,file=b.file,quote=F,sep='\t',col.names=F,row.names=F)
	out = tempfile()
	
	command = paste(path.to.bedtools,'intersectBed ',opt.string,' -a ', a.file,' -b ',b.file,' > ',out,sep='')
	cat(command,'\n')
	try(system(command))
	info=file.info(out)
	if(info$size==0){res=data.frame()}else{res = read.table(out,header=F,as.is=T)}	
	unlink(c(a.file,b.file,out))
	return(res)
}
