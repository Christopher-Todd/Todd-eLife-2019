#####  User Inputs  #####
#Set working directory
setwd("E:/CRISPRi guide design/")

#### input files ####
#set path to guide matrix from Step 2
guide.matrix.path = "all_guide_hit_mat.txt"
#set path to analysed CasOFFinder results from Step 5
casoff.output.path = "Repeat_hits.txt"


####Guide design method####
#Either as many as possible "Many"
#Or multiple guides to the same element "Multi"
design.method="Multi"

#Number of Guides wanted
number.guides=4

#If Multi method, how many guides is enough for each element?
enough.guides=3



####Preferred elements within the family?####

#Do you want guides to be exclusive to preferred elements?
#"exclusive"
#Or do you want guides to be designed only preferred towards seleceted elements?
#"preferred"
#Or do you you just want guides to be designed to the whole class?
#"no pref"
guide.pref="preferred"


#Prefered element names
prefered.names=c("Enh")


#type is either "ratio" or max "number" of non.prefered elements
exclusive.cutoff.type = "ratio"
exclusive.ratio=0.3
exclusive.number=100


####Off-target cutoffs####

#filter based on predicted exon off-targets?
exon.target.filt=F
#filter based on predicted off-targets?
off.target.filt=T
#cutoff can be either "ratio","number",or "both"
off.target.cutoff.type = "ratio"
#ratio of on-target to off-target 
on.target.ratio=0.8
#max number of off-target hits
off.target.number=100

#Is related TE family hits allowed?
#i.e. ignored from off-target numbers
related.family.tolerated=T

####Output options:####

#Output matrix file name
output.name = "RLTR13D5_guide_options"

##Plot your candidate guide hits within the family
plot.hits=T
plot.preferred.seperate=T
#option for heatmap scaled to relative size of preferred and non preferred groups
scale.preferred=F

#packages required for plotting
library(ComplexHeatmap)
library(circlize)

####  Script   ####

mat = read.delim(guide.matrix.path)
off = read.delim(casoff.output.path)

##Application of off-target filtering prior to guide picking
if(exon.target.filt){
  filt.guides=as.character(off[off$exon.hits==0,1])
  mat=mat[filt.guides,]
  rm(filt.guides)
}
if(off.target.filt){
  if(off.target.cutoff.type=="ratio"){
    if(related.family.tolerated){
      off.target.ratio=(off$correct.hits+off$related.hits)/off$total.hits
    }else{
    off.target.ratio=(off$correct.hits)/off$total.hits
    }
    filt.guides=as.character(off$guide[off.target.ratio>=on.target.ratio])
    mat=mat[filt.guides,]
    rm(filt.guides,off.target.ratio)
  }
  if(off.target.cutoff.type=="number"){
    if(related.family.tolerated){
      off.target.hits=off$total.hits-(off$correct.hits+off$related.hits)
    }else{
    off.target.hits=off$total.hits-(off$correct.hits)
    }
    filt.guides=as.character(off$guide[off.target.hits<=off.target.number])
    mat=mat[filt.guides,]
    rm(filt.guides,off.target.hits)
  }
  if(off.target.cutoff.type=="both"){
    if(related.family.tolerated){
      off.target.ratio=(off$correct.hits+off$related.hits)/off$total.hits
      off.target.hits=off$total.hits-(off$correct.hits+off$related.hits)
      }else{
      off.target.ratio=(off$correct.hits)/off$total.hits
      off.target.hits=off$total.hits-(off$correct.hits)
    }
    filt.guides=as.character(off$guide[(off.target.ratio>=on.target.ratio)&(off.target.hits<=off.target.number)])
    mat=mat[filt.guides,]
    rm(filt.guides,off.target.ratio,off.target.hits)
  }
}

##functions for many and multi method of guide picking
many.method=function(input.mat,number.guides){
  
  hit.mat=data.frame()
  for(i in 1:number.guides){
      #At beginning start with the full matrix
      if(i == 1){mat.sub=input.mat}
      #Calculate the number of elements that each guide hits
      hits=rowSums(mat.sub)
      #Select the guide with the most hits
      top.guide=row.names(mat.sub)[which.max(hits)]
      #Populate the dataframe with the selected guide
      hit.row=input.mat[top.guide,]
      hit.mat=rbind(hit.mat,hit.row)
      #Check to see if any elements within the family already have been targetted
      uncovered.elements=colnames(hit.mat)[colSums(hit.mat)<1]
      #Disregard any elements with guides hitting them for next iteration of the loop
      mat.sub=input.mat[!(row.names(input.mat)%in%row.names(hit.mat)),uncovered.elements]
    }
  return(row.names(hit.mat))
  
}

multi.method=function(input.mat,number.guides,enough.guides){
  #generate empty matrix
  hit.mat=data.frame()
  
  for(i in 1:number.guides){
    #At beginning start with the full matrix
    if(i == 1){mat.sub=input.mat}
    #Calculate the number of elements that each guide hits
    hits=rowSums(mat.sub)
    #Select the guide with the most hits
    top.guide=row.names(mat.sub)[which.max(hits)]
    #Populate the dataframe with the selected guide
    hit.row=input.mat[top.guide,]
    hit.mat=rbind(hit.mat,hit.row)
    #Check to see if any elements within the family already have enough guides targeting them
    uncovered.elements=colnames(hit.mat)[colSums(hit.mat)<enough.guides]
    #Disregard any elements with enough guides hitting them for next iteration of the loop
    mat.sub=input.mat[!(row.names(input.mat)%in%row.names(hit.mat)),uncovered.elements]
  }
  
  return(row.names(hit.mat))
}


if(guide.pref=="exclusive"){
  pref.mat=mat[,grep(paste(prefered.names,collapse = "|"),colnames(mat))]
  non.mat=mat[,grepl(paste(prefered.names,collapse = "|"),colnames(mat))==F]
  enh.hits=rowSums(pref.mat)
  non.hits=rowSums(non.mat)
  total.hits=enh.hits+non.hits
  if(exclusive.cutoff.type=="ratio"){
    exclusive.cutoff=(enh.hits/total.hits)>exclusive.ratio
  }
  if(exclusive.cutoff.type=="number"){
    exclusive.cutoff=non.hits<exclusive.number
  }
  pref.mat=pref.mat[exclusive.cutoff,]
  if(design.method=="Multi"){
    guides=multi.method(input.mat = pref.mat,number.guides,enough.guides)
  }
  if(design.method=="Many"){
    guides=many.method(input.mat = pref.mat,number.guides,enough.guides)
  }
  rm(exclusive.cutoff,non.mat,enh.hits,non.hits,total.hits,pref.mat)
  }else
if(guide.pref=="preferred"){
  pref.mat=mat[,grep(paste(prefered.names,collapse = "|"),colnames(mat))]
  if(design.method=="Multi"){
    guides=multi.method(input.mat = pref.mat,number.guides,enough.guides)
  }
  if(design.method=="Many"){
    guides=many.method(input.mat = pref.mat,number.guides,enough.guides)
  }
  rm(pref.mat)
}else
if(guide.pref=="no.pref"){
  if(design.method=="Multi"){
    guides=multi.method(input.mat = mat,number.guides,enough.guides)
  }
  if(design.method=="Many"){
    guides=many.method(input.mat = mat,number.guides,enough.guides)
  }
}
  




heat.col=colorRamp2(c(0,1),c("gray","red"))

if(plot.hits){
  if(plot.preferred.seperate){
    pref.cols=grepl(paste(prefered.names,collapse = "|"),colnames(mat))
    pref.mat=mat[guides,pref.cols==T]
    non.mat=mat[guides,pref.cols==F]
    plot.height=(number.guides+1)*0.5
    if(scale.preferred){
      heatmap.scale=(ncol(pref.mat)/(ncol(non.mat)+ncol(pref.mat)))}else{
        heatmap.scale=0.5
      }
    
    hm.pref=Heatmap(pref.mat*1,col=heat.col,
                    cluster_rows = T,cluster_columns = T,show_heatmap_legend = F,
                    show_column_names = F,show_row_names = F,show_column_dend = F,show_row_dend = F,
                    width=6*heatmap.scale,column_title = "Pref")
    hm.non=Heatmap(non.mat*1,col=heat.col,
                   cluster_rows = T,cluster_columns = T,show_heatmap_legend = F,
                   show_column_names = F,show_row_names = T,show_column_dend = F,show_row_dend = F,
                   width=6*(1-heatmap.scale),column_title = "Non.Pref")
    hm.plot=hm.pref+hm.non
  }else{
    hm.plot=Heatmap(mat[guides,]*1,col=heat.col,
                    cluster_rows = T,cluster_columns = T,show_heatmap_legend = F,
                    show_column_names = F,show_row_names = T,show_column_dend = F,show_row_dend = F)
    plot.height=(number.guides)*0.5
  }
  pdf(paste(output.name,"_heatmap.pdf",sep=""),width = 8,height=plot.height)
  draw(hm.plot)
  dev.off()
}

write.table(mat[guides,],paste(output.name,"_hit_matrix.txt",sep=""),sep="\t",col.names = T,row.names = T,quote = F)
write.table(guides,paste(output.name,"_sequences.txt",sep=""),sep="\t",col.names = F,row.names = F,quote = F)
write.table(off[off$guide%in%guides,],paste(output.name,"_CasOFF_summary.txt",sep=""),sep="\t",col.names = T,row.names = F,quote = F)

