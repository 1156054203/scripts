#!/usr/bin/Rscript
.libPaths('/online/home/chenyl/software/R')
arg <- commandArgs(T)
file <- basename(arg[1])
dir <- dirname(arg[1])
setwd(dir)

library(ggplot2)
library(reshape2)

mode <- read.table('/online/home/chenyl/dongfang/Rmode.txt',header=T,stringsAsFactors=F) ##Rmode contains 96 kinds,content is "before ref after alt context type"

myfun <- function(f1){
frame1 <- read.table(f1,header=T,stringsAsFactors=F)
frame1 <- frame1[order(frame1$ref,frame1$alt),]

for (i in 5:ncol(frame1)){

   sam <- frame1[,c(1:4,i)]
   name <- colnames(sam)[5]
   tmp1 <- sam[which(sam$ref=='C'|sam$ref=='T'),]
   tmp2 <- sam[which(sam$ref=='A'|sam$ref=='G'),]
   tmp3 <- cbind(tmp1[1:4],num=rowSums(cbind(tmp1,name=tmp2[order(tmp2$ref,
        tmp2$alt,tmp2$before,tmp2$after,decreasing=T),][,5])[,c(5,6)]))
   mono <- cbind(tmp3[,1:4],context=paste(tmp3$before,tmp3$ref,tmp3$after,
        tmp3$alt,sep=''),type=paste(tmp3$ref,'>',tmp3$alt,sep=''),tmp3[5])
   colnames(mono)[7] <- name
   mode <- cbind(mode,prop.table(mono[7]))
}

frame2 <- melt(mode,id.vars=c(1:6),variable.name='sam',value.name='num')
write.table(frame2,file=paste(strsplit(arg[1],'\\.')[[1]][1],'.frame.txt',sep=''),quote=F,sep='\t',row.names=F)

stag <- max(frame2$num)*100-0.65
tag <- max(frame2$num)*100-0.6

palette <- c('#A52A2A','#FFC125','#B03060','#76EE00','#A0522D','#BCEE68',
             '#9400D3','#8B475D','#87CEFA','#836FFF','#71C671','#6E8B3D')

ggplot(frame2,aes(factor(context,levels=context[1:96]),num*100,group=context,color=type))+
       geom_boxplot(outlier.shape=NA)+
       scale_color_manual(values=palette)+
      theme(panel.background=element_blank(),panel.grid =element_blank()
          ,plot.title=element_text(hjust=0,size=12)
          ,panel.border=element_blank(),axis.line.y =element_line(color='black')
          ,legend.position='none',axis.text.x=element_text(vjust=0.5,size=6.5)
          ,axis.title.x=element_text(vjust=11,hjust=-0.028,size=5)
          ,axis.title.y=element_text(size=10),axis.text.y=element_text(size=10))+
    annotate('segment',x=0,y=stag,xend=16,yend=stag,color='#B03060')+
    annotate('text',x=7.5,y=tag,color='black',label='C>A')+
    annotate('segment',x=17,y=stag,xend=32,yend=stag,color='#FFC125')+
    annotate('text',x=22.5,y=tag,color='black',label='C>G')+
    annotate('segment',x=33,y=stag,xend=48,yend=stag,color='#A52A2A')+
    annotate('text',x=40.5,y=tag,color='black',label='C>T')+
    annotate('segment',x=49,y=stag,xend=64,yend=stag,color='#76EE00')+
    annotate('text',x=56.5,y=tag,color='black',label='T>A')+
    annotate('segment',x=65,y=stag,xend=80,yend=stag,color='#A0522D')+
    annotate('text',x=72.5,y=tag,color='black',label='T>C')+
    annotate('segment',x=81,y=stag,xend=96,yend=stag,color='#BCEE68')+
    annotate('text',x=88.5,y=tag,color='black',label='T>G')+
    labs(title='')+
    scale_y_continuous('Proportion(%)')+
    coord_cartesian(ylim=c(0, 3))+ ## Modify axis display range,but do not remove data
    scale_x_discrete('Preceded by 5\'\nFollowed by 3\'',
                       labels=rep(c('A\nA','\n\nC','\n\n\nG','\n\n\n\nT'
                                    ,'C\nA','\n\nC','\n\n\nG','\n\n\n\nT'
                                    ,'G\nA','\n\nC','\n\n\nG','\n\n\n\nT'
                                    ,'T\nA','\n\nC','\n\n\nG','\n\n\n\nT'),times=6))
tmpname <- strsplit(arg[1],split='\\.')[[1]][1]
ggsave(paste(tmpname,'.png',sep=''),height=200,width=360,units='mm',limitsize=F)
}


myfun(file)

