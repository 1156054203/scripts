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
   mode <- cbind(mode,mono[7])
}

frame2 <- melt(mode,id.vars=c(1:6),variable.name='sam',value.name='num')
write.table(frame2,file=paste(strsplit(arg[1],'\\.')[[1]][1],'.frame.txt',sep=''),quote=F,sep='\t',row.names=F)

ggplot(frame2,aes(context,prop.table(num),fill=type))+
   geom_bar(stat='identity',width=0.5)+
   theme_bw()+facet_grid(sam~type,scales='free_x')+
   scale_x_discrete('Preceded by 5\'\nFollowed by 3\''
                                  ,labels=rep(c('A\nA','\n\nC','\n\n\nG','\n\n\n\nT'
                                  ,'C\nA','\n\nC','\n\n\nG','\n\n\n\nT'
                                  ,'G\nA','\n\nC','\n\n\nG','\n\n\n\nT'
                                  ,'T\nA','\n\nC','\n\n\nG','\n\n\n\nT'),times=6))+
   scale_y_continuous('Proportion(%)')+
   theme(legend.position='none',axis.text.x=element_text(vjust=0.5,size=5)
        ,axis.title.y=element_text(size=14),axis.title.x=element_text(vjust=9,hjust=-0.037,size=5)
        ,panel.grid=element_blank(),rect=element_blank()
        ,axis.line.y=element_line(color='black'),strip.text.y=element_text(size=9),plot.margin=unit(c(0,0,0,0.40), "cm"))
tmpname <- strsplit(arg[1],split='\\.')[[1]][1]
ggsave(paste(tmpname,'.png',sep=''),height=408,width=359,units='mm',limitsize=F)
}


myfun(file)

