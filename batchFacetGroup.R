#!/usr/bin/Rscript
.libPaths('/online/home/chenyl/software/R')
arg <- commandArgs(T)
setwd(arg[1])
library(ggplot2)
myfun <- function(sam){
tframe <- read.table(sam,header=T,stringsAsFactors=F)
ggplot(tframe,aes(transformation,ratio,fill=type))+
   geom_bar(stat='identity',width=0.5)+
   theme_bw()+facet_grid(sam~type,scales='free_x')+
   scale_x_discrete('Preceded by 5\'\nFollowed by 3\''
                                  ,labels=rep(c('A\nA','\n\nC','\n\n\nG','\n\n\n\nT'
                                  ,'C\nA','\n\nC','\n\n\nG','\n\n\n\nT'
                                  ,'G\nA','\n\nC','\n\n\nG','\n\n\n\nT'
                                  ,'T\nA','\n\nC','\n\n\nG','\n\n\n\nT'),times=12))+
   scale_y_continuous('Proportion(%)')+
   theme(legend.position='none',axis.text.x=element_text(vjust=0.5,size=5)
        ,axis.title.y=element_text(size=14),axis.title.x=element_text(vjust=9,hjust=-0.037,size=5)
        ,panel.grid=element_blank(),rect=element_blank()
        ,axis.line.y=element_line(color='black'),strip.text.y=element_text(size=9),plot.margin=unit(c(0,0,0,0.40), "cm"))
name <- paste(strsplit(sam,split='\\.')[[1]][1],'.wes','.png',sep='')
ggsave(name,height=674,width=359,units='mm',limitsize=F)
}

file1 <- list.files('.','stat')
for (x in file1){
  myfun(x)
}
