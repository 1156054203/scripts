#!/usr/bin/Rscript

#.libPaths(''c:/Users/user/Documents/R/win-library/'')
myfun <- function(file1){
tframe <- read.table(file1,header=T)
for (i in 5:(ncol(tframe)-1)){
sam <- tframe[,c(1:4,i)]
name <- colnames(sam)[5]

tmp1 <- sam[which(sam$Ref=='C'|sam$Ref=='T'),]
tmp2 <- sam[which(sam$Ref=='A'|sam$Ref=='G'),]
tmp3 <- cbind(tmp1[1:4],num=rowSums(cbind(tmp1,num=tmp2[order(tmp2$Ref,
        tmp2$Var,tmp2$Before,tmp2$After,decreasing=T),][,5])[,c(5,6)]))
mono <- cbind(tmp3[,1:4],context=paste(tmp3$Before,tmp3$Ref,tmp3$After,
        tmp3$Var,sep=''),type=paste(tmp3$Ref,'>',tmp3$Var,sep=''),tmp3[5])

#mono <- mono[order(mono$Ref,mono$Var,mono$Before),]

library(ggplot2)

stag <- max(prop.table(mono$num))+0.01
tag <- max(prop.table(mono$num))+0.03
palette <- c('#A52A2A','#FFC125','#B03060','#76EE00','#A0522D','#BCEE68',
               '#9400D3','#8B475D','#87CEFA','#836FFF','#71C671','#6E8B3D')
ggplot(mono,aes(1:96,prop.table(num),fill=type))+
    geom_bar(stat='identity',width=0.5)+
    scale_fill_manual(values=palette)+
    theme_bw()+
    theme(panel.background=element_blank(),panel.grid =element_blank()
          ,plot.title=element_text(hjust=0,size=12)
          ,panel.border=element_blank(),axis.line.y =element_line(color='black')
          ,legend.position='none',axis.text.x=element_text(vjust=0.5,size=6.5)
          ,axis.title.x=element_text(vjust=11,hjust=-0.005,size=7)
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
    #annotate('segment',x=97,y=stag,xend=112,yend=stag,color='#9400D3')+
    #annotate('text',x=104.5,y=tag,color='black',label='G>A')+
    #annotate('segment',x=113,y=stag,xend=128,yend=stag,color='#87CEFA')+
    #annotate('text',x=120.5,y=tag,color='black',label='G>T')+
    #annotate('segment',x=129,y=stag,xend=144,yend=stag,color='#8B475D')+
    #annotate('text',x=136.5,y=tag,color='black',label='G>C')+
    #annotate('segment',x=145,y=stag,xend=160,yend=stag,color='#836FFF')+
    #annotate('text',x=152.5,y=tag,color='black',label='T>A')+
    #annotate('segment',x=161,y=stag,xend=176,yend=stag,color='#6E8B3D')+
    #annotate('text',x=168.5,y=tag,color='black',label='T>G')+
    #annotate('segment',x=177,y=stag,xend=192,yend=stag,color='#71C671')+
    #annotate('text',x=185.5,y=tag,color='black',label='T>C')+
    labs(title=name)+
    scale_y_continuous('Proportion(%)')+
    scale_x_continuous('Preceded by 5\'\nFollowed by 3\'',breaks=seq(1,96),
                       labels=rep(c('A\nA','\n\nC','\n\n\nG','\n\n\n\nT'
                                    ,'C\nA','\n\nC','\n\n\nG','\n\n\n\nT'
                                    ,'G\nA','\n\nC','\n\n\nG','\n\n\n\nT'
                                    ,'T\nA','\n\nC','\n\n\nG','\n\n\n\nT'),times=6))
  ggsave(paste(name,'.png',sep=''),width=337,height=120,limitsize=F,units='mm')
}
}
myfun('taiwan.HCC.spectra.tsv')