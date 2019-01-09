#!/usr/bin/Rscript

cat('Start at',date(),'\n')

args <- commandArgs(T)
pro <- args[1]
group <- args[2]
myfun <- function(sam){

setwd(paste('/online/home/chenyl/dongfang/',pro,'/',group,sep=''))

substitution <- read.table(sam,header=F,sep='\t',stringsAsFactors=F)
#colnames(substitution) <- c('chr','pos','ref','alt','flankseq5','flankseq3','gene','strand')
tmpname <- strsplit(sam,split='\\.')[[1]][1]
name <- paste(tmpname,'.test.xls',sep='')

data1 <- data.frame('chr','pos','before','ref','after','alt','strand')
write.table(data1,name,append=T,sep='\t',quote=F,row.names=F,col.names=F)

for (i in 1:nrow(substitution)) {
   tmp <- substitution[i,]
   chr=tmp[,1]
   pos=tmp[,2]
   strand=tmp[,8]
   if (tmp[,8]==1) {
      flank5=tmp[,5]
      flank3=tmp[,6]
      if (tmp[,3]=='C'|tmp[,3]=='T') {
         ref=tmp[,3]
         alt=tmp[,4]
         before=substr(flank5,10,10)
         after=substr(flank3,1,1)

         rawdata=c(chr,pos,before,ref,after,alt,strand)
      } 
      else {
         reftmp=tmp[,3]
         if (reftmp=='A') {
            ref='T' }
         else if (reftmp=='G') {
            ref='C' }
         
         alttmp=tmp[,4]
         if (alttmp=='A') {
            alt='T' }
         else if (alttmp=='C') {
            alt='G' }
         else if (alttmp=='G') {
            alt='C' }
         else if (alttmp=='T') {
            alt='A' }
         
         betmp=substr(flank5,10,10)
         if (betmp=='A') {
            before='T' }
         else if (betmp=='C') { 
            before='G' }
         else if (betmp=='G') {
            before='C' }
         else if (betmp=='T') {
            before='A' }

         aftmp=substr(flank3,1,1)
         if (aftmp=='A') {
            after='T' }
         else if (aftmp=='C') {
            after='G' }
         else if (aftmp=='G') {
            after='C' }
         else if (aftmp=='T') {
            after='A' }

        rawdata=c(chr,pos,before,ref,after,alt,strand)
      }  
   }
   else {
      flank5=tmp[,5]
      flank3=tmp[,6]
      if (tmp[,3]=='A'|tmp[,3]=='G') {
         ref=tmp[,3]
         alt=tmp[,4]
         before=substr(flank5,10,10)
         after=substr(flank3,1,1)
         rawdata=c(chr,pos,before,ref,after,alt,strand)
      }
      else {

         reftmp=tmp[,3]
         if (reftmp=='C') {
            ref='G' }
         else if (reftmp=='T') {
            ref='A' }

         alttmp=tmp[,4]
         if (alttmp=='A') {
            alt='T' }
         else if (alttmp=='C') {
            alt='G' }
         else if (alttmp=='G') {
            alt='C' }
         else if (alttmp=='T') {
            alt='A' }

        betmp=substr(flank5,10,10)
         if (betmp=='A') {
            before='T' }
         else if (betmp=='C') {
            before='G' }
         else if (betmp=='G') {
            before='C' }
         else if (betmp=='T') {
            before='A' }

         aftmp=substr(flank3,1,1)
         if (aftmp=='A') {
            after='T' }
         else if (aftmp=='C') {
            after='G' }
         else if (aftmp=='G') {
            after='C' }
         else if (aftmp=='T') {
            after='A' }

         rawdata=c(chr,pos,before,ref,after,alt,strand)
      }
   }
data2 <- data.frame(chr=rawdata[1],pos=rawdata[2],before=rawdata[3],ref=rawdata[4],after=rawdata[5],alt=rawdata[6],strand=rawdata[7])
write.table(data2,name,append=T,sep='\t',quote=F,row.names=F,col.names=F)
}

}

files <- list.files( paste('/online/home/chenyl/dongfang/',pro,'/',group,sep=''),'ready.txt')
for (x in files){
  myfun(x)
}

cat('Ended at',date(),'\n')
