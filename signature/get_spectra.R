#!/usr/bin/Rscript

cat('Start at',date(),'\n')

args <- commandArgs(T)

myfun <- function(x){

path <- dirname(x)
sam <- basename(x)

setwd(path)

substitution <- read.table(sam,header=F,sep='\t',stringsAsFactors=F)
#colnames(substitution) <- c('chr','pos','ref','alt','flank5seq','flank3seq','gene','strand')
tmpname <- strsplit(sam,split='\\.')[[1]][1]
name <- paste(tmpname,'.ready.txt',sep='')

data1 <- data.frame('chr','pos','ref','alt','flank5seq','flank3seq','before','refcontext','after','altcontext','gene','strand')
write.table(data1,name,append=T,sep='\t',quote=F,row.names=F,col.names=F)

for (i in 1:nrow(substitution)) {
   tmp <- substitution[i,]
   chr=tmp[,1]
   pos=tmp[,2]
   ref=tmp[,3]
   alt=tmp[,4]
   flank5=tmp[,5]
   before=substr(flank5,10,10)
   flank3=tmp[,6]
   after=substr(flank3,1,1)
   gene=tmp[,7]
   strand=tmp[,8]
   if (strand==1) {
      #flank5=tmp[,5]
      #flank3=tmp[,6]
      if (ref=='C'|ref=='T') {
         refcontext=ref
         altcontext=alt

         rawdata=c(chr,pos,ref,alt,flank5,flank3,before,refcontext,after,altcontext,gene,strand)
      } 
      else {
         #reftmp=tmp[,3]
         if (ref=='A') {
            refcontext='T'}
         else if (ref=='G') {
            refcontext='C'}
         
         #alttmp=tmp[,4]
         if (alt=='A') {
            altcontext='T' }
         else if (alt=='C') {
            altcontext='G' }
         else if (alt=='G') {
            altcontext='C' }
         else if (alt=='T') {
            altcontext='A' }
         
        #betmp=substr(flank5,10,10)
        #if (betmp=='A') {
        #   before='T' }
        #else if (betmp=='C') { 
        #   before='G' }
        #else if (betmp=='G') {
        #   before='C' }
        #else if (betmp=='T') {
        #   before='A' }

        #aftmp=substr(flank3,1,1)
        #if (aftmp=='A') {
        #   after='T' }
        #else if (aftmp=='C') {
        #   after='G' }
        #else if (aftmp=='G') {
        #   after='C' }
        #else if (aftmp=='T') {
        #   after='A' }

        rawdata=c(chr,pos,ref,alt,flank5,flank3,before,refcontext,after,altcontext,gene,strand)
      }  
   }
   else {
      #flank5=tmp[,5]
      #flank3=tmp[,6]
      if (ref=='A'|ref=='G') {
         refcontext=ref
         altcontext=alt
        
         rawdata=c(chr,pos,ref,alt,flank5,flank3,before,refcontext,after,altcontext,gene,strand)
      }
      else {

         #reftmp=tmp[,3]
         if (ref=='C') {
            refcontext='G'}
         else if (ref=='T') {
            refcontext='A'}

         #alttmp=tmp[,4]
         if (alt=='A') {
            altcontext='T' }
         else if (alt=='C') {
            altcontext='G' }
         else if (alt=='G') {
            altcontext='C' }
         else if (alt=='T') {
            altcontext='A' }

         #betmp=substr(flank5,10,10)
         #if (betmp=='A') {
         #   before='T' }
         #else if (betmp=='C') {
         #   before='G' }
         #else if (betmp=='G') {
         #   before='C' }
         #else if (betmp=='T') {
         #   before='A' }

         #aftmp=substr(flank3,1,1)
         #if (aftmp=='A') {
         #   after='T' }
         #else if (aftmp=='C') {
         #   after='G' }
         #else if (aftmp=='G') {
         #   after='C' }
         #else if (aftmp=='T') {
         #   after='A' }

         rawdata=c(chr,pos,ref,alt,flank5,flank3,before,refcontext,after,altcontext,gene,strand)
      }
   }
data2 <- data.frame(chr=rawdata[1],pos=rawdata[2],ref=rawdata[3],alt=rawdata[4],flank5seq=rawdata[5],flank5seq=rawdata[6],before=rawdata[7],
refcontext=rawdata[8],after=rawdata[9],altcontext=rawdata[10],gene=rawdata[11],strand=rawdata[12])

write.table(data2,name,append=T,sep='\t',quote=F,row.names=F,col.names=F)
}

}

#files <- list.files(paste('/online/home/chenyl/dongfang/',pro,'/',group,sep=''),'ready.txt')
#for (x in files){
#  myfun(x)
#}

myfun(args[1])

cat('Ended at',date(),'\n')
