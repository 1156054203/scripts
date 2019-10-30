#!/bin/bash
#count fpkm for gene that peak annotated
geneloc=$1 #gene bed 4 columns
geneanno=$2 #peak annotate bed 9 columns
flagstat=$3 #samtools flagstat
bamfile=$4 #bam file
outpre=$5 #out prefix

if [ $# != 5 ];then
   echo -e "Usage:\n     bash $0 geneloc geneanno flagstat bamfile outpre"
   exit 0
fi

sumread=`sed -n 1p $flagstat|cut -d+ -f1`
awk 'NR==FNR{OFS=FS="\t";a[$1FS$4]=$2FS$3;next}{if($2FS$9 in a)print $2,a[$2FS$9],$9}' $geneloc $geneanno|awk '{OFS="\t";print$0,($3-$2+1)}'|bedtools multicov -bams $bamfile -bed -|awk -v sum=$sumread '{printf("%s\t%s\t%s\t%s\t%.10f\n",$1,$2,$3,$4,($6*1000000)/(sum*$5))}'|sort -V >${outpre}.fpkm.txt
