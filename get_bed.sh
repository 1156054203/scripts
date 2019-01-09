#!/bin/bash
Usage="bash $0 vcffile"

if [ $# -eq 1 ];then
   pro=$1
else
   echo $Usage
   exit 1
fi

prefix=/online/home/chenyl/dongfang

hg38_fa=/online/databases/Homo_sapiens/hg38/hg38bundle/Homo_sapiens_assembly38.fasta
hg19_fa=/online/software/speedseq/genomes/ucsc_ref_hg19/ucsc.hg19.fasta
bedtool=/online/software/bedtools-2.25.0/bin/bedtools

for var in `find $prefix/${pro}/group* -maxdepth 1 -name "*.anno.txt"`;do
    name=`echo $var|awk -F'/' '{print$NF}'|cut -c1-27`
    group=`echo $var|awk -F'/' '{print$7}'`
    out=$prefix/${pro}/${group}
    #5 flank_sequence
    cat $var|grep -v '^#'|awk 'BEGIN{OFS="\t"}{if($6=="+")print$1,$2-11,$2-1,"name",1,$6;else print$1,$2,$2+10,"name",1,$6}' >${out}/${name}.5flank.bed
    #3 flank_sequence
    cat $var|grep -v '^#'|awk 'BEGIN{OFS="\t"}{if($6=="-")print$1,$2,$2+10,"name",1,$6;else print$1,$2-11,$2-1,"name",1,$6}' >${out}/${name}.3flank.bed
done

for var in `find $prefix/${pro}/group* -name "*.bed"`;do
    out=`echo ${var%.*}.seq`
    /online/software/bedtools-2.25.0/bin/bedtools getfasta -fi $hg38_fa -bed $var -tab -s -fo $out
done

for var in `find $prefix/${pro}/group* -name "*.seq"`;do
    name=`echo $var|awk -F'/' '{print$NF}'|cut -c1-27`
    group=`echo $var|awk -F'/' '{print$7}'`
    out=$prefix/${pro}/${group}
    cat $out/${name}.anno.txt | grep -v '^#' | paste -d '\t' - $out/${name}.5flank.seq $out/${name}.3flank.seq | awk '{OFS="\t"}{print$1,$2,$3,$4,$8,$10,$5,$6}' >$out/${name}.ready.txt
    sed -i "1i #chr\tpos\tref\talt\t5flank_seq\t3flank_seq\tgene\tstrand" $out/${name}.ready.txt
    cat $out/${name}.ready.txt | grep -v '^#' | awk '{OFS="\t"}{print substr($5,10),$3,substr($6,1,1),$4}'|sort|uniq -c|awk '{OFS="\t"}{print$2,$3,$4,$5,$2$3$4$5,$1}' | join -a1 -1 5 -2 5 -o 1.1 1.2 1.3 1.4 1.6 2.6 $prefix/mode.xls - | awk '{OFS="\t"}{if($6=="")print$1,$2,$3,$4,$5;else print$1,$2,$3,$4,$6}' >$out/${name}.ready.xls
    sed -i "1i before\tref\tafter\talt\t${name}" $out/${name}.ready.xls
done

