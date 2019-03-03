#!/bin/bash
Usage="bash $0 vcffile"

#if [ $# -eq 1 ];then
#   pro=$1
#else
#   echo $Usage
#   exit 1
#fi

prefix=/online/home/chenyl/dongfang/intersection

hg38_fa=/online/databases/Homo_sapiens/hg38/hg38bundle/Homo_sapiens_assembly38.fasta
hg19_fa=/online/software/speedseq/genomes/ucsc_ref_hg19/ucsc.hg19.fasta
bedtool=/online/software/bedtools-2.25.0/bin/bedtools
Rscript=/online/software/R-3.4.1/bin/Rscript
getSpectra=/online/home/chenyl/dongfang/get_spectra.R

for var in `find $prefix/group* -maxdepth 1 -name "*.isec.anno.vcf"`;do
    name=$(basename $var '.isec.anno.vcf')
    group=`echo $var|awk -F'/' '{print$7}'`
    out=$prefix/${group}
    #5 flank_sequence
    cat $var|grep -v '^#'|awk -F '[\t|]' '{print$1,$2,$27}'|awk 'BEGIN{OFS="\t"}{if($3=="+")print$1,$2-11,$2-1,"name",1,$3;else print$1,$2,$2+10,"name",1,$3}' >${out}/${name}.5flank.bed
    #3 flank_sequence
    cat $var|grep -v '^#'|awk -F '[\t|]' '{print$1,$2,$27}'|awk 'BEGIN{OFS="\t"}{if($3=="-")print$1,$2,$2+10,"name",1,$3;else print$1,$2-11,$2-1,"name",1,$3}' >${out}/${name}.3flank.bed
done

for var in `find $prefix/group* -name "*.bed"`;do
    out=`echo ${var%.*}.seq`
    /online/software/bedtools-2.25.0/bin/bedtools getfasta -fi $hg38_fa -bed $var -tab -s -fo $out
done

for var in `find $prefix/group* -name "*.seq"`;do
    name=`echo $var|awk -F'/' '{print$NF}'|cut -c1-19`
    group=`echo $var|awk -F'/' '{print$7}'`
    out=$prefix/${group}
    cat $out/${name}.isec.anno.vcf | grep -v '^#' |awk -F '[\t|]' '{print$1,$2,$4,$5,$11,$27}'|paste -d '\t' - $out/${name}.5flank.seq $out/${name}.3flank.seq | awk '{OFS="\t"}{print$1,$2,$3,$4,$8,$10,$5,$6}' >$out/${name}.tmp.txt
    sed -i "1i #chr\tpos\tref\talt\t5flank_seq\t3flank_seq\tgene\tstrand" $out/${name}.tmp.txt
done

for var in `find $prefix/group* -name "*.tmp.txt"`;do
    name=`echo $var|awk -F'/' '{print$NF}'|cut -c1-19`
    group=`echo $var|awk -F'/' '{print$7}'`
    out=$prefix/${group}
    $Rscript $getSpectra $out/${name}.tmp.txt
    #cat $out/${name}.ready.txt | grep -v '^#' | awk '{OFS="\t"}{print substr($5,10),$3,substr($6,1,1),$4}'|sort|uniq -c|awk '{OFS="\t"}{print$2,$3,$4,$5,$2$3$4$5,$1}' | join -a1 -1 5 -2 5 -o 1.1 1.2 1.3 1.4 1.6 2.6 /online/home/chenyl/dongfang/mode.xls - | awk '{OFS="\t"}{if($6=="")print$1,$2,$3,$4,$5;else print$1,$2,$3,$4,$6}' >$out/${name}.ready.xls
    cat $out/${name}.ready.txt | grep -v '^#' | awk '{OFS="\t"}{print $7,$8,$9,$10}'|sort|uniq -c|awk '{OFS="\t"}{print$2,$3,$4,$5,$2$3$4$5,$1}' | join -a1 -1 5 -2 5 -o 1.1 1.2 1.3 1.4 1.6 2.6 /online/home/chenyl/dongfang/mode.xls - | awk '{OFS="\t"}{if($6=="")print$1,$2,$3,$4,$5;else print$1,$2,$3,$4,$6}' >$out/${name}.ready.xls
    sed -i "1i before\tref\tafter\talt\t${name}" $out/${name}.ready.xls
done

