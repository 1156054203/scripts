#!/bin/bash
Usage="bash $0 project"

prefix=/offline/Analysis/WGS
wd=`pwd`

if [ $# -eq 1 ];then
   pro=$1
else
   echo $Usage
   exit 1
fi


for var in `find ${prefix}/${pro}/CHG*/out -maxdepth 1 -name "CHG*.ready.snp.tsv.xls"`;do
    name=`echo $var|cut -d/ -f6`
    test -d ${wd}/${pro}||mkdir -p ${wd}/${pro}
    out=${wd}/${pro}
    grep '^chr' $var|awk '{if(length($5)==1&&$2>=2)print$1"\t"$2-2"\t"$2+1}' >${out}/${name}.bed
done

for var in `find ${pro} -maxdepth 1 -name "*.bed"`;do
    test -f ${var%.*}.seq&&rm ${var%.*}.seq 
    /online/software/bedtools-2.25.0/bin/bedtools getfasta -fi /online/databases/Homo_sapiens/hg38/hg38bundle/Homo_sapiens_assembly38.fasta -bed ${var} -tab -fo ${var%.*}.seq
    cut -f2 ${var%.*}.seq|awk 'BEGIN{OFS="\t";print"pre","ref","post"}{print toupper($0)}' >${var%.*}.tmp
    mv ${var%.*}.tmp ${var%.*}.seq
    sed -i -r 's/([ATGC].)/&\t/;s/([ATGC])/&\t/' ${var%.*}.seq
    name=`echo $var|cut -d/ -f2`
    grep '^chr' ${prefix}/${pro}/${name%.*}/out/${name%.*}.ready.snp.tsv.xls|awk 'BEGIN{print"alt"}{if(length($5)==1&&$2>=2)print$5}'|paste -d'\t' ${var%.*}.seq -|awk 'BEGIN{OFS="\t"}{print$2,$4,$1,$3}' >${var%.*}.info
    rm ${var%.*}.seq
    rm ${var%.*}.bed
done

for var in `find ${wd}/${pro} -maxdepth 1 -name "*.info"`;do
    total=`grep -v 'ref' $var|wc -l`
    grep -v 'ref' $var|awk '{if($1=="A"&&$2=="T")print$0}'|sort -k3,4|uniq -c|awk 'BEGIN{print"type""\t""pre""\t""post""\t""transformation""\t""ratio"}{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\n",$1/"'$total'"*100}' >${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="A"&&$2=="G")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\n",$1/"'$total'"*100}' >>${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="A"&&$2=="C")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\n",$1/"'$total'"*100}' >>${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="T"&&$2=="A")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\n",$1/"'$total'"*100}' >>${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="T"&&$2=="G")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\n",$1/"'$total'"*100}' >>${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="T"&&$2=="C")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\n",$1/"'$total'"*100}' >>${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="G"&&$2=="A")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\n",$1/"'$total'"*100}' >>${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="G"&&$2=="T")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\n",$1/"'$total'"*100}' >>${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="G"&&$2=="C")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\n",$1/"'$total'"*100}' >>${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="C"&&$2=="A")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\n",$1/"'$total'"*100}' >>${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="C"&&$2=="T")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\n",$1/"'$total'"*100}' >>${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="C"&&$2=="G")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\n",$1/"'$total'"*100}' >>${var%.*}.xls
done
