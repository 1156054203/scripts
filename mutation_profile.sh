#!/bin/bash
Usage="bash $0 project"

prefix=/online/home/chenyl/temp
wd=`pwd`

if [ $# -eq 1 ];then
   pro=$1
else
   echo $Usage
   exit 1
fi

ref_hg19=/online/databases/gatk_bundle_2.8_hg19/ucsc.hg19.fasta
ref_hg38=/online/databases/Homo_sapiens/hg38/hg38bundle/Homo_sapiens_assembly38.fasta

for var in `find $prefix/${pro}/filter/CHG* -maxdepth 1 -name "*anno.filter.vcf"`;do
    name=`echo $var|awk -F'/' '{print$NF}'|cut -d'_' -f1`
    out=$prefix/${pro}/filter/$name
    #head=`sed -n '/#CHROM/p' $var`
    #awk 'BEGIN{OFS="\t"}{if($4!="-"&&$5!="-"&&length($4)==1&&length($5)==1&&$2>=2)print$0}' $var >${out}/${name}.vcf
    #sed -i "1i$head" ${out}/${name}.vcf
    awk '{if(length($4)==1&&length($5)==1&&$2>=2)print$1"\t"$2-2"\t"$2+1}' $var >${out}/${name}.bed
done


for var in `find $prefix/${pro}/filter/CHG* -maxdepth 1 -name "*.bed"`;do
    test -f ${var%.*}.seq&&rm ${var%.*}.seq
    /online/software/bedtools-2.25.0/bin/bedtools getfasta -fi ${ref_hg38} -bed ${var} -tab -fo ${var%.*}.seq
    cut -f2 ${var%.*}.seq|awk 'BEGIN{OFS="\t";print"pre","ref","post"}{print toupper($0)}' >${var%.*}.tmp
    mv ${var%.*}.tmp ${var%.*}.seq
    sed -i -r 's/([ATGC].)/&\t/;s/([ATGC])/&\t/' ${var%.*}.seq
    #name=`echo $var|awk -F '/' '{print$NF}'|cut -d. -f1`
   grep '^chr' ${var%.*}_anno.filter.vcf|awk 'BEGIN{print"alt"}{if(length($5)==1)print$5}'|paste -d'\t' ${var%.*}.seq -|awk 'BEGIN{OFS="\t"}{print$2,$4,$1,$3}' >${var%.*}.info
    rm ${var%.*}.seq
    rm ${var%.*}.bed
done

for var in `find $prefix/${pro}/filter/CHG* -maxdepth 1 -name "*.info"`;do
    total=`grep -v 'ref' $var|wc -l`
    sample=`echo $var|awk -F/ '{print$NF}'|cut -d. -f1`
    grep -v 'ref' $var|awk '{if($1=="A"&&$2=="T")print$0}'|sort -k3,4|uniq -c|awk 'BEGIN{print"type""\t""pre""\t""post""\t""transformation""\t""ratio""\t""sample"}{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\t",$1/"'$total'"*100;print"'$sample'"}' >${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="A"&&$2=="G")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\t",$1/"'$total'"*100;print"'$sample'"}' >>${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="A"&&$2=="C")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\t",$1/"'$total'"*100;print"'$sample'"}' >>${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="T"&&$2=="A")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\t",$1/"'$total'"*100;print"'$sample'"}' >>${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="T"&&$2=="G")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\t",$1/"'$total'"*100;print"'$sample'"}' >>${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="T"&&$2=="C")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\t",$1/"'$total'"*100;print"'$sample'"}' >>${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="G"&&$2=="A")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\t",$1/"'$total'"*100;print"'$sample'"}' >>${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="G"&&$2=="T")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\t",$1/"'$total'"*100;print"'$sample'"}' >>${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="G"&&$2=="C")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\t",$1/"'$total'"*100;print"'$sample'"}' >>${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="C"&&$2=="A")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\t",$1/"'$total'"*100;print"'$sample'"}' >>${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="C"&&$2=="T")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\t",$1/"'$total'"*100;print"'$sample'"}' >>${var%.*}.xls
    grep -v 'ref' $var|awk '{if($1=="C"&&$2=="G")print$0}'|sort -k3,4|uniq -c|awk '{printf"%s\t",$2">"$3;printf"%s\t",$4;printf"%s\t",$5;printf"%s\t",$4$2$5">"$4$3$5;printf"%0.2f\t",$1/"'$total'"*100;print"'$sample'"}' >>${var%.*}.xls
done

for var in `find $prefix/${pro}/filter/CHG* -maxdepth 1 -name "*.xls"`;do
    sample=`echo $var|awk -F/ '{print$NF}'|cut -d. -f1`
    
    cat /online/home/chenyl/temp/mode.xls|sed -n '2,$p'|while read line;do 
         context=`echo $line|awk '{OFS="\t"}{print$2,$3,$4,$5,$6,"'${sample}'"}'`
         index=`echo $line|awk '{print$1}'`
         str=`echo $line|awk '{print$5}'`
         if [[ -z `grep "$str" $var` ]];then
              sed -i "${index}i${context}" $var
         fi
    done
done
