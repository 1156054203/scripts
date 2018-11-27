#!/bin/bash
Usage="bash $0 project hg19"

prefix=/online/home/chenyl/temp/others
wd=`pwd`

if [ $# -eq 2 ];then
   pro=$1
   tp=$2
else
   echo $Usage
   exit 1
fi

ref_hg19=/online/databases/gatk_bundle_2.8_hg19/ucsc.hg19.fasta
ref_hg38=/online/databases/Homo_sapiens/hg38/hg38bundle/Homo_sapiens_assembly38.fasta

WES_hg19=`cat /online/home/chenyl/temp/WES/wes-hg19|grep hg19|awk -F. '{print$1}'`
WES_hg38=`cat /online/home/chenyl/temp/WES/wes-hg38|grep hg38|awk -F. '{print$1}'`
WGS_hg19=`cat /online/home/chenyl/temp/WGS/wgs-hg19|grep hg19|awk -F. '{print$1}'`
WGS_hg38=`cat /online/home/chenyl/temp/WGS/wgs-hg38|grep hg38|awk -F. '{print$1}'`

bed=(snp.vcf norm.vcf)
WES=${bed[0]}
WGS=${bed[1]}

for var in $(eval echo \$${pro}_${tp});do
    sample=/online/home/chenyl/temp/others/${var}.$(eval echo \$${pro})
    name=${var}_${pro}
    test -d ${wd}/${pro}||mkdir -p ${wd}/${pro}
    out=${wd}/${pro}
    head=`sed -n '/#CHROM/p' $sample`
    awk 'BEGIN{OFS="\t"}{if($4!="-"&&$5!="-"&&$4!="."&&$5!="."&&length($4)==1&&length($5)==1&&$2>=2)print$0}' $sample >${out}/${name}.vcf
    sed -i "1i$head" ${out}/${name}.vcf
    awk '{if(length($4)==1&&length($5)==1&&$2>=2)print$1"\t"$2-2"\t"$2+1}' ${out}/${name}.vcf >${out}/${name}.bed
done


for sample in $(eval echo \$${pro}_${tp});do
    var=`find ${pro} -maxdepth 1 -name "${sample}*.bed"`
    test -f ${var%.*}.seq&&rm ${var%.*}.seq
    /online/software/bedtools-2.25.0/bin/bedtools getfasta -fi $(eval echo \$ref_${tp}) -bed ${var} -tab -fo ${var%.*}.seq  ## 4 Modify the reference genome
    cut -f2 ${var%.*}.seq|awk 'BEGIN{OFS="\t";print"pre","ref","post"}{print toupper($0)}' >${var%.*}.tmp
    mv ${var%.*}.tmp ${var%.*}.seq
    sed -i -r 's/([ATGC].)/&\t/;s/([ATGC])/&\t/' ${var%.*}.seq
    name=`echo $var|cut -d/ -f2|cut -d. -f1`
    grep -v '^#' ${out}/${name}.vcf|awk 'BEGIN{print"alt"}{if(length($4)==1&&length($5)==1&&$2>=2)print$5}'|paste -d'\t' ${var%.*}.seq -|awk 'BEGIN{OFS="\t"}{print$2,$4,$1,$3}' >${var%.*}.info
    sed -i '/N\|Y/d' ${var%.*}.info
    rm ${var%.*}.seq
    rm ${var%.*}.bed
done

for sample in $(eval echo \$${pro}_${tp});do
    var=`find ${pro} -maxdepth 1 -name "${sample}*.info"`
    total=`grep -v 'ref' $var|wc -l`
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

for sample in $(eval echo \$${pro}_${tp});do
    var=`find ${prefix}/${sample}_${pro} -name "${sample}*.xls"`
    cat /online/home/chenyl/temp/mode.xls|sed -n '2,$p'|while read line;do
         context=`echo $line|awk '{OFS="\t"}{print$2,$3,$4,$5,$6,"'${sample}'"}'`
         index=`echo $line|awk '{print$1}'`
         str=`echo $line|awk '{print$5}'`
         if [[ -z `grep "$str" $var` ]];then
            if [[ -z `sed -n "${index}p" $var` ]];then
                let num=$index-1
                sed -i "${num}a${context}" $var
            else
                sed -i "${index}i${context}" $var
            fi
         fi
    done
done
