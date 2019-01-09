#!/bin/bash
Usage="bash $0 project"

if [ $# -eq 1 ];then
   pro=$1
else
   echo $Usage
   exit 1
fi

prefix=/offline/Analysis/WGS
outdir=/online/home/chenyl/dongfang/$pro

for var in `find $prefix/${pro}/somatic_analysis/group*/somatic_*/results -maxdepth 1 -name "all.somatic.snvs.vcf"`;do
    name=`echo $var|cut -d/ -f8`
    group=`echo $var|cut -d/ -f7`
    test -d ${outdir}/${group}||mkdir -p ${outdir}/${group}
    out=${outdir}/${group}
    head=`sed -n '/#CHROM/p' $var`
    awk 'BEGIN{OFS="\t"}{if($4!="-"&&$4!="."&&$4!="Y"&&$4!="N"&&$5!="-"&&$5!="."&&$5!="Y"&&$5!="N"&&length($4)==1&&length($5)==1&&$2>=2)print$0}' $var >${out}/${name}.vcf
    sed -i "1i$head" ${out}/${name}.vcf
done
