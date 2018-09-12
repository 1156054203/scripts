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

if [ -d ${prefix}/${pro}/somatic_analysis ];then
  for var in `find ${prefix}/${pro}/somatic_analysis/group*/somatic_*/results -name "all.somatic.snvs.vcf"`;do
      name=`echo $var|cut -d/ -f8`
      group=`echo $var|cut -d/ -f7`
      test -d ${wd}/${pro}/${group}||mkdir -p ${wd}/${pro}/${group}
      out=${wd}/${pro}/${group}
total=`awk 'length($4)==1&&length($5)==1' $var|wc -l`
AT=`awk 'BEGIN{sum=0}{if($4=="A"&&$5=="T")sum+=1}END{print sum}' $var`
AG=`awk 'BEGIN{sum=0}{if($4=="A"&&$5=="G")sum+=1}END{print sum}' $var`
AC=`awk 'BEGIN{sum=0}{if($4=="A"&&$5=="C")sum+=1}END{print sum}' $var`
TA=`awk 'BEGIN{sum=0}{if($4=="T"&&$5=="A")sum+=1}END{print sum}' $var`
TG=`awk 'BEGIN{sum=0}{if($4=="T"&&$5=="G")sum+=1}END{print sum}' $var`
TC=`awk 'BEGIN{sum=0}{if($4=="T"&&$5=="C")sum+=1}END{print sum}' $var`
CG=`awk 'BEGIN{sum=0}{if($4=="C"&&$5=="G")sum+=1}END{print sum}' $var`
CA=`awk 'BEGIN{sum=0}{if($4=="C"&&$5=="A")sum+=1}END{print sum}' $var`
CT=`awk 'BEGIN{sum=0}{if($4=="C"&&$5=="T")sum+=1}END{print sum}' $var`
GC=`awk 'BEGIN{sum=0}{if($4=="G"&&$5=="C")sum+=1}END{print sum}' $var`
GA=`awk 'BEGIN{sum=0}{if($4=="G"&&$5=="A")sum+=1}END{print sum}' $var`
GT=`awk 'BEGIN{sum=0}{if($4=="G"&&$5=="T")sum+=1}END{print sum}' $var`
point=`awk 'BEGIN{sum=0}{if($4=="."||$5==".")sum+=1}END{print sum}' $var`

TiTv=`echo "scale=2;($AG+$GA+$CT+$TC)/($AT+$TA+$AC+$CA+$GT+$TG+$CG+$GC)"|bc`

echo "Total SNP	$total" >${out}/${name}.xls
echo "Ti/Tv=	$TiTv" >>${out}/${name}.xls
echo "A->T ratio	`echo $AT $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "A->G ratio	`echo $AG $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "A->C ratio	`echo $AC $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "T->A ratio	`echo $TA $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "T->G ratio	`echo $TG $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "T->C ratio	`echo $TC $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "C->G ratio	`echo $CG $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "C->A ratio	`echo $CA $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "C->T ratio	`echo $CT $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "G->C ratio	`echo $GC $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "G->A ratio	`echo $GA $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "G->T ratio	`echo $GT $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "point ratio	`echo $point $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
done
else
    for var in `find ${prefix}/${pro}/CHG*/out -maxdepth 1 -name "CHG*.ready.snp.tsv.xls"`;do
        name=`echo $var|cut -d/ -f6`
        test -d ${wd}/${pro}||mkdir -p ${wd}/${pro}
        out=${wd}/${pro}
total=`awk 'length($4)==1&&length($5)==1' $var|wc -l`
AT=`awk 'BEGIN{sum=0}{if($4=="A"&&$5=="T")sum+=1}END{print sum}' $var`
AG=`awk 'BEGIN{sum=0}{if($4=="A"&&$5=="G")sum+=1}END{print sum}' $var`
AC=`awk 'BEGIN{sum=0}{if($4=="A"&&$5=="C")sum+=1}END{print sum}' $var`
TA=`awk 'BEGIN{sum=0}{if($4=="T"&&$5=="A")sum+=1}END{print sum}' $var`
TG=`awk 'BEGIN{sum=0}{if($4=="T"&&$5=="G")sum+=1}END{print sum}' $var`
TC=`awk 'BEGIN{sum=0}{if($4=="T"&&$5=="C")sum+=1}END{print sum}' $var`
CG=`awk 'BEGIN{sum=0}{if($4=="C"&&$5=="G")sum+=1}END{print sum}' $var`
CA=`awk 'BEGIN{sum=0}{if($4=="C"&&$5=="A")sum+=1}END{print sum}' $var`
CT=`awk 'BEGIN{sum=0}{if($4=="C"&&$5=="T")sum+=1}END{print sum}' $var`
GC=`awk 'BEGIN{sum=0}{if($4=="G"&&$5=="C")sum+=1}END{print sum}' $var`
GA=`awk 'BEGIN{sum=0}{if($4=="G"&&$5=="A")sum+=1}END{print sum}' $var`
GT=`awk 'BEGIN{sum=0}{if($4=="G"&&$5=="T")sum+=1}END{print sum}' $var`
point=`awk 'BEGIN{sum=0}{if($4=="."||$5==".")sum+=1}END{print sum}' $var`

TiTv=`echo "scale=2;($AG+$GA+$CT+$TC)/($AT+$TA+$AC+$CA+$GT+$TG+$CG+$GC)"|bc`

echo "Total SNP	$total" >${out}/${name}.xls
echo "Ti/Tv=	$TiTv" >>${out}/${name}.xls
echo "A->T ratio	`echo $AT $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "A->G ratio	`echo $AG $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "A->C ratio	`echo $AC $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "T->A ratio	`echo $TA $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "T->G ratio	`echo $TG $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "T->C ratio	`echo $TC $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "C->G ratio	`echo $CG $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "C->A ratio	`echo $CA $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "C->T ratio	`echo $CT $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "G->C ratio	`echo $GC $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "G->A ratio	`echo $GA $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "G->T ratio	`echo $GT $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
echo "point ratio	`echo $point $total|awk '{printf"%0.2f\n",$1/$2*100}'`" >>${out}/${name}.xls
done

fi
