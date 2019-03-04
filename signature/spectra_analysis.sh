#!/bin/bash

usage="bash $0 sheet"
if [ $# -eq 2 ];then
    user=$1
    sheet=$2
else 
   echo $usage
   exit 1
fi

list=(`cat $sheet|sed -n '2,$p'|cut -f1|cut -d- -f3`)
prefix=/online/home/chenyl/dongfang/tcga-maf
bcftools=/online/home/chenyl/software/bcftools-1.9/bcftools
bgzip=/online/home/chenyl/software/htslib-develop/bgzip
tabix=/online/home/chenyl/software/htslib-develop/tabix 
hg38_fa=/online/databases/Homo_sapiens/hg38/hg38bundle/Homo_sapiens_assembly38.fasta
cache=/online/home/wanghn/data/vep
vep_path=/online/home/chenyl/software/ensembl-vep-release-92
bedtool=/online/software/bedtools-2.25.0/bin/bedtools
python=/online/software/python-3.5.1/bin/python3
getSpectra=/online/home/chenyl/dongfang/get_spectra2.py
    

for sample in ${list[@]};do
    echo
    echo "${sample}...generating pbs file!"
    test -d $prefix/intersection/${sample}||mkdir $prefix/intersection/${sample}
    outdir=$prefix/intersection/${sample}

    cat >${prefix}/intersection/run_${sample}_spectra.pbs <<stop

#PBS -N ${sample}
#PBS -o ${outdir}/${sample}.o
#PBS -e ${outdir}/${sample}.e
#PBS -l nodes=1:ppn=12
#PBS -r y
#PBS -u $user

export tabix=/online/software/tabix

time $bcftools view $prefix/LIHC-muse-362pair/tcga.lihc.muse.${sample}.vcf -o $prefix/LIHC-muse-362pair/tcga.lihc.muse.${sample}.vcf.gz -O z
time $bcftools index $prefix/LIHC-muse-362pair/tcga.lihc.muse.${sample}.vcf.gz

wait
time $bcftools view $prefix/LIHC-mutect-364pair/tcga.lihc.mutect.${sample}.vcf -o $prefix/LIHC-mutect-364pair/tcga.lihc.mutect.${sample}.vcf.gz -O z
time $bcftools index $prefix/LIHC-mutect-364pair/tcga.lihc.mutect.${sample}.vcf.gz

wait
time $bcftools view $prefix/LIHC-varscan-364pair/tcga.lihc.varscan.${sample}.vcf -o $prefix/LIHC-varscan-364pair/tcga.lihc.varscan.${sample}.vcf.gz -O z
time $bcftools index $prefix/LIHC-varscan-364pair/tcga.lihc.varscan.${sample}.vcf.gz

wait
time $bcftools isec -n +2 $prefix/LIHC-muse-362pair/tcga.lihc.muse.${sample}.vcf.gz $prefix/LIHC-mutect-364pair/tcga.lihc.mutect.${sample}.vcf.gz $prefix/LIHC-varscan-364pair/tcga.lihc.varscan.${sample}.vcf.gz -O v | awk '{OFS="\t"}{if(\$3!="-"&&\$3!="."&&\$3!="Y"&&\$3!="N"&&\$4!="-"&&\$4!="."&&\$4!="Y"&&\$4!="N"&&length(\$3)==1&&length(\$4)==1)print\$0}' >$outdir/${sample}.isec.vcf

wait
cat $outdir/${sample}.isec.vcf | awk '{OFS="\t"}{print\$1,\$2,".",\$3,\$4,".",\$5}' | $vep_path/vep -i stdin --cache --dir $cache --format vcf --a GRCh38 --offline --merged --fasta ${hg38_fa} --sift b --polyphen b --af_1kg --af_esp --af_gnomad  --fork 4 --vcf --no_stats -o stdout | $vep_path/filter_vep -filter "Consequence is not intergenic_variant" --force -o $outdir/${sample}.isec.anno.vcf

wait
#5 flank_sequence
cat $outdir/${sample}.isec.anno.vcf|grep -v '^#'|awk -F '[\t|]' '{print\$1,\$2,\$27}'|awk 'BEGIN{OFS="\t"}{if(\$3=="+")print\$1,\$2-11,\$2-1,"name",1,\$3;else print\$1,\$2,\$2+10,"name",1,\$3}' >${outdir}/${sample}.5flank.bed
#3 flank_sequence
cat $outdir/${sample}.isec.anno.vcf|grep -v '^#'|awk -F '[\t|]' '{print\$1,\$2,\$27}'|awk 'BEGIN{OFS="\t"}{if(\$3=="-")print\$1,\$2,\$2+10,"name",1,\$3;else print\$1,\$2-11,\$2-1,"name",1,\$3}' >${outdir}/${sample}.3flank.bed

wait
/online/software/bedtools-2.25.0/bin/bedtools getfasta -fi $hg38_fa -bed ${outdir}/${sample}.5flank.bed -tab -s -fo $outdir/${sample}.5flank.seq
/online/software/bedtools-2.25.0/bin/bedtools getfasta -fi $hg38_fa -bed ${outdir}/${sample}.3flank.bed -tab -s -fo $outdir/${sample}.3flank.seq

wait
cat $outdir/${sample}.isec.anno.vcf | grep -v '^#' |awk -F '[\t|]' '{print\$1,\$2,\$4,\$5,\$11,\$27}'|paste -d '\t' - $outdir/${sample}.5flank.seq $outdir/${sample}.3flank.seq | awk '{OFS="\t"}{print\$1,\$2,\$3,\$4,\$8,\$10,\$5,\$6}' >$outdir/${sample}.tmp.txt

sed -i "1i #chr\tpos\tref\talt\t5flank_seq\t3flank_seq\tgene\tstrand" $outdir/${sample}.tmp.txt

$python $getSpectra $outdir/${sample}.tmp.txt

cat $outdir/${sample}.ready.txt | grep -v '^#' | awk '{OFS="\t"}{print \$7,\$8,\$9,\$10}'|sort|uniq -c|awk '{OFS="\t"}{print\$2,\$3,\$4,\$5,\$2\$3\$4\$5,\$1}' | join -a1 -1 5 -2 5 -o 1.1 1.2 1.3 1.4 1.6 2.6 /online/home/chenyl/dongfang/mode.xls - | awk '{OFS="\t"}{if(\$6=="")print\$1,\$2,\$3,\$4,\$5;else print\$1,\$2,\$3,\$4,\$6}' >$outdir/${sample}.ready.xls

sed -i "1i before\tref\tafter\talt\t${name}" $outdir/${sample}.ready.xls
stop

chmod 755 ${prefix}/intersection/run_${sample}_spectra.pbs
echo "qsub ${prefix}/intersection/run_${sample}_spectra.pbs" >>$prefix/intersection/run_spectra.sh

done

sed -i "1i #!/bin/bash" $prefix/intersection/run_spectra.sh
chmod 755 $prefix/intersection/run_spectra.sh
