#!/bin/bash

Usage="bash $0 chenyl project"

if [ $# -eq 2 ];then
   user=$1
   pro=$2
else
   echo $Usage
   exit 1
fi

vcftools=/online/software/vcftools-v0.1.14/bin/vcftools
vep_dir=/online/home/chenyl/software/ensembl-vep-release-92
agilent_v6=/online/databases/Homo_sapiens/hg_exome/Agilent_V6/hg38_Agilent_V6-col6.bed
cache=/online/home/wanghn/data/vep
hg38_fa=/online/databases/Homo_sapiens/hg38/hg38bundle/Homo_sapiens_assembly38.fasta
hg19_fa=/online/software/speedseq/genomes/ucsc_ref_hg19/ucsc.hg19.fasta
rna_fa=/online/home/chenyl/software/databases/Homo_sapiens.GRCh38.dna.primary_assembly.fa

prefix=/online/home/chenyl/dongfang

for var in `find ${prefix}/${pro}/group* -maxdepth 1 -name "*.vcf"`;do
    sample=`echo $var|awk -F'/' '{print$NF}'|cut -d'.' -f1`
    echo
    echo "${sample}...generating pbs file!"
    outdir=$(dirname $var)
    test -d $outdir||mkdir -p $outdir
   
    cat >${outdir}/run_${sample}_anno.pbs <<stop

#PBS -N ${sample}
#PBS -o ${outdir}/${sample}.o
#PBS -e ${outdir}/${sample}.e
#PBS -l nodes=1:ppn=12
#PBS -r y
#PBS -u $user

export tabix=/online/software/tabix

time $vcftools --vcf $var --min-meanDP 4 --bed $agilent_v6 --recode --recode-INFO-all --stdout | $vep_dir/vep -i stdin --cache --dir $cache --format vcf --a GRCh38 --offline --merged --fasta ${hg38_fa} -o stdout --sift b --polyphen b --af_1kg --af_esp --af_gnomad  --fork 4 --force --vcf --no_stats | grep -v '^#'|awk '{OFS="\t";FS="[\t;|]"}{if(\$34=="1"||\$34=="-1")print\$1,\$2,\$4,\$5,\$18,\$34}' >${var%%.*}.anno.txt

sed -i "1i #chr\tpos\tref\talt\tgene\tstrand" ${var%%.*}.anno.txt

time $vcftools --vcf $var --min-meanDP 4 --bed $agilent_v6 --recode --recode-INFO-all --stdout | $vep_dir/vep -i stdin --cache --dir $cache --format vcf --a GRCh38 --offline --merged --fasta ${hg38_fa} -o stdout --sift b --polyphen b --af_1kg --af_esp --af_gnomad  --fork 4 --force --vcf --no_stats | $vep_dir/filter_vep --filter "SIFT < 0.05 or PolyPhen > 0.15 and EAS_AF < 0.05 or gnomAD_EAS_AF < 0.05" -o stdout | grep -v '^#' | awk '{OFS="\t";FS="[\t;|]"}{if(\$34=="1"||\$34=="-1")print\$1,\$2,\$4,\$5,\$18,\$34}' >${var%%.*}.filter.txt

sed -i "1i #chr\tpos\tref\talt\tgene\tstrand" ${var%%.*}.filter.txt

stop

chmod 755 ${outdir}/run_${sample}_anno.pbs
echo "qsub run_${sample}_anno.pbs" >>${outdir}/run_${pro}_anno.sh

done

sed -i "1i #!/bin/bash" $prefix/$pro/group12/run_${pro}_anno.sh
sed -i "1i #!/bin/bash" $prefix/$pro/group34/run_${pro}_anno.sh

chmod 755 $prefix/$pro/group12/run_${pro}_anno.sh
chmod 755 $prefix/$pro/group34/run_${pro}_anno.sh
