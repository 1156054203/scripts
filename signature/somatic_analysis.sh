#!/bin/bash
##somatic mutation analysis by strelka and Mutect and Mutect2

usage="bash $0 user sheet kind"
if [ $# -eq 3 ];then
   user=$1
   sheet=$2
   kind=$3
else 
   echo $usage
   exit 1
fi

group=`echo $sheet|cut -d. -f1`
tumor=(`cat $sheet|sed -n '1~2p'|cut -f1`)
normal=(`cat $sheet|sed -n '2~2p'|cut -f1`)
tpro=(`cat $sheet|sed -n '1~2p'|cut -f3`)
npro=(`cat $sheet|sed -n '2~2p'|cut -f3`)

prefix=/offline/Analysis/WGS
ref=/online/databases/Homo_sapiens/hg38/hg38bundle/Homo_sapiens_assembly38.fasta
bed=/online/databases/Homo_sapiens/hg_exome/Agilent_V6/hg38_Agilent_V6-col6.bed
bed_3col=/online/home/chenyl/software/databases/Agilent_V6-3col.bed
dbsnp=/online/databases/Homo_sapiens/hg38/hg38bundle/dbsnp_144.hg38.vcf.gz
gatk=/online/home/chenyl/software/gatk-4.0.12.0/gatk
strelkaDir=/online/software/strelka_workflow-1.0.15
config=$strelkaDir/etc/strelka_config_bwa_default.ini
java8=/online/software/jdk1.8.0_111/bin/java
muse=/online/home/chenyl/software/MuSE-master/MuSE

outdir=/online/home/chenyl/dongfang/$kind/$group
test -d $outdir||mkdir -p $outdir

if [ $kind = 'strelka' ];then
   for i in `seq 0 $((${#tumor[@]}-1))`;do
       t=${tumor[$i]}
       n=${normal[$i]}
       tp=${tpro[$i]}
       np=${npro[$i]}
       tbam=$prefix/$tp/$t/out/${t}.ready.bam
       nbam=$prefix/$np/$n/out/${n}.ready.bam
       pair=${t}-${n}
      
       #mkdir -p $outdir/$pair
       echo $pair...generating pbs file...

cat >${outdir}/map_${pair}.pbs <<stop
#PBS -N $pair
#PBS -o ${outdir}/${pair}.o
#PBS -e ${outdir}/${pair}.e
#PBS -l nodes=1:ppn=10
#PBS -r y
#PBS -u $user
#PBS -q high

cd ${outdir}
time $strelkaDir/bin/hg38_configureStrelkaWorkflow.pl --normal=$nbam --tumor=$tbam --ref=$ref --config=$config --output-dir=./$pair
cd ./$pair
make -j 8
     
stop
      
      echo $pair...Finished!
      echo "qsub ${outdir}/map_${pair}.pbs" >>${outdir}/${kind}_run.sh
      done
      sed -i '1i#!/bin/bash' ${outdir}/${kind}_run.sh
      chmod 755 ${outdir}/${kind}_run.sh


elif [ $kind = 'muse' ];then
   for i in `seq 0 $((${#tumor[@]}-1))`;do
       t=${tumor[$i]}
       n=${normal[$i]}
       tp=${tpro[$i]}
       np=${npro[$i]}
       tbam=$prefix/$tp/$t/out/${t}.ready.bam
       nbam=$prefix/$np/$n/out/${n}.ready.bam
       pair=${t}-${n}

       mkdir -p $outdir/$pair
       echo $pair...generating pbs file...

cat >${outdir}/map_${pair}.pbs <<stop

#PBS -N $pair
#PBS -o ${outdir}/$pair/${pair}.o
#PBS -e ${outdir}/$pair/${pair}.e
#PBS -l nodes=1:ppn=10
#PBS -r y
#PBS -u $user
#PBS -q high

cd $outdir/$pair
time $muse call \
$tbam \
$nbam \
-f $ref \
-l $bed_3col \
-O ${pair}

time $muse sump \
-I ${pair}.MuSE.txt \
-D $dbsnp \
-E \
-O ${pair}.filter.vcf

stop

      echo $pair...Finished!
      echo "qsub ${outdir}/map_${pair}.pbs" >>${outdir}/${kind}_run.sh
      done
      sed -i '1i#!/bin/bash' ${outdir}/${kind}_run.sh
      chmod 755 ${outdir}/${kind}_run.sh


elif [ $kind = 'gatk' ];then
   for i in `seq 0 $((${#tumor[@]}-1))`;do
       t=${tumor[$i]}
       n=${normal[$i]}
       tp=${tpro[$i]}
       np=${npro[$i]}
       tbam=$prefix/$tp/$t/out/${t}.ready.bam
       nbam=$prefix/$np/$n/out/${n}.ready.bam
       pair=${t}-${n}

       mkdir -p $outdir/$pair
       echo $pair...generating pbs file...

cat >${outdir}/map_${pair}.pbs <<stop

#PBS -N $pair
#PBS -o ${outdir}/$pair/${pair}.o
#PBS -e ${outdir}/$pair/${pair}.e
#PBS -l nodes=1:ppn=10
#PBS -r y
#PBS -u $user
#PBS -q high

cd $outdir/$pair
time $gatk --java-options "-Xmx20g" Mutect2 \
-R $ref \
-I $tbam \
-I $nbam \
-tumor $t \
-normal $n \
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
-L $bed \
-O ${pair}.vcf

time $gatk FilterMutectCalls \
-V ${pair}.vcf \
-O ${pair}.filter.vcf

stop
      echo $pair...Finished!
      echo "qsub ${outdir}/map_${pair}.pbs" >>${outdir}/${kind}_run.sh
      done
      sed -i '1i#!/bin/bash' ${outdir}/${kind}_run.sh
      chmod 755 ${outdir}/${kind}_run.sh
fi
