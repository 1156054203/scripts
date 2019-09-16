#!/bin/bash

while :
do
  ps -p 49780
  if [ $? -eq 0 ];then  
     echo "multiprocess bwamapping is running"
  else
     echo -e "multiprocess bwamapping is done\n"
     break
  fi

  sleep 600
done


trap "exec 112>&-;exec 112<&-;exit 0" 2
tmp_fifo="/tmp/fifo"
mkfifo $tmp_fifo
exec 112<>$tmp_fifo
rm $tmp_fifo

thread=8
for ((i=0;i<$thread;i++));do
echo >&112
done

cat list.txt|while read var;do
sample=`echo $var|cut -d' ' -f1`
callpeakip=`echo $var|awk '{printf("macs2 callpeak --nomodel -f BAMPE --keep-dup 1 -g 2.7e9 -t bamfile/%s.unique.bam -c bamfile/%s.unique.bam -n %sVS%s -q 0.05 -B --SPMR\n",$1,$2,$1,$2)}'`
callpeakinput=`echo $var|awk '{printf("macs2 callpeak --nomodel -f BAMPE --keep-dup 1 -g 2.7e9 -t bamfile/%s.unique.bam -n %s_input -q 0.05 -B --SPMR\n",$2,$1)}'`
casetobw=`echo $var|awk '{printf("script/bedGraphToBigWig %sVS%s_treat_pileup.bdg reflen1/hg38_gencode.len %s.treat.bw\n",$1,$2,$1)}'`
ctrltobw=`echo $var|awk '{printf("script/bedGraphToBigWig %s_input_treat_pileup.bdg reflen1/hg38_gencode.len %s.control.bw\n",$1,$1)}'`
rmbdgip=`echo $var|awk '{printf("bash script/rm_bdg.sh %sVS%s_control_lambda.bdg %sVS%s_treat_pileup.bdg\n",$1,$2,$1,$2)}'`
rmbdginput=`echo $var|awk '{printf("bash script/rm_bdg.sh %s_input_control_lambda.bdg %s_input_treat_pileup.bdg\n",$1,$1)}'`
rmbedinput=`echo $var|awk '{printf("bash script/rm_bdg.sh %s_input_peaks.narrowPeak %s_input_summits.bed\n",$1,$1)}'`
toxls=`echo $var|awk '{printf("bash script/simple.sh %sVS%s_peaks.xls %sVS%s_peaks.narrowPeak %sVS%s_summits.bed\n",$1,$2,$1,$2,$1,$2)}'`
cutbed=`echo $var|awk '{printf("bash script/cut.sh %sVS%s_peaks.narrowPeak %sVS%s.bed 9\n",$1,$2,$1,$2)}'`
addheader=`echo $var|awk '{printf("python script/addHeader.py %sVS%s_peaks.narrowPeak > exampler.xls\n",$1,$2)}'`
read -u112
echo "$sample Start at: `date`"
starttime=$(date +%s)
{
      { $callpeakip
        $callpeakinput
        $casetobw
        $ctrltobw
        $rmbdgip
        $rmbdginput
        $rmbedinput
        $toxls
        $cutbed
        $addheader
      } && {
       echo "peakcalling is finished"
      } || {
       echo "peakcalling is error"
      }
      echo >&112
} &
endtime=$(date +%s)
echo "$sample Ended at: `date`"
echo "$sample,run time: $(expr $endtime - $starttime)"
echo

done

wait
exec 112>&-
exec 112<&-

python script/countPeak.py `echo *.narrowPeak|sed 's/ /,/g'` > All_peakNum.xls

echo
echo 'All Down !!!'
exit 0
