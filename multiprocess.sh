#!/bin/bash

function a_sub {
echo "$var Start at: `date`"
starttime=$(date +%s)
sleep 3    ## Command to execute
echo "$var Ended at: `date`"
endtime=$(date +%s)
echo "$var,run time: $(expr $endtime - $starttime)"
}

tmp_fifo="/tmp/fifo"
mkfifo $tmp_fifo
exec 11<>$tmp_fifo
rm $tmp_fifo

thread=6
for ((i=0;i<$thread;i++));do
echo >&6
done

for ((i=0;i<50;i++));do
read -u6 
{
      a_sub && {
       echo "a_sub is finished"
      } || {
       echo "sub error"
      }
      echo >&6
} &

 
done

wait
exec 6>&-
exec 6<&-
exit 0
