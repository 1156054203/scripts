#!/bin/bash

usage="Usage:\n\tbash $0 -r/-s\n
      \noptions:\n\t -r\trun  tomcat
      \n          \t -s\tstop tomcat"

if [ $# -ne 1 ];then
   echo -e $usage
   exit 1
elif [[ $1 == '-r' ]];then
   bash /home/chenyl/tomcat/bin/startup.sh
   echo -e "\ntomcat is runnig"
elif [[ $1 == '-s' ]];then
   bash /home/chenyl/tomcat/bin/shutdown.sh
   echo -e "\ntomcat has been shut down"
else
   echo -e $usage
   exit 1
fi
