#!/bin/bash

usage="Usage:\n\tbash $0 -r/-s [idxdir]\n
      \noptions:\n\t -r\trun  tomcat
      \n          \t -s\tstop tomcat
      \n          \t idxdir\tthe absolute path of directory for index.html"

tomcatdir=/home/chenyl/tomcat
if [ $# -eq 2 ];then
   sed -i "s#docBase=\(\".*\"\) #docBase=\"$2\" #g" $tomcatdir/conf/server.xml
   if [[ $1 == '-r' ]];then
      bash $tomcatdir/bin/startup.sh
      echo -e "\ntomcat is runnig"
   elif [[ $1 == '-s' ]];then
      bash $tomcatdir/bin/shutdown.sh
      echo -e "\ntomcat has been shut down"
   else
      echo -e $usage
      exit 1
   fi
elif [ $# -eq 1 ];then
   if [[ $1 == '-r' ]];then
      bash $tomcatdir/bin/startup.sh
      echo -e "\ntomcat is runnig"
   elif [[ $1 == '-s' ]];then
      bash $tomcatdir/bin/shutdown.sh
      echo -e "\ntomcat has been shut down"
   else
      echo -e $usage
      exit 1
   fi
else
   echo -e $usage
   exit 1
fi
