#!/bin/bash
set -o errexit
###################################################################################
help_file="\tWe could get n pairs of simulated data through running this script and this script is a supplementary for simulation.sh\n
\tThe basic format of run this script is like below:\n
\tsh run_simulation.sh -o ${outdir} -n ${number}\n"
 
TEMP=`getopt -o o:t:n:h --long outdir:,time:,number:,help -- "$@"` 
eval set -- "$TEMP"
while true ; do
    case "$1" in
	-o|--outdir)
		outdir=$2 ; shift 2 ;;
	-n|--number)
		number=$2 ; shift 2 ;;
    -h|--help)
		echo -e ${help_file}; exit 1 ;;
	--) shift ; break ;;
    *) echo "Internal error!" ; exit 1 ;;
    esac
done

number=${number:-4}
outdir=${outdir:-"/data/sixone/lllab/DMR/BiB_final/pROC/simulation/1"}
######################################################################################
count_wait(){
    limit=$1
    count=$2
    if [ $count == $limit ];then
        count=0
        wait
    fi
    let count=count+1
}
#####################################################################################
count=1
for((time=1;time<=${number};time++));
do 
    if [ ! -f ${outdir} ];then mkdir -p ${outdir}; fi
    sh /data/sixone/lllab/DMR/BiB_final/pROC/script/multiple_simulation.sh -o ${outdir} -t ${time}&
    #limit the number of process at a time
    count_wait 2 ${count}
done 



 