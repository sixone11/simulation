#!/bin/bash
set -o errexit
########################################################################
##The workflow of this script
#Step 1: Run RRBSsim
#Step 2: Sequence stitching
#Step 3: Trim and Fastqc
#Step 4: Run Bismark
#Step 5: Convert coverage file to bed file
########################################################################
help_file="\tWe could get the 1 test simulation data and 1 control simulation data through running this script\n
\tThe basic format of run this script is like below:\n
\tsh simulation.sh -o ${outdir} -t ${time}\n"
TEMP=`getopt -o o:t:h --long outdir:,time:,help -- "$@"` 
eval set -- "$TEMP"
while true ; do
    case "$1" in
	-o|--outdir)
		outdir=$2 ; shift 2 ;;
    -t|--time)
		time=$2 ; shift 2 ;;
	-h|--help)
		echo -e ${help_file}; exit 1 ;;
	--) shift ; break ;;
    *) echo "Internal error!" ; exit 1 ;;
    esac
done

#Parameter initialization
outdir=${outdir:-"output"}
time=${time:-1}
cd ${outdir}

###############################################################################
#The assignment of path 

fa=../input
RRBSsim=../input/RRBSsim-master
#Please specify the path by yourself
GRCh37=/data/GRCh37
###############################################################################
#Create the outfile  folder
mkdir -p  out
###############################################################################

#Step 1: Run RRBSsim
mkdir -p rrbssim
#Step 1.1: The simulation of treatment group
python2 ${RRBSsim}/RRBSsim -f ${fa}/test.fa/methy_level_0.0.fa --CG_level 0.05 --SAM -R -d 30 -o test${time}_0.0 --seed ${RANDOM} --multicore 5&

for i in {1..9}
do 
python2 ${RRBSsim}/RRBSsim -f ${fa}/test.fa/methy_level_0.${i}.fa --CG_level 0.${i} --SAM -R -d 30 -o test${time}_0.${i} --seed ${RANDOM} --multicore 5&
done

python2 ${RRBSsim}/RRBSsim -f ${fa}/test.fa/methy_level_1.0.fa --CG_level 0.95 --SAM -R -d 30 -o test${time}_1.0 --seed ${RANDOM} --multicore 5&

wait 

#Step 1.2: The simulation of control group 
python2 ${RRBSsim}/RRBSsim -f ${fa}/control.fa/methy_level_0.0.fa --CG_level 0.05 --SAM -R -d 30 -o control${time}_0.0 --seed ${RANDOM} --multicore 5&

for i in {1..9}
do 
	python2 ${RRBSsim}/RRBSsim -f ${fa}/control.fa/methy_level_0.${i}.fa --CG_level 0.${i} --SAM -R -d 30 -o control${time}_0.${i} --seed ${RANDOM} --multicore 5&
done

python2 ${RRBSsim}/RRBSsim -f ${fa}/control.fa/methy_level_1.0.fa --CG_level 0.95 --SAM -R -d 30 -o control${time}_1.0 --seed ${RANDOM} --multicore 5&

wait 

#Step 2: Sequence stitching
#Step 2.1: In treatment group
mkdir -p  fq
cat rrbssim/test${time}_*.1.fq >fq/test${time}.1.fq 
cat rrbssim/test${time}_*.2.fq >fq/test${time}.2.fq 
#Step2.2: In control group
cat rrbssim/control${time}_*.1.fq >fq/control${time}.1.fq 
cat rrbssim/control${time}_*.2.fq >fq/control${time}.2.fq 
wait 

#######################################################################################
#Step 3: Trim and Fastqc

#Step 3.1: Fastqc
#Step 3.1.1: In treatment group 
mkdir -p  fastqc_initial
fastqc -o fastqc_initial --noextract -t 10 fq/test${time}.1.fq fq/test${time}.2.fq>out/fastqc_initial_test${time}.out;

#Step3.1.2: In control group 
fastqc -o fastqc_initial --noextract -t 10 fq/control${time}.1.fq fq/control${time}.2.fq>out/fastqc_initial_control${time}.out;


#Step 3.2: Trim
mkdir -p  trim_fq

cd trim_fq
trim_galore -a GATCGGAAGAGCA -a2 GATCGGAAGAGCA  --rrbs  --paired ../fq/test${time}.1.fq   ../fq/test${time}.2.fq &

trim_galore -a GATCGGAAGAGCA -a2 GATCGGAAGAGCA  --rrbs  --paired ../fq/control${time}.1.fq   ../fq/control${time}.2.fq &


wait 
cd .. 

#Step 3.3: Fqstqc again
#Step 3.3.1: In treatment group
mkdir -p  fastqc_trim
fastqc -o fastqc_trim --noextract -t 10 trim_fq/test${time}.1_val_1.fq ./trim_fq/test${time}.2_val_2.fq>out/fastqc_trim_test${time}.out;

#Step 3.3.2: In control group
fastqc -o fastqc_trim --noextract -t 10 trim_fq/control${time}.1_val_1.fq trim_fq/control${time}.2_val_2.fq>out/fastqc_trim_control${time}.out;

#######################################################################################
#Step 4: Run Bismark
#Step 4.1: generate bam file
mkdir -p  bam
cd bam
#test group 
bismark ${GRCh37} -1 ../trim_fq/test${time}.1_val_1.fq -2 ../trim_fq/test${time}.2_val_2.fq --bowtie2 --multicore  5 > ../out/bismark_bam_test${time}.out;
#control group 
bismark ${GRCh37} -1 ../trim_fq/control${time}.1_val_1.fq -2 ../trim_fq/control${time}.2_val_2.fq --bowtie2 --multicore  5 > ../out/bismark_bam_control${time}.out;

#Bismark 4.2: generate cov file
cd ..
mkdir -p  cov

#test group 
bismark_methylation_extractor -p --gzip --no_overlap --genome_folder ${GRCh37} --bedGraph -o cov --multicore 10 bam/test${time}.1_val_1_bismark_bt2_pe.bam > out/bismark_cov_test${time}.out ; 
gunzip cov/test${time}.1_val_1_bismark_bt2_pe.bismark.cov.gz;
mv cov/test${time}.1_val_1_bismark_bt2_pe.bismark.cov cov/test${time}.cov;

#control group
bismark_methylation_extractor -p --gzip --no_overlap --genome_folder ${GRCh37} --bedGraph -o cov --multicore 10 bam/control${time}.1_val_1_bismark_bt2_pe.bam > out/bismark_cov_control${time}.out ; 
gunzip cov/control${time}.1_val_1_bismark_bt2_pe.bismark.cov.gz;
mv cov/control${time}.1_val_1_bismark_bt2_pe.bismark.cov cov/control${time}.cov;
###############################################################################
#Step 5: Convert coverage file to bed file
#Cov file to bed file and basic process of bed file
mkdir -p  bed
awk 'BEGIN{FS="\t";OFS="\t"} {print $1,$2,$2+1,$4,$5,$6}' cov/test${time}.cov|sed '/chrX/d;/chrY/d;/chrM/d'|grep chr21 |sort -k1,1V -k2,2n -k3,3n > bed/test${time}.bed;

awk 'BEGIN{FS="\t";OFS="\t"} {print $1,$2,$2+1,$4,$5,$6}' cov/control${time}.cov|sed '/chrX/d;/chrY/d;/chrM/d'|grep chr21 |sort -k1,1V -k2,2n -k3,3n > bed/control${time}.bed;
#################################################################################














