#!/bin/bash
set -o errexit
###################################################################################
help_file="sh -p seq default:seq -o outdir -s sample_numbers default:(2 4 6 8)\n
 -d {mean_depths} (array) default:(1 5 10 20 30) -r ${rates[@]} default:rates=(0.05 0.1 0.15 0.2)"

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
outdir=${outdir:-/data/sixone/lllab/DMR/BiB_final/pROC/1}
time=${time:-1}
cd ${outdir}
##################################################################################
#fa=/data/sixone/lllab/DMR/BiB_final/custom_code/fa;
fa=/data/sixone/lllab/DMR/BiB_final/code/data
#RRBSsim=/data/sixone/software/DMR/RRBSsim-master;
RRBSsim=/data/sixone/lllab/DMR/BiB_final/duplicate_header/RRBSsim-master
GRCh37=/data/GRCh37;
####################################################################################
#outfile 
mkdir -p  out
###################################################################################

#The simulation of different types of DMRs

mkdir -p rrbssim
#The simulation of DMRs in treatment  
python2 ${RRBSsim}/RRBSsim -f ${fa}/test.fa/methy_level_0.0.fa --CG_level 0.05 --SAM -R -d 30 -o test${time}_0.0 --seed ${RANDOM} --multicore 5&

for i in {1..9}
do 
python2 ${RRBSsim}/RRBSsim -f ${fa}/test.fa/methy_level_0.${i}.fa --CG_level 0.${i} --SAM -R -d 30 -o test${time}_0.${i} --seed ${RANDOM} --multicore 5&
done

python2 ${RRBSsim}/RRBSsim -f ${fa}/test.fa/methy_level_1.0.fa --CG_level 0.95 --SAM -R -d 30 -o test${time}_1.0 --seed ${RANDOM} --multicore 5&

wait 

#control group
python2 ${RRBSsim}/RRBSsim -f ${fa}/control.fa/methy_level_0.0.fa --CG_level 0.05 --SAM -R -d 30 -o control${time}_0.0 --seed ${RANDOM} --multicore 5&

for i in {1..9}
do 
	python2 ${RRBSsim}/RRBSsim -f ${fa}/control.fa/methy_level_0.${i}.fa --CG_level 0.${i} --SAM -R -d 30 -o control${time}_0.${i} --seed ${RANDOM} --multicore 5&
done

python2 ${RRBSsim}/RRBSsim -f ${fa}/control.fa/methy_level_1.0.fa --CG_level 0.95 --SAM -R -d 30 -o control${time}_1.0 --seed ${RANDOM} --multicore 5&

wait 

#Sequence stitching
#test group
mkdir -p  fq
cat rrbssim/test${time}_*.1.fq >fq/test${time}.1.fq 
cat rrbssim/test${time}_*.2.fq >fq/test${time}.2.fq 
#control group
cat rrbssim/control${time}_*.1.fq >fq/control${time}.1.fq 
cat rrbssim/control${time}_*.2.fq >fq/control${time}.2.fq 
wait 

#######################################################################################
#Trim and Fastqc

#Fastqc
#test group 
mkdir -p  fastqc_initial
fastqc -o fastqc_initial --noextract -t 10 fq/test${time}.1.fq fq/test${time}.2.fq>out/fastqc_initial_test${time}.out;

#control group 
fastqc -o fastqc_initial --noextract -t 10 fq/control${time}.1.fq fq/control${time}.2.fq>out/fastqc_initial_control${time}.out;


#Trim
mkdir -p  trim_fq

cd trim_fq
trim_galore -a GATCGGAAGAGCA -a2 GATCGGAAGAGCA  --rrbs  --paired ../fq/test${time}.1.fq   ../fq/test${time}.2.fq &

trim_galore -a GATCGGAAGAGCA -a2 GATCGGAAGAGCA  --rrbs  --paired ../fq/control${time}.1.fq   ../fq/control${time}.2.fq &

wait 
cd .. 

#Fqstqc again

#test group
mkdir -p  fastqc_trim
fastqc -o fastqc_trim --noextract -t 10 trim_fq/test${time}.1_val_1.fq ./trim_fq/test${time}.2_val_2.fq>out/fastqc_trim_test${time}.out;

#control group 
fastqc -o fastqc_trim --noextract -t 10 trim_fq/control${time}.1_val_1.fq trim_fq/control${time}.2_val_2.fq>out/fastqc_trim_control${time}.out;

#######################################################################################
#Bismark 1:generate bam file
mkdir -p  bam
cd bam
#test group 
bismark ${GRCh37} -1 ../trim_fq/test${time}.1_val_1.fq -2 ../trim_fq/test${time}.2_val_2.fq --bowtie2 --multicore  5 > ../out/bismark_bam_test${time}.out;
#control group 
bismark ${GRCh37} -1 ../trim_fq/control${time}.1_val_1.fq -2 ../trim_fq/control${time}.2_val_2.fq --bowtie2 --multicore  5 > ../out/bismark_bam_control${time}.out;

#Bismark 2:generate cov file
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
#######################################################################################
#Cov file to bed file and basic process of bed file
mkdir -p  bed
awk 'BEGIN{FS="\t";OFS="\t"} {print $1,$2,$2+1,$4,$5,$6}' cov/test${time}.cov|sed '/chrX/d;/chrY/d;/chrM/d'|grep chr21 |sort -k1,1V -k2,2n -k3,3n > bed/test${time}.bed;

awk 'BEGIN{FS="\t";OFS="\t"} {print $1,$2,$2+1,$4,$5,$6}' cov/control${time}.cov|sed '/chrX/d;/chrY/d;/chrM/d'|grep chr21 |sort -k1,1V -k2,2n -k3,3n > bed/control${time}.bed;
#######################################################################################

















