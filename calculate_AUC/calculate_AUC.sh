#!/usr/bin/env bash
. aa
dir_t=${1:-/data/sixone/lllab/DMR/BiB_final/pROC}
time=${2:-"1"}
dir=${dir_t}/${time}/bed

cd ${dir}
mkdir process

#awk 'BEGIN{FS="\t";OFS="\t";} {if(NR>1 && $7<=0.05 && $8<=0.05 && $6>10) print $1,$2,$3,$6}' metilene_DMR>process/metilene_DMR

cp HMM-DM/HMM_DM_DMR ./
cp HMM-Fisher/HMM_Fisher_DMR ./

awk 'BEGIN{FS="\t";OFS="\t";} {if(NR>1) print $1,$2,$3}' metilene_DMR>process/metilene_DMR

awk 'BEGIN{FS="\t";OFS="\t";} {if(NR>1) print $1,$2,$3}' DMRfinder_DMR>process/DMRfinder_DMR

awk 'BEGIN{FS="\t";OFS="\t";} {if(NR>1) print $1,$2,$3}' DSS_DMR>process/DSS_DMR

awk 'BEGIN{FS=" ";OFS="\t";} {print $1,$2,$3}' BiSeq_DMR>process/BiSeq_DMR

awk 'BEGIN{FS="\t";OFS="\t";} {if(NR>1) print $1,$2,$3}' HMM_Fisher_DMR>process/HMM_Fisher_DMR

awk 'BEGIN{FS="\t";OFS="\t";} {if(NR>1) print $1,$2,$3}' HMM_DM_DMR>process/HMM_DM_DMR

awk 'BEGIN{FS="\t";OFS="\t";} {if(NR>1) print $1,$2,$3}' methylSig_DMR>process/methylSig_DMR

awk 'BEGIN{FS="\t";OFS="\t";} {if(NR>1) print $1,$2,$3}' methylKit_DMR>process/methylKit_DMR

for i in metilene DMRfinder HMM_DM HMM_Fisher DSS methylSig methylKit BiSeq
do
    bedtools intersect -a process/${i}_DMR -b merge_total.bed -c > process/${i}_DMR_t
    mv process/${i}_DMR_t process/${i}_DMR
done 

####################################################################################
#intersection
#all_CpG_site=`cat merge_total.bed|wc -l`
#union
#all_CpG_site=`cat merge_union.bed|wc -l`

#要说明一下这个文件是怎么来的！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！/data/sixone/lllab/DMR/BiB_final/pROC/DMRs_unDMRs_signal.bed
#bedtools intersect -a /data/sixone/lllab/DMR/BiB_final/pROC/DMRs_unDMRs_signal.bed -b merge_total.bed -c >DMRs_unDMRs_signal_CpGcount_intersect.bed
#bedtools intersect -a /data/sixone/lllab/DMR/BiB_final/pROC/DMRs_unDMRs_signal.bed -b merge_union.bed -c >DMRs_unDMRs_signal_CpGcount_union.bed

#fileter DMRs
python3 -c '
import pandas as pd;
df = pd.read_csv("merge_total.bed", sep="\t",header=None);
df_t=df[[3,4,5,6,7,8,9,10]]
df=pd.concat([df[[0,1,2]], df_t.mean(1)],axis=1);
df.to_csv("merge_total_mean.bed",index=False,sep="\t",header=None);' 


bedtools intersect -a /data/sixone/lllab/DMR/BiB_final/code/data/total_DMRs -b merge_total.bed -c|awk 'BEGIN{FS="\t";OFS="\t";} {if($4>=10) print $1,$2,$3}' >DMRs_1

bedtools intersect -a DMRs_1 -b merge_total_mean.bed -loj > DMRs_2
python3 -c '
import pandas as pd
df = pd.read_csv("DMRs_2", sep="\t",header=None)
out=pd.DataFrame(df.groupby([0,1,2])[6].mean()).reset_index()
out=out[out[6]>=5]
out=out[[0,1,2]]
out.to_csv("DMRs.bed",sep="\t",header=None,index=False);' 


rm DMRs_1 DMRs_2

bedtools subtract -a /data/sixone/lllab/DMR/BiB_final/pROC/chr21.bed -b DMRs.bed>non_DMRs.bed

awk 'BEGIN{FS="\t";OFS="\t";} {print $0,1}' DMRs.bed>DMRs.bed_t
awk 'BEGIN{FS="\t";OFS="\t";} {print $0,0}' non_DMRs.bed>non_DMRs.bed_t
cat DMRs.bed_t non_DMRs.bed_t> DMRs_unDMRs_signal.bed_t
sortBed -i DMRs_unDMRs_signal.bed_t >DMRs_unDMRs_signal.bed
rm *_t 

bedtools intersect -a DMRs_unDMRs_signal.bed -b merge_total_mean.bed  -c >DMRs_unDMRs_signal_CpGcount.bed
###########################################################################


###########################################################################
cd process
cp ../merge_total.bed ./
cp ../merge_union.bed ./

#DMRfinder
function calculate_rate(){
    name=$1
    #bedtools intersect -a ${name}_DMR -b ../merge_total_mean.bed  -c >${name}_DMR_CpGcount.bed
    cp ${name}_DMR ${name}_DMR_CpGcount.bed

    bedtools intersect -a ../DMRs_unDMRs_signal_CpGcount.bed -b ${name}_DMR_CpGcount.bed -loj > ${name}_rate_t

    python3 /data/sixone/lllab/DMR/BiB_final/multiple_simulation/script/simulationdata_overlap.py -i ${name}_rate_t -o ${name}_rate.bed
}


for i in metilene DMRfinder HMM_DM HMM_Fisher DSS methylSig methylKit BiSeq
do 
    calculate_rate ${i}
done 


#debug
echo ${dir_t}/${time} >>/data/sixone/lllab/DMR/BiB_final/multiple_simulation/results/${time}.out
#calculate AUC

Rscript /data/sixone/lllab/DMR/BiB_final/pROC/script/calculate_AUC_all.R ${dir_t}/results/${time}.out

Rscript /data/sixone/lllab/DMR/BiB_final/pROC/script/calculate_PR.R ${dir_t}/results/${time}.out

rm *_t
