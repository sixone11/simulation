#!/usr/bin/env bash

dir=${1:-input}
cd ${dir}
###########################################################################
# calculate the number of CpG sites in 
for i in metilene DMRfinder HMM_DM HMM_Fisher DSS methylSig methylKit BiSeq
do
    bedtools intersect -a ${i}_DMR -b merge_total.bed -c > ${i}_DMR_t
    mv ${i}_DMR_t ${i}_DMR
done 

################################################################################
# Step 1: fileter DMRs
# Step 1.1: Calculate the mean depth of each CpG site
python3 -c '
import pandas as pd;
df = pd.read_csv("merge_total.bed", sep="\t",header=None);
df_t=df[[3,4,5,6,7,8,9,10]]
df=pd.concat([df[[0,1,2]], df_t.mean(1)],axis=1);
df.to_csv("merge_total_mean.bed",index=False,sep="\t",header=None);' 
# Step 1.2: Filter DMRs by CpG covered in the region
bedtools intersect -a total_DMRs -b merge_total.bed -c|awk 'BEGIN{FS="\t";OFS="\t";} {if($4>=10) print $1,$2,$3}' >DMRs_1
# Step 1.3 : Filter DMRs by mean depth in the region
bedtools intersect -a DMRs_1 -b merge_total_mean.bed -loj > DMRs_2
python3 -c '
import pandas as pd
df = pd.read_csv("DMRs_2", sep="\t",header=None)
out=pd.DataFrame(df.groupby([0,1,2])[6].mean()).reset_index()
out=out[out[6]>=5]
out=out[[0,1,2]]
out.to_csv("DMRs.bed",sep="\t",header=None,index=False);' 

rm DMRs_1 DMRs_2

bedtools subtract -a chr21.bed -b DMRs.bed>non_DMRs.bed

# Add signal to regions: 1 represents the DMR and 0 non-DMRs
awk 'BEGIN{FS="\t";OFS="\t";} {print $0,1}' DMRs.bed>DMRs.bed_t
awk 'BEGIN{FS="\t";OFS="\t";} {print $0,0}' non_DMRs.bed>non_DMRs.bed_t
cat DMRs.bed_t non_DMRs.bed_t> DMRs_unDMRs_signal.bed_t
sortBed -i DMRs_unDMRs_signal.bed_t >DMRs_unDMRs_signal.bed

# Calculate the number of CpG sites in each region
bedtools intersect -a DMRs_unDMRs_signal.bed -b merge_total_mean.bed  -c >DMRs_unDMRs_signal_CpGcount.bed
###########################################################################


# Step 2: calculate the overlap rate
function calculate_rate(){
    name=$1
    
    bedtools intersect -a DMRs_unDMRs_signal_CpGcount.bed -b ${name}_DMR -loj > ${name}_rate_t

    python3 ../simulationdata_overlap.py -i ${name}_rate_t -o ../output/${name}_rate.bed
}


for i in metilene DMRfinder HMM_DM HMM_Fisher DSS methylSig methylKit BiSeq
do 
    calculate_rate ${i}
done 

# Delete the temporary variables
rm ../input/*_t
rm ../input/*_t_*