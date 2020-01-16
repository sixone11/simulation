#!/bin/bash
. aa 
dir_t=${1:-/data/sixone/lllab/DMR/BiB_final/pROC}
time=${2:-1}
dir=${dir_t}/${time}/bed
cd ${dir}
mkdir MDS
cp merge_total.bed  MDS

awk 'BEGIN{FS="\t";OFS="\t";} {if(NR>1) print $1,$2,$3}' metilene_DMR>MDS/metilene_DMR
awk 'BEGIN{FS="\t";OFS="\t";} {if(NR>1 && $4>3) print $1,$2,$3}' DMRfinder_DMR>MDS/DMRfinder_DMR


cd MDS
cat metilene_DMR DMRfinder_DMR >metilene_DMRfinder_merge_DMR
sortBed -i metilene_DMRfinder_merge_DMR>metilene_DMRfinder_merge_DMR_t
bedtools merge -i metilene_DMRfinder_merge_DMR_t >metilene_DMRfinder_merge_DMR
bedtools intersect -a DMRfinder_DMR -b merge_total.bed -c>DMRfinder_DMR_CpGcount
cat metilene_DMR DMRfinder_DMR >metilene_DMRfinder_merge_DMR

bedtools intersect -a metilene_DMR -b merge_total.bed -c>metilene_DMR_CpGcount

sortBed -i metilene_DMRfinder_merge_DMR>metilene_DMRfinder_merge_DMR_t
bedtools merge -i metilene_DMRfinder_merge_DMR_t >metilene_DMRfinder_merge_DMR
bedtools intersect -a metilene_DMRfinder_merge_DMR -b merge_total.bed -c>metilene_DMRfinder_merge_CpGcount
#head metilene_DMRfinder_merge_CpGcount
#head metilene_DMR_CpGcount 
#head DMRfinder_DMR_CpGcount 

bedtools intersect -a metilene_DMRfinder_merge_CpGcount -b metilene_DMR_CpGcount -loj >mDmerge_metilene_CpGcount

python3 /data/sixone/lllab/DMR/BiB_final/real_overlap/script/heatmap_overlap.py -i mDmerge_metilene_CpGcount -o mDmerge_metilene_CpGrate

bedtools intersect -a metilene_DMRfinder_merge_CpGcount -b DMRfinder_DMR_CpGcount -loj >mDmerge_DMRfinder_CpGcount

python3 /data/sixone/lllab/DMR/BiB_final/real_overlap/script/heatmap_overlap.py -i mDmerge_DMRfinder_CpGcount -o mDmerge_DMRfinder_CpGrate

rm *_t
rm *_t_*

####################################
#暂时没法评价 后续拼接改变阈值那个脚本