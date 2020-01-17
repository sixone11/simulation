#!/bin/bash
#Initialization of variable
dir=${1:-input}
#######################################################################
cd ${dir}


cat metilene_DMR DMRfinder_DMR >metilene_DMRfinder_merge_DMR
sortBed -i metilene_DMRfinder_merge_DMR>metilene_DMRfinder_merge_DMR_t
# Get the union of 2 tools' DMRs
bedtools merge -i metilene_DMRfinder_merge_DMR_t >metilene_DMRfinder_merge_DMR
#sort
sortBed -i metilene_DMRfinder_merge_DMR > metilene_DMRfinder_merge_DMR_t
mv metilene_DMRfinder_merge_DMR_t metilene_DMRfinder_merge_DMR

# Calculate the number of CpG sites that covered in each DMRfinder's DMR
bedtools intersect -a DMRfinder_DMR -b merge_total.bed -c>DMRfinder_DMR_CpGcount
# Calculate the number of CpG sites that covered in each metilene's DMR
bedtools intersect -a metilene_DMR -b merge_total.bed -c>metilene_DMR_CpGcount
# Calculate the number of CpG sites that covered in the union of 2 tools' DMRs
bedtools intersect -a metilene_DMRfinder_merge_DMR -b merge_total.bed -c>metilene_DMRfinder_merge_CpGcount

# Calculate the overlap between the union of 2 tools and metilene
bedtools intersect -a metilene_DMRfinder_merge_CpGcount -b metilene_DMR_CpGcount -loj >mDmerge_metilene_CpGcount

#Calculate the CpG rate that metilene covered in every union of region.
python3 ../calculate_overlap.py -i mDmerge_metilene_CpGcount -o ../output/mDmerge_metilene_CpGrate

# Calculate the overlap between the union of 2 tools and DMRfinder
bedtools intersect -a metilene_DMRfinder_merge_CpGcount -b DMRfinder_DMR_CpGcount -loj >mDmerge_DMRfinder_CpGcount

#Calculate the CpG rate that DMRfinder covered in every union of region.
python3 ../calculate_overlap.py -i mDmerge_DMRfinder_CpGcount -o ../output/mDmerge_DMRfinder_CpGrate

#delete the temporary variables
rm *_t
rm *_t_*