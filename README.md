#  PART 1:  preparation before simulating

1. Input
   - hg19.fa.gz : The reference genome file
   - DMRs_unDMRs_signal.bed : The first three columns are location information and the fourth column represents the regions' classification(0: non-DMRs 2:hyper-DMRs 3:hypo-DMRs)
   - real_data_mean_chr21.bed: The mean CpG level of CpG sites from real data.

2. Output
   -  test.fa: A folder that contains the fasta file for treatment group.
   -  control.fa : A folder that contains the fasta file for control group.
   -  test.bed: A folder that contains the fasta file for treatment group.
   -  control.bed: A folder that contains the fasta file for control group.
   -  total DMRs: The first three columns are the location information and the fourth column represents the regions' classification(0: non-DMRs 1:DMRs).
3. Script
   - preparation.sh

# PART 2: The process of simulation

1. Input
   - control.fa : A folder that contains the fasta file for control group.
   - test.fa : A folder that contains the fasta file for treatment group.

2. Output
   -  rrbssim: A folder that contains all simulated files.
   -  fastqc_initial:  A folder that contains the html files before trimming.
   -  fq:  A folder that contains the simulated fastq files.
   -  fastqc_trim： A folder that contains the file after trimming.
   -  bam:  A folder that contains the bam files which were processed by bismark.
   -  cov:  A folder that contains the coverage files which were processed by bismark
   -  bed:   A folder that contains the bed files which generated from cov file. 
3. Script
   - simulation.sh : We could get the 1 test simulation data and 1 control simulation data through running this script
     - The basic format of run this script is like: sh simulation.sh -o ${outdir} -t ${time}
   - run_simulation.sh : We could get n pairs of simulated data through running this script and this script is a supplementary for simulation.sh
     - The basic format of run this script is like: sh run_simulation.sh -o ${outdir} -n ${number}

#  PART 3: The process of calculating the overlap in simulated data

1. Input

   - *_DMR: The results of different tools about DMRs.
   - chr21.bed: The bed file of the start and end of chromosome21.
   - total DMRs: The first three columns are the location information and the fourth column represents the regions' classification(0: non-DMRs 1:DMRs).

2. Output

   -  *_rate.bed: The overlap rate between the DMRs/non-DMRs and the tool.

3. Script

   - simulateddata_overlap.sh: The main script of calculating overlap rate in simulated data.

   - simulattedata_overlap.py: The sub script of calculating overlap rate in simulated data.



# PART 4: The process of calculating the overlap in real data

1. Input
   - *_DMR: The results of different tools about DMRs.
   - merge_total.bed: The bed file that contains the depth information of the intersection of CpG sites of all samples.

2. Output
   -  mDmerge\_*_CpGrate: The overlap rate of  2 different tools.
3. Script
   - calculate_realdata_overlaprate.sh: The main script of calculating the overlap in real data.
   - calculate_overlap.py: The sub script that used to calculate the AUC in real data.



# PART 5: The process of calculating  the AUC

1. Input
   - *_DMR: The results of different tools about DMRs.
   - merge_total.bed: The bed file that contains the depth information of the intersection of CpG sites of all samples.
   - total DMRs: The first three columns are the location information and the fourth column represents the regions' classification(0: non-DMRs 1:DMRs).
   - chr21.bed: The bed file of the start and end of chromosome21.

2. Output
   -  *_rate.bed: The overlap rate between the DMRs/non-DMRs and the tool.
   -  AUC_results: The file that contains the value of AUC.
3. Script
   - calculate_AUC.sh: The main script of calculating AUC
   - simulationdata_overlap.py: The sub script that used to calculated the overlapped rate in simulated data.
   - calculate_AUC_all.R The sub script that used to calculate the AUC.
