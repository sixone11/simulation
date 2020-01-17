
#Necessary file 1: real_data_mean_chr21.bed : The file is the mean value of each CpG site based on real data
#Necessary file 2: DMRs_unDMRs_signal.bed : The file is the position that includes non-DMRs region and DMRs region, and the fourth column represents the region type 0: non-DMRs 2:hyper DMRs 3:hypo-DMRs 
#Necessary file 3: hg19.fa :the reference sequence of hg19


#step 1:get the mean methylation level of regions
#1.1
bedtools intersect -a DMRs_unDMRs_signal.bed -b real_data_mean_chr21.bed -loj >position_CpG_mean.bed

#1.2:calculate the mean
python -c'
import pandas as pd
import numpy as np
def main():
    df = pd.read_csv("position_CpG_mean.bed", sep="\t",header=None);
    df[df[7]=="."]=np.nan
    df[7]=df[7].astype('float64')
    out=pd.DataFrame(df.groupby([0,1,2,3])[7].mean()).reset_index()
    out[7]=round(out[7],1)
    out[1]=out[1].astype('int64')
    out[2]=out[2].astype('int64')
    out[3]=out[3].astype('int64')
    out.to_csv("position_mean.bed",index=False,sep="\t",header=None);
if __name__=="__main__":
    main()
'


#stpe 2:get the bed file. 
# In a bed file the mean methylation level of each region is the same, and we keep one decimal place of mean methylation level. 
# We can get the bed file of control group which the mean methylation level wasn't changed anymore.
#At the same time, we changed the methylation level of DMRs according to the DMRs type(hyper or hypo), and the variation is random and larger than 0.1(It is OK if you want to set the different variation to a certain ratio). 

python -c'
import pandas as pd
import numpy as np
import random
import os 

def main():
    df = pd.read_csv("position_mean.bed", sep="\t",header=None)
    out=list()
    for i in [round(i,1) for i in np.arange(0,1.1,0.1)]:
        out.append(df[df[4]==i].reset_index(drop=True))
    for i in [k for k in range(0,len(out))]:
        data=out[i]
        if i==0:
            for j in [k for k in range(0,data.shape[0])]:
                if data[3][j]==3:
                    data[3][j]=2
        elif i==10:
            for j in [k for k in range(0,data.shape[0])]:
                if data[3][j]==2: data[3][j]=3
        out[i]=data
    for i in [k for k in range(0,len(out))]:
        data=out[i]
        for j in [k for k in range(0,data.shape[0])]:
            if data[3][j]==2:
                limit=(1-round(data[4][0],1))*10
                if limit<1: limit=1
                data[4][j]=round(random.randint(1,limit)*0.1,1)
            if data[3][j]==3:
                data[4][j]=round(random.randint(1,data[4][0]*10)*0.1,1) 
            out[i]=data
    output=pd.concat(out,axis=0)
    output.dtypes[4]="float"
    os.system("mkdir -p test.bed")
    for i in [i*0.1 for i in range(0,11)]:
        df_out=output[output[4]==round(i,1)].sort_index(by = [0,1])
        name="test.bed/methy_level_"+str(round(i,1))+".bed"
        df_out.to_csv(name,index=False,sep="\t",header=None)
    os.system("mkdir -p control.bed")
    for i in [i*0.1 for i in range(0,11)]:
        df_out=df[df[4]==round(i,1)].sort_index(by = [0,1])  
        name="control.bed/methy_level_"+str(round(i,1))+".bed"
        df_out.to_csv(name,index=False,sep="\t",header=None)
if __name__=="__main__":
    main()
'

#step 3:get the reference sequence according to the bed file above.
mkdir -p test.fa
mkdir -p control.fa
for i in {0..9}
do 
	bedtools getfasta -fi hg19.fa -bed test.bed/methy_level_0.${i}.bed -fo test.fa/methy_level_0.${i}.fa&
done
bedtools getfasta -fi hg19.fa -bed test.bed/methy_level_1.0.bed -fo test.fa/methy_level_1.0.fa&

for i in {0..9}
do 
	bedtools getfasta -fi hg19.fa -bed control.bed/methy_level_0.${i}.bed -fo control.fa/methy_level_0.${i}.fa&
done

bedtools getfasta -fi hg19.fa -bed test.bed/methy_level_1.0.bed -fo control.fa/methy_level_1.0.fa&


#step 4:get the DMRs' position
awk 'BEGIN{FS="\t";OFS="\t";} {if($4!=0) print $1,$2,$3}' DMRs_unDMRs_signal.bed>total_DMRs
