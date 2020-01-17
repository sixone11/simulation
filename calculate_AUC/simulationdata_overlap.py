#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import pandas as pd
import sys
import os 


class simulationdata_overlap(object):
    '''
   
    '''

    def run(self,args):
        self.name = os.path.basename(args[0])
        parser = self.create_parser()
        self.opts = parser.parse_args(args[1:])
        self.main_func(self.opts)
        

    def main_func(self,opts):
        df=pd.read_csv(opts.i,sep="\t",header=None)
        df[8]=df[8].replace(".",0)
        df[8]=df[8].astype(int)
        no_overlap=df[df[5]=="."].reset_index(drop=True)
        overlap_regions=df[df[5]!="."].reset_index(drop=True)

        out=[]
        for i in range(0,overlap_regions.shape[0]):
            out.append(['chr21',max(overlap_regions[1][i],overlap_regions[6][i]),min(overlap_regions[2][i],overlap_regions[7][i])])
        outname=str(opts.i)+'_t' 
        with open(outname, 'w') as f:
            for items in out:
                f.write('\t'.join(str(i) for i in items) + '\n')
        #command="bedtools intersect -a "+ str(outname)+ " -b merge_union.bed -c |awk  '{print $4}' >"+str(outname)+"_CGIcount"
        command="bedtools intersect -a "+ str(outname)+ " -b merge_total.bed -c |awk  '{print $4}' >"+str(outname)+"_CGIcount"
        os.system(command)
        CGIcount_name=str(outname)+"_CGIcount"
        CGIcount=pd.read_csv(CGIcount_name, sep="\t",header=None)
        overlap_regions[8]=CGIcount
        out=pd.DataFrame(overlap_regions.groupby([0,1,2,3,4])[8].sum()).reset_index()
        out=pd.concat([out[[0,1,2,3]],out[8]/out[4]],axis=1)
        out.columns=[0,1,2,3,8]
        out=pd.concat([out,no_overlap[[0,1,2,3,8]]],axis=0).reset_index(drop=True)
        out.sort_values(by=1,ascending=True,inplace=True)
        out.reset_index(drop=True,inplace=True)
        #out[out[7]>=0.5].shape[0]/out.shape[0]
        out[8]=round(out[8],2)
        out.to_csv(opts.o,sep="\t",index=False,na_rep='NaN',header=None)


    def create_parser(self):
        parser = argparse.ArgumentParser(description = "Calculate the rate for heatmap plot")
        parser.add_argument('-i', metavar = 'in-file', type = str, help = "input file that has been sort and merge before",default="")
        parser.add_argument('-o', metavar = 'out-file', type = str, help = "outfile",default="")
        return parser



if __name__=='__main__':
    simulationdata_overlap().run(sys.argv)


