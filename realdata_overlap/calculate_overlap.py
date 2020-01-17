#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import pandas as pd
import sys
import os 


class calculate_overlap(object):
    '''
    Calculate the overlap rate based on the position of the first three columns.

    '''

    def run(self,args):
        self.name = os.path.basename(args[0])
        parser = self.create_parser()
        self.opts = parser.parse_args(args[1:])
        self.main_func(self.opts)
        

    def main_func(self,opts):
        df=pd.read_csv(opts.i,sep="\t",header=None)
        df[7]=df[7].replace(".",0)
        df[7]=df[7].astype(int)
        #First: select the areas that do not have overlap with the union of tools
        no_overlap=df[df[4]=="."].reset_index(drop=True)
        overlap_regions=df[df[4]!="."].reset_index(drop=True)
        out=[]
        #Second: Redefine the starting position of the region that have overlap with the union of tools.
        #For example: the union of region is chr21 30 50 and the region in DMRfinder is chr21 40 60. 
        # Then we should convert the position to chr21 max(30,40) min(50,60) i.e. chr21 40 50. 
        for i in range(0,overlap_regions.shape[0]):
            out.append(['chr21',max(overlap_regions[1][i],overlap_regions[5][i]),min(overlap_regions[2][i],overlap_regions[6][i])])
        outname=str(opts.i)+'_t' 
        with open(outname, 'w') as f:
            for items in out:
                f.write('\t'.join(str(i) for i in items) + '\n')
        # Third: recalculate the number of CpG site in regions
        command="bedtools intersect -a "+ str(outname)+ " -b merge_total.bed -c |awk  '{print $4}' >"+str(outname)+"_CGIcount"
        os.system(command)
        CGIcount_name=str(outname)+"_CGIcount"
        CGIcount=pd.read_csv(CGIcount_name, sep="\t",header=None)
        overlap_regions[7]=CGIcount
        # Fourth: Calculate the total number of CpG sites that overlapped with a region according to the first three columns.
        out=pd.DataFrame(overlap_regions.groupby([0,1,2,3])[7].sum()).reset_index()
        out=pd.concat([out[[0,1,2]],out[7]/out[3]],axis=1)
        out.columns=[0,1,2,7]
        out=pd.concat([out,no_overlap[[0,1,2,7]]],axis=0).reset_index(drop=True)
        out.sort_values(by=1,ascending=True,inplace=True)
        out.reset_index(drop=True,inplace=True)
        #out[out[7]>=0.5].shape[0]/out.shape[0]
        out=out.fillna(0)
        out.to_csv(opts.o,sep="\t",index=False,na_rep='NaN',header=None)

    def create_parser(self):
        parser = argparse.ArgumentParser(description = "Calculate the rate for heatmap plot")
        parser.add_argument('-i', metavar = 'in-file', type = str, help = "input file that has been sort and merge before",default="")
        parser.add_argument('-o', metavar = 'out-file', type = str, help = "outfile",default="")
        return parser



if __name__=='__main__':
    calculate_overlap().run(sys.argv)


