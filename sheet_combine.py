#!/bin/env python3
#Merge two tables according to common columns
import numpy as np
import pandas as pd
import os,sys,glob

if len(sys.argv) < 4 :
    sys.exit('Usage:\n     python3 sys.argv[0] inpath1 inpath2 outpath')

def merge(inpath1,inpath2,outpath):
    files=glob.glob(inpath1+'/'+'*annotation_simple.txt')
    for anno in files:
        frename=os.path.basename(anno).split('_annotation')[0]
        go=frename+'_simple.GO-Analysis_BP_All.xlsx'
        outname=frename+'_mergeAnnotation.xlsx'
        frame1=pd.read_csv(anno,sep='\t')
        frame2=pd.read_excel(inpath2+'/'+'GOAnalysis_result'+'/'+go,sheet_name=2,usecols=[0,2,3,4])
        frame3=pd.merge(frame1,frame2,left_on='Nearest PromoterID',right_on='QueryID',how='outer').drop('QueryID',axis=1)
        frame3.to_excel(outpath+'/'+outname,header=True,index=False,na_rep='NA')

if __name__ == "__main__":
    merge(sys.argv[1],sys.argv[2],sys.argv[3])
