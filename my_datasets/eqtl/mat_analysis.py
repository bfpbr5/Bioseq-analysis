import os
import random
from pprint import pprint

import pandas as pd


def merge_mat(cpoint=0, end=56200, batch_size=50):
    df_res = pd.read_table("result/result_"+str(cpoint), sep='\t', header=None, index_col=None, skip_blank_lines=False)
    for i in range(cpoint+batch_size,end, batch_size):
        # with open("result/result_"+str(i), 'r') as f:
        #     if(f.readline() == ''):

        df = pd.read_table("result/result_"+str(i), sep='\t', header=None, index_col=None, skip_blank_lines=False)
        df_res = pd.concat([df_res,df],axis=0)
    df_res.fillna('',inplace=True)
    # print(df_res.iloc[51:80,:])
    return df_res

def chr_label():

    pass

def drop_data(data):
    # if data == 0:
    #     return ''
    data.strip('N')
    a = random.random()
    if a <= 0.25:
        data.replace('N','A')
    elif a <= 0.5:
        data.replace('N','T')
    elif a <- 0.75:
        data.replace('N','C')
    else:
        data.replace('N','G')
    return data


def main():
    cpoint = 3350
    end = 56200
    df_seq = merge_mat(cpoint, end)
    df_seq[1].apply(drop_data)
    df_exp = pd.read_table("/Users/bfpbr5/Documents/GitHub/genomic_benchmarks/data/eQTL/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", sep="\t", header=0)
    df_exp.insert(1,'seq',0)
    df_exp.insert(2,'chr',0)
    # print(df_seq.iloc[:,1])
    df_exp.iloc[cpoint:end,1] = pd.Series(df_seq.iloc[:,1])
    df_exp.iloc[cpoint:end,2] = pd.Series(df_seq.iloc[:,0])
    df_exp.drop(df_exp[(df_exp['seq']=='')].index, inplace=True)
    df_exp.drop(df_exp[(df_exp['seq']=='\n')].index, inplace=True)
    df_exp.drop(df_exp[(df_exp['seq']==0)].index, inplace=True)
    df_exp.to_json("eqtl_data.json")

if __name__ == '__main__':
    main()
