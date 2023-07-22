import pandas as pd
import os
import gzip

import gentools

import os

import json
import pandas as pd
from sklearn.model_selection import train_test_split



# read files from folder, combine them into a df
def name_roll(filePath):
    for i,j,k in os.walk(filePath):
        return k

def rewrite(filename):
    with open(filename, mode="r", encoding="utf8") as f:
        data = f.readlines()
        data = data[2:]
        n = open(filename,'w+')
        n.writelines(data)
        n.close()

def read_mat(filename):
    data = pd.read_table(filename, sep="\t", header=0)
    return data

def main():
    # data_list = name_roll("/Users/bfpbr5/Documents/GitHub/DESSO-master-whole/data/TfbsUniform_hg19_ENCODE")
    negfile = "/Users/bfpbr5/Documents/GitHub/DESSO-master-whole/data/encode_101_background/posi_back.txt"
    posfiles = "/Users/bfpbr5/Documents/GitHub/DESSO-master-whole/data/encode_101"
    outputs = []
    i = 0
    for gz in name_roll(posfiles):
        try:
            g = gzip.open(posfiles+'/'+gz,"r")
            glines = g.readlines()
        except:
            print(g)
            continue
        data = [x.decode('utf-8').split('\t') for x in glines]
        for peak in data[1:]:
            output = {}
            seq = peak[2]
            output['seq'] = seq
            output['tag'] = 1
            outputs.append(output)
                # pos_seqs.append(seq)
        i += 1
        if(i%100 == 0):
            print("Processed Gzip File: %d"%(i))
    genome = gentools.readGenome("/Users/bfpbr5/Documents/GitHub/DESSO-master-whole/data/GRCh37.p13.genome.fa.gz")
    with open(negfile,"r") as f:
        data = [x.split('\t') for x in f.readlines()]
        data = data[1:]
        i = 0
        for item in data:
            output = {}
            loci = (item[0],item[1],item[2][:-1])
            seq = gentools.loc2seq(genome,loci=loci)
            output['seq'] = seq
            output['tag'] = 0
            outputs.append(output)
            i += 1
            if(i%10000 == 0):
                print("Processed Back Line: %d"%(i))
    
    output_name = "DESSO.json"
    fp = open(output_name,"w")
    json.dump(outputs,fp)
    fp.close()

if __name__ == "__main__":
    main()
