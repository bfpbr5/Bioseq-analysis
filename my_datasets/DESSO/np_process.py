# import pandas as pd
import os
import gzip

# import gentools

# import os

# import json
# import pandas as pd
# from sklearn.model_selection import train_test_split



# read files from folder, combine them into a df
def name_roll(filePath):
    for i,j,k in os.walk(filePath):
        return k

def rewrite(file, do_open=True, filename=None):
    if do_open:
        f = open(file, mode="r", encoding="utf8")
        filename = file
    else:
        f = file
    data = f.readlines()
    data = data[1:]
    n = open(filename+'.bed','w+')
    for line in data:
        line = line.decode('utf-8')
        line = line.split('\t')
        line = line[0:3]
        line = '\t'.join(line)
        line = line+'\n'
        n.write(line)
    n.close()
    if do_open:
        f.close()

# def read_mat(filename):
#     data = pd.read_table(filename, sep="\t", header=0)
#     return data

def main():
    peakFolder = "DESSO-master-whole/data/TfbsUniform_hg19_ENCODE"
    for gz in name_roll(peakFolder):
        try:
            g = gzip.open(peakFolder+'/'+gz,"rb")
            rewrite(g,False,gz)
            g.close()
        except:
            print(gz)
            continue


if __name__ == '__main__':
    main()