import os
from xml.dom.minidom import parseString

import Bio
import pandas as pd
from Bio import Entrez

import search_NCBI


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

# expr data to gene-expr pair
def mat_extract(mat):
    pass

def geneLoci(gene):
    parseString

def main():
    ''' foldername = "/Users/bfpbr5/Documents/GitHub/genomic_benchmarks/data/eQTL/tissue"
    # namelist = name_roll(foldername)
    # vecs = []
    # name = namelist[0]
    # rewrite(name)
    # tissuename = name[23:-4]
    # data = read_mat(name)
    # for name in namelist:
    #     rewrite(name)
    #     tissuename = name[23:-4]
    #     data = read_mat(name)
    #     vec = mat_extract(data) # [0.2,0.4,...] 56200 genes
    #     vecs.append(vec)
    # print(vecs)
    # id_list = []
    # Entrez.email = "biff19@mails.tsinghua.edu.cn"     # Always tell NCBI who you are
    # search_results = Entrez.read(Entrez.epost("pubmed", id=",".join(id_list)))
    # webenv = search_results["WebEnv"]
    # query_key = search_results["QueryKey"]
    # batch_size = 3
    # out_handle = []
    # count = 56200
    # for start in range(0,count,batch_size):
    #     end = min(count, start+batch_size)
    #     print "Going to download record %i to %i" % (start+1, end)
    #     fetch_handle = Entrez.efetch(db="gene", rettype="fasta", retmode="text",
    #                                 retstart=start, retmax=batch_size,
    #                                 webenv=webenv, query_key=query_key)
    #     data = fetch_handle.read()
    #     fetch_handle.close()
    #     out_handle.append(data)
    '''
    data = read_mat("/Users/bfpbr5/Documents/GitHub/genomic_benchmarks/data/eQTL/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct")
    errors1, errors2 = search_NCBI.getSeq(list(data['Description']), 12400)
    print(errors1)
    print(errors2)
    with open("result/erros1","w") as f:
        f.writelines(errors1)
    with open("result/erros2","w") as f:
        f.writelines(errors2)

if __name__ == "__main__":
    main()
