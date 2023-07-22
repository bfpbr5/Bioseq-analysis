# -*- encoding: utf-8 -
# env: python3.7

import json
import os
import sqlite3
import time
from pprint import pprint

import requests
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

# from PyEntrezId import Conversion


def read():
    # 读取包含单个序列 Fasta 格式文件
    fa_seq = SeqIO.read("res/sequence1.fasta", "fasta")
    # print fa_seq
    # 读取包含多个序列的 fasta 格式文件
    fa_list = []
    res_multi_data = SeqIO.parse("res/multi.fasta", "fasta")
    for fa in res_multi_data:
        fa_list.append(fa.seq)
    # 一个多序列文件中的所有序列
    seqs = [fa.seq for fa in res_multi_data]
    print (seqs)
    # 如果不想要seq对象中的字母表，可以用str()来强制类型转换
    seqs = [str(fa.seq) for fa in res_multi_data]
    print (seqs)

    # 读取包含单个序列的 gb 格式文件
    gb_seq = SeqIO.read("res/sequence1.gb", "genbank")
    print (gb_seq)

def show_data():
    # 读取包含单个序列 Fasta 格式文件
    fa_seq = SeqIO.read("res/sequence1.fasta", "fasta")

    # =====获取详细的信息=====
    # 提取基因ID，name
    # Fasta 文件中序列名所在行的第一个词被作为 id 和 name
    print ("id: ", fa_seq.id)
    print ("name: ", fa_seq.name)
    # 基因 Description 是fasta文件格式中的第一行
    print ("description: ", fa_seq.description)
    # 序列
    print ("seq: ", fa_seq.seq)
    # 序列来源库信息（NCBI的数据库信息会包括数据库交叉引用）
    print ("dbxrefs: ", fa_seq.dbxrefs)
    # 全部序列的注释信息
    print ("annotations: ", fa_seq.annotations)
    # 序列中每个字母的注释信息
    print ("letter_annotations: ", fa_seq.letter_annotations)
    # 部分序列的注释信息
    print ("features: ", fa_seq.features)

def show_genebank():
    '''浏览 genebank 序列文件内容'''
    # 读取包含单个序列的 gb 格式文件
    gb_seq = SeqIO.read("res/sequence1.gb", "genbank")
    print (gb_seq)

    # =====获取详细的信息=====
    # 提取基因ID，name
    # gb文件中序列名包含比fasta更加详细的序列信息，下面分别是 id 和 name
    print ("id: ", gb_seq.id)
    print ("name: ", gb_seq.name)
    # 基因 Description 是fasta文件格式中的第一行
    print ("description: ",  gb_seq.description)
    # 序列信息, 这里的序列信息是以 bioPython 中的seq对象存储
    print ("seq: ", gb_seq.seq)
    # 序列来源库信息（NCBI的数据库信息会包括数据库交叉引用）
    print ("dbxrefs: ", gb_seq.dbxrefs)
    # 全部序列的注释信息
    print ("annotations: ", gb_seq.annotations)
    # 序列中每个字母的注释信息
    print ("letter_annotations: ", gb_seq.letter_annotations)
    # 部分序列的注释信息,SeqFeature 对象的形式保存了features table中的所有entries（如genes和CDS等）
    print ("features: ", gb_seq.features)
    # 该基因的物种信息
    print ("organism: ", gb_seq.annotations["organism"])
    # 关于序列的注释信息，相关数据库的交叉引用号
    print ("comment: ", gb_seq.annotations["comment"])
    # 序列来源的物种名
    print ("source: ", gb_seq.annotations["source"])
    # 该基因的分类学信息
    print ("taxonomy: ", gb_seq.annotations["taxonomy"])
    # 该基因的整理后的注释信息
    print ("structured_comment: ", gb_seq.annotations["structured_comment"])
    # 该基因序列相关的关键词
    print ("keywords: ", gb_seq.annotations["keywords"])
    # 该基因的相关文献编号，或递交序列的注册信息
    print ("references: ", gb_seq.annotations["references"])
    # 该基因的入库时，给的基因编号，以及在染色体上的位点信息
    print ("accessions: ", gb_seq.annotations["accessions"])
    # 该基因的分子类型，一般为 DNA
    print ("molecule_type: ", gb_seq.annotations["molecule_type"])
    # 该基因的数据文件划分方式
    print ("data_file_division: ", gb_seq.annotations["data_file_division"])
    # 基因发布时间
    print ("date: ", gb_seq.annotations["date"])
    # 该基因的更新版本
    print ("sequence_version: ", gb_seq.annotations["sequence_version"])
    # 该基因的拓扑结构
    print ("topology: ", gb_seq.annotations["topology"])


def create_seq_file():
    # 新建一个DNA序列对象
    dna_seq = Seq("GGATGGTTGTCTATTAACTTGTTCAAAAAAGTATCAGGAGTTGTCAAGGCAGAGAAGAGAGTGTTTGCA", IUPAC.unambiguous_dna)
    pprint(dna_seq)
    # 新建一个RNA序列对象
    rna_seq = Seq("GGATGGTTGTCTATTAACTTGTTCAAAAAAGTATCAGGAGTTGTCAAGGCAGAGAAGAGAGTGTTTGCA", IUPAC.unambiguous_rna)
    pprint(rna_seq)
    # # 新建一个蛋白质序列对象
    protein_seq = Seq("GGATGGTTGTCTATTAACTTGTTCAAAAAAGTATCAGGAGTTGTCAAGGCAGAGAAGAGAGTGTTTGCA", IUPAC.protein)
    pprint(protein_seq)


# def get_gene_info():
#     EnsemblId = 'ENSG00000141510'
#     Id = Conversion('A.N.Other@example.com')
#     EntrezId = Id.convert_ensembl_to_entrez(EnsemblId)
#     print("EntrezId", EntrezId)
#     handle = Entrez.esummary(db="gene", id=EntrezId)
#     record = Entrez.read(handle)
#     DocSSet = record['DocumentSummarySet']['DocumentSummary'][0]
#     # with open('aaa.json', 'w', encoding='utf8') as fp:
#     #     json.dump(DocSSet, fp)
#     result = [
#         DocSSet['ChrStart'],
#         DocSSet['GenomicInfo'][0]
#     ]
#     pprint(result)


def batch_down(db_name= 'nucleotide', term_type= "Opuntia[orgn] and rpl16", ret_type= 'fasta', ret_mode= 'text', gene_id = ''):
    Entrez.email = "biff19@mails.tsinghua.edu.cn"
    batch_size = 3  # 批次
    if db_name == 'nucleotide':
        search_handle = Entrez.esearch(db=db_name,term=term_type, usehistory="y")
        search_results = Entrez.read(search_handle)
        search_handle.close()
        # gi_list = search_results["IdList"]
        count = int(search_results["Count"])
    elif db_name == 'gene':
        search_handle = Entrez.esearch(db=db_name,term=term_type, usehistory="y")
        search_results = Entrez.read(search_handle)
        search_handle.close()
        gene_id = search_results["IdList"][0]
        webenv=search_results["WebEnv"]
        query_key=search_results["QueryKey"]
        search_handle = Entrez.esummary(db=db_name,id=gene_id, webenv=webenv, query_key=query_key)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        result = search_results['DocumentSummarySet']['DocumentSummary'][0]
        print(result)
        chr = result["GenomicInfo"][0]["ChrLoc"]
        loci_t = result["GenomicInfo"][0]["ChrStart"]
        loci_p = result["GenomicInfo"][0]["ChrStop"]
    return chr,loci_t,loci_p

    # out_handle = open(out_file, "a+")
    # for start in range(0, count, batch_size):
    #     end = min(count, start+batch_size)
    #     print("Going to download record %i to %i" % (start+1, end))
        # fetch_handle = Entrez.efetch(db=db_name, rettype=ret_type, retmode=ret_mode,
        #                             retstart=start, retmax=batch_size,
        #                             webenv=search_results["WebEnv"], query_key=search_results["QueryKey"])
    #     data = fetch_handle.read()
    #     fetch_handle.close()
    #     out_handle.write(data)
    # out_handle.close()


def main():
    print(batch_down(db_name="gene", term_type="WASH7P"))


if __name__ == '__main__':
    main()

