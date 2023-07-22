import json
import random
import tarfile
from pathlib import Path

import pandas as pd
import yaml
from genomic_benchmarks.seq2loc import fasta2loc
from genomic_benchmarks.utils.datasets import _fastagz2dict
from sklearn.model_selection import train_test_split
from tqdm.notebook import tqdm


def readGenome(filename='/Users/bfpbr5/Documents/GitHub/genomic_benchmarks/data/hg38.analysisSet.fa.gz'):
    genome = _fastagz2dict(filename)
    return genome
# print(genome.keys())

def romanToInt(s):
        """
        :type s: str
        :rtype: int
        """
        result = 0
        temp_dict = {"I":1,"V":5,"X":10,"L":50,"C":100,"D":500,"M":1000}
        for i in range(len(s)-1):
            if temp_dict[s[i]]<temp_dict[s[i+1]]:
                result -= temp_dict[s[i]]
            else:
                result += temp_dict[s[i]]
        result += temp_dict[s[-1]]
        return result

def norm_chr(chr):
    if chr[0] != 'c': # 是数字
        result = "chr"+chr
        return result
    if chr[2] == 'r': # 是chrN
        if (chr[3].isdecimal()) or (chr[3] in ['X', 'Y', 'M']): # 是Chr1-n或ChrX、Y
            return chr
        result = 'chr'+str(romanToInt(chr[3:])) # chr+罗马数字
        return result
    elif chr[1] == 'h': # chN
        result = 'chr'+str(int(chr[2:]))
        return result
    else:
        return 'chr0'

def longer_loc(chrstart, chrend, start = 0, end = 0, var=500, loci=(0,0,0)):
    # random.seed(seed)
    # rd = random.random() * var + 100
    rd = var
    start = int(start)
    end = int(end)
    start = int(start - rd)
    end = int(end + rd)
    if start < int(chrstart):
        start = chrstart
    if end > int(chrend):
        end = chrend
    return int(start), int(end)

def is_in_seq(seq, tfs):
    if seq[0] != tfs[0]:
        return False
    if (int(seq[2]) - int(tfs[1])) * (int(seq[1])-int(tfs[2])) > 0:
        return False
    return True

def rand_loc(genome):
    chr = random.randint(1,24)
    if chr == 23:
        chr = 'X'
    if chr == 24:
        chr = 'Y'
    rand_chr = "chr"+str(chr)
    rd = random.randint(1,len(genome[rand_chr])-502)
    rd_end = rd+random.randint(100,500)
    return rand_chr, rd, rd_end

def rev(seq, strand='+'):
    if strand == '-':
        return str(Seq(seq).reverse_complement())
    else:
        return seq

def snv2seq(chr, start, end, strand, allele, allele_loc_start, genome):
    seq = rev(genome[chr][int(start):int(allele_loc_start)],strand)
    if (allele in ['A','T','C','G']) or len(allele) > 1:
        seq = seq + allele + rev(genome[chr][(int(allele_loc_start)+len(allele)):int(end)],strand)
    elif allele == '-':
        seq = seq + rev(genome[chr][(int(allele_loc_start)+1):int(end)],strand)
    else:
        seq = 'ERROR'
    return seq

from Bio.Seq import Seq


def rev(seq, strand='+'):
    if strand == '-':
        return str(Seq(seq).reverse_complement())
    else:
        return seq

def readData(filename):
    with open(filename, "r", encoding="utf8") as f:
        lines = f.readlines()
    return lines

def loc2seq(genome, strand='+', chr=0, start=0, end=0, loci=(0,0,0)): #loci:chr, start, end
    if(chr == 0):
        chr = norm_chr(loci[0])
        start = int(loci[1])
        end = int(loci[2])
    if start > end:
        strand = '-'
    seq = rev(genome[chr][int(start):int(end)],strand)
    seq = seq.upper()
    return seq


def seq2json(filename, lines, genome, is_lib=False):
    data_list = []
    for i in range(len(lines)):
        if not lines[i]:
            continue
        if is_lib:
            data = lines[i].split(',')
            chrom = data[1]
            if i == 0:
                continue
            seq = loc2seq(chrom,int(data[2]),int(data[3]),data[4],genome)
            data_list.append({'chrom':data[1], 'start': data[2], 'end': data[3], 'seq': seq})
        else:
            data = lines[i].split('\t')
            chrom = data[0]
            seq = loc2seq(chrom,int(data[1]),int(data[2]),data[3],genome)
            data_list.append({'chrom':data[0], 'start': data[1], 'end': data[2], 'seq': seq})

    with open(filename, 'w', encoding='utf8') as f:
        json.dump(data_list, f)

