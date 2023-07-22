import json
import os
import random

import Bio
import pandas as pd

import to_json

# from pyjaspar import jaspardb



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
    if chr[0] != 'c':
        result = "chr"+chr
        return result
    if chr[2] == 'r':
        if (chr[3].isdecimal()) or (chr[3] in ['X', 'Y']):
            return chr
        result = 'chr'+str(romanToInt(chr[3:]))
        return result
    elif chr[1] == 'h':
        result = 'chr'+str(int(chr[2:]))
        return result
    else:
        return 'chr0'

def longer_loc(start, end, chrstart, chrend, var=500):
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




genome = to_json.readGenome('/Users/bfpbr5/Downloads/hg38.analysisSet.fa.gz')
# jdb_obj = jaspardb()
refs = {}
pos1s = {}
pos2s = {}
errors = []
df = pd.read_csv("~/maf_test.csv")
df1 = df.loc[:,["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Strand", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"]]
for idx,data in df1.iterrows():
    gene = data[0]
    chr = data[1]
    start_pos = data[2]
    end = data[3]
    strand = data[4]
    ref_allele = data[5]
    pos_allele1 = data[6]
    pos_allele2 = data[7]
    if (chr[-1] == 'L') or (chr[-1] == 'R') :
            chr = chr[:-1]
    chr = norm_chr(chr)
    try:
        chrend = len(genome[chr]) - 1
    except (KeyError):
        errors.append({'id': gene, 'chr':chr})
        continue
    start, end = longer_loc(start_pos, end, 1, chrend, var=400)
    ref_seq = snv2seq(chr, start, end, strand, ref_allele, start_pos, genome)
    pos1_seq = snv2seq(chr, start, end, strand, pos_allele1, start_pos, genome)
    pos2_seq = snv2seq(chr, start, end, strand, pos_allele2, start_pos, genome)
    ref = {"gene":gene, "seq":ref_seq}
    pos1 = {"gene":gene, "seq":pos1_seq}
    pos2 = {"gene":gene, "seq":pos2_seq}
    refs.append(ref)
    pos1s.append(pos1)
    pos2s.append(pos2)

# j = 0
# # class_sample = ["MA1273.1"]
# for key in class_sample.keys():
#     motifs[key] = []
#     for mid in class_sample[key]:
#         # if j < 3:
#         #     j += 1
#         #     continue
#         # motif = jdb_obj.fetch_motif_by_id(mid)
#         motif = 1
#         if motif:
#             one_motif = {}
#             one_motif['id'] = mid
#             one_motif['name'] = motif.name
#             one_motif['class'] = motif.tf_class
#             one_motif['family'] = motif.tf_family

#         seqs = []
#         bedname = bedpath+mid+'.bed'
#         with open(bedname,"r",encoding="utf8") as f:
#             line = f.readline()
#             i = 0
#             while line:
#                 if line == '\n':
#                     line = f.readline()
#                     continue
#                 data = line.split('\t')
#                 chr = data[0]
#                 start = data[1]
#                 end = data[2]
#                 if (chr[-1] == 'L') or (chr[-1] == 'R') :
#                         chr = chr[:-1]
#                 chr = norm_chr(chr)
#                 try:
#                     chrend = len(genome[chr]) - 1
#                 except (KeyError):
#                     errors.append({'id': mid, 'chr':chr})
#                     i += 1
#                     line = f.readline()
#                     continue
#                 start, end = longer_loc(start, end, 1, chrend, var=400)
#                 strand = data[-1][0]
#                 seq = to_json.loc2seq(chr, start, end, strand, genome)
#                 if ('N' not in seq) and (seq != ''):
#                     seqs.append(seq)
#                 line = f.readline()
#                 print("line %d done\n"%(i))
#                 i += 1
#                 if i > 100:
#                     break
#         one_motif['seq'] = seqs
#         motifs[key].append(one_motif)

    # j += 1
    # print("="*38+str(j)+"="*38)
outname = "tcga_lung_pos1.json"
with open(outname, 'w', encoding='utf8') as f:
    json.dump(pos1s, f)

outname = "tcga_lung_pos2.json"
with open(outname, 'w', encoding='utf8') as f:
    json.dump(pos2s, f)

outname = "tcga_lung_neg.json"
with open(outname, 'w', encoding='utf8') as f:
    json.dump(refs, f)
