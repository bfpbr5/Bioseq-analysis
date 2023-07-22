import json
import os
import random
from pyjaspar import jaspardb
import pandas as pd
import to_json

# Convert Roman numerals to integers
def romanToInt(s):
    result = 0
    temp_dict = {"I":1,"V":5,"X":10,"L":50,"C":100,"D":500,"M":1000}
    for i in range(len(s)-1):
        if temp_dict[s[i]]<temp_dict[s[i+1]]:
            result -= temp_dict[s[i]]
        else:
            result += temp_dict[s[i]]
    result += temp_dict[s[-1]]
    return result

# Normalize chromosome notation
def norm_chr(chr):
    if chr[0] != 'c':
        result = "chr"+chr
        return result
    if chr[2] == 'r':
        if chr[3].isdecimal():
            return chr
        result = 'chr'+str(romanToInt(chr[3:]))
        return result
    elif chr[1] == 'h':
        result = 'chr'+str(int(chr[2:]))
        return result
    else:
        return 'chr0'

# Randomize the location by adding a random offset
def randomize_loc(start, end, chrstart, chrend, var=400):
    rd = random.random() * var + 100
    start = int(start - rd)
    end = int(end + rd)
    if start < int(chrstart):
        start = chrstart
    if end > int(chrend):
        end = chrend
    return int(start), int(end)

# Check if a sequence is within transcription factor binding sites
def is_in_seq(seq, tfs):
    if seq[0] != tfs[0]:
        return False
    if (int(seq[2]) - int(tfs[1])) * (int(seq[1])-int(tfs[2])) > 0:
        return False
    return True

# Randomly generate a location in the genome
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

# Generate negative samples by finding sequences not in transcription factor binding sites
def neg_generator(genome, randn=360):
    neg_result = []
    n = 0
    while n < randn:
        rand_seq = rand_loc(genome)
        rand_in_tfs = False
        with open("/Users/bfpbr5/Documents/GitHub/genomic_benchmarks/data/bed/all.bed", 'r', encoding='utf8') as f:
            line = f.readline()
            while line:
                if line[0] in ['\n', '\\']:
                    line = f.readline()
                    continue
                data = line.split('\t')
                loc = (data[0], int(data[1]), int(data[2]))
                rand_in_tfs = is_in_seq(rand_seq, loc)
                line = f.readline()
        if not rand_in_tfs:
            rand_seq = to_json.loc2seq(rand_seq[0], rand_seq[1], rand_seq[2], '+', genome)
            if 'N' not in rand_seq:
                neg_result.append(rand_seq)
                n += 1
    return neg_result

genome = to_json.readGenome('/Users/bfpbr5/Downloads/hg38.analysisSet.fa.gz')
bedpath = "/Users/bfpbr5/Documents/GitHub/genomic_benchmarks/data/bed/"
motif_ids = [x[:-4] for x in os.listdir(bedpath) if x.endswith('.bed')]

noSeqFilename = "/Users/bfpbr5/Documents/GitHub/jaspar_no_seq.json"
df = pd.read_json(noSeqFilename,encoding="utf-8", orient="records")
df['class'] = df['class'].apply(lambda row: row[0] if row else '')
df['family'] = df['family'].apply(lambda row: row[0] if row else '')
classes = df['class'].unique()

# Generate a sample of transcription factors for each class
class_sample = {}
for item in list(classes):
    sample_list = list(df['id'].loc[df['class'] == item])
    sn_max = 4 # maximum TFs per class
    class_sample[item] = sample_list[:sn_max]

# Fetch motifs from JASPAR database
jdb_obj = jaspardb()
motifs = {}
errors = []
for key in class_sample.keys():
    motifs[key] = []
    for mid in class_sample[key]:
        motif = jdb_obj.fetch_motif_by_id(mid)
        if motif:
            one_motif = {'id': mid, 'name': motif.name, 'class': motif.tf_class, 'family': motif.tf_family}

        seqs = []
        bedname = bedpath+mid+'.bed'
        with open(bedname,"r",encoding="utf8") as f:
            line = f.readline()
            while line:
                data = line.split('\t')
                chr = norm_chr(data[0])
                try:
                    chrend = len(genome[chr]) - 1
                except KeyError:
                    errors.append({'id': mid, 'chr':chr})
                    line = f.readline()
                    continue
                start, end = randomize_loc(data[1], data[2], 1, chrend)
                strand = data[-1][0]
                seq = to_json.loc2seq(chr, start, end, strand, genome)
                if 'N' not in seq and seq != '':
                    seqs.append(seq)
                line = f.readline()
        one_motif['seq'] = seqs
        motifs[key].append(one_motif)

# Save positive samples
outname = "jaspar_pos.json"
with open(outname, 'w', encoding='utf8') as f:
    json.dump(motifs, f)

# Generate and save negative samples
neg = neg_generator(genome)
neg_result = {'negative_seq':neg}
with open('neg_seq.json', 'w', encoding='utf8') as f:
    json.dump(neg_result, f)

with open('jaspar_pos.json', 'r', encoding='utf8') as f:
    result = json.load(f)
