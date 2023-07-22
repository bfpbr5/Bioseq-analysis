

loci_file = "new_expr_loci2.txt"
expr_file = "genomic_benchmarks/data/eQTL/gtex_expr.gct"
metadata = "genomic_benchmarks/my_datasets/eqtl/metadata.txt"
train = "genomic_benchmarks/my_datasets/eqtl/train.txt"
valid = "genomic_benchmarks/my_datasets/eqtl/valid.txt"
test = "genomic_benchmarks/my_datasets/eqtl/test.txt"
validchr = ["chr6","chr9"]
testchr = ["chr7","chr8"]
with open(loci_file,"r") as f1, open(expr_file, 'r') as f2, open(train, 'w') as tr, open(valid, 'w') as vl, open(test, 'w') as tt:
    names = f2.readline().strip('\n').split('\t')
    names = names[2:]
    with open(metadata,'w') as mt:
        mt.write(str(names))

    while True:
        line = {}
        loci = f1.readline().strip('\n').split('\t')
        expr = f2.readline().strip('\n').split('\t')
        if loci[0] == 'n':
            continue
        if not loci[0]:
            break

        loci[0] = 'chr'+str(loci[0])
        loci.append(loci[2])
        loci[2] = int(loci[1])
        loci[1] = int(loci[1])-1
        line['index'] = loci
        line['label'] = expr[2:]

        if(loci[0] in validchr):
            vl.write(str(line)+'\n')
        elif(loci[0] in testchr):
            tt.write(str(line)+'\n')
        else:
            tr.write(str(line)+'\n')
