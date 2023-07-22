


filename = "/Users/bfpbr5/Documents/GitHub/new_expr_loci2.txt"
outname = "/Users/bfpbr5/Documents/GitHub/expr_loci_nones"
nones = []
i = 0
with open(filename, 'r') as f:
    while True:
        line = f.readline()
        if not line:
            break
        if line[0] == 'n':
            nones.append(str(i))
        i += 1

with open (outname, 'w') as f:
    for none in nones:
        f.write(none+'\n')