import pandas as pd

# Read a tab-delimited file
def read_mat(filename, header=0):
    data = pd.read_table(filename, sep="\t", header=header)
    return data

# Extract a list of gene IDs from the first column of the DataFrame
def get_idlist(gtf_df):
    result = [row[1][:15] for row in gtf_df.itertuples()]
    return result

# Search for a gene's transcription start site and strand direction in the GTF DataFrame
def search_gtf(enid, gtf):
    for row in gtf.itertuples():
        attrs = row[-1].split(';')
        if(enid == attrs[0][9:-1]):
            tss = (row[1],row[4])
            strand = row[7]
            return tss, strand
    return 'none', 'none'

if __name__=='__main__':
    # File paths
    gtf_name = "/Users/bfpbr5/Documents/GitHub/genomic_benchmarks/data/Homo_sapiens.GRCh38.107.chr.gtf"
    gct_name = "/Users/bfpbr5/Documents/GitHub/genomic_benchmarks/data/eQTL/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"
    
    # Read input files
    gct = read_mat(gct_name)
    gtf = read_mat(gtf_name, None)

    # Get list of gene IDs from GCT file
    idlist = get_idlist(gct)
    
    # Open output file
    outname = "expr_loci"
    f = open(outname, 'w')
    
    # Search for each gene in the GTF file and write results to output file
    for i, gene in enumerate(idlist, start=1):
        tss, strand = search_gtf(gene, gtf)
        f.write(f'{tss[0]}\t{tss[1]-1}\t{strand}\n')
        print(f"{i} genes found.")
    
    # Close output file
    f.close()
