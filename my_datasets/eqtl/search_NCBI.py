from Bio import Entrez
import gentools

def batch_down(db_name= 'nucleotide', term_type= "WASH7P"):
    """
    Fetch the genomic information of a gene from NCBI.
    """
    # Set email for NCBI Entrez
    Entrez.email = "biff19@mails.tsinghua.edu.cn"
    
    # Search for gene in the specified NCBI database
    search_handle = Entrez.esearch(db=db_name, term=term_type, usehistory="y")
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    if db_name == 'gene':
        gene_id = search_results["IdList"][0]
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]
        
        # Get gene summary
        search_handle = Entrez.esummary(db=db_name, id=gene_id, webenv=webenv, query_key=query_key)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        result = search_results['DocumentSummarySet']['DocumentSummary'][0]
        chr = result["GenomicInfo"][0]["ChrLoc"]
        loci_t = result["GenomicInfo"][0]["ChrStart"]
        loci_p = result["GenomicInfo"][0]["ChrStop"]

    return chr, loci_t, loci_p

def getSeq(data, cpoint):
    """
    Fetch the genomic sequence of a list of genes and write to a file.
    """
    error1 = []  # Errors in searching
    error2 = []  # Errors in chromosome numbering
    genome = gentools.readGenome()
    batch_size = 50

    for start in range(cpoint, len(data), batch_size):
        end = min(len(data), start + batch_size)
        with open("/Users/bfpbr5/Documents/GitHub/genomic_benchmarks/result/result_"+str(start),"w") as f:
            for j in range(start, end):
                try:
                    loci = batch_down("gene", data[j])
                    chr = gentools.norm_chr(loci[0])
                    chrlen = len(genome[chr])
                    
                    if(loci[1] < loci[2]):
                        strand = '+'
                        loci_t = max(1, int(loci[1]) - 20000)
                        loci_p = min(chrlen, int(loci[2]) + 20000)
                    else:
                        strand = '-'
                        loci_t = max(1, int(loci[2]) - 20000)
                        loci_p = min(chrlen, int(loci[1]) + 20000)
                    
                    seq = gentools.loc2seq(genome, strand, chr=chr, start=loci_t, end=loci_p)
                except:
                    error1.append(j) if db_name == 'gene' else error2.append(j)
                    f.write('chr\t\n')
                    continue
                
                f.write(f"{chr}\t{seq}\n")
            
            print(f"Current Searching {start} records.\n")
    
    return error1, error2
