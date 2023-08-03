# Import libraries
import pyBigWig
import numpy as np 
import pandas as pd

# Import data
total_rna_1 = pyBigWig.open("LB04_totalRNA1_forward.bigWig") # 'chr': 4,641,652
total_rna_2 = pyBigWig.open("LB04_totalRNA2_forward.bigWig")
termseq_1 = pyBigWig.open("LB_04_termseq_1_forward.bigWig")
termseq_2 = pyBigWig.open("LB_04_termseq_2_forward.bigWig")
tex_minus = pyBigWig.open("LB04_texminus_forward.bigWig")
tex_plus = pyBigWig.open("LB_0_4_dRNAseq_TEX_PLUS_FWD.bigWig")
tss = pyBigWig.open("LB04_TSSmarks_forward.bigBed")
gene_annotation = pyBigWig.open("Gene_annotations.bigBed")

# Define Functions

""" TOTAL RNA EXPRESSION
    Param: reading frame [start, stop)
    Return: region of total RNA expression [start, stop, peak, sample #]"""
def totalRNAExpression(frame):
    read_frame = frame
    df_total1 = total_rna_1.values('chr', read_frame[0], read_frame[1])
    df_total2 = total_rna_2.values('chr', read_frame[0], read_frame[1])

    total_exp = []
    start_index = 0
    stop_index = read_frame[1] - read_frame[0] -1

    """SAMPLE 1 for LB 0.4 Total RNA"""
    # Case 0: start or stop frame splits total RNA
    in_region = df_total1[start_index] != 0
    out_region = df_total1[stop_index] != 0
    if (in_region):
        while(total_rna_1.values('chr', start_index + read_frame[0], start_index + read_frame[0] + 1)[0] != 0):
            start_index -= 1
        stop_index -= start_index
        df_total1 = total_rna_1.values('chr', read_frame[0] + start_index, read_frame[1])
        read_frame = [read_frame[0] + start_index, read_frame[1]]
        start_index = 0
        in_region = False
    
    if (out_region):
        df_total1 = total_rna_1.values('chr', start_index + read_frame[0], stop_index + read_frame[0] + 1000) # note: boundary of 500 --> LATER
        while (df_total1[stop_index] != 0):
            stop_index += 1
        stop_index += 1
        # Case 0.5: circular chromosome returns to 0 --> LATER

    for i in range(stop_index-start_index):
        # Case 1: total rna expression begin + end within frame
        if(df_total1[i] != 0 and not in_region):
            total_exp.append(i)
            in_region = True
        
        if(df_total1[i] == 0 and in_region):
            total_exp.append(i-1)
            in_region = False
    
    #print(total_exp)

    total_rna = []
    for x in range(0, len(total_exp), 2):
        peak = max(total_rna_1.values('chr', total_exp[x] + read_frame[0], total_exp[x+1] + read_frame[0]))
        total_rna.append([total_exp[x] + read_frame[0], total_exp[x+1] + read_frame[0], peak, 1, "+"])
    
    #print(total_rna)
    
    """SAMPLE 2 for LB 0.4 Total RNA"""
    read_frame = frame
    total_exp = []
    start_index2 = 0
    stop_index2 = read_frame[1] - read_frame[0] -1

    # Case 0: start or stop frame splits total RNA
    in_region = df_total2[start_index2] != 0
    out_region = df_total2[stop_index2] != 0
    
    if (in_region):
        while(total_rna_2.values('chr', start_index2 + read_frame[0], start_index2 + read_frame[0] + 1)[0] != 0):
            start_index2 -= 1
        stop_index2 -= start_index2
        df_total2 = total_rna_2.values('chr', read_frame[0] + start_index, read_frame[1])
        read_frame = [read_frame[0] + start_index2, read_frame[1]]
        start_index2 = 0
        in_region = False
    
    if (out_region):
        df_total2 = total_rna_2.values('chr', start_index2 + read_frame[0], stop_index2 + read_frame[0] + 1000) # note: boundary of 1000 --> LATER********
        while (df_total2[stop_index2] != 0):
            stop_index2 += 1
        stop_index2 += 1

    for i in range(stop_index2-start_index2):
        # Case 1: total rna expression begin + end within frame
        if(df_total2[i] != 0 and not in_region):
            total_exp.append(i)
            in_region = True
        
        if(df_total2[i] == 0 and in_region):
            total_exp.append(i-1)
            in_region = False

    #print(total_exp)

    for x in range(0, len(total_exp), 2):
        peak = max(total_rna_2.values('chr', total_exp[x] + read_frame[0], total_exp[x+1] + read_frame[0]))
        total_rna.append([total_exp[x] + read_frame[0], total_exp[x+1] + read_frame[0], peak, 2, "+"])
    
    total_rna.sort()
    return total_rna, read_frame[0] - min(start_index, start_index2), read_frame[0] + max(stop_index, stop_index2)

""" KNOWN GENE Function
    Param: reading frame
    Return: Genes in region [start, end, gene name, strand]
"""

def knownGene(read_frame):
    region = gene_annotation.entries('chr', read_frame[0], read_frame[1])
    gene_list = []
    if not region: return gene_list

    for gene in region:
        gene_start = gene[0]
        gene_end = gene[1]
        gene_info = gene[2].split('\t')
        gene_name = gene_info[0]
        gene_strand = gene_info[2]

        gene_list.append([gene_start, gene_end, gene_name, gene_strand])

    return gene_list

""" FIVE END (TEX) Function
    Param: reading frame
    Return: best candidate for 5' in region [location, expression]
"""
def findFiveEnds(start, stop):
    df_texplus = tex_plus.values('chr', start, stop)
    df_texminus = tex_minus.values('chr', start, stop)
    df_tss = tss.entries('chr', start, stop)

    five_found = []

    tex_plus_found = list(map(lambda x: x[0], filter(lambda x: not np.isnan(x[1]), enumerate(df_texplus))))
    tex_minus_found = list(map(lambda x: x[0], filter(lambda x: not np.isnan(x[1]), enumerate(df_texminus))))

    # check for tex-minus peak in same bp-- later tolerance? *********************************************************************
    for tex_loc in tex_plus_found:
        if tex_loc in tex_minus_found:
            #if df_texplus[tex_loc] > df_texminus[tex_loc]:
            five_found.append([tex_loc + start, df_texplus[tex_loc]])
        else:
            five_found.append([tex_loc + start, df_texplus[tex_loc]])

    if five_found:     
        return max(five_found, key=lambda x:x[1])
    return None

    # use tss in qc ************************************************************

""" THREE END (TERMSEQ) Function
    Param: reading frame
    Return: best candidate for 3' in region [location, expression]
"""
def findThreeEnds(start, stop):
    df_three1 = termseq_1.values('chr', start, stop)
    df_three2 = termseq_2.values('chr', start, stop)

    three_found = []

    three_found1 = list(map(lambda x: x[0], filter(lambda x: x[1] != 0, enumerate(df_three1))))
    three_found2 = list(map(lambda x2: x2[0], filter(lambda x2: x2[1] != 0, enumerate(df_three2))))
    for t1 in three_found1: three_found.append([t1, df_three1[t1]])
    for t2 in three_found2: three_found.append([t2, df_three2[t2]])
   
    if three_found: 
        best = max(three_found, key=lambda x:x[1])
        return [best[0] + start, best[1]]

    return None

""" COMPARE TO KNOWN GENES Function
    Param: genes in region, candidate region
    Return: string with relationship to gene
"""
def compGenes(gene_list, region): #********************************************possible to be both antisense and 5'/3', add to an array to be returned later
    region_range = set(range(region[0], region[1]))
    #print(region_range)
    for gene in gene_list:
        gene_fwd = gene[3] == "+"
        gene_range = range(gene[0], gene[1])
        intersection = list(region_range.intersection(gene_range))
        
        if intersection:
            # Case 4: Antisense ************************************************************fix later when negative strand
            if not gene_fwd: return "antisense: " + str(gene[2])

            # Case 0: Is the gene ********************************************* create tolerance
            if gene[0] in intersection and gene[1] in intersection: return "is_gene"

            # Case 1: 5' end
            elif gene[0] in intersection: return 'five_prime: ' + str(gene[2])

            # Case 2: 3' end
            elif gene[1] in intersection: return 'three_prime: ' + str(gene[2])
            
            # Case 3: ORF-Internal
            else: return "orf_internal: " + str(gene[2])
    
    # Case 4: no intersection in any genes, intergenic
    else: return "intergenic"

def toTable(final_genes):
    df = pd.DataFrame(final_genes, columns=['Gene starts', 'Gene stops', 'Strand', '5\'', 'TSS expression', '3\'', 'Termseq expression', 'Gene classification', 'Length'])
    print(df.to_markdown())
    df.to_csv('found_genes.csv')

# Main Function: *only forward strand for now***
def main():
    # Define Reading Frame
    read_frame = [1500000,1645500 ] 

    # Find Regions of continuous Total Expression
    total_regions, start_frame, stop_frame = totalRNAExpression(read_frame) # Regions with continuous total RNA expression, updated indices

    # Final Genes With Information
    final_genes = []

    # Search for 5' and 3' near Regions of Total Expression
    for candidate in total_regions[:]:
        # check strand later *************************
        print("Region: " + str(candidate[0]) + ", " + str(candidate[1]))

        # Compare to Known Gene Annotations********************************************************
        gene_list = knownGene([candidate[0]-500, candidate[1]+500])
        if gene_list: loc = compGenes(gene_list, candidate)
        else: loc = "intergenic"
        print("found in: " + str(loc))

        # Find 5' Regions with TEX Plus
        five_end = findFiveEnds(candidate[0]-15, candidate[1]) # tolerance for distance between 5' and start of total expression ******************************
        print("5' end: " + str(five_end))

        # Find 3' Regions with Termseq
        three_end = findThreeEnds(candidate[0], candidate [1] + 15) # tolerance for distance between 3' and end of total expression ********************************
        print("3' end: " + str(three_end))

        # Remove Regions with no 5' nor 3'
        if (not three_end) or (not five_end) or (loc == "is_gene"):
            total_regions.remove(candidate)
            print("region removed")
        
        else:
            final_genes.append([candidate[0], candidate[1], candidate[4], five_end[0], five_end[1], three_end[0], three_end[1], loc, abs(three_end[0]-five_end[0])])

        print()

    print("-------------------------------------------------------------------------------------------------------------------------------------------------")

    #print(final_genes)
    
    toTable(final_genes)
    

main()