#!/usr/bin/python

import tables
import numpy as np
import sys
import argparse


def main():
    sys.stderr.write("command: %s\n" % " ".join(sys.argv))
    parser = argparse.ArgumentParser(description = "get SNPs that are +/- 100kb away from exons")
    parser.add_argument("snp_genotype")
    parser.add_argument("snp_index")
    parser.add_argument("f_overlappedSNPs")
    parser.add_argument("f_counts")
    parser.add_argument("output")
    
    options = parser.parse_args()
    f_input = options.f_counts
    chrom = f_input.split("/")[-1].split("_")[0]
    counts_matrix = open(f_input, 'r')
    SNPs_matrix = open(options.f_overlappedSNPs, 'r')
    f_genotype = tables.open_file(options.snp_genotype)
    f_snp_index = tables.open_file(options.snp_index)
    output = options.output
    f_output = open(output, 'w')

       
    h_counts = counts_matrix.readline().split()
    rows = counts_matrix.readlines()   
    # open snp_tab.h5 and index file and define nodes for data
    node = f_genotype.get_node("/%s" % chrom)
    index_node = f_snp_index.get_node("/%s" % chrom)
    index_array = index_node[:]
    
    # #set up a dictionary to quickly match SNP_id with the corrsponding row in the genotype matrix
    # d = {}
    # SNPs_header = SNPs_matrix.readline().split()[2:]
    # #Create a list which containts columns
    # col_num = []
    # newHeader_SNPs = []
    # for i in range(5, len(h_counts)):
    #     if h_counts[i] in array:
    #         selected_col = np.where(array == h_counts[i])[0][0] 
    #         col_num.append(selected_col)
    #         newHeader_SNPs.append(h_counts[i])
    #generate an array of indexes for selected columns corresponding to list of individuals
    
    SNPs_header = SNPs_matrix.readline().split()[2:]
    d = {}
    for line in SNPs_matrix:
        col = line.split()
        snp_chr = col[0]
        if snp_chr == chrom:
            genotype = col[2:]
            SNP_id = col[1]
            if(genotype.count("0") > 3) and (genotype.count("1") > 3) and (genotype.count("2") > 3):
                d[SNP_id] = col[2:]
    # #find cutoff of count matrix
    # percentile = []
    # for row in total_counts:
    #     col = np.array(row.split()[5:])
    #     num = np.percentile(col, 75)
    #     if num  == 0:
    #         continue
    #     else:
    #         percentile.append(num)
    # cutoff = np.log10(percentile)

    count = 0
    new_header = ["GENE", "EXONS.START", "EXONS.END", "LENGTH", "SNP_ID"] + SNPs_header + h_counts 
    #Add header to the output file
    f_output.write("\t".join(new_header) + '\n')
    #loop over counts matrix amnd use the coordinates of peaks to find SNPs within +/- 100kb away from the peaks
    for i in range(0, len(rows)):
        col = rows[i].split()
        gene_name = col[0]
        exon_start = int(col[2])
        exon_end = int(col[3])
        exon_min = exon_start - 100000
        if exon_min <= 0:
            exon_min = 1
        exon_max = exon_end + 100000
        t_length = int(col[4])
        index = index_array[(exon_min-1):exon_max]
       #sys.stderr.write("index_length: %d" % len(np.where(index != -1)[0]))
        count += 1
        if count > 50: #stand out the dots which represent 50 loops have been run for each dot
            sys.stderr.write(".")
            count = 0 
        #loop over the SNP index to quickly look up for the SNP information from node
        for j in np.where(index != -1)[0]:  #np.where: return elements eithxer from [0](x) or [1](y) based on the given conditions
            SNP = node[index[j]]
            SNP_name = SNP['name']
            if SNP_name in d:
                genotype = d[SNP_name]
                SNP_overlaps = [str(x) for x in genotype]
                norm_counts = [str(y) for y in col[0:]]
                line =  SNP_overlaps + norm_counts
                f_output.write("%s\t%d\t%d\t%d\t%s\t%s\n" % (gene_name, exon_start, exon_end, t_length, SNP_name, "\t".join(line)))

    f_output.close()
    counts_matrix.close()
    SNPs_matrix.close()

main()


    
        
        


