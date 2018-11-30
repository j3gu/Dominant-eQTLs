#!/usr/bin/env python
import sys
import tables
import operator
import gzip
import time
import numpy as np
import argparse
from operator import itemgetter


def main():
        sys.stderr.write("command: %s\n" % " ".join(sys.argv))

        parser = argparse.ArgumentParser(description = "get overlapped SNPs for each individual")
        parser.add_argument("chrom")
        parser.add_argument("snp_tab")
        parser.add_argument("snp_hapl")
        parser.add_argument("f_genotype")
        parser.add_argument("f_count")
        parser.add_argument("f_output")

        options = parser.parse_args()                             
       
        #Get subject IDs from phenotypes matrix
        f_subjectID = open(options.f_count,'r').readline()
        pre_subjectID = f_subjectID.split()[5:]
        subjectID = [x.split(".")[1] for x in pre_subjectID] #only 4 digits after GTEx
        target = open(options.f_output, 'w')


        #obtain all the subject IDs from genotype file
        with gzip.open(options.f_genotype, 'r') as f:
                for row in f:
                        if row.startswith('#CHROM'):
                                genotypeHeader = row.split("\t")
                                break

        allIDs = [x.split("-")[1] for x in genotypeHeader[9:]] #only 4 digits after GTEx 


        #get columns of individual genotypes 
        columns = []
        for i in range(0, len(subjectID)):
                for j in range(0, len(allIDs)):
                        if subjectID[i] == allIDs[j]:
                                columns.append(j)

        output_header = [genotypeHeader[x+9] for x in columns]
        new_header = [i.split("\n")[0] for i in output_header]

        chr_Name = options.chrom
        snps_h5 = tables.open_file(options.snp_tab)
        #contains sub-folders which represent the chromosomes
        haplo_h5 = tables.open_file(options.snp_hapl)        
        snchr = snps_h5.get_node("/%s" % chr_Name)
        hnchr = haplo_h5.get_node("/%s" % chr_Name)

        #Obtain the haplotype info for each SNP
        def getHaplotype(inMatr, colPos,rowPos):
                haplTyp = 0
                if inMatr[rowPos,colPos]< 0 or inMatr[rowPos,colPos+1]< 0:
                        return("NA")
                else:
                        haplTyp += inMatr[rowPos,colPos]
                        haplTyp += inMatr[rowPos,colPos+1]
                        return (haplTyp)

        snpHaplTypeInfo = 'Chr' + "\t" + "id" +  "\t" + "\t".join(new_header) + '\n'
        target.write(snpHaplTypeInfo)
        #loop through all the SNPs
        count = 0
        #start_time = time.time()
        #sys.stderr.write("loop num:%d" % hnchr.shape[0])
        snp_lists = []
        for i in range(0,hnchr.shape[0]):
                snp_name = snchr[i]['name']
                temp_snpHapl = chr_Name + "\t" + snp_name

                count +=1
                if count > 500:
                        sys.stderr.write(".")
                        count = 0
                        #elapsed_time = time.time() - start_time
                        #sys.stderr.write("--- %s seconds ---" % elapsed_time)
                        #start_time = time.time()
                #remove duplicated SNPs        
                if snp_name in snp_lists:
                        continue
                else:
                        snp_lists.append(snp_name)

                        #adding up the number of non-reference alleles from each haplotype
                        overlapped_snps = []
                        for j in range(0,len(columns)):
                                overlapped_snps.append(getHaplotype(hnchr,2*columns[j],i))

                        #remove NAs to compute for minor allele frequency
                        indiv_snps = [x for x in overlapped_snps if x is not "NA"]
                        minor_allele = min((2*len(indiv_snps)-sum(indiv_snps)), sum(indiv_snps))
                        maf = float(minor_allele)/(2*len(indiv_snps))

                        #filter by having at least four individuals per genotype group
                        min_num = min(sum([i == 1 for i in indiv_snps]), sum([i == 0 for i in indiv_snps]), sum([i == 2 for i in indiv_snps]))
                        if maf >= 0.05 and min_num > 3:#filter out MAF < 0.05 and snps with less than 4 indivs per genotype group
                                temp_hapl = '\t' + "\t".join([str(i) for i in overlapped_snps])
                                temp_snpHapl += temp_hapl
                                temp_snpHapl += '\n'
                                target.write(temp_snpHapl)

        target.close()
main()


        
