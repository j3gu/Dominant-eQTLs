#!/usr/bin/python
import sys
import argparse

def main():
   # sys.stderr.write("arguments:%s\n" % "".join(sys.argv))
    parser = argparse.ArgumentParser(description = "Divide into blocks of genes per file")
    parser.add_argument("f_expr")
    parser.add_argument("fixed_num")
    parser.add_argument("f_output", nargs = '+')

    options = parser.parse_args()
    
    fnames = options.f_output
    f = open(options.f_expr, 'r')
    #header
    header = f.readline()
    split_index = open(options.fixed_num).readlines()
    
    s = {}
    start = ""
    for i in range(0, len(fnames)):
        term = "_".join(fnames[i].split('/')[-1].split("_")[0:2])
        s[term] = open(fnames[i], 'w')

    for line in split_index:
        col = line.split()
        chrom_name = col[0]
        index_num = col[1]
        s[chrom_name] = index_num
        
        
    d = {}
    initial = ""
    for line in f.readlines():
        col = line.split()
        chrom = col[1]
        if chrom == initial:
            d[chrom].append(line)
        elif chrom in ["chrX", "chrY", "chrMT"]:
            continue
        else:
            d[chrom] = []
            d[chrom].append(line)
        initial = chrom

    for key in sorted(d.keys()):
        mat = d[key]
        num_files = int(s[key])
        num = len(mat)/float(num_files)
        #deal with integer operation problem
        if num > int(num):
            num = int(num) + 1
        else:
            num = int(num)
        z = 0
        #sys.stderr.write("fixed_nums:%d" % num)
        last_num = range(0, len(mat), num)[-1]
        for i in range(0, len(mat), num):
            name = key + "_subset" + str(100*z)
            s[name].write("%s" % header)
            if i < last_num:
                for j in range(i, (i+num)):
                    s[name].write("%s" % mat[j])
            else:
                for j in range(i, len(mat)):
                    s[name].write("%s" % mat[j])       
            s[name].close()
            z += 1
    f.close()

main()

 
    
