#!/usr/bin/python
import sys
import argparse

def main():
    sys.stderr.write("arguments:%s\n" % "".join(sys.argv))
    parser = argparse.ArgumentParser(description = "merge genotype files into one")
    parser.add_argument("f_output")
    parser.add_argument("f_genotypes", nargs = '+')

    options = parser.parse_args()
    
    fnames = options.f_genotypes
    output = open(options.f_output, 'w')

    d = {}
    for name in fnames:
        chrom = name.split("/")[-1].split("_")[0]
        d[chrom] = open(name, 'r')
        
    header = d["chr22"].readline()
    output.write(header)
    sys.stderr.write("arguments:%s\n" % "".join(d.keys()))
    for key in sorted(d.keys()):
        header = d[key].readline()
        for line in d[key].readlines():
            output.write(line)

    output.close()

main()

 
    
