'''
this script computes a file of chromosome sizes for 
'''

import sys

# 0. user define variables
inputFile='/Volumes/omics4tb/alomana/projects/csp.jgi/data/genome/Creinhardtii_281_v5.0.fa'
outputFile='/Volumes/omics4tb/alomana/projects/csp.jgi/data/genome/chrom.sizes.txt'

# 1. computing the sizes
g=open(outputFile,'w')

newChromosome=False
firstChromosome=True
with open(inputFile,'r') as f:
    for line in f:
        formatted=line.replace('\n','')
        
        if formatted[0] == '>' and firstChromosome == True:
            name=formatted[1:]
            size=0
            firstChromosome=False
        elif formatted[0] == '>' and firstChromosome == False:
            g.write(name+'\t'+str(size)+'\n')
            name=formatted[1:]
            size=0
        else:
            size=size+len(formatted)

g.close()
