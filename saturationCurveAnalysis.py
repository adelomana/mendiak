'''
this script creates the FASTQ files for the saturation curve
'''

import os,sys,random


def butcher(sample):

    '''
    this function creates crumbles of FASTQ files for curve saturation analysis
    '''

    # f.0. read reads (NPI)
    possibleFiles=os.listdir(readsFilesDir)
    inputFile=readsFilesDir+[element for element in possibleFiles if correspondance[sample] in element][0]

    allReads={}    
    with open(inputFile,'r') as f:
        count=0
        for line in f:
            vector=line[:-2]
            count=count+1
            if count == 1:
                readName=vector
            if count == 2:
                sequence=vector
            if count == 4:
                quality=vector
                allReads[readName]=[sequence,quality]
                count=0

    # f.2. go over resolution decrements, saving FASTQ files of appropriate size
    allKeys=list(allReads.keys())
    numberOfReads=len(allKeys)
    mr=int(numberOfReads/1e6) 

    
    while mr >= 1:

        print 'selecting %s mr...'%mr

        # select n random reads
        k=int(mr*1e6) 
        selectedReads=random.sample(allKeys,k)
       
        # write a file
        outputFile=piecesDir+sample+'.resolution.%s.fastq'%str(mr)
        g=open(outputFile,'w')
        for read in selectedReads:
            g.write('%s\n'%read)
            g.write('%s\n'%allReads[read][0])
            g.write('+\n')
            g.write('%s\n'%allReads[read][1])
        g.close()

        # next iteration
        mr=mr-resolution

    return None

# 0. user defined variables
readsFilesDir='/proj/omics4tb/alomana/projects/csp.jgi/data/fastq/clean/'
piecesDir='/proj/omics4tb/alomana/projects/csp.jgi/data/fastq.crumbles/'
resolution=2

correspondance={}
correspondance['0hA']='ASCAO'
#correspondance['0hB']='ASCAP'
correspondance['24hA']='ASCAS'
correspondance['24hB']='ASCAT'
#correspondance['48hA']='ASCAU'
correspondance['48hB']='ASCAW'

# 1. iterate over samples
for sample in correspondance.keys():
    
    # 2. create FASTQ crumbles
    print 'working with sample'+sample+'...'
    butcher(sample)

# 3. figure maker
