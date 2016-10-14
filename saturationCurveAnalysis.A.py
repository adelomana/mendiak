'''
this script creates the FASTQ files for the saturation curve
'''

import os,sys,shutil


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
        print(allReads)

    # f.2. go over resolution decrements, saving FASTQ files of appropriate size
    allKeys=list(allReads.keys())
    numberOfReads=len(allKeys)
    mr=int(numberOfReads/1e6)
    while mr > 1:

        print('selecting %s mr...'%mr)

        # select n random reads
        k=int(mr*1e3) ### THIS SHOULD BE 1e6!!!!
        print(k)
        selectedReads=random.sample(allKeys,k)
        print(selectedReads[:10]+len(selectedReads))

        # write a file
        outputFile=piecesDir+sample+'resolution.%s.fastq'%str(mr)
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
readsFilesDir='/Users/alomana/scratch/clean/'
piecesDir='/Users/alomana/scratch/pieces/'
resolution=3

correspondance={}
correspondance['0hA']='ASCAO'
correspondance['0hB']='ASCAP'
correspondance['24hA']='ASCAS'
correspondance['24hB']='ASCAT'
correspondance['48hA']='ASCAU'
correspondance['48hB']='ASCAW'
correspondance['test']='test'

# 1. iterate over samples
print('iterating over samples...')
for sample in correspondance.keys():

    sample='test'
    
    # 2. build a pipeline and submit it to the cluster
    print('working with sample',sample,'...')
    butcher(sample)
