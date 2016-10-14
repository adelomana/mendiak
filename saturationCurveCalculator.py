'''
this script builds saturation curves for number of reads for available samples.
we are using pooled ICs 
'''

import os,sys

def figuresMaker():

    '''
    function to creates the figures for the saturation analysis
    '''

    return None

def pipeline(sample):

    '''
    computational steps for the saturation curve
    '''

    # f.0. define the number of jobs: mr!
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
    print allReads

    # f.1. select random reads (also consider the original sample)

    # f.2. create a sender file that map the reads and call peaks

    # f.3. submit job

    return None

# 0. user defined variables
readsFilesDir='/proj/omics4tb/alomana/projects/csp.jgi/data/fastq/clean/'
STARexecutable='/proj/omics4tb/alomana/software/STAR-master/bin/Linux_x86_64/STAR'
scratchDir='/proj/omics4tb/alomana/scratch/csp/'
numberOfThreads=1
resolution=5

correspondance={}
correspondance['0hA']='ASCAO'
correspondance['0hB']='ASCAP'
correspondance['24hA']='ASCAS'
correspondance['24hB']='ASCAT'
correspondance['48hA']='ASCAU'
correspondance['48hB']='ASCAW'
correspondance['test']='test'

# 1. iterate over samples
print 'iterating over samples...'
for sample in correspondance.keys():

    sample='test'
    
    # 2. build a pipeline and submit it to the cluster
    print 'working with sample',sample,'...'
    pipeline(sample)

# 3. build figures from tables of peaks
figuresMaker()
