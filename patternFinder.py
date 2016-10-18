import os

# 0. user defined variables
peaksDir='/Volumes/omics4tb/alomana/projects/csp.jgi/data/macs2.run3/'

correspondance={}
correspondance['0hA']='ASCAO'
correspondance['0hB']='ASCAP'
correspondance['24hA']='ASCAS'
correspondance['24hB']='ASCAT'
correspondance['48hA']='ASCAU'
correspondance['48hB']='ASCAW'

# 1. selecting the samples
allFiles=os.listdir(peaksDir)
peaksFileNames=[element for element in allFiles if '_peaks.xls' in element if 'callerC' in element]

# 2. filter peaks: at least 5-fold and no longer than 1 kb
for peaksFileName in peaksFileNames:
    peaksFile=peaksDir+peaksFileName
    print(peaksFile)

# 2.1. filtering peaks

# 2.2. plot the distribution of peaks before and after filtering peaks


# 3. define all genes that have matching patterns
# 3.1. 101 signature

# 3.1. 010 signature

# 4. define significance
