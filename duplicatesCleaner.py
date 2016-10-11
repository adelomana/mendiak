import sys,os

# 0. user defined variables
bamFilesDir='/Volumes/omics4tb/alomana/projects/csp.jgi/data/bam/'
scratchDir='/Volumes/omics4tb/alomana/scratch/'

# 1. select files
inputDirs=os.listdir(bamFilesDir)

# 2. calling piccard
samToolsExecutable='time samtools rmdup -s'

for inputDir in inputDirs:

    label=inputDir.split('-')[0]

    flag1='%s%s/Aligned.sortedByCoord.out.bam'%(bamFilesDir,inputDir)
    flag2='%s%s/removedDuplicates.bam'%(bamFilesDir,inputDir)

    pieces=[samToolsExecutable,flag1,flag2]
    fullCmd=' '.join(pieces)
    print()
    print(fullCmd)
    print()
    os.system(fullCmd)


    
