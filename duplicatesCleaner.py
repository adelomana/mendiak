import sys,os

# 0. user defined variables
bamFilesDir='/Volumes/omics4tb/alomana/projects/csp.jgi/data/bam/'
piccardPath='/Users/alomana/software/picard-tools-1.135/picard.jar'
javaPath='/usr/bin/java'
scratchDir='/Volumes/omics4tb/alomana/scratch/'

# 1. select files
inputDirs=os.listdir(bamFilesDir)

# 2. calling piccard
piccardExecutable='%s -Xmx2g -jar %s MarkDuplicates'%(javaPath,piccardPath)

for inputDir in inputDirs:

    label=inputDir.split('-')[0]

    flag1='I=%s%s/Aligned.sortedByCoord.out.bam'%(bamFilesDir,inputDir)
    flag2='O=%s%s/noDuplicates.bam'%(bamFilesDir,inputDir)
    flag3='METRICS_FILE=%spiccard.%s.duplicates.txt'%(scratchDir,label)    
    flag4='REMOVE_DUPLICATES=true'

    pieces=[piccardExecutable,flag1,flag2,flag3,flag4]
    fullCmd=' '.join(pieces)
    print()
    print(fullCmd)
    print()
    os.system(fullCmd)


    
