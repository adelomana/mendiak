import os,sys

'''
this script finds the clean FASTQ files and calls STAR for the reads alignment.
'''

def genomeIndexer():

    '''
    this function creates the genome index. Should be run only once.
    '''

    flag1=' --runMode genomeGenerate'
    flag2=' --runThreadN %s'%numberOfThreads
    flag3=' --genomeDir %s'%genomeIndexDir
    flag4=' --genomeFastaFiles %s'%genomeFastaFile
    flag5=' --sjdbGTFfile %s'%genomeAnnotationFile
    flag6=' --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 99'

    cmd=STARexecutable+flag1+flag2+flag3+flag4+flag5+flag6

    print()
    print(cmd)
    print()
    os.system(cmd)

    return None

def STARcalling(inputFile):

    '''
    this function calls STAR
    '''
    
    finalDir=bamFilesDir+inputFile+'/'
    if os.path.exists(finalDir) == False:
        os.mkdir(finalDir)

    readF1=readsFilesDir+inputFile

    flag0='time '
    flag1=' --genomeDir %s'%genomeIndexDir
    flag2=' --runThreadN %s'%numberOfThreads
    flag3=' --readFilesIn %s'%readF1
    flag4=' --outFileNamePrefix %s'%finalDir
    flag5=' --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 9976258636 --outFilterMismatchNoverLmax 0.0'

    cmd=flag0+STARexecutable+flag1+flag2+flag3+flag4+flag5
    
    # writing a sender file
    flag=inputFile.split('-')[0]

    senderFile=sendersDir+flag+'.sh'
    f=open(senderFile,'w')
    f.write('#!/bin/bash\n\n')
    f.write('#$ -N %s\n'%flag)
    f.write('#$ -o /proj/omics4tb/alomana/scratch/csp/messages.%s.o.txt\n'%flag)
    f.write('#$ -e /proj/omics4tb/alomana/scratch/csp/messages.%s.e.txt\n'%flag)
    f.write('#$ -P Bal_alomana\n')
    f.write('#$ -pe serial %s\n'%str(numberOfThreads))
    f.write('#$ -q baliga\n')
    f.write('#$ -S /bin/bash\n\n')
    f.write('cd /users/alomana\n')
    f.write('source .bash_profile\n\n')
    f.write('%s\n'%cmd)
    f.close()

    # submitting a sender file
    os.system('qsub %s'%senderFile)
    
    return None

# 0. defining several input/output paths
readsFilesDir='/proj/omics4tb/alomana/projects/csp.jgi/data/fastq/clean/'
bamFilesDir='/proj/omics4tb/alomana/projects/csp.jgi/data/bam/'
STARexecutable='/proj/omics4tb/alomana/software/STAR-master/bin/Linux_x86_64/STAR'
genomeIndexDir='/proj/omics4tb/alomana/projects/csp.jgi/data/genomeIndex'
genomeFastaFile='/proj/omics4tb/alomana/projects/csp.jgi/data/genome/Creinhardtii_281_v5.0.fa'
genomeAnnotationFile='/proj/omics4tb/alomana/projects/csp.jgi/data/genome/Creinhardtii_281_v5.5.gene_exons.gff3'
numberOfThreads=40
sendersDir='senders/'

# 1. recover the clean FASTQ files
print('reading FASTQ files...')
allTags=[]
allFiles=os.listdir(readsFilesDir)
inputFiles=list(set(allFiles))

# 2. making genome indexes
#print('making genome index...')
#genomeIndexer()
       
# 3. calling STAR through SGE...
print('calling STAR through SGE...')
for inputFile in inputFiles:
    STARcalling(inputFile)

