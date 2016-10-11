'''
this script calls MACS 2.0 to call peaks from aligned bam files.
'''

import os,sys

def callerA():

    '''
    calling MACS2.0 using their corresponding ICs
    '''

    for experiment in experiments:

        inputFile=bamDir+[element for element in bamFiles if experiment[0] in element][0]+'/removedDuplicates.bam'
        controlFiles=[bamDir+[element for element in bamFiles if experiment[1] in element][0]+'/removedDuplicates.bam']
        name='callerA.'+experiment[2]

        universalCaller(inputFile,controlFiles,name)

    return None

def callerB():

    '''
    calling MACS2.0 using their sample ICs
    '''

    for experiment in experiments:

        inputFile=bamDir+[element for element in bamFiles if experiment[0] in element][0]+'/removedDuplicates.bam'

        tag=experiment[2][:-1]
        labels=[element[1] for element in experiments if tag in element[2]]
        controlFiles=[]
        for bam in bamFiles:
            for label in labels:
                if label in bam:
                    controlFiles.append(bamDir+bam+'/removedDuplicates.bam')
        name='callerB.'+experiment[2]

        universalCaller(inputFile,controlFiles,name)

    return None

def callerC():

    '''
    calling MACS2.0 using pooled ICs
    '''

    for experiment in experiments:

        inputFile=bamDir+[element for element in bamFiles if experiment[0] in element][0]+'/removedDuplicates.bam'

        labels=[element[1] for element in experiments]
        controlFiles=[]
        for bam in bamFiles:
            for label in labels:
                if label in bam:
                    controlFiles.append(bamDir+bam+'/removedDuplicates.bam')
        name='callerC.'+experiment[2]

        universalCaller(inputFile,controlFiles,name)

    return None

def callerD():

    '''
    calling MACS2.0 using exhaustively all single ICs
    '''

    for iter1 in experiments:
        for iter2 in experiments:
            inputFile=bamDir+[element for element in bamFiles if iter1[0] in element][0]+'/removedDuplicates.bam'
            controlFiles=[bamDir+[element for element in bamFiles if iter2[1] in element][0]+'/removedDuplicates.bam']
            name='callerD.'+iter1[2]+'.'+iter2[2]

            universalCaller(inputFile,controlFiles,name)

    return None

def SGEsubmitter(commands,jobName):

    '''
    this function submits a job to the cluster under the SGE environment
    '''

    # writing a sender file
    senderFile=sendersDir+jobName+'.sh'
    f=open(senderFile,'w')
    f.write('#!/bin/bash\n\n')
    f.write('#$ -N %s\n'%jobName)
    f.write('#$ -o %s/messages.%s.o.txt\n'%(scratchDir,jobName))
    f.write('#$ -e %s/messages.%s.e.txt\n'%(scratchDir,jobName))
    f.write('#$ -P Bal_alomana\n')
    f.write('#$ -pe serial %s\n'%str(numberOfThreads))
    f.write('#$ -q baliga\n')
    f.write('#$ -S /bin/bash\n\n')
    f.write('cd /users/alomana\n')
    f.write('source .bash_profile\n\n')
    for command in commands:
        f.write('%s\n\n'%command)
    f.close()

    # submitting a sender file
    os.system('qsub %s'%senderFile)
    sys.exit()

    return None

def universalCaller(inputFile,controlFiles,name):

    '''
    this function does the actual call, with specified input and output files
    '''

    # f.1. peak calling
    formattedControlFiles=' '.join(element for element in controlFiles)
    
    flag0='time'
    flag1='callpeak'
    flag2='-t '+inputFile
    flag3='-c '+formattedControlFiles
    flag4='-n '+name
    flag5='--outdir '+outputDir
    flag6='-g 100e6'
    flag7='-q 0.05'
    flag8='--call-summits'
    flag9='--verbose 4'
    flag10='--bw 150'
    flag11='-m 2 50'
    flag12='--bdg'
    flag13='> %smessages.%s.txt'%(scratchDir,name)
    
    pieces=[flag0,macs2Executable,flag1,flag2,flag3,flag4,flag5,flag6,flag7,flag8,flag9,flag10,flag11,flag12,flag13]
    fullCmd1=' '.join(pieces)

    # f.2. working with the generation of FE and logLR files
    inputTrack=outputDir+name+'_treat_pileup.bdg'
    controlTrack=outputDir+name+'_control_lambda.bdg'
    flag0='time'
    flag1='bdgcmp'
    flag2='-t '+inputTrack
    flag3='-c '+controlTrack
    flag4='-o %s%s_FE.bdg'%(outputDir,name)
    flag5='-m FE'
    
    pieces=[flag0,macs2Executable,flag1,flag2,flag3,flag4,flag5]
    fullCmd2=' '.join(pieces)

    flag6='-o %s%s_logLR.bdg'%(outputDir,name)
    flag7='-m logLR -p 0.00001'

    pieces=[flag0,macs2Executable,flag1,flag2,flag3,flag6,flag7]
    fullCmd3=' '.join(pieces)

    # f.3. sorting bedgraph file
    allFiles=os.listdir(outputDir)
    bedGraphFiles=[name+'_control_lambda.bdg',name+'_treat_pileup.bdg',name+'_FE.bdg',name+'_logLR.bdg']
    sortingCommands=[]
    for element in bedGraphFiles:
        sortedFile=element.replace('.bdg','.sorted.bdg')
        flag0='time'
        executable='sort'
        flag1='-k1,1 -k2,2n'
        flag2='%s%s'%(outputDir,element)
        flag3='> %s%s'%(outputDir,sortedFile)
        pieces=[flag0,executable,flag1,flag2,flag3]
        instance=' '.join(pieces)
        sortingCommands.append(instance)
        
    # f.4. converting bedgraph into bigwig files
    allFiles=os.listdir(outputDir)
    sortedFiles=[element for element in allFiles if '.sorted.bdg' in element if name in element]
    sortedFiles=[name+'_control_lambda.sorted.bdg',name+'_treat_pileup.sorted.bdg',name+'_FE.sorted.bdg',name+'_logLR.sorted.bdg']
    conversionCommands=[]
    for element in sortedFiles:
        bwFile=element.replace('.sorted.bdg','.bw')
        flag0='time'
        flag1='%s%s'%(outputDir,element)
        flag2='%s%s'%(outputDir,bwFile)
        pieces=[flag0,bigWigExecutable,flag1,chromoSizeFile,flag2]
        instance=' '.join(pieces)
        conversionCommands.append(instance)
        
    # f.5. submitting jobs
    allCommands=[fullCmd1,fullCmd2,fullCmd3]

    for element in sortingCommands:
        allCommands.append(element)

    for element in conversionCommands:
        allCommands.append(element)
    
    SGEsubmitter(allCommands,name)

    return None

# 0. user defined variables
macs2Executable='/users/alomana/anaconda2/bin/macs2'
bamDir='/proj/omics4tb/alomana/projects/csp.jgi/data/bam/'
outputDir='/proj/omics4tb/alomana/projects/csp.jgi/data/macs2.run3/'
scratchDir='/proj/omics4tb/alomana/scratch/macs2/'
sendersDir=scratchDir+'senders/'
chromoSizeFile='/proj/omics4tb/alomana/projects/csp.jgi/data/genome/chrom.sizes.txt'
bigWigExecutable='/proj/omics4tb/alomana/software/bedGraphToBigWig'

numberOfThreads=1

# 1. selecting pair of files to run
experiments=[
    ['ASCAO','ASCAX','0hA',141], # FAIRE_0A, IC_0A
    ['ASCAP','ASCAY','0hB',182], # FAIRE_0B, IC_0B
    ['ASCAS','ASCAZ','24hA',167], # FAIRE_24A, IC_24A
    ['ASCAT','ASCAG','24hB',163], # FAIRE_24B, IC_24B
    ['ASCAU','ASCAH','48hA',162], # FAIRE_48A, IC_48A
    ['ASCAW','ASCAN','48hB',128], # FAIRE_48B, IC_48B
    ]

bamFiles=os.listdir(bamDir)

# 2. actual calling

# 2.1. calling MACS2.0 using their corresponding ICs
print '### STARTING CALLER A...'
callerA()

# 2.2. calling MACS2.0 using their sample ICs
print '### STARTING CALLER B...'
callerB()

# 2.3. calling MACS2.0 using pooled ICs
print '### STARTING CALLER C...'
callerC()

# 2.4. calling MACS2.0 using exhaustively all single ICs
print '### STARTING CALLER D...'
callerD()

