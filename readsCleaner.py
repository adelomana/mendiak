###
### this script calls Trimmomatic to clean raw reads
###

import os,sys

def trimmomaticCaller(instance):
    
    '''
    this function deals with the trimming of the reads using Trimmomatic. Recommended options, ILLUMINACLIP:path2AdaptersFile:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    '''

    print('working with file'+instance)

    logFile=logFilesDir+instance+'.messagesFromTrimming.txt'

    inputFile=rawReadsDir+instance+'.fastq'

    outputFile=cleanReadsDir+instance+'.trimmed.fastq'

    cmd=javaPath+' -jar '+trimmomaticPath+' SE -threads 4 -phred33 -trimlog %s %s %s ILLUMINACLIP:%s:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'%(logFile,inputFile,outputFile,path2Adapter)  

    print(cmd)
    os.system(cmd)
    print()
    
    return None

### MAIN

# 0. defining user variables
rawReadsDir='/Users/alomana/scratch/csp/raw/'
cleanReadsDir='/Users/alomana/scratch/csp/clean/'
logFilesDir='/Users/alomana/scratch/csp/'
path2Adapter='/Users/alomana/software/Trimmomatic-0.36/adapters/TruSeq3-SE.fa'

javaPath='/Library/Internet\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java'
trimmomaticPath='/Users/alomana/software/Trimmomatic-0.36/trimmomatic-0.36.jar'

# 1. reading the files
tag='ASCA'
fastqFiles=os.listdir(rawReadsDir)
readFiles=[element for element in fastqFiles if tag in element]

#readFiles=['test.fastq']

# 2. calling Trimmomatic
for readFile in readFiles:
    instance=readFile.split('.fastq')[0]
    trimmomaticCaller(instance)
