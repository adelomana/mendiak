'''
this script creates the FASTQ files for the saturation curve
'''

import os,sys,random,numpy
import scipy,scipy.interpolate
import matplotlib,matplotlib.pyplot

def butcher(sample,iteration):

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

        print('selecting %s mr...'%mr)

        # select n random reads
        k=int(mr*1e6) 
        selectedReads=random.sample(allKeys,k)
       
        # write a file
        outputFile=piecesDir+sample+'.resolution.%s.iteration.%s.fastq'%(str(mr),str(iteration))
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
readsFilesDir='/Volumes/omics4tb/alomana/projects/csp.jgi/data/fastq/clean/'
piecesDir='/Volumes/omics4tb/alomana/projects/csp.jgi/data/fastq.crumbles.n5/'
peaksDir='/Volumes/omics4tb/alomana/projects/csp.jgi/data/macs2.crumbles.n5/'
figuresDir='/Users/alomana/gDrive2/tmp/'
resolution=2
iterations=5

correspondance={}
correspondance['0hA']='ASCAO'
correspondance['0hB']='ASCAP'
correspondance['24hA']='ASCAS'
correspondance['24hB']='ASCAT'
correspondance['48hA']='ASCAU'
correspondance['48hB']='ASCAW'

finalSamples={}
finalSamples['0hA']=[33.25,30028]
finalSamples['0hB']=[10.67,9134]
finalSamples['24hA']=[21.15,35467]
finalSamples['24hB']=[15.2, 32323]
finalSamples['48hA']=[19.24,13712]
finalSamples['48hB']=[33.73,18778]

'''
# 1. iterate over samples
for sample in correspondance.keys():    
    # 2. create FASTQ crumbles
    print('working with sample'+sample+'...')
    for iteration in range(iterations):
        iteration=iteration
        print('working with iteration %s...'%str(iteration))
        butcher(sample,iteration)
'''

# 3. figure maker
allFiles=os.listdir(peaksDir)
peakFiles=[element for element in allFiles if "_peaks.xls" in element]

# 3.1. reading data
print('reading peak data...')
saturationValues={}
for peakFileName in peakFiles:

    pieces=peakFileName.split('_peaks')[0].split('.')
    sampleName=pieces[1]
    mr=int(pieces[3])
    if sampleName not in saturationValues.keys():
        saturationValues[sampleName]={}
    if mr not in saturationValues[sampleName].keys():
        saturationValues[sampleName][mr]=[]
    
    peakFile=peaksDir+peakFileName
    with open(peakFile,'r') as f:
        lastLine=f.readlines()[-1]
    vector=lastLine.split('\t')
    peaksTextForm=vector[-1].split('_')[-1]
    peaksTextForm=peaksTextForm.replace('\n','')

    if peaksTextForm[0] == '#':
        formatted=0
    elif peaksTextForm[-1].isdigit() == False:
        formatted=int(peaksTextForm[:-1])
    else:
        formatted=int(peaksTextForm)
    
    saturationValues[sampleName][mr].append(formatted)

# 3.2. building graph
print('creating figure...')

colorSchema={}
colorSchema['0hA']='darkgreen'
colorSchema['0hB']='green'
colorSchema['24hA']='darkred'
colorSchema['24hB']='red'
colorSchema['48hA']='blue'
colorSchema['48hB']='darkblue'

allSamples=list(saturationValues.keys())
allSamples.sort()

for sample in allSamples:
    
    x=list(saturationValues[sample].keys())
    x.sort()
    y=[];err=[]
    for element in x:
        y.append(numpy.mean(saturationValues[sample][element]))
        err.append(numpy.std(saturationValues[sample][element]))

    # setting up the color
    theColor=colorSchema[sample]

    # plotting the rest of the trajectory
    matplotlib.pyplot.errorbar(x,y,yerr=err,fmt='o',color=theColor,alpha=0.5,mew=0.)

    # adding final sample
    finalx=finalSamples[sample][0]
    finaly=finalSamples[sample][1]
    matplotlib.pyplot.plot(finalx,finaly,'s',color=theColor,mew=0.)
    
    # interpolating lines
    interpolation=scipy.interpolate.PchipInterpolator(x,y)
    xnew=numpy.linspace(min(x),max(x),num=500,endpoint=True)
    ynew=interpolation(xnew)
    matplotlib.pyplot.plot(xnew,ynew,'-',lw=2,color=theColor,label=sample)
    

matplotlib.pyplot.xlim([0,35])
matplotlib.pyplot.ylim([-1000,40000])

matplotlib.pyplot.xlabel('sequencing depth (mr)',fontsize=28)
matplotlib.pyplot.ylabel('peaks (x1e3)',fontsize=28)

matplotlib.pyplot.legend(loc=2)

matplotlib.pyplot.yticks([5000,10000,15000,20000,25000,30000,35000,40000],['5','10','15','20','25','30','35','40'],fontsize=20)
matplotlib.pyplot.xticks(fontsize=20)

matplotlib.pyplot.tight_layout()

matplotlib.pyplot.savefig(figuresDir+'saturationCurve.png')
