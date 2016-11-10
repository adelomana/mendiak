###
### this script explore the distribution of patterns between expresion and peak locations
###


import os,sys,numpy,pickle
import matplotlib,matplotlib.pyplot
import multiprocessing, multiprocessing.pool

def consistentPeakFinder(task):

    '''
    this function iterates over the peaks of another sample to find consistency
    '''

    consistentPeaks=[]
    
    # recovering task arguments
    workingPeakName=task[0]
    labelA=task[1]
    labelB=task[2]

    # searching the other sample for consistency
    founds={}
    peakA=workingPeaks[labelA][workingPeakName] 
    for peakNameB in workingPeaks[labelB]: 
        peakB=workingPeaks[labelB][peakNameB] 

        # check if they are in the same contig
        if peakA[0] == peakB[0]:
            flag,overlap=isConsistent([peakA,peakB])
            if flag == True:
                pairA=labelA+'.'+workingPeakName
                pairB=labelB+'.'+peakNameB
                consistentPeak=[pairA,pairB]
                founds[overlap]=consistentPeak

    if founds != {}:
        overlaps=list(founds.keys())
        if len(overlaps) > 1:
            overlaps.sort()
        consistentPeaks=founds[overlaps[0]]

    return consistentPeaks

def generalConsistency(workingPeaks):

    '''
    this function checks the consistency of peaks over replicates
    '''

    hydra=multiprocessing.pool.Pool(numberOfThreads)
    
    # working with all comparisons
    consistentPeaks={}
    samples=list(workingPeaks.keys())
    samples.sort()

    for i in range(len(samples)):
        for j in range(len(samples)):
            if i <= j:
                comparison=[]
                theKey=(samples[i],samples[j])
                tasks=tasks=[[peakName,samples[i],samples[j]] for peakName in workingPeaks[samples[i]]]
                print('\tcomparing %s peaks between samples %s and %s...'%(len(tasks),theKey[0],theKey[1]))
                output=hydra.map(consistentPeakFinder,tasks)
                for element in output:
                    if element != []:
                        comparison.append(element[1])
                consistentPeaks[theKey]=comparison

    # creating variable for graphical representation
    M=[]
    for i in range(len(samples)):
        v=[]
        for j in range(len(samples)):
            localKey=(samples[i],samples[j])
            inverseKey=(samples[j],samples[i])
            if localKey in consistentPeaks.keys():
                workingKey=localKey
            else:
                workingKey=inverseKey
            v.append(len(consistentPeaks[workingKey]))
        M.append(v)
        
    return consistentPeaks,M

def genomeReader():

    '''
    this function returns the position of genes
    '''

    genePositions={} # genePositions[geneName]=[chr,start,stop]

    with open(gff3File,'r') as f:
        next(f)
        next(f)
        for line in f:
            vector=line.split('\t')
            if vector[2] == 'gene':
                geneName=vector[8].split('.v5.5')[0].split('ID=')[1]
                chromosome=vector[0]
                start=int(vector[3])
                stop=int(vector[4])
                genePositions[geneName]=[chromosome,start,stop]

    return genePositions

def getSignature10(peakNameA):

    '''
    this function returns the gene name for a peak with signature 10
    '''
    
    peakA=workingPeaks['0hA'][peakNameA]
    currentChromosome=peakA[0]

    # check that is not found at 24hA
    founds=[]
    peaksB=[workingPeaks['24hA'][peakName] for peakName in workingPeaks['24hA'] if workingPeaks['24hA'][peakName][0] == currentChromosome]
    for peakB in peaksB:
        flag,overlap=isConsistent([peakA,peakB])
        founds.append(flag)
    foundA=sum(founds)

    # check that is not found at 24hB
    founds=[]
    peaksB=[workingPeaks['24hB'][peakName] for peakName in workingPeaks['24hB'] if workingPeaks['24hB'][peakName][0] == currentChromosome]
    for peakB in peaksB:
        flag,overlap=isConsistent([peakA,peakB])
        founds.append(flag)
    foundB=sum(founds)

    if foundA == 0 and foundB == 0:
        geneName=peakLocator(peakA)
    else:
        geneName=None

    return geneName
    
def isConsistent(pair):

    '''
    this function compares 2 peaks: they are considered consistent if they share 50% of their average length.
    '''

    flag=False

    peakA=pair[0]
    peakB=pair[1]
    # compute overlap
    min1=peakA[1]
    max1=peakA[2]
    min2=peakB[1]
    max2=peakB[2]
    overlap=max(0, min(max1, max2) - max(min1, min2))
    if overlap > 0:
        # compute that overlap is at least 50% of average length
        averageLength=numpy.mean([peakA[3],peakB[3]])
        threshold=averageLength*0.5
        if overlap >= threshold:
            flag=True                
        
    return flag,overlap

def matrixGrapher(M,labels):

    '''
    this figure builds the correlation matrix
    '''

    figureName=figuresDir+'peakConsistency.png'

    matplotlib.pyplot.imshow(M,interpolation='none',cmap='viridis')
    cb=matplotlib.pyplot.colorbar(label='Consistent Peaks',orientation='vertical',fraction=0.01) 
    cb.ax.tick_params(labelsize=10)
    matplotlib.pyplot.grid(False)

    # setting the numbers
    x=0.
    y=0.
    deltax=1.
    deltay=1.
    for i in range(len(M)):
        for j in range(len(M)):
            stringValue=str(M[i][j])
            matplotlib.pyplot.text(x+deltax*i,y+deltay*j,stringValue,fontsize=10,color='white',horizontalalignment='center',verticalalignment='center')

    matplotlib.pyplot.xticks(range(len(labels)),labels,size=18,rotation=90)
    matplotlib.pyplot.yticks(range(len(labels)),labels,size=18)
       
    matplotlib.pyplot.tick_params(axis='x',which='both',bottom='off',top='off')
    matplotlib.pyplot.tick_params(axis='y',which='both',right='off',left='off')
    matplotlib.pyplot.axes().set_aspect('equal')
    matplotlib.pyplot.title('consistent peaks')
    matplotlib.pyplot.tight_layout(0.5)
    matplotlib.pyplot.savefig(figureName)

    matplotlib.pyplot.clf()

    return None

def peakLocator(peak):

    '''
    this function returns the gene closer to a peak:
    1) it searches what gene has closer midpoint to peak midpoint
    2) it returns a valid gene if it is within 20% of the length of the gene
    '''
    
    gene4Peak=None

    # 1. computing the gene with shortest distance
    distance=float('Inf')
    peakChromosome=peak[0]
    peakCenter=peak[1]+(peak[2]-peak[1])/2
    for geneName in genePositions.keys():
        workingChromosome=genePositions[geneName][0]
        if workingChromosome == peakChromosome:
            start=genePositions[geneName][1]
            stop=genePositions[geneName][2]
            interval=stop-start
            geneCenter=start+interval/2
            workingDistance=abs(geneCenter-peakCenter)
            if workingDistance < distance:
                distance=workingDistance
                gene4Peak=geneName
                
    # 2. check that the peak is within 20% of the length of the gene
    start=genePositions[gene4Peak][1]
    stop=genePositions[gene4Peak][2]
    interval=stop-start
    bottom=start-0.2*interval
    top=stop+0.2*interval

    if peak[1] > bottom and peak[2] < top:
        pass
    else:
        gene4Peak=None
    
    return gene4Peak

def peakReader():

    '''
    this function reads specific information from the peaks file generated by MACS2.0
    '''

    peaks={} # a dictionary with the following structure: peaks[name]=[chro,start,end,length,fe]
    peaksFile=peaksDir+peaksFileName
    with open(peaksFile,'r') as f:
        for line in f:
            vector=line.split('\t')
            if len(vector) == 10:
                if 'peak' in vector[9]:

                    # name
                    brokenName=vector[9].split('_')
                    name=brokenName[1]+'.'+brokenName[2].replace('\n','')
                    # chr
                    chro=vector[0]
                    # start
                    start=int(vector[1])
                    # end
                    end=int(vector[2])
                    # length
                    length=int(vector[3])
                    # fold-enrichment
                    fe=float(vector[7])

                    peaks[name]=[chro,start,end,length,fe]

    lastLetter=name[-1]
    if lastLetter.isdigit():
        numberOfPeaks=int(name.split('.')[1])
    else:
        numberOfPeaks=int(name.split('.')[1][:-1])
    numberOfSummits=len(peaks.keys())

    print('\t%s peaks found; %s summits found.'%(numberOfPeaks,numberOfSummits))

    return peaks

def peaksDistributionPlotter(peaks,flag):

    '''
    this function build 2D histograms of peaks fe and size
    '''

    x=[];y=[]
    for name in peaks.keys():
        fe=peaks[name][-1]
        size=peaks[name][-2]
        
        x.append(fe)
        y.append(numpy.log10(size))

    feRange=[1,5]
    sizeRange=[2,4.05]

    h,xedges,yedges,tempo=matplotlib.pyplot.hist2d(x,y,bins=100,range=[feRange,sizeRange])
    z=numpy.log10(h+1).T
    zm=numpy.ma.masked_where(z == 0,z)

    
    newViridis=matplotlib.cm.viridis
    newViridis.set_bad('white')
    matplotlib.pyplot.imshow(zm,extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]],cmap=newViridis,interpolation='nearest',origin='lower',aspect='auto',vmin=0,vmax=2)
    cb=matplotlib.pyplot.colorbar(label='log10 Peak Count',fraction=0.05)
    cb.ax.tick_params(labelsize=10)

    # highlightling area of best peaks
    matplotlib.pyplot.plot([2,5],[3,3],'-k',color='red')
    matplotlib.pyplot.plot([2,2],[2,3],'-k',color='red')
    
    matplotlib.pyplot.xlabel('Fold Enrichment')
    matplotlib.pyplot.ylabel('Peak Size (bp)')

    positions=numpy.log10(numpy.array([100,200,300,500,750,1000,2000,4000,8000]))
    names=['100','200','300','500','750','1,000','2,000','4,000','8,000']
    matplotlib.pyplot.yticks(positions,names)

    matplotlib.pyplot.title(flag)
    
    matplotlib.pyplot.savefig(figuresDir+'figure.%s.png'%flag)
    matplotlib.pyplot.clf()

    return None

def peaksFilter():

    '''
    this function removes any peak that is lower than 5-fold and extends for longer than 1 kb. 
    it also plots the distribution of filtred and non-filtred peaks.
    '''

    filteredPeaks={}

    Dfe=[]
    Dsize=[]
    Efe=[]
    Esize=[]

    for name in peaks.keys():
        fe=peaks[name][-1]
        size=peaks[name][-2]

        Dfe.append(fe)
        Dsize.append(size)
        
        if fe >= peakFEThreshold and size <= peakLengthThreshold:
            filteredPeaks[name]=peaks[name]
            
            Efe.append(fe)
            Esize.append(size)

    # printing the number of filtered peaks and summits
    allKeys=filteredPeaks.keys()
    numberOfSummits=len(allKeys)
    uniquePeaks=[]
    for element in allKeys:
        lastLetter=element[-1]
        if lastLetter.isdigit():
            value=int(element.split('.')[1])
        else:
            value=int(element.split('.')[1][:-1])
        uniquePeaks.append(value)
    uniquePeaks=list(set(uniquePeaks))
    numberOfPeaks=len(uniquePeaks)

    print('\t%s filtered peaks found; %s filtered summits found.'%(numberOfPeaks,numberOfSummits))    
    
    return filteredPeaks

##
## MAIN
##

# 0. user defined variables
peaksDir='/proj/omics4tb/alomana/projects/csp.jgi/data/macs2.test/'
figuresDir='/proj/omics4tb/alomana/scratch/'
jarFile='/proj/omics4tb/alomana/projects/csp.jgi/data/jars/jarFile_full.pckl'
gff3File='/proj/omics4tb/alomana/projects/csp.jgi/data/genome/Creinhardtii_281_v5.5.gene.gff3'

correspondance={}
correspondance['0hA']='ASCAO'
correspondance['0hB']='ASCAP'
correspondance['24hA']='ASCAS'
correspondance['24hB']='ASCAT'
correspondance['48hA']='ASCAU'
correspondance['48hB']='ASCAW'

peakFEThreshold=2
peakLengthThreshold=1000
genomeSize=111098438

numberOfThreads=16

# 0.1. reading gene locations
genePositions=genomeReader()

# 1. selecting the samples
print('selecting samples...')
allFiles=os.listdir(peaksDir)
peaksFileNames=[element for element in allFiles if '_peaks.xls' in element if 'callerE' in element]
peaksFileNames.sort()

# 2. filter peaks: at least 2-fold and no longer than 1 kb
print('filtering samples...')

rawPeaks={}
selectedPeaks={}
for peaksFileName in peaksFileNames:

    label=peaksFileName.split('_')[0].split('.')[1]
    print('filtering sample %s...'%label)
    # 2.1. reading peaks
    peaks=peakReader()
    rawPeaks[label]=peaks

    # 2.2. filtering peaks
    filteredPeaks=peaksFilter()
    selectedPeaks[label]=filteredPeaks

    # 2.3. plot the distribution of peaks before and after filtering peaks
    #flag=peaksFileName.split('_peaks')[0]
    #peaksDistributionPlotter(peaks,flag)

# 3. define all genes that have matching patterns.
print('finding consistency among peaks...')
workingPeaks=rawPeaks

# 3.1. finding consistent peaks
#consistentPeaks,M=generalConsistency(workingPeaks)
#f=open(jarFile,'wb')
#pickle.dump(M,f)
#f.close()
f=open(jarFile,'rb')
M=pickle.load(f)
f.close()

# 3.3. plotting consistency numbers
#sampleNames=list(workingPeaks.keys())
#sampleNames.sort()
#matrixGrapher(M,sampleNames)

# 4. locating genes with variable signature
print('locating genes with specific signatures...')
hydra=multiprocessing.pool.Pool(numberOfThreads)

# 4.1. 10 signature
genesSignature10=[]
peakNames=list(workingPeaks['0hA'].keys())
output=hydra.map(getSignature10,peakNames)
genesSignature10=list(set(output))
print(genesSignature10)
sys.exit()



sys.exit()
# 4.2. 01 signature
genesSignature10=[]

# define consistent peaks at 24 h
for peakNameA in workingPeaks['24hA']:
    peakA=workingPeaks['0hA'][peakNameA]
    currentChromosome=peakA[0]
    peaksB=[workingPeaks['24hA'][peakName] for peakName in workingPeaks['24hA'] if workingPeaks['24hA'][peakName][0] == currentChromosome]
    tasks=[[peakA,element] for element in peaksB]
    output=hydra.map(isConsistent,tasks)
# check that they do not exist at 0hA

# 4.3. distribution of fold changes
print('finding expression values...')

# 4.4. building figure
print('building figure...')


# 4. define significance
#print('assessing signficance of found patterns...')
