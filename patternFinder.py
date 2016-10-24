import os,sys,numpy
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
    peakA=rawPeaks[labelA][workingPeakName] #! rawPeaks
    found=[]
    for peakNameB in rawPeaks[labelB]: #! rawPeaks
        peakB=rawPeaks[labelB][peakNameB] #! rawPeaks

        # check if they are in the same contig
        if peakA[0] == peakB[0]:
            flag=isConsistent(peakA,peakB)
            if flag == True:
                pairA=labelA+'.'+workingPeakName
                pairB=labelB+'.'+peakNameB
                consistentPeak=[pairA,pairB]
                consistentPeaks.append(consistentPeak)
                
    return consistentPeaks

def generalConsistency():

    '''
    this function checks the consistency of peaks over replicates
    '''

    hydra=multiprocessing.pool.Pool(numberOfThreads)
    
    # working with all comparisons
    consistentPeaks={}
    samples=list(rawPeaks.keys()) #! selectedPeaks
    samples.sort()

    for i in range(len(samples)):
        for j in range(len(samples)):
            if i <= j:
                comparison=[]
                theKey=(samples[i],samples[j])
                tasks=tasks=[[peakName,samples[i],samples[j]] for peakName in rawPeaks[samples[i]]] # selectedPeaks
                output=hydra.map(consistentPeakFinder,tasks)
                for element in output:
                    if element != []:
                        for subelement in element:
                            comparison.append(subelement)
                consistentPeaks[theKey]=comparison

    # making a graphical representation
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
            v.append(consistentPeaks[workingKey])
        M.append(v)
        
    matrixGrapher(M,samples)

    return consistentPeaks
    
def isConsistent(peakA,peakB):

    '''
    this function compares 2 peaks: they are considered consistent if they share 50% of their average length.
    '''

    flag=False

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
        
    return flag

def matrixGrapher(M,labels):

    '''
    this figure builds the correlation matrix
    '''

    figureName=figuresDir+'peakConsistency.png'

    print(M)

    matplotlib.pyplot.imshow(M,interpolation='none',cmap='viridis')
    cb=matplotlib.pyplot.colorbar(label='Consistent Peaks',orientation='vertical',fraction=0.01) 
    cb.ax.tick_params(labelsize=10)
    matplotlib.pyplot.grid(False)

    # setting the numbers
    x=-0.3
    y=0.1
    deltax=1.
    deltay=1.
    for i in range(len(M)):
        for j in range(len(M)):
            if i != j:
                value=str(M[i][j])
                matplotlib.pyplot.text(x+deltax*i,y+deltay*j,value,fontsize=10,color='white')

    matplotlib.pyplot.xticks(range(len(labels)),labels,size=18,rotation=90)
    matplotlib.pyplot.yticks(range(len(labels)),labels,size=18)
    
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.tick_params(axis='x',which='both',bottom='off',top='off')
    matplotlib.pyplot.tick_params(axis='y',which='both',right='off',left='off')
    matplotlib.pyplot.axes().set_aspect('equal')
    matplotlib.pyplot.savefig(figureName)

    matplotlib.pyplot.clf()

    return None

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

# 0. user defined variables
#peaksDir='/Volumes/omics4tb/alomana/projects/csp.jgi/data/macs2.run3/'
#peaksDir='/Volumes/omics4tb/alomana/projects/csp.jgi/data/macs2.test/'
peaksDir='/Users/adriandelomana/tempo/macs2.test/'
#figuresDir='/Users/alomana/gDrive2/tmp/'
figuresDir='/Users/adriandelomana/gDrive/tmp/'


correspondance={}
correspondance['0hA']='ASCAO'
correspondance['0hB']='ASCAP'
correspondance['24hA']='ASCAS'
correspondance['24hB']='ASCAT'
correspondance['48hA']='ASCAU'
correspondance['48hB']='ASCAW'

peakFEThreshold=2
peakLengthThreshold=1000

numberOfThreads=4

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
print('finding matching pattern peaks...')

# 3.1. finding consistent peaks
consistentPeaks=generalConsistency()

for key in consistentPeaks.keys():
    print(key,len(consistentPeaks[key]))


sys.exit()

# 3.1. 101 signature


# 3.1. 010 signature

# checking consistency at 24 h samples
#@ consider making first loop parallel, also, make a test folder
#@ make this a function, check for all samples
consistentPeaks=[]


# check the absence of peak for samples 0 h and 48 h make this function parallel too
#for peak in consistentPeaks:
#    print()
#    occupiedArea=None



sys.exit()

# checking absence at 0 h and 48 h. 

# 4. define significance
#print('assessing signficance of found patterns...')
