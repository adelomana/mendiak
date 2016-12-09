###
### this script explore the distribution of patterns between expresion and peak locations
###

import os,sys,numpy,pickle,copy
import multiprocessing, multiprocessing.pool

import matplotlib
matplotlib.use('Agg') # necessary for saving figures at remove machine
import matplotlib.pyplot

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
    founds=[]
    peakA=rawPeaks[labelA][workingPeakName] 
    for peakNameB in rawPeaks[labelB]: 
        peakB=rawPeaks[labelB][peakNameB] 

        # check if they are in the same contig
        if peakA[0] == peakB[0]:
            flag,overlap=isConsistent([peakA,peakB])
            if flag == True:
                pairA=labelA+'.'+workingPeakName
                pairB=labelB+'.'+peakNameB
                consistentPair=[pairA,pairB]
                founds.append([consistentPair,overlap])
                
    # dealing with multiple hits
    if founds != []:
        if len(founds) == 1:
            consistentPeaks=founds[0][0]
        else:
            overlaps=[element[1] for element in founds]
            sortedOverlaps=copy.deepcopy(overlaps)
            sortedOverlaps.sort(reverse=True)
            if sortedOverlaps[0] != sortedOverlaps[1]:
                theIndex=overlaps.index(sortedOverlaps[0])
                consistentPeaks=founds[theIndex][0]
            else:
                feRanks=[]
                putativePeaks=[element[0][1] for element in founds if element[1] == sortedOverlaps[0]]
                feRanks=[rawPeaks[labelB][element.split(labelB+'.')[1]][-1] for element in putativePeaks]
                sortedFEranks=copy.deepcopy(feRanks)
                sortedFEranks.sort(reverse=True)
                selected=putativePeaks[feRanks.index(sortedFEranks[0])]
                consistentPeaks=[founds[0][0][0],selected]
                
    return consistentPeaks

def generalConsistency():

    '''
    this function checks the consistency of peaks over replicates
    '''

    hydra=multiprocessing.pool.Pool(numberOfThreads)
    
    consistentAllPeaks={}
    consistentFilteredPeaks={}
    
    samples=list(rawPeaks.keys())
    samples.sort()

    for i in range(len(samples)):
        filteredPeakNamesA=list(selectedPeaks[samples[i]].keys())
        for j in range(len(samples)):
            if i < j:
                filteredPeakNamesB=list(selectedPeaks[samples[j]].keys())
                comparisonA=[]
                comparisonF=[]
                theKey=(samples[i],samples[j])

                workNames=list(rawPeaks[samples[i]].keys())
                tasks=[[peakName,samples[i],samples[j]] for peakName in workNames]
                print('\t comparing %s peaks between samples %s and %s...'%(len(tasks),theKey[0],theKey[1]))

                output=hydra.map(consistentPeakFinder,tasks)                
                for element in output:
                    if element != []:
                        comparisonA.append(element)

                        peakNameA=element[0].split(theKey[0]+'.')[1]
                        peakNameB=element[1].split(theKey[1]+'.')[1]
                        if peakNameA in filteredPeakNamesA and peakNameB in filteredPeakNamesB:
                            comparisonF.append(element)
                        
                consistentAllPeaks[theKey]=comparisonA
                consistentFilteredPeaks[theKey]=comparisonF
                print('\t found {} all and {} filtered consistent peaks.'.format(len(comparisonA),len(comparisonF)))
                print()
                
    # creating variable for graphical representation
    MA=[]; MF=[]
    for i in range(len(samples)):
        va=[]; vf=[]
        for j in range(len(samples)):
            if i == j:
                valueA=len(rawPeaks[samples[i]])
                valueF=len(selectedPeaks[samples[i]])
            else:
                localKey=(samples[i],samples[j])
                inverseKey=(samples[j],samples[i])
                if localKey in consistentAllPeaks.keys():
                    workingKey=localKey
                else:
                    workingKey=inverseKey
                valueA=len(consistentAllPeaks[workingKey])
                valueF=len(consistentFilteredPeaks[workingKey])
            va.append(valueA)
            vf.append(valueF)
        MA.append(va)
        MF.append(vf)

    return consistentAllPeaks,consistentFilteredPeaks,MA,MF

def genomeOccupancyCalculator(assessingPeaks,label):

    '''
    this function computes the size of the genome peaks occupy
    '''

    sumLength=0
    allLabels=[]
    for peakName in assessingPeaks.keys():
        currentLabel=assessingPeaks[peakName][0]+'.'+str(assessingPeaks[peakName][1])+'.'+str(assessingPeaks[peakName][2])
        if currentLabel not in allLabels:
            allLabels.append(currentLabel)
            peakSize=assessingPeaks[peakName][3]
            sumLength=sumLength+peakSize
    percentage=sumLength/genomeSize
    print('\tsum of',label,'peak lengths is',sumLength,'relative size',percentage)

    return None

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

def getSignature10(task):

    '''
    this function returns the gene name for a peak with signature 10
    '''
    
    peakPair=task[0]
    flatCondition=task[1]
    extension=task[2]

    if extension == 'all':
        workingPeaks=rawPeaks
    elif extension == 'filtered':
        workingPeaks=selectedPeaks
    else:
        print('error from getSignature10. exiting...')
        sys.exit()

    currentSampleLabel=peakPair[0].split('.')[0]
    currentPeakName=peakPair[0].split(currentSampleLabel+'.')[1]
    currentChromosome=rawPeaks[currentSampleLabel][currentPeakName][0]
    
    # check that none of the two peaks of the pair is at the flatten condition
    founds=[]
    flattenPeaksA=[workingPeaks[flatCondition[0]][peakName] for peakName in workingPeaks[flatCondition[0]] if workingPeaks[flatCondition[0]][peakName][0] == currentChromosome]
    flattenPeaksB=[workingPeaks[flatCondition[1]][peakName] for peakName in workingPeaks[flatCondition[1]] if workingPeaks[flatCondition[1]][peakName][0] == currentChromosome]
    flattenPeaks=flattenPeaksA+flattenPeaksB
    #print('# of flatten peaks to probe',len(flattenPeaks),'chr',currentChromosome)
    for single in peakPair:
        peakA=rawPeaks[currentSampleLabel][currentPeakName]
        for peakB in flattenPeaks:
            flag,overlap=isConsistent([peakA,peakB])
            founds.append(flag)
    #print('# of answered flags',len(founds))
    #print('# of founds in flatten',sum(founds))
    if sum(founds) == 0:
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

def matrixGrapher(M,labels,figureTitle):

    '''
    this figure builds the correlation matrix
    '''

    figureName=figuresDir+'peakConsistency.{}.{}.png'.format(figureTitle.split(' ')[0],figureTitle.split(' ')[1])

    matplotlib.pyplot.imshow(M,interpolation='none',cmap='viridis')
    cb=matplotlib.pyplot.colorbar(orientation='vertical',fraction=0.05) 
    cb.set_label(label='Consistent Peaks',size=20)
    cb.ax.tick_params(labelsize=16)
    matplotlib.pyplot.grid(False)

    # setting the numbers
    x=0.
    y=0.
    deltax=1.
    deltay=1.
    for i in range(len(M)):
        for j in range(len(M)):
            stringValue=str(M[i][j])
            matplotlib.pyplot.text(x+deltax*i,y+deltay*j,stringValue,fontsize=16,color='white',horizontalalignment='center',verticalalignment='center',fontweight='bold')

    matplotlib.pyplot.xticks(range(len(labels)),labels,size=20,rotation=90)
    matplotlib.pyplot.yticks(range(len(labels)),labels,size=20)
       
    matplotlib.pyplot.tick_params(axis='x',which='both',bottom='off',top='off')
    matplotlib.pyplot.tick_params(axis='y',which='both',right='off',left='off')
    matplotlib.pyplot.axes().set_aspect('equal')
    matplotlib.pyplot.title(figureTitle,fontsize=24)
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
                
    # 2. check that the peak is within 20% of the length of the gene #### this function needs to be improved to avoid overlaps and search for at least 0.5 kb if gene is very small
    #! consider doing 33%. 
    if gene4Peak != None:
        start=genePositions[gene4Peak][1]
        stop=genePositions[gene4Peak][2]
        interval=stop-start
        bottom=start-int(0.2*interval)
        top=stop+int(0.2*interval)

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
    cb=matplotlib.pyplot.colorbar(fraction=0.05)
    cb.set_label(label='log$_{10}$ Peak Count',size=20)
    cb.ax.tick_params(labelsize=18)

    # highlightling area of best peaks
    matplotlib.pyplot.plot([2,5],[3,3],'-k',color='red',lw=2)
    matplotlib.pyplot.plot([2,2],[2,3],'-k',color='red',lw=2)
    
    matplotlib.pyplot.xlabel('Fold Enrichment',fontsize=28)
    matplotlib.pyplot.ylabel('Site Length (bp)',fontsize=28)

    positions=numpy.log10(numpy.array([100,200,300,500,750,1000,2000,4000,8000]))
    names=['100','200','300','500','750','1,000','2,000','4,000','8,000']
    matplotlib.pyplot.yticks(positions,names,fontsize=20)
    matplotlib.pyplot.xticks(fontsize=20)

    theTitle='sample '+flag.split('.')[1]
    matplotlib.pyplot.title(theTitle,fontsize=36)
    
    matplotlib.pyplot.tight_layout()
    
    matplotlib.pyplot.savefig(figuresDir+'figure.%s.png'%flag)
    matplotlib.pyplot.clf()

    return None

def peaksFilter():

    '''
    this function removes any peak that is lower than 2-fold and extends for longer than 1 kb. 
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
#peaksDir='/proj/omics4tb/alomana/projects/csp.jgi/data/macs2.test/'
peaksDir='/proj/omics4tb/alomana/projects/csp.jgi/data/macs2.run3/'
jarDir='/proj/omics4tb/alomana/projects/csp.jgi/results/jars.wk/'
gff3File='/proj/omics4tb/alomana/projects/csp.jgi/data/genome/Creinhardtii_281_v5.5.gene.gff3'
figuresDir='/proj/omics4tb/alomana/projects/csp.jgi/results/figures/'

numberOfThreads=40

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

# 0.1. reading gene locations
genePositions=genomeReader()

# 1. selecting the samples
print('selecting samples...')
allFiles=os.listdir(peaksDir)
peaksFileNames=[element for element in allFiles if '_peaks.xls' in element if 'callerC' in element]
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
    filteredPeaks=peaksFilter()
    selectedPeaks[label]=filteredPeaks

    # 2.2. plot the distribution of peaks
    flag=peaksFileName.split('_peaks')[0]
    peaksDistributionPlotter(peaks,flag)

    # 2.3. computing the size of genome that peaks occupy
    genomeOccupancyCalculator(peaks,'all')
    genomeOccupancyCalculator(filteredPeaks,'filtered')

# 3. define all genes that have matching patterns
print()
print('finding consistency among peaks...')

# 3.1. finding consistent peaks among all peaks
consistentAllPeaks,consistentFilteredPeaks,MA,MF=generalConsistency()

"""
# 3.2. saving consistent peaks and reading consistent peaks
jarFile=jarDir+'consistentPeaks.pickle'
f=open(jarFile,'wb')
pickle.dump([consistentAllPeaks,consistentFilteredPeaks,MA,MF],f)
f.close()
"""

# 3.3. plotting heat map of consistency counts for all and filtered peaks
sampleNames=list(rawPeaks.keys())
sampleNames.sort()
matrixGrapher(MA,sampleNames,'full set')
matrixGrapher(MF,sampleNames,'filtered set')

# 4. locating genes with variable signature
print()
print('defining genes with specific signatures...')
hydra=multiprocessing.pool.Pool(numberOfThreads)

# 4.1. 10 signature
print('\t working with 10 signature...')

genesAllSignature10=[]
genesFilteredSignature10=[]

positiveCondition=('0hA','0hB')
flatCondition=('24hA','24hB')

tasks=[[element,flatCondition,'all'] for element in consistentAllPeaks[positiveCondition]]
print('\t interrogating',len(tasks),'consistent pairs of peaks...')
output=hydra.map(getSignature10,tasks)
genesAllSignature10=list(set(output))
genesAllSignature10.remove(None)

tasks=[[element,flatCondition,'filtered'] for element in consistentFilteredPeaks[positiveCondition]]
print('\t interrogating',len(tasks),'consistent and filtered pairs of peaks...')
output=hydra.map(getSignature10,tasks)
genesFilteredSignature10=list(set(output))
genesFilteredSignature10.remove(None)

print('%s genes found with broad 10 signature.'%len(genesAllSignature10))
print('%s genes found with filtered 10 signature.'%len(genesFilteredSignature10))

jarFile=jarDir+'signature10.pickle'
f=open(jarFile,'wb')
pickle.dump([genesAllSignature10,genesFilteredSignature10],f)
f.close()

# 4.2. 01 signature
print('\t working with 01 signature...')

genesAllSignature01=[]
genesFilteredSignature01=[]

flatCondition=('0hA','0hB')
positiveCondition=('24hA','24hB')

tasks=[[element,flatCondition,'all'] for element in consistentAllPeaks[positiveCondition]]
print('\t interrogating',len(tasks),'consistent pairs of peaks...')
output=hydra.map(getSignature10,tasks)
genesAllSignature01=list(set(output))
genesAllSignature01.remove(None)

tasks=[[element,flatCondition,'filtered'] for element in consistentFilteredPeaks[positiveCondition]]
print('\t interrogating',len(tasks),'consistent and filtered pairs of peaks...')
output=hydra.map(getSignature10,tasks)
genesFilteredSignature01=list(set(output))
genesFilteredSignature01.remove(None)

print('%s genes found with broad 01 signature.'%len(genesAllSignature01))
print('%s genes found with filtered 01 signature.'%len(genesFilteredSignature01))

jarFile=jarDir+'signature01.pickle'
f=open(jarFile,'wb')
pickle.dump([genesAllSignature01,genesFilteredSignature01],f)
f.close()

# 5. final message
print()
print('... all done.')
