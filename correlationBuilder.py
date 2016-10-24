# this script builds correlation among samples

import os,sys,matplotlib
import matplotlib.pyplot

def correlationsComputer(inputFiles,resultsFile):

    '''
    this function computes the correlation between bw files using wigCorrelate from http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigCorrelate
    '''

    fullPaths=[bigWigDir+element for element in inputFiles]
    
    command='wigCorrelate '+' '.join(fullPaths)+' > '+resultsFile

    #print command
    #os.system(command)
    #print

    correlations={}
    with open(resultsFile,'r') as f:
        for line in f:
            vector=line.split()
            sample1=vector[0].split('/')[-1].split('.bw')[0]
            sample2=vector[1].split('/')[-1].split('.bw')[0]
            value=float(vector[2])
            fullLabel=sample1+'.'+sample2
            correlations[fullLabel]=value

    return correlations

def matrixGrapher(M,labels,resultsFile):

    '''
    this figure builds the correlation matrix
    '''

    figureName='figures/correlation.%s.png'%resultsFile.split('.txt')[0]

    matplotlib.pyplot.imshow(M,interpolation='none',cmap='viridis',vmin=0.,vmax=1.)
    cb=matplotlib.pyplot.colorbar(label='Pearson correlation coefficient (PCC)',orientation='vertical',fraction=0.01) # need to improve label font size
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
                value='%.2f'%M[i][j]
                shortValue=value.replace('0.','.')
                matplotlib.pyplot.text(x+deltax*i,y+deltay*j,shortValue,fontsize=10,color='white')

    matplotlib.pyplot.xticks(range(len(labels)),labels,size=18,rotation=90)
    matplotlib.pyplot.yticks(range(len(labels)),labels,size=18)
    
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.tick_params(axis='x',which='both',bottom='off',top='off')
    matplotlib.pyplot.tick_params(axis='y',which='both',right='off',left='off')
    matplotlib.pyplot.axes().set_aspect('equal')
    matplotlib.pyplot.savefig(figureName)

    matplotlib.pyplot.clf()

    return None
    
def metadataReader():

    '''
    this function builds the metadata
    '''

    metadata={}

    inputFile='../data/metadata.txt'
    with open(inputFile,'r') as f:
        for line in f:
            vector=line.split()
            sampleName=vector[5]
            sampleCode=vector[6]
            metadata[sampleCode]=sampleName
            
    return metadata

# 0. user defined variables
bigWigDir='../data/faireData/'

# 1. reading metadata
metadata=metadataReader()

## # 2. building correlations among FAIRE and IC
resultsFile='faireAndIC.txt'
someFiles=os.listdir(bigWigDir)
workingFiles=[element for element in someFiles if '.bw' in element and 'fe_' not in element and 'diff_' not in element]

# 2.1. computing the correlations
correlations=correlationsComputer(workingFiles,resultsFile)

# 2.2. building the figure
orderedListOfIDs=sorted(metadata,key=metadata.__getitem__)
correlationLabels=correlations.keys()
M=[]
for element1 in orderedListOfIDs:
    v=[]
    for element2 in orderedListOfIDs:
        if element1 == element2:
            value=1.
        else:
            matches=[label for label in correlationLabels if element1 in label and element2 in label]
            if len(matches) != 1:
                print matches
                print 'error a'
                sys.exit()
            value=correlations[matches[0]]
            
        v.append(value)
    M.append(v)

figureLabels=[]
for label in orderedListOfIDs:
    figureLabels.append(metadata[label])
matrixGrapher(M,figureLabels,resultsFile)
        
# 3. building correlations among fold-enrichment signal samples
resultsFile='fe.txt'
someFiles=os.listdir(bigWigDir)
workingFiles=[element for element in someFiles if '.bw' in element and 'fe_' in element]

# 3.1. computing the correlations
correlations=correlationsComputer(workingFiles,resultsFile)

# 3.2. building the figure
fullListOfIDs=[]
names=correlations.keys()
for name in names:
    splitSamples=name.split('.')
    for element in splitSamples:
        if element not in fullListOfIDs:
            fullListOfIDs.append(element)
fullListOfIDs.sort()

M=[]
for element1 in fullListOfIDs:
    v=[]
    for element2 in fullListOfIDs:
        if element1 == element2:
            value=1.
        else:
            matches=[element for element in names if element1 in element and element2 in element]
            if len(matches) != 1:
                print 'error selecting the samples, exiting...'
                sys.exit()
            value=correlations[matches[0]]
        v.append(value)
    M.append(v)

figureLabels=[]
for element in fullListOfIDs:
    elementa=element.split('_')[0]
    elementb=element.split('-')[1]
    elementc=element.split('-')[2][2]
    label=elementa+'.'+elementb+'.'+elementc
    figureLabels.append(label)

matrixGrapher(M,figureLabels,resultsFile)
