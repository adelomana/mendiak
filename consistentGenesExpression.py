###
### this script reads genes with 10 or 01 signature and looks for relationships with expression
###

import pickle,sys,numpy
import matplotlib,matplotlib.pyplot
import scipy,scipy.stats,scipy.interpolate
import multiprocessing,multiprocessing.pool

def backgroundDistributionFinder(tempo):

    '''
    this function computes the background distribution for genes without FSSs
    '''

    sampledGenes=numpy.random.choice(otherSideGenes,size=len(expression24h))
    v=[]
    for gene in sampledGenes:
        fpkm=expression[gene]
        v.append(numpy.log10(fpkm+1))
    v.sort()

    return v

def expressionReader():
    
    '''
    this function reads the expression in the form of a dictionary
    '''

    expression={}
    with open(expressionDataFile,'r') as f:
        next(f)
        for line in f:
            vector=line.split('\t')
            transcriptName=vector[0].split('.t')[0]
            fpkm=float(vector[-2])
            expression[transcriptName]=fpkm
            
    return expression

###
### MAIN
###

# 0. user defined variables
jarDir='/Volumes/omics4tb/alomana/projects/csp.jgi/results/jars/'
consistentGenesJar=jarDir+'consistentGenes.filtered.pickle'
expressionDataFile='/Volumes/omics4tb/alomana/projects/csp.jgi/data/expression/matrix.19506.v.5.5.FPKM.txt'
figuresDir='/Volumes/omics4tb/alomana/projects/csp.jgi/results/figures/'
numberOfThreads=4
iterations=100

# 1. reading data
print('reading data...')

# 1.1. reading signature genes
f=open(consistentGenesJar,'rb')
consistentGenes=pickle.load(f)
f.close()

# 1.2. reading expression data at 24 h
expression=expressionReader()

# 2. compute distribution of expression for signature genes
print('retrieving distribution of signature genes...')

genes24h=consistentGenes[('24hA', '24hB')]
expression24h=[]
for gene in genes24h:
    fpkm=expression[gene]
    expression24h.append(numpy.log10(fpkm+1))

print('working with {} consistent positive transcripts at 24 h'.format(len(genes24h)))

# 3. computing statistical difference
print('computing statistical analysis...')

# working with signature genes
probS,binEdges=numpy.histogram(expression24h,bins='auto',density=True)
cumsumSignature=numpy.cumsum(probS)/sum(probS)
halfBin=binEdges[1]/2
x=[element+halfBin for element in binEdges]
x.pop()

matplotlib.pyplot.plot(x,cumsumSignature,'or',mew=0,alpha=0.5)
model=scipy.interpolate.PchipInterpolator(x,cumsumSignature)
newTime=numpy.arange(min(x),max(x),(max(x)-min(x))/100)
inter=model(newTime)
matplotlib.pyplot.plot(newTime,inter,'-r',lw=2,label='genes with FSSs')

# working with background genes. Most time consuming block of the entire script
allGenes=list(expression.keys())
otherSideGenes=[gene for gene in allGenes if gene not in expression24h]

hydra=multiprocessing.pool.Pool(numberOfThreads)
tasks=[i for i in range(100)]
backgroundDist=hydra.map(backgroundDistributionFinder,tasks)

background=numpy.array(backgroundDist)
meanBackground=numpy.mean(background,axis=0)
ksF,pvalueF=scipy.stats.ks_2samp(expression24h,meanBackground)
print('against mean background distribution: ks={}, p-value={}'.format(ksF,pvalueF))

probB,tempo=numpy.histogram(meanBackground,bins=binEdges,density=True)
cumsumBackground=numpy.cumsum(probB)/sum(probB)

matplotlib.pyplot.plot(x,cumsumBackground,'ok',mew=0,alpha=0.5)
model=scipy.interpolate.PchipInterpolator(x,cumsumBackground)
newTime=numpy.arange(min(x),max(x),(max(x)-min(x))/100)
inter=model(newTime)
matplotlib.pyplot.plot(newTime,inter,'-k',lw=2,label='genes without FSSs')

matplotlib.pyplot.xlim([-0.1,4.1])
matplotlib.pyplot.ylim([-0.05,1.05])

matplotlib.pyplot.xlabel('expression, $e$ (log$_{10}$(FPKM+1))')
matplotlib.pyplot.ylabel('CDF (P($E \leq e$))')

matplotlib.pyplot.legend(loc=4)

figureFile=figuresDir+'cdf.24h.png'
matplotlib.pyplot.savefig(figureFile)

# 4. final statement
print('... all done.')
