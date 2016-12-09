###
### this script reads genes with 10 or 01 signature and looks for relationships with expression
###

import pickle,sys,numpy
import matplotlib,matplotlib.pyplot
import scipy,scipy.stats

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
jarDir='/Volumes/omics4tb/alomana/projects/csp.jgi/results/jars.wk/'
consistentGenesJar=jarDir+'consistentGenes.filtered.pickle'
expressionDataFile='/Volumes/omics4tb/alomana/projects/csp.jgi/data/expression/matrix.19506.v.5.5.FPKM.txt'
iterations=100

# 1. reading data

# 1.1. reading signature genes
f=open(consistentGenesJar,'rb')
consistentGenes=pickle.load(f)
f.close()

for key in consistentGenes.keys():
    print(key)
    print(len(consistentGenes[key]))
    print(consistentGenes[key][:10])
    print()
sys.exit()

# 1.2. reading expression data at 24 h
expression=expressionReader()

# 2. compute distribution of expression for signature genes
signatureExpression=[]
for gene in genesAllSignature01:
    fpkm=expression[gene]
    signatureExpression.append(numpy.log10(fpkm+1))

# 4. plot figure
print('computing statistical analysis...')

# working with signature genes
probS,binEdges=numpy.histogram(signatureExpression,bins='auto',density=True)
cumsumSignature=numpy.cumsum(probS)/sum(probS)
halfBin=binEdges[1]/2
x=[element+halfBin for element in binEdges]
x.pop()
matplotlib.pyplot.plot(x,cumsumSignature,'or',mew=0)

# working with background genes
allGenes=list(expression.keys())
otherSideGenes=[gene for gene in allGenes if gene not in genesAllSignature01]

backgroundHistogram=[]

for i in range(iterations):
    sampledGenes=numpy.random.choice(otherSideGenes,size=len(genesAllSignature01))
    v=[]
    for gene in sampledGenes:
        fpkm=expression[gene]
        v.append(numpy.log10(fpkm+1))
    prob,binEdges=numpy.histogram(v,bins=binEdges,density=True)
    backgroundHistogram.append(prob)

    print(scipy.stats.ks_2samp(signatureExpression,v))

print()

sys.exit()

allProb=numpy.array(backgroundHistogram)
meanProb=numpy.mean(allProb,axis=0)

cumsumBackground=numpy.cumsum(meanProb)/sum(meanProb)
matplotlib.pyplot.plot(x,cumsumBackground,'ok',mew=0)

# statistical difference
print(scipy.stats.ks_2samp(probS,meanProb))

matplotlib.pyplot.savefig('thefigure.pdf')
