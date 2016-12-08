###
### this script reads genes with 10 or 01 signature and looks for relationships with expression
###

import pickle,sys,numpy
import matplotlib,matplotlib.pyplot

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
signature01Pickle=jarDir+'signature01.pickle'
expressionDataFile='/Volumes/omics4tb/alomana/projects/csp.jgi/data/expression/matrix.19506.v.5.5.FPKM.txt'

# 1. reading data

# 1.1. reading signature genes
f=open(signature01Pickle,'rb')
[genesAllSignature01,genesFilteredSignature01]=pickle.load(f)
f.close()

# 1.2. reading expression data at 24 h
expression=expressionReader()

# 2. compute distribution of expression for signature genes
signatureExpression=[]
for gene in genesAllSignature01:
    fpkm=expression[gene]
    signatureExpression.append(numpy.log10(fpkm+1))

# 3. compute distribution of n samples of random genes

# 4. plot figure

# working with signature genes
prob,binEdges=numpy.histogram(signatureExpression,bins='auto',density=True)
halfBin=binEdges[1]/2
x=[element+halfBin for element in binEdges]
x.pop()

print(len(signatureExpression))
print(prob,prob.shape)
print(binEdges,binEdges.shape)
print(x)

matplotlib.pyplot.plot(x,prob,'ok')

matplotlib.pyplot.savefig('thefigure.pdf')
