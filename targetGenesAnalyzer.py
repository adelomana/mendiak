### this script finds transcripts that behave binary and investigates their chromatin profile
import sys,numpy

# 0. user defined variables and paths
expressionDataFile='/Volumes/omics4tb/alomana/projects/csp.jgi/data/expression/matrix.19506.v.5.5.FPKM.txt'
lowerThreshold=1.
upperThreshold=2.
minDifference=5.

# 1. finding transcripts of binary expression
print('finding transcripts of binary expression...')

openingGenes={}
closingGenes={}

with open(expressionDataFile,'r') as f:
    next(f)
    for line in f:
        vector=line.split('\t')

        # 1.1. finding transcripts that open: from <1. for all points before 2m	4m	8m	12m	18m and >2. at 12h, 24h and 48h. Also a difference of >5. FPKM
        name=vector[0]
        earlyValues=numpy.array([float(element) for element in vector[1:6]])
        lateValues=numpy.array([float(element) for element in vector[-3:]])
        difference=abs(numpy.mean(lateValues)-numpy.mean(earlyValues))
        flag3=False

        # 1.1. finding transcripts that open: from <1. for all points before 2m	4m	8m	12m	18m and >2. at 12h, 24h and 48h
        flag1=((earlyValues < lowerThreshold).sum() == earlyValues.size).astype(numpy.int)
        flag2=((lateValues > upperThreshold).sum() == lateValues.size).astype(numpy.int)
        if flag1 == 1 and flag2 == 1 and difference > minDifference:
            openingGenes[name]=[earlyValues,lateValues,difference]
            
        # 1.2. finding transcripts that close: from > 2. for all points at 2m	4m	8m	12m	18m and <1. at 12h, 24h and 48h
        flag1=((earlyValues > upperThreshold).sum() == earlyValues.size).astype(numpy.int)
        flag2=((lateValues < lowerThreshold).sum() == lateValues.size).astype(numpy.int)
        if flag1 == 1 and flag2 == 1 and difference > minDifference:
            closingGenes[name]=[earlyValues,lateValues,difference]

print(str(len(openingGenes.keys()))+' opening genes detected.')
print(str(len(closingGenes.keys()))+' closing genes detected.')

# 2. sorting genes based on expression differences
print('sorting genes based on expression differences...')

sortedOpeningGenes=[y[1] for y in sorted([(openingGenes[x][2], x) for x in openingGenes.keys()],reverse=True)]
sortedClosingGenes=[y[1] for y in sorted([(closingGenes[x][2], x) for x in closingGenes.keys()],reverse=True)]

####### PRINTING
for opportunity in sortedOpeningGenes[:8]:
    print(opportunity+'\t'+str(openingGenes[opportunity][2])+'\t'+str(openingGenes[opportunity][0])+'\t'+str(openingGenes[opportunity][1]))
print()
for challenge in sortedClosingGenes[:8]:
    print(challenge+'\t'+str(closingGenes[challenge][2])+'\t'+str(closingGenes[challenge][0])+'\t'+str(closingGenes[challenge][1]))
####### END PRINTING

# 2. display chromatin profile for target genes
print(len(sortedOpeningGenes))
print(len(sortedClosingGenes))

# 4. statistical analysis of overlap of genes ON/OFF with OPEN/CLOSE

# 4.1. hypergeometric test

# 4.2. permutation analysis

