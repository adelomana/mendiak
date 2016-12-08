###
### this script reads genes with 10 or 01 signature and looks for relationships with expression
###

import pickle

# 0. user defined variables
jarDir='/Volumes/omics4tb/alomana/projects/csp.jgi/results/jars/'
signature01Pickle=jarDir+'signature01.pickle'

# 1. reading signature genes
f=open(signature01Pickle,'rb')
[genesAllSignature01,genesFilteredSignature01]=pickle.load(f)
f.close()

# 2. compute distribution of expression for signature genes
signatureExpression=

# 3. compute distribution of n samples of random genes

# 4. plot figure

print(len(genesAllSignature01),len(genesFilteredSignature01))
