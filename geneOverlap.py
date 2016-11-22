###
### this script reads genes with 10 or 01 signature and looks for relationships with expression
###

import pickle

# 0. user defined variables
jarDir='/Volumes/omics4tb/alomana/projects/csp.jgi/results/jars/'
signature10Pickle=jarDir+'signature10.pickle'

# 1. reading signature genes
f=open(signature10Pickle,'rb')
[genesAllSignature10,genesFilteredSignature10]=pickle.load(f)
f.close()

print(genesAllSignature10)
