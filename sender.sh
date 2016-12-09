#!/bin/bash

#$ -N deployment
#$ -o /proj/omics4tb/alomana/scratch/messages.deployment.o.txt
#$ -e /proj/omics4tb/alomana/scratch/messages.deployment.e.txt
#$ -P Bal_alomana
#$ -pe serial 40
#$ -q baliga
#$ -S /bin/bash

cd /users/alomana
source .bash_profile

cd /proj/omics4tb/alomana/projects/csp.jgi/src
time python patternFinder.py > messages.deployment.txt
