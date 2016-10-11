# mendiak
Quantitative tools for the analysis of FAIRE-seq data.  

**readsMapper.py**  
Tool to map reads to genome using STAR.

duplicatesCleaner.py  
Tool to remove duplicates from BAM files using samtools.

chromosomeSizeProvider.py  
Tool to obtain chromosome sizes from a FASTA file. Required by bedGraphToBigWig.

peakCaller.macs2.py  
Tool to call peaks using MACS2.0 in a SGE environment.
