# mendiak
Quantitative tools for the analysis of FAIRE-seq data.  

**readsCleaner.py**  
Tool to clean reads using Trimmomatic.

**readsMapper.py**  
Tool to map reads to genome using STAR.

**duplicatesCleaner.py**  
Tool to remove duplicates from BAM files using samtools.

**chromosomeSizeProvider.py**  
Tool to obtain chromosome sizes from a FASTA file. Required by bedGraphToBigWig.

**correlationBuilder.py**  
Tool to obtain a figure of correlation between samples based on wigCorrelate

**peakCaller.macs2.py**  
Tool to call FSSs using MACS2.0 in a SGE environment.

**saturationCurveAnalysis.py**  
Tool to generate subsampled FASTQ files for the generation of saturation curves. It also create the final figures. Intermediate steps of mapping reads and calling peaks is done by other tools, i.e., readsMapper.py and peakCaller.macs2.py respectively.

**patternFinder.py**  
Tool to find genes that have 01 or 10 signature of FSSs.
