# TKGWV2
# An ancient DNA relatedness pipeline for ultra-low coverage whole genome shotgun data

## Description
TKGWV2 is a pipeline to estimate biological relatedness between individuals specifically aimed at ultra-low coverage ancient DNA data obtained from whole genome sequencing.
It is a massive update to the original TKrelated method published in 2017 (Fernandes et al.), with ease of use and efficiency as major concerns.
From the initial suggested minimum coverage threshold of 0.1X, TKGWV2 can be applied to BAM files of individuals with 0.025X, or, in some cases, as little as 0.005X when the other individual has above 0.03X. These characteristics have the potential to offer relatedness estimation during screening sequencing steps at early stages of an ancient DNA project, and can therefore be very useful for project planning.
Although mainly designed for shotgun data, TKGWV2 also works with, for example, the widely-used 1240K capture dataset, as long as allele frequencies are provided.

## Requirements
TKGWV2 was developed for Linux. The user is required to have the following software available as a system-wide installation:
- Python 3
- PLINK 1.9
- samtools (tested on version 1.7)
- R (tested on version 3+)
- R package: data.table (install.packages("data.table"))

## Installation and usage
There is no installation required for TKGWV2. All you need to do is download and unzip this package. Keep the folder structure and make sure that all Python and R scripts are executable.




You will need two sets of files to run TKGWV2. First, download and unzip this package to a folder of your choice. Then, download the set of support files from:
https://drive.google.com/drive/folders/1Aw-0v_7CUorHJOLpCJ0QVCwEdH43HlO4?usp=sharing

The set of support files were tested in the original publication of TKGWV2 and shown to work on any ancient population from Europe






- PLINK frequency file was generated from 1000 Genomes Phase 3 CEU population. The only QC step required is to exclude fixed SNPs, reducing from 77 to 22 million variants
(maybe have default EUR files, where both the BED file and FRQ file have no fixed SNPs, will make things much faster)






### Converting BAMs to individual text-PLINKs
# Programs with system-wide installation required # 
# 1. samtools
# 2. plink
# 3. Python3
#
# Description of arguments #
# referenceGenome = Path to index fasta file used to generate BAM files
# gwvList = List of positions to generate pileup file from in BED format
# gwvPlink = PLINK dataset with same positions as gwvList, to generate MAP files. A dummy dataset with a single individual is provided
# bamExtension = Default ".bam". Specific extension/suffix of the BAM files to be used. Anything before this will be considered as the sample ID
# minMQ = Default 30. Minimum mapping quality
# minBQ = Default 30. Minimum base quality
# excludeTerminalReadBases = Default FALSE. Exclude terminal read positions from pileup file
#
# Input:
# - individual BAM
# Output:
# - individual text-PLINK




# Programs with system-wide installation required # 
# 1. plink
# R package required # 
# 1. data.table  
#
# Description of arguments #
# dyads = Default TRUE will analyse all possible dyads between all individuals in current folder.
#         Otherwise, user can input 'dyads' as a tab-spaced text file with each line containing a pair to be analysed.
#         Sample1 Sample2
#         Sample1 Sample3
#         (...)
# freqFile = Allele frequencies file in PLINK format (FRQ) for the same SNP set used in bam2plink() (or for example for the 1240K dataset)
# ignoreThresh = Default 1. Threshold for the minimum number of SNPs allowed to estimate relatedness.
#
# Input:
# - individual text-PLINK
# - pairwise FRQ
# Output:
# - pairwise relatedness coefficient result text file  




## Citation
Fernandes, DM, Cheronet, O, Gelabert, P, Pinhasi, R. TKGWV2: An ancient DNA relatedness pipeline for ultra-low coverage whole genome shotgun data (2021). Preprint.
