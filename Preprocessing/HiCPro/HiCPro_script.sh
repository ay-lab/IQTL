#!/bin/bash

##===============
## although I have written this script in simple bash format
## in practice, please use some form of job submission enviroment
## like qsub or SLURM
## since these HiCPro computation are time consuming
##===============

## sample name - user can specify custom sample / donor names 
prefix='DONOR1_XXXX'

# input directory containing all the samples and fastq files
# we recommend keeping all the samples (and their fastq files) under a common directory
# which is referred here as the "inpbasedir"
inpbasedir='/path/to/HiChIP/data/'

# input directory for this particular data (for this sample)
# we assume that the fastq files are placed under the "rawdata" folder for this particular sample
inpdir=$inpbasedir$prefix'/rawdata/'

# configuration file containing the execution parameters of HiCPro
# see the "configfile" provided along with this GitHub repository
configfile='configfile'

# directory which will contain the output results for this sample
outdir=$inpbasedir$prefix'/HiCPro'

# execute the HiCPro pipeline
HiC-Pro -i $inpdir -o $outdir -c $configfile


