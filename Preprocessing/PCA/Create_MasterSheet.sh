#!/bin/bash

##===============
## although I have written this script in simple bash format
## in practice, please use some form of job submission enviroment
## like qsub or SLURM
##===============

CodeExec='Create_MasterSheet_Donor_Specific.R'

# corresponds to the value K as in top-K loops
TOPLOOPCOUNT=10000	#5000

## lower distance threshold in FitHiChIP loops
DistLow=10000	#10000

## upper distance threshold in FitHiChIP loops
DistHigh=3000000

## peak-to-peak background employed in FitHiChIP loops
## 0: loose background, 1: stringent background
p2p=0

## bin size employed in FitHiChIP loops
binsize=5000

##=== base directory containing all the sample specific FitHiChIP loops
##=== If the sample names are *sample1*, *sample2*, etc.. then the output directory structure would be *FitHiChIPDir/Sample1*, *FitHiChIPDir/Sample2*, ... where *FitHiChIPDir* is the base directory name.
FitHiChIPBaseDir='/path/to/FitHiChIPDir'

## list of donors along with their annotations (sample-wise information)
## provided in the "Data" folder
DonorListFile='../Data/DonorList_Annotated.txt'

## FitHiChIP directory name (within each sample specific folders)
## basically the output directory name provided in the FitHiChIP configuration file
## we assume that all samples have the same output folder name
dirprefix='FitHiChIP_10Kb_3Mb_Resolution_5Kb'

##=== base output directory
OutDir='MasterSheet_top_'$TOPLOOPCOUNT'_loops'
mkdir -p $OutDir

Rscript ${CodeExec} ${p2p} ${FitHiChIPBaseDir} ${OutDir} ${binsize} ${dirprefix} ${DistLow} ${DistHigh} ${TOPLOOPCOUNT} ${DonorListFile} 1

