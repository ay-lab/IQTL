#!/bin/bash 

##===============
## although I have written this script in simple bash format
## in practice, please use some form of job submission enviroment
## like qsub or SLURM
##===============

# corresponds to the value K as in top-K loops
# same as used in the script 
TOPLOOPCOUNT=10000

CodeExec='Perform_PCA.R'

## Input directory containing the mastersheet
BaseInpDir='MasterSheet_top_'${TOPLOOPCOUNT}'_loops'

## input file to be used for PCA
inploopfile=$BaseInpDir'/MasterSheet_loops.bed'

## sample wise annotation file
DonorListFile='../Data/DonorList_Annotated.txt'

## --colDonor indicates the column in "DonorListFile" containing sample names
## --colListAnnot indicates the list of columns containing categorical information for the donors (comma or colon separated list)
## (which we want to use for PCA plot and dividing the donors in those specific categories)
Rscript $CodeExec --InpFile ${inploopfile} --AnnotFile ${DonorListFile} --colDonor 1 --colListAnnot 5:6:7:8 


