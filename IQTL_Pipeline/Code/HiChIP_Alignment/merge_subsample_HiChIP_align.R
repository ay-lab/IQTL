#!/usr/bin/env Rscript

#===========================
## 1) Subsamples merged HiChIP alignment files to include only the CIS valid pair reads (from HiCPro) within specified distance range
#===========================

library(data.table)
options(scipen = 10)
options(datatable.fread.datatable=FALSE)

args <- commandArgs(trailingOnly = TRUE)

DonorListFile <- args[1]
HiChIPDataDir <- args[2]
IQTLResDir <- args[3]
AutosomalChr <- as.integer(args[4])
TargetDonor <- as.character(args[5])
Targetchr <- as.character(args[6])

##==== read the complete donor list (annotated)
Donor_DF <- read.table(DonorListFile, header=T, sep="\t", stringsAsFactors=F)
DonorList <- unique(Donor_DF[,1])

inpvalidpairfile <- paste0(IQTLResDir, '/HiChIP_CIS_Reads/', TargetDonor, '/Merged_CIS_Validpairs.txt.gz')
mergedalignfilename <- paste0(IQTLResDir, '/HiChIP_Complete_Alignments/', TargetDonor, '/chrwise/merged_HiChIP_', Targetchr, '_sorted.bam')
SubsampleAlignDir <- paste0(IQTLResDir, '/HiChIP_Alignments_Subsample_CIS_Reads/', TargetDonor, '/chrwise')
system(paste("mkdir -p", SubsampleAlignDir))
subsamplealignfile <- paste0(SubsampleAlignDir, '/merged_HiChIP_', Targetchr, '_subsampled_CIS_sorted.bam')
	
## first extract the read names from the HiChIP valid pairs file (subsampled)
readnamefile <- paste0(SubsampleAlignDir, '/target_readnames_', Targetchr, '.txt')
system(paste("zcat", inpvalidpairfile, "| awk \'{if ($2==\"", Targetchr, "\") {print $1}}' - >", readnamefile))	

## then extract these reads (corresponding to the read names) from the input alignment file
system(paste("samtools view -h -N", readnamefile, mergedalignfilename, "| samtools sort -o", subsamplealignfile, "-@ 8 -"))	

## then index the bam file
system(paste("samtools index", subsamplealignfile))	

## then remove the temporary read name file
system(paste("rm", readnamefile))	



