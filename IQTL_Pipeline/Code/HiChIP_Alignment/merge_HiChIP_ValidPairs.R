#!/usr/bin/env Rscript

#===========================
## 1) Merges donor-wise HiChIP valid pairs (one donor may have multiple samples, depending on the sequencing runs)
## We merge CIS reads subject to the specified distance thresholds
#===========================

library(data.table)
options(scipen = 10)
options(datatable.fread.datatable=FALSE)

args <- commandArgs(trailingOnly = TRUE)

DonorListFile <- args[1]
HiChIPDataDir <- args[2]
IQTLResDir <- args[3]
LowDistThr <- as.integer(args[4])
HighDistThr <- as.integer(args[5])
AutosomalChr <- as.integer(args[6])
TargetDonor <- as.character(args[7])

##==== read the complete donor list (annotated)
Donor_DF <- read.table(DonorListFile, header=T, sep="\t", stringsAsFactors=F)
DonorList <- unique(Donor_DF[,1])

## this directory will store the merged HiChIP CIS reads for the specific donor
HiChIPReadDir <- paste0(IQTLResDir, '/HiChIP_CIS_Reads/', TargetDonor)
system(paste("mkdir -p", HiChIPReadDir))

mergedCISLongFile <- paste0(HiChIPReadDir, '/Merged_CIS_Validpairs.txt')

samplelist <- Donor_DF[which(Donor_DF[,1] == TargetDonor), 2]
cat(sprintf("\n\n ==>> TargetDonor : %s ", TargetDonor))
cat(sprintf("\n\n ==>> samplelist : %s ", paste(samplelist, collapse=" ")))

##=========================
## merge the HiChIP valid pairs for the target donor
## and for the constituent samples
## CIS reads with respect to the specified distance thresholds are merged
##=========================

for (i in 1:length(samplelist)) {
	cat(sprintf("\n Processing sample : %s ", samplelist[i]))
	inpfile1 <- paste0(HiChIPDataDir, '/', samplelist[i], '/HiCPro/hic_results/data/rawdata/rawdata.allValidPairs')
	inpfile2 <- paste0(HiChIPDataDir, '/', samplelist[i], '/HiCPro/hic_results/data/rawdata/rawdata_allValidPairs')

	if (file.exists(inpfile1)) {
		inpfile <- inpfile1
	} else {
		inpfile <- inpfile2
	}

	if (i == 1) {
		if (AutosomalChr == 0) {
			## reads for all chromosomes, including X and Y
			system(paste0("awk -v l=", LowDistThr, " -v h=", HighDistThr, " \'(($2==$5) && (($6-$3)>=l) && (($6-$3)<=h))\' ", inpfile, " > ", mergedCISLongFile))
		} else {
			system(paste0("awk -v l=", LowDistThr, " -v h=", HighDistThr, " \'(($2==$5) && ($2 ~ /^chr([1-9]|2[0-2]|1[0-9])$/ ) && (($6-$3)>=l) && (($6-$3)<=h))\' ", inpfile, " > ", mergedCISLongFile))
		}
	} else {
		if (AutosomalChr == 0) {
			## reads for all chromosomes, including X and Y
			system(paste0("awk -v l=", LowDistThr, " -v h=", HighDistThr, " \'(($2==$5) && (($6-$3)>=l) && (($6-$3)<=h))\' ", inpfile, " >> ", mergedCISLongFile))
		} else {
			system(paste0("awk -v l=", LowDistThr, " -v h=", HighDistThr, " \'(($2==$5) && ($2 ~ /^chr([1-9]|2[0-2]|1[0-9])$/ ) && (($6-$3)>=l) && (($6-$3)<=h))\' ", inpfile, " >> ", mergedCISLongFile))
		}
	}
}

## compress the merged valid pairs file
system(paste0("gzip ", mergedCISLongFile))

