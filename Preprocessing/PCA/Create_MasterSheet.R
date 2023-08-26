#!/usr/bin/env Rscript

library(data.table)
library(UtilRPckg)
library(dplyr)

options(scipen = 10)
options(datatable.fread.datatable=FALSE)

args <- commandArgs(trailingOnly = TRUE)

p2p <- as.integer(args[1])
FitHiChIPBaseDir <- args[2]
MasterSheetDir <- args[3]
Binsize_Val <- as.integer(args[4])
FitHiChIP_Main_OutDir_Name <- args[5]
Dist_Low_Range <- as.integer(args[6])
Dist_High_Range <- as.integer(args[7])
TOPLOOPCOUNT <- as.integer(args[8])
DonorListFile <- args[9]
if (length(args) > 9) {
	DonorListFileHeader <- as.integer(args[10])
} else {
	DonorListFileHeader <- 0
}
if (length(args) > 10) {
	DonorCol <- as.integer(args[11])
} else {
	DonorCol <- 1
}

# FDR threshold of FitHiChIP loops
FDR_THR <- 0.01

CHRLIST_NAMENUM <- paste("chr", seq(1,22), sep="")

Parameter_File <- paste0(MasterSheetDir, '/Parameters.txt')
fp_Param <- file(Parameter_File, "w")
outtext <- paste0("\n *** Parameters involved in master sheet creation *** \n p2p : ", p2p, "\n FitHiChIPBaseDir: ", FitHiChIPBaseDir, "\n MasterSheetDir : ", MasterSheetDir, "\n DonorListFile : ", DonorListFile, "\n Binsize_Val: ", Binsize_Val, "\n FitHiChIP_Main_OutDir_Name: ", FitHiChIP_Main_OutDir_Name, "\n Dist_Low_Range : ", Dist_Low_Range, " \n Dist_High_Range : ", Dist_High_Range)
writeLines(outtext, con=fp_Param, sep="\n")
close(fp_Param)

# list of IQTL donors
if (DonorListFileHeader == 0) {
	DonorListData <- read.table(DonorListFile, header=F, sep="\t", stringsAsFactors=F, check.names=F)
} else {
	DonorListData <- read.table(DonorListFile, header=T, sep="\t", stringsAsFactors=F, check.names=F)
}
Complete_DonorList <- as.vector(DonorListData[,DonorCol])

ExtractChrData <- function(InpFile, chrName, OutFile=NULL, header=TRUE, dist=c(-1,-1), mid=FALSE) {
	if (is.null(OutFile)) {
		OutFile <- paste0(dirname(InpFile), "/temp_Chr_data.bed")
	}

	# process the distance thresholds
	# and insert in two variables
	if ((dist[1] > 0) & (dist[2] > 0) & (dist[2] > dist[1])) {
		distthrlow <- dist[1]
		distthrhigh <- dist[2]		
	} else {
		distthrlow <- -1
		distthrhigh <- -1
	}

	if (file_ext(InpFile) == "gz") {
		if (header == TRUE) {
			if (mid == TRUE) {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("zcat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if (NR>1 && $1==\"", chrName, "\" && $3==\"", chrName, "\" && (abs($2-$4)>=dt) && (abs($2-$4)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("zcat ", InpFile, " | awk \' {if (NR>1 && $1==\"", chrName, "\" && $3==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			} else {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("zcat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if (NR>1 && $1==\"", chrName, "\" && $4==\"", chrName, "\" && (abs($2-$5)>=dt) && (abs($2-$5)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("zcat ", InpFile, " | awk \' {if (NR>1 && $1==\"", chrName, "\" && $4==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			}
		} else {
			if (mid == TRUE) {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("zcat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if ($1==\"", chrName, "\" && $3==\"", chrName, "\" && (abs($2-$4)>=dt) && (abs($2-$4)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("zcat ", InpFile, " | awk \' {if ($1==\"", chrName, "\" && $3==\"", chrName, "\") {print $0}}\' -  > ", OutFile))						
				}
			} else {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("zcat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if ($1==\"", chrName, "\" && $4==\"", chrName, "\" && (abs($2-$5)>=dt) && (abs($2-$5)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("zcat ", InpFile, " | awk \' {if ($1==\"", chrName, "\" && $4==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			}
		}
	} else {
		if (header == TRUE) {
			if (mid == TRUE) {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("cat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if (NR>1 && $1==\"", chrName, "\" && $3==\"", chrName, "\" && (abs($2-$4)>=dt) && (abs($2-$4)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("cat ", InpFile, " | awk \' {if (NR>1 && $1==\"", chrName, "\" && $3==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			} else {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("cat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if (NR>1 && $1==\"", chrName, "\" && $4==\"", chrName, "\" && (abs($2-$5)>=dt) && (abs($2-$5)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("cat ", InpFile, " | awk \' {if (NR>1 && $1==\"", chrName, "\" && $4==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			}
		} else {
			if (mid == TRUE) {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("cat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if ($1==\"", chrName, "\" && $3==\"", chrName, "\" && (abs($2-$4)>=dt) && (abs($2-$4)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("cat ", InpFile, " | awk \' {if ($1==\"", chrName, "\" && $3==\"", chrName, "\") {print $0}}\' -  > ", OutFile))						
				}
			} else {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("cat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if ($1==\"", chrName, "\" && $4==\"", chrName, "\" && (abs($2-$5)>=dt) && (abs($2-$5)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("cat ", InpFile, " | awk \' {if ($1==\"", chrName, "\" && $4==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			}
		}
	}

}	# end function

GetNumLines <- function(inpfile) {
	nline <- as.integer(system(paste("cat", inpfile, "| wc -l"), intern = TRUE))
	return(nline)
}

##===============
## function
##===============
FillFeatureValues <- function(UnionLoopFile, UnionLoopFeatureFile, P2ALoopFileList, CHRLIST_NAMENUM, BinSize, Complete_DonorList) {

	valid_chr_count <- 0
	UnionLoopTempFile1 <- paste0(dirname(UnionLoopFile), '/temp_CurrChr_Merged_Loops.bed')

	for (chr_idx in (1:length(CHRLIST_NAMENUM))) {
		chrName <- CHRLIST_NAMENUM[chr_idx]
		cat(sprintf("\n ===>>> Within function FillFeatureValues --- processing chromosome : %s ", chrName))

		ExtractChrData(UnionLoopFile, chrName, UnionLoopTempFile1, header=F)
		nreadCurr <- GetNumLines(UnionLoopTempFile1)
		if (nreadCurr == 0) {
			next
		}
		MergedIntTempData <- data.table::fread(UnionLoopTempFile1, header=F)

		AllLoop_BinDF <- cbind.data.frame((MergedIntTempData[,2] / BinSize), (MergedIntTempData[,5] / BinSize))
		colnames(AllLoop_BinDF) <- c('B1', 'B2')

		RawCC_Categ <- matrix(0, nrow=nrow(MergedIntTempData), ncol=length(P2ALoopFileList))
		ExpCC_Categ <- matrix(0, nrow=nrow(MergedIntTempData), ncol=length(P2ALoopFileList))
		QVal_Categ <- matrix(1, nrow=nrow(MergedIntTempData), ncol=length(P2ALoopFileList))

		InpTempFitHiChIPLoopFile <- paste0(dirname(UnionLoopFile), '/temp_CurrChr_InpFitHiChIPLoopFile.bed')
			
		for (i in (1:length(P2ALoopFileList))) {
			inpfile <- P2ALoopFileList[i]
			
			system(paste0("awk -F\'[\t]\' -v b=", BinSize, " \'{if ((NR>1) && ($1==\"", chrName, "\")) {print ($2/b)\"\t\"($5/b)\"\t\"$7\"\t\"$(NF-4)\"\t\"$NF}}\' ", inpfile, " > ", InpTempFitHiChIPLoopFile))
			nreadInp <- UtilRPckg::GetNumLines(InpTempFitHiChIPLoopFile)
			
			if (nreadInp > 0) {
				
				InpTempData <- data.table::fread(InpTempFitHiChIPLoopFile, header=F)
				colnames(InpTempData) <- c('B1', 'B2', 'RawCC', 'ExpCC', 'qval')
				CN <- colnames(InpTempData)				
				cat(sprintf("\n ***** Computing overlap of merged loops with the FitHiChIP loop file index: %s name: %s for the chromosome : %s ***** \n", i, inpfile, chrName))

				mergeDF <- dplyr::left_join(AllLoop_BinDF, InpTempData)

				idx <- which(is.na(mergeDF$RawCC))
				mergeDF$RawCC[idx] <- 0
				mergeDF$ExpCC[idx] <- 0
				mergeDF$qval[idx] <- 1
				cat(sprintf("\n ----> number of master loops : %s  number of overlapping loops in input file : %s ", nrow(AllLoop_BinDF), (nrow(AllLoop_BinDF) - length(idx))))

				RawCC_Categ[, i] <- mergeDF$RawCC
				ExpCC_Categ[, i] <- mergeDF$ExpCC
				QVal_Categ[, i] <- mergeDF$qval
			
			}

			cat(sprintf("\n ***** Assigned raw contact count and q-values for the FitHiChIP loop file index: %s for the chromosome : %s ***** \n", i, chrName))

		}

		if (file.exists(InpTempFitHiChIPLoopFile) == TRUE) {
			system(paste("rm", InpTempFitHiChIPLoopFile))
		}

		valid_chr_count <- valid_chr_count + 1
		namesvec <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
		for (i in (1:length(P2ALoopFileList))) { 
			MergedIntTempData <- cbind.data.frame(MergedIntTempData, RawCC_Categ[, i], ExpCC_Categ[, i], QVal_Categ[, i])
		}
		appendnamevec <- c()
		for (i in (1:length(Complete_DonorList))) {
			for (v in c('_RawCC', '_ExpCC', '_QVal')) {
				appendnamevec <- c(appendnamevec, paste0(Complete_DonorList[i], v))
			}
		}
		colnames(MergedIntTempData) <- c(namesvec, appendnamevec)

		if (valid_chr_count == 1) {
			write.table(MergedIntTempData, UnionLoopFeatureFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)
		} else {
			write.table(MergedIntTempData, UnionLoopFeatureFile, row.names=F, col.names=F, sep="\t", quote=F, append=T)	
		}
	
	}

	if (file.exists(UnionLoopTempFile1) == TRUE) {
		system(paste("rm", UnionLoopTempFile1))
	}

}

#================================
# main code
#================================

Complete_P2A_FitHiChIP_Loop_FileList <- c()

UnionLoopFile <- paste0(MasterSheetDir, '/union_loops.bed')
UnionLoopTempFile <- paste0(MasterSheetDir, '/union_loops_temp.bed')
UnionLoopFeatureFile <- paste0(MasterSheetDir, '/MasterSheet_loops.bed')
if (file.exists(UnionLoopFile)) {
	system(paste("rm", UnionLoopFile))
}
if (file.exists(UnionLoopTempFile)) {
	system(paste("rm", UnionLoopTempFile))
}
if (file.exists(UnionLoopFeatureFile)) {
	system(paste("rm", UnionLoopFeatureFile))
}

for (i in 1:length(Complete_DonorList)) {
	cat(sprintf("\n ==>> creating union of top-K significant loops (autosomal chr) - donor idx : %s ", i))
	
	P2A_loopfile <- paste0(FitHiChIPBaseDir, '/', Complete_DonorList[i], '/', FitHiChIP_Main_OutDir_Name, '/FitHiChIP_Peak2ALL_b', Binsize_Val, '_L', Dist_Low_Range, '_U', Dist_High_Range, '/P2PBckgr_', p2p, '/Coverage_Bias/FitHiC_BiasCorr/FitHiChIP.interactions_FitHiC.bed')

	Complete_P2A_FitHiChIP_Loop_FileList <- c(Complete_P2A_FitHiChIP_Loop_FileList, P2A_loopfile)

	P2A_sigloopfile <- paste0(FitHiChIPBaseDir, '/', Complete_DonorList[i], '/', FitHiChIP_Main_OutDir_Name, '/FitHiChIP_Peak2ALL_b', Binsize_Val, '_L', Dist_Low_Range, '_U', Dist_High_Range, '/P2PBckgr_', p2p, '/Coverage_Bias/FitHiC_BiasCorr/FitHiChIP.interactions_FitHiC_Q', FDR_THR, '.bed')

	cat(sprintf("\n ==>> P2A_sigloopfile : %s ", P2A_sigloopfile))
	
	numsigloop <- as.integer(system(paste0("awk \'((NR>1) && ($1!=\"chrX\") && ($1!=\"chrY\"))\' ", P2A_sigloopfile, " | wc -l"), intern = TRUE))
	if (numsigloop >= TOPLOOPCOUNT) {
		system(paste0("awk \'{if ((NR>1) && ($1!=\"chrX\") && ($1!=\"chrY\")) {print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$7\"\t\"$NF}}\' ", P2A_sigloopfile, " | sort -k8,8g | head -n ", TOPLOOPCOUNT, " | cut -f1-6 >> ", UnionLoopTempFile))
	} else {
		system(paste0("awk \'{if ((NR>1) && ($1!=\"chrX\") && ($1!=\"chrY\")) {print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$7\"\t\"$NF}}\' ", P2A_loopfile, " | sort -k8,8g | head -n ", TOPLOOPCOUNT, " | cut -f1-6 >> ", UnionLoopTempFile))
	}
}
system(paste("sort -k1,1 -k2,2n -k5,5n", UnionLoopTempFile, "| uniq >", UnionLoopFile))
FillFeatureValues(UnionLoopFile, UnionLoopFeatureFile, Complete_P2A_FitHiChIP_Loop_FileList, CHRLIST_NAMENUM, Binsize_Val, Complete_DonorList)


