#!/usr/bin/env Rscript

#====================================
# script for iQTL project
# creating donor specific master sheet containing union of FitHiChIP loops
# and their contact count and significance values per donor

## Note: in IQTL, all the donors employ the same peak information for FitHiChIP loops
#====================================

library(optparse)
library(data.table)
library(dplyr)
options(scipen = 10)
options(datatable.fread.datatable=FALSE)

#====================================================
option_list = list(
	make_option(c("--P2P"), type="integer", action="store", default=0, help="FitHiChIP background model. 0: loose, 1: stringent."),
	make_option(c("--FitHiChIPBaseDir"), type="character", default=NULL, help="Base directory containing FitHiChIP loops for all samples. Mandatory parameter."),
	make_option(c("--MasterSheetDir"), type="character", default=NULL, help="Output directory to contain the mastersheet. Mandatory parameter."),
	make_option(c("--DonorFile"), type="character", default=NULL, help="File containing donor list. One column file, where each sample represents donor name. Mandatory parameter."),
	make_option(c("--BinSize"), type="integer", action="store", default=5000, help="Bin size employed. Default 5000 (5 Kb)."),
	make_option(c("--FitHiChIPDirPrefix"), type="character", default=NULL, help="Prefix of FitHiChP output folder. Mandatory parameter."),
	make_option(c("--DistLow"), type="integer", action="store", default=20000, help="Lower distance threshold for FitHiChIP loops. Default 20000 (20 Kb)."),
	make_option(c("--DistHigh"), type="integer", action="store", default=2000000, help="Higher distance threshold for FitHiChIP loops. Default 2000000 (2 Mb)."),
	make_option(c("--UseAllLoop"), type="integer", action="store", default=0, help="If 1, master sheet of all loops (significant or not). Default = 0"),
	make_option(c("--RefLoopFile"), type="character", default=NULL, help="If provided, contains the reference loops for mastersheet. Otherwise, by default, union of donor wise loops is used."),
	make_option(c("--Auto"), type="integer", action="store", default=1, help="Binary (1/0) value, indicates if the autosomal chromosomes are only used. Default 1."),
	make_option(c("--RefGenome"), type="character", default=NULL, help="Reference genome. Available options: hg19, hg38, mm9, mm10, mm39. Mandatory parameter."),
	make_option(c("--FDRThr"), type="numeric", action="store", default=0.01, help="FDR Threshold of FitHiChIP loops. Default 0.01.")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# FitHiChIP background model
p2p <- as.integer(opt$P2P)
FitHiChIPBaseDir <- opt$FitHiChIPBaseDir
MasterSheetDir <- opt$MasterSheetDir
DonorListFile <- opt$DonorFile
Binsize_Val <- as.integer(opt$BinSize)
FitHiChIP_Main_OutDir_Name <- opt$FitHiChIPDirPrefix
Dist_Low_Range <- as.integer(opt$DistLow)
Dist_High_Range <- as.integer(opt$DistHigh)
if (!is.null(opt$RefLoopFile)) {
	RefLoopFile <- opt$RefLoopFile
}

# set of chromosomes which would be individually analyzed for testing overlap
if ((opt$RefGenome == "hg19") | (opt$RefGenome == "hg38")) {
	if (opt$Auto == 0) {
		CHRLIST_NAMENUM <- c(paste("chr", seq(1,22), sep=""), "chrX", "chrY")
	} else {
		CHRLIST_NAMENUM <- c(paste("chr", seq(1,22), sep=""))
	}
} else if ((opt$RefGenome == "mm9") | (opt$RefGenome == "mm10") | (opt$RefGenome == "mm39")) {
	if (opt$Auto == 0) {
		CHRLIST_NAMENUM <- c(paste("chr", seq(1,19), sep=""), "chrX", "chrY")
	} else {
		CHRLIST_NAMENUM <- c(paste("chr", seq(1,19), sep=""))
	}
}

system(paste("mkdir -p", MasterSheetDir))

# FDR threshold of FitHiChIP loops
FDR_THR <- as.numeric(opt$FDRThr)	# 0.01

Parameter_File <- paste0(MasterSheetDir, '/Parameters.txt')
fp_Param <- file(Parameter_File, "w")
outtext <- paste0("\n *** Parameters involved in master sheet creation *** \n p2p : ", p2p, "\n FitHiChIPBaseDir: ", FitHiChIPBaseDir, "\n MasterSheetDir : ", MasterSheetDir, "\n DonorListFile : ", DonorListFile, "\n Binsize_Val: ", Binsize_Val, "\n FitHiChIP_Main_OutDir_Name: ", FitHiChIP_Main_OutDir_Name, "\n Dist_Low_Range : ", Dist_Low_Range, " \n Dist_High_Range : ", Dist_High_Range)
writeLines(outtext, con=fp_Param, sep="\n")
close(fp_Param)

# list of IQTL donors
DonorListData <- read.table(DonorListFile, header=F, sep="\t", stringsAsFactors=F, check.names=F)
Complete_DonorList <- as.vector(DonorListData[,1])

#================================
# this function annotates individual loops in the merged union set of loops
# according to the feature vectors of individual loops 
# parameters:
# UnionLoopFile: file containing merged set of loops from all input replicates and all categories
# UnionLoopFeatureFile: file which will contain all features corresponding to all loops in UnionLoopFile
# AllLoopList: Loops with FitHiChIP significance for different input replicates and categories
# CHRLIST_NAMENUM: list of chromosome names
# BinSize: bin size of input loops
# Complete_DonorList: Label of donors 
#================================
FillFeatureValues <- function(UnionLoopFile, UnionLoopFeatureFile, P2ALoopFileList, CHRLIST_NAMENUM, BinSize, Complete_DonorList) {

	##====== valid chromosome counter for which there is at least one loop in the complete set of loops
	valid_chr_count <- 0

	##======= temporary file
	UnionLoopTempFile1 <- paste0(dirname(UnionLoopFile), '/temp_CurrChr_Merged_Loops.bed')

	for (chr_idx in (1:length(CHRLIST_NAMENUM))) {
		chrName <- CHRLIST_NAMENUM[chr_idx]
		cat(sprintf("\n ===>>> Within function FillFeatureValues --- processing chromosome : %s ", chrName))

		##=========== first extract the loops involving current chromosome		
		UtilRPckg::ExtractChrData(UnionLoopFile, chrName, UnionLoopTempFile1, header=F)
		nreadCurr <- UtilRPckg::GetNumLines(UnionLoopTempFile1)
		if (nreadCurr == 0) {
			next
		}
		MergedIntTempData <- data.table::fread(UnionLoopTempFile1, header=F)

		##======== get the interacting bins (start position / BinSize)
		AllLoop_BinDF <- cbind.data.frame((MergedIntTempData[,2] / BinSize), (MergedIntTempData[,5] / BinSize))
		colnames(AllLoop_BinDF) <- c('B1', 'B2')

		##======== extract CC, q-values for loops of the current chromosome
		RawCC_Categ <- matrix(0, nrow=nrow(MergedIntTempData), ncol=length(P2ALoopFileList))
		QVal_Categ <- matrix(1, nrow=nrow(MergedIntTempData), ncol=length(P2ALoopFileList))
		ExpCC_Categ <- matrix(0, nrow=nrow(MergedIntTempData), ncol=length(P2ALoopFileList))

		# stores temporary FitHiChIP loops for the current chromosome
		InpTempFitHiChIPLoopFile <- paste0(dirname(UnionLoopFile), '/temp_CurrChr_InpFitHiChIPLoopFile.bed')

		##====== first process the P2ALoopFileList (peak to all contacts), and get the raw contact counts, expected contact counts and q-values
		for (i in (1:length(P2ALoopFileList))) {
			inpfile <- P2ALoopFileList[i]
			##==== extract the interacting bins (fields 2 and 5), raw contact count (field 7)
			system(paste0("awk -v b=", BinSize, " \'{if ((NR>1) && ($1==\"", chrName, "\")) {print ($2/b)\"\t\"($5/b)\"\t\"$7\"\t\"$(NF-4)\"\t\"$NF}}\' ", inpfile, " > ", InpTempFitHiChIPLoopFile))
			nreadInp <- UtilRPckg::GetNumLines(InpTempFitHiChIPLoopFile)
			
			if (nreadInp > 0) {
				InpTempData <- data.table::fread(InpTempFitHiChIPLoopFile, header=F)
				colnames(InpTempData) <- c('B1', 'B2', 'RawCC', 'ExpCC', 'qval')
				CN <- colnames(InpTempData)			
				cat(sprintf("\n ***** Computing overlap of merged loops with the FitHiChIP loop file index: %s name: %s for the chromosome : %s ***** \n", i, inpfile, chrName))

				## merge - keep all elements of "AllLoop_BinDF"
				mergeDF <- dplyr::left_join(AllLoop_BinDF, InpTempData)
				
				## fill the NA entries
				idx <- which(is.na(mergeDF$RawCC))
				mergeDF$RawCC[idx] <- 0
				mergeDF$ExpCC[idx] <- 0
				mergeDF$qval[idx] <- 1
				cat(sprintf("\n ----> number of master loops : %s  number of overlapping loops in input file : %s ", nrow(AllLoop_BinDF), (nrow(AllLoop_BinDF) - length(idx))))

				##===== store the rawCC and q-value
				RawCC_Categ[, i] <- mergeDF$RawCC
				ExpCC_Categ[, i] <- mergeDF$ExpCC
				QVal_Categ[, i] <- mergeDF$qval

			}	# end number of reads condition

			cat(sprintf("\n ***** Assigned raw + expected contact counts, and q-values for the FitHiChIP loop file index: %s for the chromosome : %s ***** \n", i, chrName))

		}	# end processing peak-to-all loop list

		# delete the temporary file
		if (file.exists(InpTempFitHiChIPLoopFile) == TRUE) {
			system(paste("rm", InpTempFitHiChIPLoopFile))
		}

		##======= for the current chromosome, merge the interacting regions along with the features accumulated
		##======= and construct the final data frame
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

		# write the data frame sequentially
		if (valid_chr_count == 1) {
			write.table(MergedIntTempData, UnionLoopFeatureFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)
		} else {
			write.table(MergedIntTempData, UnionLoopFeatureFile, row.names=F, col.names=F, sep="\t", quote=F, append=T)	
		}
	
	} 	# end chromosome index loop

	# delete the temporary file
	if (file.exists(UnionLoopTempFile1) == TRUE) {
		system(paste("rm", UnionLoopTempFile1))
	}

}	# end function

#================================
# main code
#================================

# list of files storing peak-to-all interactions (contact counts + q-values are stored) for different donors
Complete_P2A_FitHiChIP_Loop_FileList <- c()

# list of files storing peak-to-all significant interactions (contact counts + q-values are stored) for different donors
Complete_P2A_FitHiChIP_SigLoop_FileList <- c()

UnionLoopFile <- paste0(MasterSheetDir, '/union_loops.bed')
UnionLoopTempFile <- paste0(MasterSheetDir, '/union_loops_temp.bed')
UnionLoopFeatureFile <- paste0(MasterSheetDir, '/MasterSheet_loops.bed')

##========= note down the donor specific interaction files
for (i in 1:length(Complete_DonorList)) {
	# peak-to-all interactions
	P2A_loopfile <- paste0(FitHiChIPBaseDir, '/', Complete_DonorList[i], '/', FitHiChIP_Main_OutDir_Name, '/FitHiChIP_Peak2ALL_b', Binsize_Val, '_L', Dist_Low_Range, '_U', Dist_High_Range, '/P2PBckgr_', p2p, '/Coverage_Bias/FitHiC_BiasCorr/FitHiChIP.interactions_FitHiC.bed')
	# peak-to-all significant interactions
	P2A_sigloopfile <- paste0(FitHiChIPBaseDir, '/', Complete_DonorList[i], '/', FitHiChIP_Main_OutDir_Name, '/FitHiChIP_Peak2ALL_b', Binsize_Val, '_L', Dist_Low_Range, '_U', Dist_High_Range, '/P2PBckgr_', p2p, '/Coverage_Bias/FitHiC_BiasCorr/FitHiChIP.interactions_FitHiC_Q', FDR_THR, '.bed')

	Complete_P2A_FitHiChIP_Loop_FileList <- c(Complete_P2A_FitHiChIP_Loop_FileList, P2A_loopfile)
	Complete_P2A_FitHiChIP_SigLoop_FileList <- c(Complete_P2A_FitHiChIP_SigLoop_FileList, P2A_sigloopfile)
}

#=========================
## first create the union set of loops
## option 1: if opt$RefLoopFile is not NULL, use those loops
## option 2: otherwise, construct the union of significant FitHiChIP loops for all donors
#=========================
if (!is.null(opt$RefLoopFile)) {
	## option 1
	## use the reference loops
	if (file.exists(UnionLoopFile) == FALSE) {
		system(paste0("awk \'{if (NR>1) {print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}}\' ", RefLoopFile, " > ", UnionLoopFile))
	}
} else {
	## option 2
	## create the union of donor specific FitHiChIP significant files
	if (file.exists(UnionLoopFile) == FALSE) {
		## chromosome wise processing
		for (chrIdx in 1:length(CHRLIST_NAMENUM)) {
			currchr <- CHRLIST_NAMENUM[chrIdx]
			##========= create the union set of loops 
			cat(sprintf("\n ===> constructing union set of loops - processing chromosome : %s ", currchr))
			if (opt$UseAllLoop == 1) {
				for (i in 1:length(Complete_P2A_FitHiChIP_Loop_FileList)) {
					P2A_loopfile <- Complete_P2A_FitHiChIP_Loop_FileList[i]
					cat(sprintf("\n processing P2A loop file idx : %s ", i))
					if (i == 1) {
						system(paste0("awk -F\'[\t]\' \'((NR>1) && ($1==\"", currchr, "\") && ($4==\"", currchr, "\"))\' ", P2A_loopfile, " | cut -f1-6 > ", UnionLoopTempFile))
					} else {
						system(paste0("awk -F\'[\t]\' \'((NR>1) && ($1==\"", currchr, "\") && ($4==\"", currchr, "\"))\' ", P2A_loopfile, " | cut -f1-6 >> ", UnionLoopTempFile))
					}
				}
			} else {
				for (i in 1:length(Complete_P2A_FitHiChIP_SigLoop_FileList)) {
					P2A_sigloopfile <- Complete_P2A_FitHiChIP_SigLoop_FileList[i]
					cat(sprintf("\n processing P2A sigloop file idx : %s ", i))
					if (i == 1) {
						system(paste0("awk -F\'[\t]\' \'((NR>1) && ($1==\"", currchr, "\") && ($4==\"", currchr, "\"))\' ", P2A_sigloopfile, " | cut -f1-6 > ", UnionLoopTempFile))
					} else {
						system(paste0("awk -F\'[\t]\' \'((NR>1) && ($1==\"", currchr, "\") && ($4==\"", currchr, "\"))\' ", P2A_sigloopfile, " | cut -f1-6 >> ", UnionLoopTempFile))
					}
				}
			}

			##====== get the unique set of loops (from the All-to-All interactions, combining all donors)
			if (chrIdx == 1) {
				system(paste("sort -k2,2n -k5,5n", UnionLoopTempFile, "| uniq >", UnionLoopFile))
			} else {
				system(paste("sort -k2,2n -k5,5n", UnionLoopTempFile, "| uniq >>", UnionLoopFile))
			}
		}	# end chromosome loop
	}
}	# end file exist condition

##====== then fill the contact count, and significance 
if (file.exists(UnionLoopFeatureFile) == FALSE) {
	FillFeatureValues(UnionLoopFile, UnionLoopFeatureFile, Complete_P2A_FitHiChIP_Loop_FileList, CHRLIST_NAMENUM, Binsize_Val, Complete_DonorList)
}

##====== append two columns: 1) LoopCnt: number of donors having this contact (significant or not)
##====== 2) SigLoopCnt: number of donors having this contact as a significant contact
UnionLoopFeatureCntFile <- paste0(MasterSheetDir, '/MasterSheet_loops_with_Count.bed')
system(paste0("awk \'{if (NR==1) {print $0\"\tLoopCnt\tSigLoopCnt\"} else {l=0;s=0;for (i=7;i<=NF;i=i+3) {if ($i>0) {l=l+1}; if (($i>0) && ($(i+2)< ", FDR_THR, ")) {s=s+1}}; print $0\"\t\"l\"\t\"s}}\' ", UnionLoopFeatureFile, " > ", UnionLoopFeatureCntFile))


