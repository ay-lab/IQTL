#!/usr/bin/env Rscript

##============
## script to create FitHiChIP configuration file
## and run FitHiChIP
##============

options(scipen = 10)
# options(datatable.fread.datatable=FALSE)

args <- commandArgs(trailingOnly = TRUE)

InpConfigFile <- as.character(args[1])
OutDir <- as.character(args[2])
ValidPairFile <- as.character(args[3])
PeakFile <- as.character(args[4])
FitHiChIPDirPrefix <- as.character(args[5])
BinSize <- as.integer(args[6])
LowDistThr <- as.integer(args[7])
HighDistThr <- as.integer(args[8])
P2P <- as.integer(args[9])
FDRThr <- as.numeric(args[10])
ChrSizeFile <- as.character(args[11])
FitHiChIPExec <- as.character(args[12])

cat(sprintf("\n ==>> InpConfigFile : %s \n OutDir : %s \n ValidPairFile : %s \n PeakFile : %s \n FitHiChIPDirPrefix : %s \n BinSize : %s \n LowDistThr : %s \n HighDistThr : %s \n P2P : %s \n FDRThr : %s \n ChrSizeFile : %s \n FitHiChIPExec : %s ", InpConfigFile, OutDir, ValidPairFile, PeakFile, FitHiChIPDirPrefix, BinSize, LowDistThr, HighDistThr, P2P, FDRThr, ChrSizeFile, FitHiChIPExec))

system(paste("mkdir -p", OutDir))
OutConfigFile <- paste0(OutDir, "/", basename(InpConfigFile))
system(paste0("cp ", InpConfigFile, " ", OutConfigFile))
cat(sprintf("\n ==>> OutConfigFile : %s ", OutConfigFile))

## insert escape character "\" before every "/" before using the sed operation
ValidPairFile_Mod <- gsub("/", "\\\\/", ValidPairFile)
system(paste0("sed -i \"s/ValidPairs=/ValidPairs=", ValidPairFile_Mod, "/g\" ", OutConfigFile))
PeakFile_Mod <- gsub("/", "\\\\/", PeakFile)
system(paste0("sed -i \"s/PeakFile=/PeakFile=", PeakFile_Mod, "/g\" ", OutConfigFile))
FitHiChIPOutDir <- paste0(OutDir, "/", FitHiChIPDirPrefix)
FitHiChIPOutDir_Mod <- gsub("/", "\\\\/", FitHiChIPOutDir)
system(paste0("sed -i \"s/OutDir=/OutDir=", FitHiChIPOutDir_Mod, "/g\" ", OutConfigFile))
system(paste0("sed -i \"s/BINSIZE=/BINSIZE=", BinSize, "/g\" ", OutConfigFile))
system(paste0("sed -i \"s/LowDistThr=/LowDistThr=", LowDistThr, "/g\" ", OutConfigFile))
system(paste0("sed -i \"s/UppDistThr=/UppDistThr=", HighDistThr, "/g\" ", OutConfigFile))
system(paste0("sed -i \"s/UseP2PBackgrnd=/UseP2PBackgrnd=", P2P, "/g\" ", OutConfigFile))
system(paste0("sed -i \"s/QVALUE=/QVALUE=", FDRThr, "/g\" ", OutConfigFile))
ChrSizeFile_Mod <- gsub("/", "\\\\/", ChrSizeFile)
system(paste0("sed -i \"s/ChrSizeFile=/ChrSizeFile=", ChrSizeFile_Mod, "/g\" ", OutConfigFile))

system(paste(FitHiChIPExec, "-C", OutConfigFile))


