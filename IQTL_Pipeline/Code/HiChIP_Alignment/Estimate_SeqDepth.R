#!/usr/bin/env Rscript

library(data.table)
options(scipen = 10)
options(datatable.fread.datatable=FALSE)

args <- commandArgs(trailingOnly = TRUE)

OldSampleFile <- args[1]
NewSampleFile <- args[2]
DataDir <- args[3]

SampleDF <- read.table(OldSampleFile, header=T, sep="\t", stringsAsFactors=F)
SampleDF <- unique(SampleDF[,c(1,3:ncol(SampleDF))])

SeqDepthVec <- c()
for (i in 1:nrow(SampleDF)) {
	samplefile <- paste0(DataDir, '/', SampleDF[i, 1], '/Merged_CIS_Validpairs.txt.gz')
	val <- as.integer(system(paste("zcat", samplefile, "| wc -l"), intern = TRUE))
	SeqDepthVec <- c(SeqDepthVec, val)
}

CN <- colnames(SampleDF)
SampleDF <- cbind.data.frame(data.frame(Sample=SampleDF[,1]), data.frame(seqdepth=SeqDepthVec), SampleDF[, 2:ncol(SampleDF)])
colnames(SampleDF) <- c(CN[1], "seqdepth", CN[2:length(CN)])

write.table(SampleDF, NewSampleFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)

