#!/usr/bin/env Rscript

#===========================
## 1) Merges donor-wise HiChIP alignments (one donor may have multiple samples, depending on the sequencing runs)
## from HiC-pro output
#===========================

library(data.table)
options(scipen = 10)
options(datatable.fread.datatable=FALSE)

args <- commandArgs(trailingOnly = TRUE)

DonorListFile <- args[1]
HiChIPDataDir <- args[2]
IQTLResDir <- args[3]
RefGenome <- args[4]
TargetDonor <- as.character(args[5])
Targetchr <- as.character(args[6])

##==== read the complete donor list (annotated)
Donor_DF <- read.table(DonorListFile, header=T, sep="\t", stringsAsFactors=F)
DonorList <- unique(Donor_DF[,1])

AlignDir <- paste0(IQTLResDir, '/HiChIP_Complete_Alignments/', TargetDonor, '/chrwise')
system(paste("mkdir -p", AlignDir))

outalignfile <- paste0(AlignDir, '/merged_HiChIP_', Targetchr, '.bam')
sortalignfile <- paste0(AlignDir, '/merged_HiChIP_', Targetchr, '_sorted.bam')

samplelist <- Donor_DF[which(Donor_DF[,1] == TargetDonor), 2]

for (i in 1:length(samplelist)) {
	cat(sprintf("\n Processing sample : %s ", samplelist[i]))
	HiCProAlignFile <- paste0(HiChIPDataDir, '/', samplelist[i], '/HiCPro/bowtie_results/bwt2/rawdata/inpdata_', RefGenome, '.bwt2pairs.bam')
	
	## temporary BAM file
	## storing the contents for the current chromsome
	tempfile <- paste0(AlignDir, '/temp_chr_', Targetchr, '_sample_', i, '.bam')
	if (0) {
		## this command is most efficient to retrieve BAM reads for a given chromosome
		## however, the input HiCPro generated files are not sorted and indexed
		## so this commmand does not work
		system(paste("samtools view -h -b -o", tempfile, HiCProAlignFile, Targetchr))
	} else {
		## we use custom AWK script to extract the CIS reads for the input chromosome
		system(paste0("samtools view -h ", HiCProAlignFile, " | awk \'((substr($1,1,1) == \"@\") || (($3 == \"", Targetchr, "\") && (($3 == $7) || ($7 == \"=\"))))\' - | samtools view -bh -o ", tempfile, " - "))
	}

	## merge these extracted bam files
	if ( i == 1 ) {
		cmdstr <- paste0('samtools merge ', outalignfile, ' ', tempfile)
	} else {
		cmdstr <- paste0(cmdstr, ' ', tempfile)
	}
}
## merge the alignments (only for the specified chromosome)
system(paste(cmdstr))

system(paste("samtools sort -o", sortalignfile, " -@ 8 ", outalignfile))
system(paste("samtools index", sortalignfile))
system(paste("rm", outalignfile))

## delete the temporary files
for (i in 1:length(samplelist)) {
	tempfile <- paste0(AlignDir, '/temp_chr_', Targetchr, '_sample_', i, '.bam')
	system(paste("rm", tempfile))
}

