#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(factoextra)
library(optparse)

options(scipen = 10)
options(datatable.fread.datatable=FALSE)

Color_Types <- c('violet', 'deepskyblue', 'brown', 'orange', 'red', 'green', 'black', 'blue', 'cyan', 'yellow')

#================================================
option_list = list(
	make_option(c("--InpFile"), type="character", default=NULL, help="Input loop file (mastersheet of loops). Mandatory parameter."),
	make_option(c("--AnnotFile"), type="character", default=NULL, help="Reference Annotation file for the given set of donors. Mandatory parameter."),
  	make_option(c("--colDonor"), type="integer", action="store", default=1, help="Column number in the reference annotation file containing donor ID. Default = 1 means the first column contains the donor ID."),
  	make_option(c("--colListAnnot"), type="character", default=NULL, help="Comma or colon separated list of integers which are the column numbers containing specific annotations (categories) of samples with respect to which the PCA plot will be generated. Default = NULL. Mandatory parameter."),
  	make_option(c("--onlyDonorNum"), type="integer", action="store", default=0, help="1 means only the donor numbers need to be exracted from the input loop file, without their run ID etc. Default = 0"),
  	make_option(c("--Outfileprefix"), type="character", default=NULL, help="Output file name prefix string")	
);

parser <- OptionParser(option_list=option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
args <- arguments$args

#============================
# this function generates various plots and statistics for PCA
#============================
Dump_Plot_PCA <- function(inpPCARes, outprefix, GroupVec, ColorVec, plotWidth=7, plotHeight=6) {

	cat(sprintf("\n\n ============= within function Dump_Plot_PCA ======= \n\n"))

	# eigen values
	# three columns: 1) eigenvalues 2) % of variances 3) cumulative % of variance
	eig <- factoextra::get_eig(inpPCARes)

	# contribution of individuals (replicates)
	ind <- factoextra::get_pca_ind(inpPCARes)

	cat(sprintf("\n\n\n Dimension of inpPCARes -- %s  X  %s \n\n\n", nrow(inpPCARes), ncol(inpPCARes)))	
	plotfile <- paste0(outprefix,'_viz_pca_ind_Plot_LATEST.pdf')
	
	factoextra::fviz_pca_ind(inpPCARes, pointsize = 2, col.ind = GroupVec, repel = TRUE, invisible="quali", select.ind=rownames(inpPCARes)) + scale_fill_manual(values=ColorVec) + scale_shape_manual(values=rep(1, length(ColorVec))) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	ggsave(plotfile, width=plotWidth, height=plotHeight)

	plotfile <- paste0(outprefix,'_viz_pca_ind_Plot_LATEST_2.pdf')
	factoextra::fviz_pca_ind(inpPCARes, pointsize = 2, geom=c("point"), col.ind = GroupVec, repel = TRUE, invisible="quali", select.ind=rownames(inpPCARes)) + scale_fill_manual(values=ColorVec) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	ggsave(plotfile, width=plotWidth, height=plotHeight)

}	# end function

##===================
## main code
##===================

## parse input parameters
inploopfile <- opt$InpFile
if (is.null(inploopfile)) {
	cat(sprintf("\n No input file is provided - check the option --InpFile - quit !! "))
	return()
}
annotfile <- opt$AnnotFile
if (is.null(annotfile)) {
	cat(sprintf("\n No input file is provided - check the option --AnnotFile - quit !! "))
	return()
}
colDonor <- as.integer(opt$colDonor)
if (is.null(opt$colListAnnot)) {
	cat(sprintf("\n List of columns containing the donor specific annotations (categories) is not provided - check the option --colListAnnot - quit !! "))
} else {
	colListAnnot <- as.integer(unlist(strsplit(opt$colListAnnot, "[,:]")))
}
outprefixfilename <- opt$Outfileprefix
onlyDonorNum <- as.integer(opt$onlyDonorNum)

cat(sprintf("\n\n Parameters --- \n inploopfile : %s \n annotfile : %s \n colDonor : %s \n colListAnnot : %s \n outprefixfilename : %s \n onlyDonorNum : %s ", inploopfile, annotfile, colDonor, paste(colListAnnot, collapse=" "), outprefixfilename, onlyDonorNum))

## input directory containing the loops
BaseInpDir <- dirname(inploopfile)

## read input loops
inploopdata <- data.table::fread(inploopfile, header=T)

## read the reference annotations
if (tools::file_ext(annotfile) == "csv") {
	annotdata <- read.csv(annotfile)
} else {
	annotdata <- data.table::fread(annotfile, header=T)
}

##==== get the donor names
##==== and also extract the feature specific columns
Lbls <- colnames(inploopdata)[7:ncol(inploopdata)]
if (onlyDonorNum == 0) {
	DonorNames <- unique(as.vector(sapply(Lbls, function(s) paste(head(strsplit(s, "_")[[1]],-1), collapse="_"))))
} else {
	DonorNames <- unique(as.vector(sapply(Lbls, function(s) strsplit(s, "_")[[1]][1])))
}
numfield_per_donor <- (length(Lbls) / length(DonorNames))
cat(sprintf("\n Number of donors : %s Fields per donor : %s ", length(DonorNames), numfield_per_donor))
cat(sprintf("\n\n Donor names : %s ", paste(DonorNames, collapse=" ")))

rawCC_ColList <- seq(7, ncol(inploopdata), by=numfield_per_donor)
qval_ColList <- rawCC_ColList + (numfield_per_donor - 1)
if (numfield_per_donor == 3) {
	expCC_ColList <- rawCC_ColList + 1 
}

m <- match(DonorNames, annotdata[, colDonor])
idx_s <- which(!is.na(m))
idx_t <- m[!is.na(m)]
cat(sprintf("\n\n Matched annotation -- number of donors : %s number of matched donors : %s ", length(DonorNames), length(idx_s)))

#===============
# PCA using raw contact counts
#===============
outdir <- paste0(BaseInpDir, '/out_PCA_RawCC')
system(paste("mkdir -p", outdir))

OutDF <- inploopdata[, rawCC_ColList]
colnames(OutDF) <- DonorNames	

## transform the feature vector so that rows indicate individual replicates and columns indicate different features
t_OutDF <- t(OutDF)
rownames(t_OutDF) <- colnames(OutDF)

## compute PCA
res.pca <- prcomp(t_OutDF, scale = FALSE)

## display the PCA according to different annotations
for (i in 1:length(colListAnnot)) {
	annotcol <- colListAnnot[i]
	annotname <- colnames(annotdata)[annotcol]
	cat(sprintf("\n PCA using rawCC - annotation : %s ", annotname))
	annotvec <- annotdata[idx_t, annotcol]
	if (!is.null(outprefixfilename)) {
		Dump_Plot_PCA(res.pca, paste0(outdir, '/', outprefixfilename, '_', annotname), as.factor(annotvec), Color_Types[1:length(unique(annotvec))])
	} else {
		Dump_Plot_PCA(res.pca, paste0(outdir, '/PCA_', annotname), as.factor(annotvec), Color_Types[1:length(unique(annotvec))])
	}
}

#===============
# PCA using q-value
#===============	
outdir <- paste0(BaseInpDir, '/out_PCA_QVal')	
system(paste("mkdir -p", outdir))

OutDF <- inploopdata[, qval_ColList]
colnames(OutDF) <- DonorNames	

## create data frame of q-value (-log10 scale)
for (i in (1:ncol(OutDF))) {
	currdf <- UtilRPckg::CustomLog(OutDF[,i], base=10, mult=-1)
	if (i == 1) {
		LogQDF <- currdf
	} else {
		LogQDF <- cbind.data.frame(LogQDF, currdf)
	}
}
colnames(LogQDF) <- DonorNames

## transform the feature vector so that rows indicate individual replicates and columns indicate different features
t_OutDF <- t(LogQDF)
rownames(t_OutDF) <- colnames(LogQDF)

## compute PCA
res.pca <- prcomp(t_OutDF, scale = FALSE)

## display the PCA according to different annotations
for (i in 1:length(colListAnnot)) {
	annotcol <- colListAnnot[i]
	annotname <- colnames(annotdata)[annotcol]
	cat(sprintf("\n PCA using q-value - annotation : %s ", annotname))
	annotvec <- annotdata[idx_t, annotcol]
	if (!is.null(outprefixfilename)) {
		Dump_Plot_PCA(res.pca, paste0(outdir, '/', outprefixfilename, '_', annotname), as.factor(annotvec), Color_Types[1:length(unique(annotvec))])
	} else {
		Dump_Plot_PCA(res.pca, paste0(outdir, '/PCA_', annotname), as.factor(annotvec), Color_Types[1:length(unique(annotvec))])
	}
}

#===============
# PCA using RawCC / expCC
#===============	
if (numfield_per_donor == 3) {
	outdir <- paste0(BaseInpDir, '/out_PCA_RawCC_vs_expCC')
	system(paste("mkdir -p", outdir))
	for (i in 1:length(rawCC_ColList)) {
		if (i == 1) {
			RawVSExpCCDF <- UtilRPckg::DivideVectors(inploopdata[, rawCC_ColList[i]], inploopdata[, expCC_ColList[i]])
		} else {
			RawVSExpCCDF <- cbind.data.frame(RawVSExpCCDF, UtilRPckg::DivideVectors(inploopdata[, rawCC_ColList[i]], inploopdata[, expCC_ColList[i]]))
		}
	}
	colnames(RawVSExpCCDF) <- DonorNames

	## transform the feature vector so that rows indicate individual replicates and columns indicate different features
	t_OutDF <- t(RawVSExpCCDF)
	rownames(t_OutDF) <- colnames(RawVSExpCCDF)

	## compute PCA
	res.pca <- prcomp(t_OutDF, scale = FALSE)

	## display the PCA according to different annotations
	for (i in 1:length(colListAnnot)) {
		annotcol <- colListAnnot[i]
		annotname <- colnames(annotdata)[annotcol]
		cat(sprintf("\n PCA using RawCC / expCC - annotation : %s ", annotname))
		annotvec <- annotdata[idx_t, annotcol]
		if (!is.null(outprefixfilename)) {
			Dump_Plot_PCA(res.pca, paste0(outdir, '/', outprefixfilename, '_', annotname), as.factor(annotvec), Color_Types[1:length(unique(annotvec))])
		} else {
			Dump_Plot_PCA(res.pca, paste0(outdir, '/PCA_', annotname), as.factor(annotvec), Color_Types[1:length(unique(annotvec))])
		}
	}
}

