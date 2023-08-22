#!/usr/bin/env Rscript

##=========================
# IQTL code - sourya - for RASQUAL
# script to create offset and covariates from merged input Y.txt file 
##=========================
library(data.table)
library(ggplot2)

options(scipen = 10)
options(datatable.fread.datatable=FALSE)

args <- commandArgs(trailingOnly = TRUE)

RESULTBASEDIR <- as.character(args[1])
RASQUAL_src_dir <- as.character(args[2])
MasterSheetFile <- as.character(args[3])
DonorAnnotFile <- as.character(args[4])
CovariateCols <- as.integer(unlist(strsplit(args[5], "[,:]")))
if (length(args) > 5) {
    nchunk_per_chr <- as.integer(args[6])
} else {
    nchunk_per_chr <- 15
}

## include RASQUAL package specific randomization source code
source(paste0(RASQUAL_src_dir, "/R/randomize.R"))

## output directory
OUTDIR <- paste0(RESULTBASEDIR, '/RASQUAL_INPUTS')
system(paste("mkdir -p", OUTDIR))

Parameter_File <- paste0(OUTDIR, '/parameters.txt')
fp_Param <- file(Parameter_File, "w")
outtext <- paste0("\n\n ==>> in function calc_offset_covariates.R -- list of parameters \n Sample specific Annotation File : ", DonorAnnotFile, "\n Covariate columns : ", paste(CovariateCols, collapse=","), "\n FitHiChIP master sheet file : ", MasterSheetFile, "\n nchunk_per_chr : ", nchunk_per_chr)
writeLines(outtext, con=fp_Param, sep="\n")
close(fp_Param)

raw_count_file <- paste0(OUTDIR, '/Y.txt')
rasq_offset_file <- paste0(OUTDIR, '/K.txt')
rasq_covs_file <- paste0(OUTDIR, '/X.txt')
chunkfile <- paste0(OUTDIR, '/chunk_info.txt')

## read the donor annotation data
if (tools::file_ext(DonorAnnotFile) == "csv") {
    AnnotData <- read.csv(DonorAnnotFile)
} else {
    AnnotData <- read.table(DonorAnnotFile, header=T, sep="\t", stringsAsFactors=F)
}

##============
## generate the raw counts from the input FitHiChIP master sheet file
##============
## region information is put as chr:start1-end2 to identify both interacting bins
system(paste0("awk \'{if (NR>1) {printf \"%s:%s-%s\",$1,$2,$6; for (i=7; i<=(NF-2); i=i+3) {printf \"\t%s\",$i}; printf \"\\n\"}}\' ", MasterSheetFile, " > ", raw_count_file))

##============
##===== if UseNormCC = 1, normalize the contact count according to the sequencing depth
##============
UseNormCC <- 0

if (UseNormCC == 1) {
    YData <- data.table::fread(raw_count_file, header=F, sep="\t", stringsAsFactors=F)
    seqdepthvec <- as.integer(AnnotData[, SeqDepthCol])
    if (0) {
        ## this was the first approach
        ## unfortunately, normalizing like this does not produce significant output in the population specific analysis
        ## normalize the raw counts to the count per million reads
        for (colidx in 2:ncol(YData)) {        
            YData[, colidx] <- as.numeric((as.numeric(YData[, colidx]) * 1000000) / seqdepthvec[(colidx-1)])
        }
    }
    if (1) {
        ## second approach
        ## use the maximum sequencing depth
        ## and normalize each sample (and their contact counts) according to the maximum sequencing depth
        max_seqdepth <- max(seqdepthvec)
        for (colidx in 2:ncol(YData)) { 
            curr_seqdepth <- seqdepthvec[(colidx-1)]
            YData[, colidx] <- as.numeric((as.numeric(YData[, colidx]) * max_seqdepth) / curr_seqdepth)
        }
    }
    ## write back the normalized data
    write.table(YData, raw_count_file, row.names=F, col.names=F, sep="\t", quote=F, append=F)  
}

##============
##======= compute Sample specific offset terms (K.txt) from the count table. 
##======= use the script makeOffset.R in the R directory.
##============
setwd(RASQUAL_src_dir)
system(paste0("R --vanilla --quiet --args ", raw_count_file, " < ", RASQUAL_src_dir, "/R/makeOffset.R > ", OUTDIR, "/log1"))
## output file is stored in paste0(RASQUAL_src_dir, "/data/your.K.txt")
## so copy that in the target rasq_offset_file
system(paste0("cp ", RASQUAL_src_dir, "/data/your.K.txt ", rasq_offset_file))

##============
## compute covariates (output file: X.txt)
##============

##======= step 1: use default covariates from RASQUAL

## read the raw count (Y.txt) and the generated offset file (K.txt)
Y <- read.table(raw_count_file, as.is=T); 
fid=Y[[1]]; 
Y=as.matrix(Y[,-1]); 
n <- ncol(Y)

K <- as.matrix(read.table(rasq_offset_file, as.is=T)[,-1])

## fpkm calculation
fpkm <- t(t(Y/K+1)/apply(Y/K,2,sum))*1e6 #  /len*1e9

## Singular value decomposition
fpkm.svd <- svd((log(fpkm)-apply(log(fpkm),1,mean))/apply(log(fpkm),1,sd))
fpkm.svd.r <- svd(randomize((log(fpkm)-apply(log(fpkm),1,mean))/apply(log(fpkm),1,sd)))

## Covariate selection
sf <- log(apply(Y,2,sum))
covs <- fpkm.svd$v[,1:sum(fpkm.svd$d[-n]>fpkm.svd.r$d[-n])]
df.covs <- data.frame(covs)
if (cor(sf,covs[,1])^2 < 0.9) {
    covs <- cbind(sf, covs)
}
names(df.covs) <- paste0("PC", 1:ncol(df.covs))

##======= step 2: add input covariates - if provided
if (length(CovariateCols) > 0) {
    Inp_AnnotData <- data.frame(AnnotData[, CovariateCols])
    ## append the input covariates in the main list
    df.covs <- cbind(df.covs, Inp_AnnotData)
    ## Plotting
    size.text <- 10
    var.plot1 <- "PC1"
    var.plot2 <- "PC2"
    for (var.color in colnames(Inp_AnnotData)) {
        currplotfile <- paste0(OUTDIR, "/PCA_", var.color,".pdf")
        ggd <- ggplot(df.covs, aes_string( x=var.plot1, y=var.plot2, color=var.color)) + geom_point() + theme_classic() + theme(axis.text.x=element_text(size=size.text, angle=45, hjust=1), axis.line=element_line(size=0.85, colour="black"), axis.ticks=element_line(colour="black", size=0.85), axis.ticks.length=unit(.25, "cm"), axis.title.x=element_text(size=size.text), legend.position="right", plot.title=element_text(hjust=0.5, size=size.text), axis.text.y=element_text(size=size.text), axis.title.y=element_text(size=size.text)) + labs(x=var.plot1,y=var.plot2)
        ggsave(filename = currplotfile, plot = ggd, width = 6, height = 6)
    }
    ## convert the input covariates to numeric 
    ## before creating the final list of covariates
    for (i in 1:ncol(Inp_AnnotData)) {
        ## if there is a problem in converting the categories to numeric
        ## first force them to be a factor, and then convert to numeric
        z <- as.numeric(Inp_AnnotData[, i])
        if (is.na(z[1])) {
            z <- as.numeric(as.factor(Inp_AnnotData[, i]))
        }
        Inp_AnnotData[, i] <- z
    }
    ## prepare the final list of covariates
    covs <- cbind(Inp_AnnotData, covs)    
}
# Write covariates
write.table(covs, rasq_covs_file, row.names=F, col.names=F, sep="\t", quote=F, append=F)

##============
## convert the covariate files into binary format for RASQUAL by using the same R script
##============
setwd(RASQUAL_src_dir)
system(paste0("R --vanilla --quiet --args ", raw_count_file, " ", rasq_offset_file, " ", rasq_covs_file," < ", RASQUAL_src_dir, "/R/txt2bin.R > ", OUTDIR, "/log3"))

##============
## create the chunks per chromosome
##============
cntdata <- data.table::fread(raw_count_file, header=F)
vec.chr <- sapply(cntdata[,1], function(X) unlist(strsplit(X, ":"))[1])
chrList <- unique(vec.chr)
bool_DF <- FALSE
for (i in 1:length(chrList)) {
    currchr <- chrList[i]
    cat(sprintf("\n creating chunks - considering chromosome : %s ", currchr))
    x <- sort(which(vec.chr == currchr))
    y <- split(x, cut(x,nchunk_per_chr, labels=FALSE))
    for (j in 1:nchunk_per_chr) {
        currDF <- data.frame(chr=currchr, chunk=j, start=min(y[[j]]), end=max(y[[j]]))
        if (bool_DF == FALSE) {
            finalDF <- currDF
            bool_DF <- TRUE
        } else {
            finalDF <- rbind.data.frame(finalDF, currDF)
        }
    }
}
write.table(finalDF, chunkfile, row.names=F, col.names=T, sep="\t", quote=F, append=F)

