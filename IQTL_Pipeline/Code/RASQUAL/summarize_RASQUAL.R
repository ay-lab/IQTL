#!/usr/bin/env Rscript

##=========================
## IQTL code - sourya - for RASQUAL
## Script for summarizing RASQUAL results
## uses default model + population specific and allele specific model results
## implements paired t-test of the allele specific reads
## finally, filters the genotype specific contact counts (and significant results) 
## according to the strict ascending or descending trends
##=========================

library(data.table)
options(scipen = 10)
options(datatable.fread.datatable=FALSE)

args <- commandArgs(TRUE)

## Input arguments
currchr <- as.character(args[1])
base.out.path <- as.character(args[2])
MasterSheetFile <- as.character(args[3])
DonorAnnotFile <- as.character(args[4])
FDR_FitHiChIP <- as.numeric(args[5])
FDR_IQTL <- as.numeric(args[6])
chunk_ID <- as.integer(args[7])

## sequencing depth column - in the "Sample.txt" file
SeqDepthCol <- 2

## read the donor annotation data
if (tools::file_ext(DonorAnnotFile) == "csv") {
    AnnotData <- read.csv(DonorAnnotFile)
} else {
    AnnotData <- read.table(DonorAnnotFile, header=T, sep="\t", stringsAsFactors=F)
}    
SeqDepthVec <- AnnotData[, SeqDepthCol]    

## output path for three different RASQUAL models
out.path.default <- paste0(base.out.path, '/RASQUAL_Default/', currchr)
out.path.pop <- paste0(base.out.path, '/RASQUAL_Pop/', currchr)
out.path.AS <- paste0(base.out.path, '/RASQUAL_AS/', currchr)

## output directory to contain the combined statistics from these three models
summary.out.path.combined <- paste0(base.out.path, '/RASQUAL_Combined/', currchr, '/', chunk_ID)
system(paste("mkdir -p", summary.out.path.combined))

## output log file
outlogfile <- paste0(summary.out.path.combined, '/out.log')
sink(outlogfile)

cat("Chromosome: ", currchr, "\n")
cat("Base Output path: ", base.out.path, "\n")
cat("FDR threshold for IQTL: ", FDR_IQTL, "\n")
cat("Master sheet file (all FitHiChIP contacts) : ", MasterSheetFile, "\n")
cat("FDR threshold for FitHiChIP : ", FDR_FitHiChIP, "\n")
cat("chunk ID : ", chunk_ID, "\n")
cat("Donor list annotation file : ", DonorAnnotFile, "\n")

## column names of RASQUAL output text files
names.dat <- c("BinID", "SNPID", "Chromosome", "SNP.position", "Ref", "Alt", "Allele.Freq",
               "HWE.Chi.square", "Imputation.quality", "log10.BH.Qvalue", "Chi.Square.LR", "Effect.size",
               "Delta", "Phi", "Overdispersion", "SNPID.region", "N.Feature.SNPS", "N.Tested.SNPs", "N.iterations.null",
               "N.iterations.alt","Random.ties", "Log.likelihood", "Convergence", "Sq.corr.fSNPs", "Sq.corr.rSNPs")

##============
## navigate through individual RASQUAL models and corresponding results
##============

## temporary file to contain combined population outputs, per chunk and per loop
tempfile_RASQUAL_pop <- paste0(summary.out.path.combined, '/temp_RASQUAL_out_pop.txt')

## temporary file to contain combined allele specific outputs, per chunk and per loop
tempfile_RASQUAL_AS <- paste0(summary.out.path.combined, '/temp_RASQUAL_out_AS.txt')

## temporary file to contain default RASQUAL outputs, per chunk and per loop
tempfile_RASQUAL_def <- paste0(summary.out.path.combined, '/temp_RASQUAL_out_default.txt')

##==== RASQUAL significant output + loops - 
## best SNP along with the genotype specific contact count + allele specifc reads information
RASQUAL_out_sig_loops_donor_genotype_file <- paste0(summary.out.path.combined, '/RASQUAL_out.txt')
bool_sig_loops_donor_genotype_RASQUAL <- FALSE

RASQUAL_out_sig_loops_donor_genotype_NormCC_file <- paste0(summary.out.path.combined, '/RASQUAL_out_NormCC.txt')

## check this condition - sourya
## if this file exists, then the computation is already done 
if (file.exists(RASQUAL_out_sig_loops_donor_genotype_file) == FALSE) {

    ##========= extract mastersheet contacts for the current chromosome
    ##========= extract bin ID, similar to the RASQUAL format ($1":"$2"-"$6)
    ##========= then extract interacting bins, rawCC and q-values (in alternate columns)
    ##========= also get the final two fields, namely LoopCnt and SigLoopCnt
    dumpedMasterSheetFile <- paste0(summary.out.path.combined, '/dumped_mastersheet_', currchr, '.bed')
    if (file.exists(dumpedMasterSheetFile) == FALSE) {
        system(paste0("awk -F\'\t\' \'{if ((NR==1) || (($1==\"", currchr, "\") && ($4==\"", currchr, "\"))) {printf \"%s:%s-%s\",$1,$2,$6; for (i=1; i<=6; i=i+1) {printf \"\\t%s\",$i}; for (i=7; i<=(NF-2); i=i+3) {printf \"\\t%s\\t%s\",$i,$(i+2)}; printf \"\\t%s\\t%s\\n\",$(NF-1),$NF}}\' ", MasterSheetFile, " > ", dumpedMasterSheetFile))
    }
    loopDF <- data.table::fread(dumpedMasterSheetFile, header=T)
    CN <- colnames(loopDF)
    colnames(loopDF) <- c('BinID', 'chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', CN[8:length(CN)])
    cat(sprintf("\n --- nrow loopDF : %s ", nrow(loopDF)))

    ## path of different chunk outputs for different RASQUAL models
    curr_chunk_path.default <- paste0(out.path.default, '/', chunk_ID)
    curr_chunk_path.pop <- paste0(out.path.pop, '/', chunk_ID)
    curr_chunk_path.AS <- paste0(out.path.AS, '/', chunk_ID)

    ## if current chunk is not processed, continue
    if ((dir.exists(curr_chunk_path.default) == FALSE) | (dir.exists(curr_chunk_path.pop) == FALSE) | (dir.exists(curr_chunk_path.AS) == FALSE)) {
        next
    }

    ## re-estimated allele specific reads are always within the default output folder
    reEstASReadPath <- paste0(curr_chunk_path.default, '/SNPs_reEst_ASReads')

    ## get the index of interaction - minimum and maximum - within this folder
    ## the RASQUAL outputs are stored from t_minidx_1.txt to t_maxidx_2.txt files
    min_loop_num <- as.integer(system(paste0("ls -l ", reEstASReadPath, "/ | awk \'{print $NF}\' - | grep \"re_estimated_AS_Reads.vcf.gz\" | awk -F\'[_]\' \'{print $3}\' - | sort -k1,1n | uniq | head -n 1 "), intern = TRUE))
    max_loop_num <- as.integer(system(paste0("ls -l ", reEstASReadPath, "/ | awk \'{print $NF}\' - | grep \"re_estimated_AS_Reads.vcf.gz\" | awk -F\'[_]\' \'{print $3}\' - | sort -k1,1nr | uniq | head -n 1 "), intern = TRUE))
    cat(sprintf("\n\n\n ===>>>>>> processing chunk : %s min_loop_num : %s max_loop_num : %s ", chunk_ID, min_loop_num, max_loop_num))

    ## process individual loops within each chunk
    for (loopidx in seq(min_loop_num, max_loop_num)) {
        cat(sprintf("\n --- processing loop index : %s ", loopidx))

        ## genotype file for the current chunk and current interaction
        ## containing re-estimated AS reads
        vcf.path <- paste0(reEstASReadPath, '/out_SNPs_', loopidx, '_re_estimated_AS_Reads.vcf.gz')

        ## get the number of donors from the genotype file (vcf file)
        num_donors <- as.integer(system(paste("zcat", vcf.path, "| tail -n 1 | awk \'{print NF}\' -"), intern = TRUE)) - 9
        cat(sprintf("\n\n **** very first iteration - number of donors : %s ", num_donors))

        ##============== extract the default RASQUAL model outputs
        deffile1 <- paste0(curr_chunk_path.default, '/t_', loopidx, '_1.txt')
        deffile2 <- paste0(curr_chunk_path.default, '/t_', loopidx, '_2.txt')
        system(paste("cat", deffile1, deffile2, ">", tempfile_RASQUAL_def))            
        fin.dat.def <- data.table::fread(tempfile_RASQUAL_def, header=F)
        if (nrow(fin.dat.def) == 0) {
            next
        }
        names(fin.dat.def) <- names.dat
        cat(sprintf("\n nrow fin.dat.def (considering both output files _1.txt and _2.txt) : %s ", nrow(fin.dat.def)))
        ## compute the p-value from "HWE.Chi.square" generated by RASQUAL
        ## assume degree of freedom = 1 
        ## (df = 1 is generally used for HW even when there are three phenotypes 
        ## because the expected can be calculated starting with one of the two alleles*), the critical value is 3.84.
        ## check https://www.biologysimulations.com/post/how-to-use-chi-squared-to-test-for-hardy-weinberg-equilibrium
        fin.dat.def$Chi.Square.LR.pValue.Def <- pchisq(fin.dat.def$Chi.Square.LR, df=1, lower.tail=FALSE)
        ## adjust q-values -  add additional fields at the end
        fin.dat.def$BH.Qvalue.Def <- 10**fin.dat.def$log10.BH.Qvalue
        ## extract selective fields
        fin.dat.def <- fin.dat.def[, c(1:8,11:12,22,(ncol(fin.dat.def)-1),ncol(fin.dat.def))]
        colnames(fin.dat.def) <- c('BinID', 'SNPID', 'SNPchr', 'SNPloc', 'ref', 'alt', 'Allele.Freq', 'HWE.Chi.square', 'Chi.Square.LR', 'Effect.size', 'Log.likelihood', 'Chi.Square.LR.pValue.Def', 'BH.Qvalue.Def')
        cat(sprintf("\n --- nrow fin.dat.def : %s ", nrow(fin.dat.def)))

        ##============== extract the population specific RASQUAL model outputs
        popfile1 <- paste0(curr_chunk_path.pop, '/t_', loopidx, '_1.txt')
        popfile2 <- paste0(curr_chunk_path.pop, '/t_', loopidx, '_2.txt')
        system(paste("cat", popfile1, popfile2, ">", tempfile_RASQUAL_pop))
        fin.dat.pop <- data.table::fread(tempfile_RASQUAL_pop, header=F)
        if (nrow(fin.dat.pop) > 0) {
            names(fin.dat.pop) <- names.dat
            cat(sprintf("\n nrow fin.dat.pop (considering both output files _1.txt and _2.txt) : %s ", nrow(fin.dat.pop)))
            fin.dat.pop$Chi.Square.LR.pValue.pop <- pchisq(fin.dat.pop$Chi.Square.LR, df=1, lower.tail=FALSE)
            fin.dat.pop$BH.Qvalue.pop <- 10**fin.dat.pop$log10.BH.Qvalue
            ## extract only the following fields: BinID, SNPID, and BH.Qvalue.pop
            fin.dat.pop <- fin.dat.pop[, c(1,2,ncol(fin.dat.pop))]
            colnames(fin.dat.pop) <- c('BinID', 'SNPID', 'BH.Qvalue.pop')
        }

        ##============== extract the allele specific model output
        ASfile1 <- paste0(curr_chunk_path.AS, '/t_', loopidx, '_1.txt')
        ASfile2 <- paste0(curr_chunk_path.AS, '/t_', loopidx, '_2.txt')            
        system(paste("cat", ASfile1, ASfile2, ">", tempfile_RASQUAL_AS))
        fin.dat.AS <- data.table::fread(tempfile_RASQUAL_AS, header=F)
        if (nrow(fin.dat.AS) > 0) {
            names(fin.dat.AS) <- names.dat
            cat(sprintf("\n nrow fin.dat.AS (considering both output files _1.txt and _2.txt) : %s ", nrow(fin.dat.AS)))
            fin.dat.AS$Chi.Square.LR.pValue.AS <- pchisq(fin.dat.AS$Chi.Square.LR, df=1, lower.tail=FALSE)
            fin.dat.AS$BH.Qvalue.AS <- 10**fin.dat.AS$log10.BH.Qvalue        
            ## extract only the following fields: BinID, SNPID, and BH.Qvalue.AS
            fin.dat.AS <- fin.dat.AS[, c(1,2,ncol(fin.dat.AS))]
            colnames(fin.dat.AS) <- c('BinID', 'SNPID', 'BH.Qvalue.AS')       
        }

        ##============== now select only the significant entries from the default model
        ##============== subject to the specified "FDR_IQTL"
        idx <- which(fin.dat.def$BH.Qvalue.Def < FDR_IQTL)
        if (length(idx) == 0) {
            cat(sprintf("\n !!! No significant entries (FDR < %s ) in the default model - proceed to the next loop ", FDR_IQTL))
            next
        }
        fin.dat.def.SIG <- fin.dat.def[idx, ]
        cat(sprintf("\n nrow fin.dat.def.SIG (significant entries with respect to the specified FDR threshold %s is ) : %s ", FDR_IQTL, nrow(fin.dat.def.SIG)))

        ##=========== merge the population and allele specific outputs
        ## population model
        if (nrow(fin.dat.pop) > 0) {
            fin.dat.def.SIG <- dplyr::left_join(fin.dat.def.SIG, fin.dat.pop)
            idx <- which(is.na(fin.dat.def.SIG$BH.Qvalue.pop))
            if (length(idx) > 0) {
                fin.dat.def.SIG[idx, ncol(fin.dat.def.SIG)] <- 1
            }
        } else {
            fin.dat.def.SIG$BH.Qvalue.pop <- rep(1, nrow(fin.dat.def.SIG))
        }
        ## allele model
        if (nrow(fin.dat.AS) > 0) {
            fin.dat.def.SIG <- dplyr::left_join(fin.dat.def.SIG, fin.dat.AS)
            idx <- which(is.na(fin.dat.def.SIG$BH.Qvalue.AS))
            if (length(idx) > 0) {
                fin.dat.def.SIG[idx, ncol(fin.dat.def.SIG)] <- 1
            }
        } else {
            fin.dat.def.SIG$BH.Qvalue.AS <- rep(1, nrow(fin.dat.def.SIG))
        }
        cat(sprintf("\n nrow fin.dat.def.SIG (after merging with population and allele specific outputs) : %s ", nrow(fin.dat.def.SIG)))

        ##====== merge with loop specific contact count information
        fin.dat.def.SIG.loopDF <- dplyr::inner_join(fin.dat.def.SIG, loopDF, by='BinID')
        cat(sprintf("\n --- nrow fin.dat.def.SIG.loopDF : %s ", nrow(fin.dat.def.SIG.loopDF)))

        ##========= now insert the loop-wise contact counts, FDR, and the allele specific reads
        ## column numbers storing the raw contact counts in the "fin.dat.def.SIG.loopDF" structure
        rawcc_cols <- seq(from=(ncol(fin.dat.def.SIG)+6+1), to=(ncol(fin.dat.def.SIG)+6+1 + (num_donors-1)*2), by=2)
        ## column numbers storing the q-values (FitHiChIP)
        qval_cols <- rawcc_cols + 1

        ##=== output matrix format - genotype specific number of donors, CC and FDR + allele specific counts
        SNPLoopGenoDF <- matrix('0', nrow=nrow(fin.dat.def.SIG.loopDF), ncol=18)
        colnames(SNPLoopGenoDF) <- c('Geno_0_numDonor', 'Geno_0_SigDonor', 'Geno_0_CC', 'Geno_0_FDR', 'Geno_0_MeanCC', 'Geno_1_numDonor', 'Geno_1_SigDonor', 'Geno_1_CC', 'Geno_1_FDR', 'Geno_1_MeanCC', 'Geno_2_numDonor', 'Geno_2_SigDonor', 'Geno_2_CC', 'Geno_2_FDR', 'Geno_2_MeanCC', 'AS_hetero_ref_Cnt', 'AS_hetero_alt_Cnt', 'ASRead_t_Test_pval')

        ##=== if sequencing depth information is provided, we'll also create a file 
        ##=== with sequencing depth normalized contact count
        SNPLoopGenoDF_NormCC <- matrix('0', nrow=nrow(fin.dat.def.SIG.loopDF), ncol=18)
        colnames(SNPLoopGenoDF_NormCC) <- colnames(SNPLoopGenoDF)

        ## also read the allele specific read file (vcf)
        ## sourya - note the skip argument
        ## it reads from the #CHR line, and skips the other headers        
        vcfdata <- data.table::fread(vcf.path, header=T, skip="#CHR")

        ## sequential processing
        for (i in 1:nrow(fin.dat.def.SIG.loopDF)) {

            snpID <- fin.dat.def.SIG.loopDF[i, 2]
            binID <- fin.dat.def.SIG.loopDF[i, 1]    
            snpLoc <- fin.dat.def.SIG.loopDF[i, 4]
            cat(sprintf("\n ===>> filling genotype + Allele specific contacts - processing entry index : %s inpchr : %s snpID : %s binID : %s snpLoc : %s ", i, currchr, snpID, binID, snpLoc))
    
            ## extract current SNP and corresponding allele specific read information
            tempGenoData <- vcfdata[which(((vcfdata[,1] == currchr) & (vcfdata[, 2] == snpLoc)) | (vcfdata[,3] == snpID)), ]
            refallele <- tempGenoData[1, 4]
            altallele <- tempGenoData[1, 5]
            cat(sprintf(" -- reference allele : %s alternate allele : %s ", refallele, altallele))

            ## donor numbers (1 to N) which are heterozygous
            homo_donor_seq <- c(); hetero_donor_seq <- c(); alt_homo_donor_seq <- c()
            ## stores the corresponding allele specific counts (w.r.t HiChIP data)
            hetero_ref_allele_AS_counts <- c(); hetero_alt_allele_AS_counts <- c()

            ## starting from the "startcol", the tempGenoData stores genotype information of individual donors
            ## fill the donor sequence and allele specific count vectors
            startcol <- 10
            for (colidx in startcol:ncol(tempGenoData)) {
                currgeno <- tempGenoData[1, colidx]
                ## currgeno has the format x|y:z,w
                ## where, x|y can be 0|0, 1|1, 0|1 or 1|0 - represents the genotype
                ## 0: reference, 1: alternate
                ## z is the read count for the left side allele 
                ## y is the read count for the right side allele (nonzero only for heterozygous)
                genotypeval <- strsplit(currgeno, ":")[[1]][1]
                v <- strsplit(currgeno, ":")[[1]][2]
                allele1_val <- as.integer(strsplit(v, ",")[[1]][1])
                allele2_val <- as.integer(strsplit(v, ",")[[1]][2])

                ## important - check the escape characters \\
                if (grepl("0\\|0", currgeno)) {
                    homo_donor_seq <- c(homo_donor_seq, (colidx-startcol+1))
                } else if (grepl("0\\|1", currgeno)) {
                    hetero_donor_seq <- c(hetero_donor_seq, (colidx-startcol+1))
                    ## 0|1: first (second) allele specific count corresponds to the reference (alternate) allele
                    hetero_ref_allele_AS_counts <- c(hetero_ref_allele_AS_counts, allele1_val)
                    hetero_alt_allele_AS_counts <- c(hetero_alt_allele_AS_counts, allele2_val)
                } else if (grepl("1\\|0", currgeno)) {
                    hetero_donor_seq <- c(hetero_donor_seq, (colidx-startcol+1))
                    ## 1|0: second (first) allele specific count corresponds to the reference (alternate) allele
                    hetero_ref_allele_AS_counts <- c(hetero_ref_allele_AS_counts, allele2_val)
                    hetero_alt_allele_AS_counts <- c(hetero_alt_allele_AS_counts, allele1_val)
                } else if (grepl("1\\|1", currgeno)) {
                    alt_homo_donor_seq <- c(alt_homo_donor_seq, (colidx-startcol+1))
                }
            }

            ###============== group values according to genotype
            ## genotype = 0 (reference homozygous)
            SNPLoopGenoDF[i, 1] <- length(homo_donor_seq)
            SNPLoopGenoDF_NormCC[i, 1] <- length(homo_donor_seq)
            if (length(homo_donor_seq) > 0) {
                ccval <- fin.dat.def.SIG.loopDF[i, rawcc_cols[homo_donor_seq]]
                qval <- fin.dat.def.SIG.loopDF[i, qval_cols[homo_donor_seq]]
                SNPLoopGenoDF[i, 3] <- as.character(paste(ccval, collapse=":"))
                SNPLoopGenoDF[i, 4] <- as.character(paste(qval, collapse=":"))
                SNPLoopGenoDF[i, 2] <- length(which(qval < FDR_FitHiChIP))
                SNPLoopGenoDF[i, 5] <- mean(as.numeric(ccval))
                normccval <- ((fin.dat.def.SIG.loopDF[i, rawcc_cols[homo_donor_seq]] * 1000000) / SeqDepthVec[homo_donor_seq])
                SNPLoopGenoDF_NormCC[i, 3] <- as.character(paste(normccval, collapse=":"))
                SNPLoopGenoDF_NormCC[i, 4] <- as.character(paste(qval, collapse=":"))
                SNPLoopGenoDF_NormCC[i, 2] <- length(which(qval < FDR_FitHiChIP))
                SNPLoopGenoDF_NormCC[i, 5] <- mean(as.numeric(normccval))                        
            } else {
                SNPLoopGenoDF[i, 3] <- '-'
                SNPLoopGenoDF[i, 4] <- '-'
                SNPLoopGenoDF_NormCC[i, 3] <- '-'
                SNPLoopGenoDF_NormCC[i, 4] <- '-'
            }

            ## genotype = 1 (heterozygous)
            SNPLoopGenoDF[i, 6] <- length(hetero_donor_seq)
            SNPLoopGenoDF_NormCC[i, 6] <- length(hetero_donor_seq)
            if (length(hetero_donor_seq) > 0) {
                ccval <- fin.dat.def.SIG.loopDF[i, rawcc_cols[hetero_donor_seq]]
                qval <- fin.dat.def.SIG.loopDF[i, qval_cols[hetero_donor_seq]]
                SNPLoopGenoDF[i, 8] <- as.character(paste(ccval, collapse=":"))
                SNPLoopGenoDF[i, 9] <- as.character(paste(qval, collapse=":"))
                SNPLoopGenoDF[i, 7] <- length(which(qval < FDR_FitHiChIP))
                SNPLoopGenoDF[i, 10] <- mean(as.numeric(ccval))
                SNPLoopGenoDF[i, 16] <- as.character(paste(hetero_ref_allele_AS_counts, collapse=":"))
                SNPLoopGenoDF[i, 17] <- as.character(paste(hetero_alt_allele_AS_counts, collapse=":"))
                normccval <- ((fin.dat.def.SIG.loopDF[i, rawcc_cols[hetero_donor_seq]] * 1000000) / SeqDepthVec[hetero_donor_seq])
                SNPLoopGenoDF_NormCC[i, 8] <- as.character(paste(normccval, collapse=":"))
                SNPLoopGenoDF_NormCC[i, 9] <- as.character(paste(qval, collapse=":"))
                SNPLoopGenoDF_NormCC[i, 7] <- length(which(qval < FDR_FitHiChIP))
                SNPLoopGenoDF_NormCC[i, 10] <- mean(as.numeric(normccval))
                SNPLoopGenoDF_NormCC[i, 16] <- as.character(paste(hetero_ref_allele_AS_counts, collapse=":"))
                SNPLoopGenoDF_NormCC[i, 17] <- as.character(paste(hetero_alt_allele_AS_counts, collapse=":"))
            } else {
                SNPLoopGenoDF[i, 8] <- '-'
                SNPLoopGenoDF[i, 9] <- '-'  
                SNPLoopGenoDF[i, 16] <- '-'
                SNPLoopGenoDF[i, 17] <- '-'
                SNPLoopGenoDF_NormCC[i, 8] <- '-'
                SNPLoopGenoDF_NormCC[i, 9] <- '-'  
                SNPLoopGenoDF_NormCC[i, 16] <- '-'
                SNPLoopGenoDF_NormCC[i, 17] <- '-'
            }    

            ## genotype = 2 (alternate homozygous)
            SNPLoopGenoDF[i, 11] <- length(alt_homo_donor_seq)
            SNPLoopGenoDF_NormCC[i, 11] <- length(alt_homo_donor_seq)
            if (length(alt_homo_donor_seq) > 0) {
                ccval <- fin.dat.def.SIG.loopDF[i, rawcc_cols[alt_homo_donor_seq]]
                qval <- fin.dat.def.SIG.loopDF[i, qval_cols[alt_homo_donor_seq]]
                SNPLoopGenoDF[i, 13] <- as.character(paste(ccval, collapse=":"))
                SNPLoopGenoDF[i, 14] <- as.character(paste(qval, collapse=":"))
                SNPLoopGenoDF[i, 15] <- mean(as.numeric(ccval))
                SNPLoopGenoDF[i, 12] <- length(which(qval < FDR_FitHiChIP))
                normccval <- ((fin.dat.def.SIG.loopDF[i, rawcc_cols[alt_homo_donor_seq]] * 1000000) / SeqDepthVec[alt_homo_donor_seq])
                SNPLoopGenoDF_NormCC[i, 13] <- as.character(paste(normccval, collapse=":"))
                SNPLoopGenoDF_NormCC[i, 14] <- as.character(paste(qval, collapse=":"))
                SNPLoopGenoDF_NormCC[i, 15] <- mean(as.numeric(normccval))
                SNPLoopGenoDF_NormCC[i, 12] <- length(which(qval < FDR_FitHiChIP))
            } else {
                SNPLoopGenoDF[i, 13] <- '-'
                SNPLoopGenoDF[i, 14] <- '-'
                SNPLoopGenoDF_NormCC[i, 13] <- '-'
                SNPLoopGenoDF_NormCC[i, 14] <- '-'
            }

            ## process the allele specific read counts
            ## apply the paired t-Test for the allele specific reads (for the heterozygous donors)
            ## check for the error conditions
            if (length(hetero_ref_allele_AS_counts) <= 1) {
                ## case 1: number of heterozygous donors <= 1
                ## can't apply t-test
                SNPLoopGenoDF[i, 18] <- 'NA'
                SNPLoopGenoDF_NormCC[i, 18] <- 'NA'
            } else {
                ## case 2: if the number of unique elements in both vectors is 1
                ## or if the standard deviation of the difference of these vectors is 0
                ## then put a little trick
                ## increase the first element of first vector by some threshold
                ## increase the last element of second vector by some threshold
                ## to break the uniformity of these vectors (and thereby 0 standard deviation)
                if ((length(unique(hetero_ref_allele_AS_counts)) == 1) & (length(unique(hetero_alt_allele_AS_counts)) == 1)) {
                    hetero_ref_allele_AS_counts[1] <- hetero_ref_allele_AS_counts[1] + 0.1  # some threshold
                    hetero_alt_allele_AS_counts[length(hetero_alt_allele_AS_counts)] <- hetero_alt_allele_AS_counts[length(hetero_alt_allele_AS_counts)] + 0.1  # some threshold
                }
                if (sd(hetero_ref_allele_AS_counts - hetero_alt_allele_AS_counts) == 0) {
                    hetero_ref_allele_AS_counts[1] <- hetero_ref_allele_AS_counts[1] + 0.1  # some threshold
                    hetero_alt_allele_AS_counts[length(hetero_alt_allele_AS_counts)] <- hetero_alt_allele_AS_counts[length(hetero_alt_allele_AS_counts)] + 0.1  # some threshold
                }
                # pval_two_sided <- t.test(hetero_ref_allele_AS_counts, hetero_alt_allele_AS_counts, paired = TRUE, alternative = "two.sided")$p.value
                pval_greater <- t.test(hetero_ref_allele_AS_counts, hetero_alt_allele_AS_counts, paired = TRUE, alternative = "greater")$p.value
                pval_less <- t.test(hetero_ref_allele_AS_counts, hetero_alt_allele_AS_counts, paired = TRUE, alternative = "less")$p.value
                if (is.na(pval_greater) | is.na(pval_less)) {
                    SNPLoopGenoDF[i, 18] <- 'NA'
                    SNPLoopGenoDF_NormCC[i, 18] <- 'NA'
                } else {
                    t_test_pval <- min(as.numeric(c(pval_greater, pval_less)))
                    SNPLoopGenoDF[i, 18] <- t_test_pval
                    SNPLoopGenoDF_NormCC[i, 18] <- t_test_pval
                }                
            }
            ## print the vector
            cat(sprintf("\n ---> genotype processing - entry idx : %s SNP ID: %s bin / loop ID : %s complete vector - %s ", i, snpID, binID, paste(SNPLoopGenoDF[i, ], collapse=" | ")))

        }   # end loop - sequential processing

        ## construct the output matrix
        CN <- colnames(fin.dat.def.SIG.loopDF)
        outDF <- cbind.data.frame(fin.dat.def.SIG.loopDF[, c(1:(ncol(fin.dat.def.SIG)+6), ncol(fin.dat.def.SIG.loopDF)-1, ncol(fin.dat.def.SIG.loopDF))], SNPLoopGenoDF)
        colnames(outDF) <- c(CN[c(1:(ncol(fin.dat.def.SIG)+6), (length(CN)-1), length(CN))], colnames(SNPLoopGenoDF))
        outDF_NormCC <- cbind.data.frame(fin.dat.def.SIG.loopDF[, c(1:(ncol(fin.dat.def.SIG)+6), ncol(fin.dat.def.SIG.loopDF)-1, ncol(fin.dat.def.SIG.loopDF))], SNPLoopGenoDF_NormCC)
        colnames(outDF_NormCC) <- c(CN[c(1:(ncol(fin.dat.def.SIG)+6), (length(CN)-1), length(CN))], colnames(SNPLoopGenoDF_NormCC))

        if (bool_sig_loops_donor_genotype_RASQUAL == FALSE) {
            write.table(outDF, RASQUAL_out_sig_loops_donor_genotype_file, row.names=F, col.names=T, sep="\t", quote=F, append=F)
            write.table(outDF_NormCC, RASQUAL_out_sig_loops_donor_genotype_NormCC_file, row.names=F, col.names=T, sep="\t", quote=F, append=F)
            bool_sig_loops_donor_genotype_RASQUAL <- TRUE
        } else {
            write.table(outDF, RASQUAL_out_sig_loops_donor_genotype_file, row.names=F, col.names=F, sep="\t", quote=F, append=T)
            write.table(outDF_NormCC, RASQUAL_out_sig_loops_donor_genotype_NormCC_file, row.names=F, col.names=F, sep="\t", quote=F, append=T)
        }

    }   # end processing individual loops

    ## remove the dumped mastersheet file
    if (file.exists(dumpedMasterSheetFile) == TRUE) {
         system(paste("rm", dumpedMasterSheetFile))
    }

    if (file.exists(tempfile_RASQUAL_def) == TRUE) {
         system(paste("rm", tempfile_RASQUAL_def))
    }

    if (file.exists(tempfile_RASQUAL_pop) == TRUE) {
         system(paste("rm", tempfile_RASQUAL_pop))
    }

    if (file.exists(tempfile_RASQUAL_AS) == TRUE) {
         system(paste("rm", tempfile_RASQUAL_AS))
    }

}   # end file exist condition

##========================
## dump a few specific categories (based on significance from different models)
## we may use them for debugging later
##========================
if (0) {

    if (file.exists(RASQUAL_out_sig_loops_donor_genotype_file)) {
        DebugDir <- paste0(dirname(RASQUAL_out_sig_loops_donor_genotype_file), '/Debug')
        system(paste("mkdir -p", DebugDir))

        ## get entries which are significant in RASQUAL allele specific model
        ## not in the population model
        ## and also the p-value of t-test from allele specific read is not significant
        outfile1 <- paste0(DebugDir, '/Only_RASQUAL_AS_no_tTest.bed')
        system(paste0("awk \'(NR==1) || (($14>", FDR_IQTL, ") && ($15<", FDR_IQTL, ") && ($NF>", FDR_IQTL, "))\' ", RASQUAL_out_sig_loops_donor_genotype_file, " > ", outfile1))

        outfile2 <- paste0(DebugDir, '/Both_RASQUAL_AS_tTest.bed')
        system(paste0("awk \'(NR==1) || (($14>", FDR_IQTL, ") && ($15<", FDR_IQTL, ") && ($NF<", FDR_IQTL, "))\' ", RASQUAL_out_sig_loops_donor_genotype_file, " > ", outfile2))

        outfile3 <- paste0(DebugDir, '/Only_tTest_no_RASQUAL_AS.bed')
        system(paste0("awk \'(NR==1) || (($14>", FDR_IQTL, ") && ($15>", FDR_IQTL, ") && ($NF<", FDR_IQTL, "))\' ", RASQUAL_out_sig_loops_donor_genotype_file, " > ", outfile3))

        outfile4 <- paste0(DebugDir, '/Pop_Sig.bed')
        system(paste0("awk \'(NR==1) || ($14<", FDR_IQTL, ")\' ", RASQUAL_out_sig_loops_donor_genotype_file, " > ", outfile4))

        if (file.exists(RASQUAL_out_sig_loops_donor_genotype_NormCC_file)) {
            outfile1_norm <- paste0(DebugDir, '/Only_RASQUAL_AS_no_tTest_NormCC.bed')
            system(paste0("awk \'(NR==1) || (($14>", FDR_IQTL, ") && ($15<", FDR_IQTL, ") && ($NF>", FDR_IQTL, "))\' ", RASQUAL_out_sig_loops_donor_genotype_NormCC_file, " > ", outfile1_norm))

            outfile2_norm <- paste0(DebugDir, '/Both_RASQUAL_AS_tTest_NormCC.bed')
            system(paste0("awk \'(NR==1) || (($14>", FDR_IQTL, ") && ($15<", FDR_IQTL, ") && ($NF<", FDR_IQTL, "))\' ", RASQUAL_out_sig_loops_donor_genotype_NormCC_file, " > ", outfile2_norm))

            outfile3_norm <- paste0(DebugDir, '/Only_tTest_no_RASQUAL_AS_NormCC.bed')
            system(paste0("awk \'(NR==1) || (($14>", FDR_IQTL, ") && ($15>", FDR_IQTL, ") && ($NF<", FDR_IQTL, "))\' ", RASQUAL_out_sig_loops_donor_genotype_NormCC_file, " > ", outfile3_norm))

            outfile4_norm <- paste0(DebugDir, '/Pop_Sig_NormCC.bed')
            system(paste0("awk \'(NR==1) || ($14<", FDR_IQTL, ")\' ", RASQUAL_out_sig_loops_donor_genotype_NormCC_file, " > ", outfile4_norm))
        }
    }

}   # end dummy if

##====================
## now we filter the above file
## conditions:
## 1) all three genotypes should be present in at least 2 donors
## 2) for the population model significant entries, strict ascending or strict descending trend should be present (with respect to genotype specific changes in the contact count)
## 3) for the allele specific significant entries, it should be significant in RASQUAL AS model, plus ASRead trend should have significant p-value, plus the allele specific trend should be similar to the population specific trends (i.e. allele specific read variation should be similar to the genotype specific read variation)
##====================
if (file.exists(RASQUAL_out_sig_loops_donor_genotype_file)) {

    SNP_Loop_Data_RawCC <- data.table::fread(RASQUAL_out_sig_loops_donor_genotype_file, header=T)
    if (nrow(SNP_Loop_Data_RawCC) == 0) {
        return()
    }
    CN <- colnames(SNP_Loop_Data_RawCC)    
    if (file.exists(RASQUAL_out_sig_loops_donor_genotype_NormCC_file)) {
        SNP_Loop_Data_NormCC <- data.table::fread(RASQUAL_out_sig_loops_donor_genotype_NormCC_file, header=T)
        SNP_Loop_Data <- SNP_Loop_Data_NormCC
    } else {
        SNP_Loop_Data <- SNP_Loop_Data_RawCC
    }
    cat(sprintf("\n *** initial number of entries : %s ", nrow(SNP_Loop_Data)))

    ## columns storing the fields
    colidx_geno_0_numdonor <- which(grepl('Geno_0_numDonor', CN) == TRUE)
    colidx_geno_1_numdonor <- which(grepl('Geno_1_numDonor', CN) == TRUE)
    colidx_geno_2_numdonor <- which(grepl('Geno_2_numDonor', CN) == TRUE)
    colidx_geno_0_CC <- colidx_geno_0_numdonor + 2
    colidx_geno_1_CC <- colidx_geno_1_numdonor + 2
    colidx_geno_2_CC <- colidx_geno_2_numdonor + 2
    
    ## FDR for the RASQUAL default model
    FDR_Default_Model <- SNP_Loop_Data_RawCC$BH.Qvalue.Def

    ## FDR for the population specific model
    FDR_Pop_Model <- SNP_Loop_Data_RawCC$BH.Qvalue.pop

    ## FDR for the allele specific model
    FDR_AS_Model <- SNP_Loop_Data_RawCC$BH.Qvalue.AS

    ## ASRead specific paired t-test outcome
    t_test_pval <- SNP_Loop_Data_RawCC[, ncol(SNP_Loop_Data_RawCC)]

    ## initialize the following structures: 
    ## 1) mean and median contact counts and allele specific read counts
    MeanCC_DF <- matrix(0, nrow=nrow(SNP_Loop_Data), ncol=3)
    colnames(MeanCC_DF) <- c('Geno_0_MeanCC', 'Geno_1_MeanCC', 'Geno_2_MeanCC')    
    MedianCC_DF <- matrix(0, nrow=nrow(SNP_Loop_Data), ncol=3)
    colnames(MedianCC_DF) <- c('Geno_0_MedianCC', 'Geno_1_MedianCC', 'Geno_2_MedianCC')
    ## 2) mean AS read for the reference and alternate alleles
    Mean_ASRead_DF <- matrix(0, nrow=nrow(SNP_Loop_Data), ncol=2)
    colnames(Mean_ASRead_DF) <- c('Mean_ASRead_RefAllele', 'Mean_ASRead_AltAllele')
    ## 3) median AS read for the reference and alternate alleles
    Median_ASRead_DF <- matrix(0, nrow=nrow(SNP_Loop_Data), ncol=2)
    colnames(Median_ASRead_DF) <- c('Median_ASRead_RefAllele', 'Median_ASRead_AltAllele')    
    ## 4) mean deviation of AS read counts between the reference and alternate alleles
    ## that is, mean of (ASRead_reference_allele - ASRead_alternate_allele)
    Mean_Trend_ASRead_Vec <- rep(0, nrow(SNP_Loop_Data))

    ## fill these matrices for all individual entries
    for (i in 1:nrow(SNP_Loop_Data)) {
        numDonor_Homo <- SNP_Loop_Data[i, colidx_geno_0_numdonor]
        numDonor_Hetero <- SNP_Loop_Data[i, colidx_geno_1_numdonor]
        numDonor_Alt_Homo <- SNP_Loop_Data[i, colidx_geno_2_numdonor]        
        
        ## genotype specific contact counts - mean and median
        if (numDonor_Homo > 1) {        
            ccvec_geno_0 <- as.numeric(unlist(strsplit(SNP_Loop_Data[i, colidx_geno_0_CC], ":")))
            MedianCC_DF[i, 1] <- median(ccvec_geno_0)
            MeanCC_DF[i, 1] <- mean(ccvec_geno_0)
        } else if (numDonor_Homo > 0) {
            MedianCC_DF[i, 1] <- as.numeric(SNP_Loop_Data[i, colidx_geno_0_CC])
            MeanCC_DF[i, 1] <- as.numeric(SNP_Loop_Data[i, colidx_geno_0_CC])
        }
        if (numDonor_Hetero > 1) {        
            ccvec_geno_1 <- as.numeric(unlist(strsplit(SNP_Loop_Data[i, colidx_geno_1_CC], ":")))
            MedianCC_DF[i, 2] <- median(ccvec_geno_1)
            MeanCC_DF[i, 2] <- mean(ccvec_geno_1)
        } else if (numDonor_Hetero > 0) {
            MedianCC_DF[i, 2] <- as.numeric(SNP_Loop_Data[i, colidx_geno_1_CC])
            MeanCC_DF[i, 2] <- as.numeric(SNP_Loop_Data[i, colidx_geno_1_CC])
        }
        if (numDonor_Alt_Homo > 1) {        
            ccvec_geno_2 <- as.numeric(unlist(strsplit(SNP_Loop_Data[i, colidx_geno_2_CC], ":")))
            MedianCC_DF[i, 3] <- median(ccvec_geno_2)
            MeanCC_DF[i, 3] <- mean(ccvec_geno_2)
        } else if (numDonor_Alt_Homo > 0) {
            MedianCC_DF[i, 3] <- as.numeric(SNP_Loop_Data[i, colidx_geno_2_CC])
            MeanCC_DF[i, 3] <- as.numeric(SNP_Loop_Data[i, colidx_geno_2_CC])
        }

        ## allele specific reads for heterozygous donors
        if (numDonor_Hetero > 1) {
            RefAllele_ReadVec <- as.numeric(unlist(strsplit(SNP_Loop_Data[i, (ncol(SNP_Loop_Data)-2)], ":")))
            AltAllele_ReadVec <- as.numeric(unlist(strsplit(SNP_Loop_Data[i, (ncol(SNP_Loop_Data)-1)], ":")))
            Mean_ASRead_DF[i, 1] <- mean(RefAllele_ReadVec)
            Mean_ASRead_DF[i, 2] <- mean(AltAllele_ReadVec)
            Median_ASRead_DF[i, 1] <- median(RefAllele_ReadVec)
            Median_ASRead_DF[i, 2] <- median(AltAllele_ReadVec)            
            Mean_Trend_ASRead_Vec[i] <- mean(RefAllele_ReadVec - AltAllele_ReadVec)
        } else if (numDonor_Hetero > 0) {
            RefAllele_ReadCnt <- as.numeric(SNP_Loop_Data[i, (ncol(SNP_Loop_Data)-2)])
            AltAllele_ReadCnt <- as.numeric(SNP_Loop_Data[i, (ncol(SNP_Loop_Data)-1)])
            Mean_ASRead_DF[i, 1] <- RefAllele_ReadCnt
            Mean_ASRead_DF[i, 2] <- AltAllele_ReadCnt
            Median_ASRead_DF[i, 1] <- RefAllele_ReadCnt
            Median_ASRead_DF[i, 2] <- AltAllele_ReadCnt            
            Mean_Trend_ASRead_Vec[i] <- (RefAllele_ReadCnt - AltAllele_ReadCnt)
        }
    }

    ## if both raw and normalized matrices are provided
    ## following variables are used to define more stringent conditions
    if (file.exists(RASQUAL_out_sig_loops_donor_genotype_NormCC_file)) {
        MeanCC_DF_2 <- matrix(0, nrow=nrow(SNP_Loop_Data_RawCC), ncol=3)
        colnames(MeanCC_DF_2) <- c('Geno_0_MeanCC', 'Geno_1_MeanCC', 'Geno_2_MeanCC')
    
        MedianCC_DF_2 <- matrix(0, nrow=nrow(SNP_Loop_Data_RawCC), ncol=3)
        colnames(MedianCC_DF_2) <- c('Geno_0_MedianCC', 'Geno_1_MedianCC', 'Geno_2_MedianCC')

        ## fill these matrices for all individual entries
        for (i in 1:nrow(SNP_Loop_Data_RawCC)) {
            numDonor_Homo <- SNP_Loop_Data_RawCC[i, colidx_geno_0_numdonor]
            numDonor_Hetero <- SNP_Loop_Data_RawCC[i, colidx_geno_1_numdonor]
            numDonor_Alt_Homo <- SNP_Loop_Data_RawCC[i, colidx_geno_2_numdonor]        
            
            ## population specific contact counts - mean and median
            if (numDonor_Homo > 1) {        
                rawccvec_geno_0 <- as.numeric(unlist(strsplit(SNP_Loop_Data_RawCC[i, colidx_geno_0_CC], ":")))
                MedianCC_DF_2[i, 1] <- median(rawccvec_geno_0)
                MeanCC_DF_2[i, 1] <- mean(rawccvec_geno_0)
            } else if (numDonor_Homo > 0) {
                MedianCC_DF_2[i, 1] <- as.numeric(SNP_Loop_Data_RawCC[i, colidx_geno_0_CC])
                MeanCC_DF_2[i, 1] <- as.numeric(SNP_Loop_Data_RawCC[i, colidx_geno_0_CC])
            }
            if (numDonor_Hetero > 1) {        
                rawccvec_geno_1 <- as.numeric(unlist(strsplit(SNP_Loop_Data_RawCC[i, colidx_geno_1_CC], ":")))
                MedianCC_DF_2[i, 2] <- median(rawccvec_geno_1)
                MeanCC_DF_2[i, 2] <- mean(rawccvec_geno_1)
            } else if (numDonor_Hetero > 0) {
                MedianCC_DF_2[i, 2] <- as.numeric(SNP_Loop_Data_RawCC[i, colidx_geno_1_CC])
                MeanCC_DF_2[i, 2] <- as.numeric(SNP_Loop_Data_RawCC[i, colidx_geno_1_CC])
            }
            if (numDonor_Alt_Homo > 1) {        
                rawccvec_geno_2 <- as.numeric(unlist(strsplit(SNP_Loop_Data_RawCC[i, colidx_geno_2_CC], ":")))
                MedianCC_DF_2[i, 3] <- median(rawccvec_geno_2)
                MeanCC_DF_2[i, 3] <- mean(rawccvec_geno_2)
            } else if (numDonor_Alt_Homo > 0) {
                MedianCC_DF_2[i, 3] <- as.numeric(SNP_Loop_Data_RawCC[i, colidx_geno_2_CC])
                MeanCC_DF_2[i, 3] <- as.numeric(SNP_Loop_Data_RawCC[i, colidx_geno_2_CC])
            }            
        }
    }

    ## entries which have at least a donor for all three genotypes
    idx_All_Geno <- which((SNP_Loop_Data_RawCC[, colidx_geno_0_numdonor] > 0) & (SNP_Loop_Data_RawCC[, colidx_geno_1_numdonor] > 0) & (SNP_Loop_Data_RawCC[, colidx_geno_2_numdonor] > 0))
    cat(sprintf("\n *** With respect to the complete data - number of entries having donors for all three genotypes : %s ", length(idx_All_Geno)))

    ## entries having at least two donors for all three genotypes
    idx_All_Geno_2 <- which((SNP_Loop_Data_RawCC[, colidx_geno_0_numdonor] > 1) & (SNP_Loop_Data_RawCC[, colidx_geno_1_numdonor] > 1) & (SNP_Loop_Data_RawCC[, colidx_geno_2_numdonor] > 1))
    cat(sprintf("\n *** With respect to the complete data - number of entries having at least two donors for all three genotypes : %s ", length(idx_All_Geno_2)))

    ## loops significant in the default model
    ## here, all the loops are significant in the default model - but still generalizing the condition
    def_Sig_Idx <- which(!is.na(FDR_Default_Model) & (FDR_Default_Model < FDR_IQTL))

    ## loops significant in the population specific model
    pop_Sig_Idx <- which(!is.na(FDR_Pop_Model) & (FDR_Pop_Model < FDR_IQTL))

    ## loops having ascending or descending population specific trends
    ## (with respect to both mean and median contact counts) 
    ## considers only one CC data - either normalized contact count (if available) or raw contact count    
    pop_Trend_idx <- which(((MeanCC_DF[,1] >= MeanCC_DF[,2]) & (MeanCC_DF[,2] >= MeanCC_DF[,3]) & (MedianCC_DF[,1] >= MedianCC_DF[,2]) & (MedianCC_DF[,2] >= MedianCC_DF[,3])) | ((MeanCC_DF[,1] <= MeanCC_DF[,2]) & (MeanCC_DF[,2] <= MeanCC_DF[,3]) & (MedianCC_DF[,1] <= MedianCC_DF[,2]) & (MedianCC_DF[,2] <= MedianCC_DF[,3])))
    cat(sprintf("\n *** With respect to the complete data - number of entries having ascending or descending trend in population specific model : %s ***** \n", length(pop_Trend_idx)))

    ## loops whose Allele specific read based t-test and p-values are significant
    ASRead_tTest_Sig_Idx <- which(!is.na(t_test_pval) & (t_test_pval < FDR_IQTL))
    cat(sprintf("\n *** With respect to the complete data - number of entries having allele specific read t-test p-value significant : %s ***** \n", length(ASRead_tTest_Sig_Idx)))

    ## loops having allele specific reads and genotype specific contact counts in the similar direction
    ## that is, suppose reference allele is A and alternate allele is G
    ## further, read count for A > read count for G
    ## if contact count for A|A > contact count for G|G + contact count for A|A > contact count for A|G, 
    ## then we'll say that the both trends are in similar direction
    ## considers only one CC data - either normalized contact count (if available) or raw contact count   
    ## checks both mean and median    
    ASRead_Population_same_Trend_idx <- which(((MeanCC_DF[,1] >= MeanCC_DF[,3]) & (MeanCC_DF[,1] >= MeanCC_DF[,2]) & (MedianCC_DF[,1] >= MedianCC_DF[,3]) & (MedianCC_DF[,1] >= MedianCC_DF[,2]) & (Mean_Trend_ASRead_Vec >= 0)) | ((MeanCC_DF[,1] <= MeanCC_DF[,3]) & (MeanCC_DF[,1] <= MeanCC_DF[,2]) & (MedianCC_DF[,1] <= MedianCC_DF[,3]) & (MedianCC_DF[,1] <= MedianCC_DF[,2]) & (Mean_Trend_ASRead_Vec <= 0)))
    cat(sprintf("\n\n *** number of entries having similar trend between population specific CC and allele specific reads : %s ***** \n", length(ASRead_Population_same_Trend_idx)))

    # ## further stringent condition
    # ## the population trend should be ascending / descending for both mean and median CC
    # ## ASRead trend should be similar
    # if (0) {
    #     ASRead_Population_same_Trend_stringent_idx <- which(((MeanCC_DF[,1] >= MeanCC_DF[,2]) & (MeanCC_DF[,2] >= MeanCC_DF[,3]) & (MedianCC_DF[,1] >= MedianCC_DF[,2]) & (MedianCC_DF[,2] >= MedianCC_DF[,3]) & (Mean_Trend_ASRead_Vec >= 0) & (Mean_ASRead_DF[,1] >= Mean_ASRead_DF[,2]) & (Median_ASRead_DF[,1] >= Median_ASRead_DF[,2])) | ((MeanCC_DF[,1] <= MeanCC_DF[,2]) & (MeanCC_DF[,2] <= MeanCC_DF[,3]) & (MedianCC_DF[,1] <= MedianCC_DF[,2]) & (MedianCC_DF[,2] <= MedianCC_DF[,3]) & (Mean_Trend_ASRead_Vec <= 0) & (Mean_ASRead_DF[,1] <= Mean_ASRead_DF[,2]) & (Median_ASRead_DF[,1] <= Median_ASRead_DF[,2])))
    #     cat(sprintf("\n\n *** number of entries having similar trend between population specific CC and allele specific reads - stringent condition : %s ***** \n", length(ASRead_Population_same_Trend_stringent_idx)))        
    # }

    ## following metrics are used, when the normalized contact count is available
    if (file.exists(RASQUAL_out_sig_loops_donor_genotype_NormCC_file)) {

        ## loops having ascending or descending population specific trends
        ## (both mean and median contact counts should have ascending or descending population specific trends)
        ## the trend should be satisfied by both raw and normalized CC values (more stringent condition)        
        pop_Trend_idx_2 <- which(((MeanCC_DF[,1] >= MeanCC_DF[,2]) & (MeanCC_DF[,2] >= MeanCC_DF[,3]) & (MedianCC_DF[,1] >= MedianCC_DF[,2]) & (MedianCC_DF[,2] >= MedianCC_DF[,3]) & (MeanCC_DF_2[,1] >= MeanCC_DF_2[,2]) & (MeanCC_DF_2[,2] >= MeanCC_DF_2[,3]) & (MedianCC_DF_2[,1] >= MedianCC_DF_2[,2]) & (MedianCC_DF_2[,2] >= MedianCC_DF_2[,3])) | ((MeanCC_DF[,1] <= MeanCC_DF[,2]) & (MeanCC_DF[,2] <= MeanCC_DF[,3]) & (MedianCC_DF[,1] <= MedianCC_DF[,2]) & (MedianCC_DF[,2] <= MedianCC_DF[,3]) & (MeanCC_DF_2[,1] <= MeanCC_DF_2[,2]) & (MeanCC_DF_2[,2] <= MeanCC_DF_2[,3]) & (MedianCC_DF_2[,1] <= MedianCC_DF_2[,2]) & (MedianCC_DF_2[,2] <= MedianCC_DF_2[,3])))
        cat(sprintf("\n\n *** With respect to the complete data - number of entries having ascending or descending trend in population specific model (both raw and normalized CC) : %s ***** \n", length(pop_Trend_idx_2)))

        ## loops having allele specific reads and genotype specific contact counts in the similar direction
        ## that is, suppose reference allele is A and alternate allele is G
        ## further, read count for A > read count for G
        ## if contact count for AA > contact count for GG then we'll say that the both trends are in similar direction
        ## both raw and normalized CC are accounted for (more stringent condition)
        ## new code 
        ASRead_Population_same_Trend_idx_2 <- which(((MeanCC_DF[,1] >= MeanCC_DF[,3]) & (MeanCC_DF[,1] >= MeanCC_DF[,2]) & (MedianCC_DF[,1] >= MedianCC_DF[,3]) & (MedianCC_DF[,1] >= MedianCC_DF[,2]) & (MeanCC_DF_2[,1] >= MeanCC_DF_2[,3]) & (MeanCC_DF_2[,1] >= MeanCC_DF_2[,2]) & (MedianCC_DF_2[,1] >= MedianCC_DF_2[,3]) & (MedianCC_DF_2[,1] >= MedianCC_DF_2[,2]) & (Mean_Trend_ASRead_Vec >= 0)) | ((MeanCC_DF[,1] <= MeanCC_DF[,3]) & (MeanCC_DF[,1] <= MeanCC_DF[,2]) & (MedianCC_DF[,1] <= MedianCC_DF[,3]) & (MedianCC_DF[,1] <= MedianCC_DF[,2]) & (MeanCC_DF_2[,1] <= MeanCC_DF_2[,3]) & (MeanCC_DF_2[,1] <= MeanCC_DF_2[,2]) & (MedianCC_DF_2[,1] <= MedianCC_DF_2[,3]) & (MedianCC_DF_2[,1] <= MedianCC_DF_2[,2]) & (Mean_Trend_ASRead_Vec <= 0)))
        cat(sprintf("\n\n *** With respect to the complete data - number of entries having similar trend between population specific CC and allele specific reads (both raw and normalized CC) : %s ***** \n", length(ASRead_Population_same_Trend_idx_2)))

        # ## further stringent condition
        # ## the population trend should be ascending / descending for both mean and median CC
        # ## ASRead trend should be similar
        # ## the trend should be satisfied by both raw and normalized CC values (more stringent condition)
        # if (0) {
        #     ASRead_Population_same_Trend_stringent_idx_2 <- which(((MeanCC_DF[,1] >= MeanCC_DF[,2]) & (MeanCC_DF[,2] >= MeanCC_DF[,3]) & (MedianCC_DF[,1] >= MedianCC_DF[,2]) & (MedianCC_DF[,2] >= MedianCC_DF[,3]) & (MeanCC_DF_2[,1] >= MeanCC_DF_2[,2]) & (MeanCC_DF_2[,2] >= MeanCC_DF_2[,3]) & (MedianCC_DF_2[,1] >= MedianCC_DF_2[,2]) & (MedianCC_DF_2[,2] >= MedianCC_DF_2[,3]) & (Mean_Trend_ASRead_Vec >= 0) & (Mean_ASRead_DF[,1] >= Mean_ASRead_DF[,2]) & (Median_ASRead_DF[,1] >= Median_ASRead_DF[,2])) | ((MeanCC_DF[,1] <= MeanCC_DF[,2]) & (MeanCC_DF[,2] <= MeanCC_DF[,3]) & (MedianCC_DF[,1] <= MedianCC_DF[,2]) & (MedianCC_DF[,2] <= MedianCC_DF[,3]) & (MeanCC_DF_2[,1] <= MeanCC_DF_2[,2]) & (MeanCC_DF_2[,2] <= MeanCC_DF_2[,3]) & (MedianCC_DF_2[,1] <= MedianCC_DF_2[,2]) & (MedianCC_DF_2[,2] <= MedianCC_DF_2[,3]) & (Mean_Trend_ASRead_Vec <= 0) & (Mean_ASRead_DF[,1] <= Mean_ASRead_DF[,2]) & (Median_ASRead_DF[,1] <= Median_ASRead_DF[,2])))
        #     cat(sprintf("\n\n *** With respect to the complete data - number of entries having similar trend between population specific CC and allele specific reads (both raw and normalized CC) - stringent condition : %s ***** \n", length(ASRead_Population_same_Trend_stringent_idx_2)))        
        # }

    }   # end normalized contact count condition


    ##========================
    ## filtered entries
    ## using population and allele specific trends
    ##========================

    FiltOutDir <- paste0(summary.out.path.combined, '/FILTERED')
    system(paste("mkdir -p", FiltOutDir))

    ## significance condition
    ## 1. significant in default model
    ## 2. present in all three genotypes in at least 2 donors
    ## 3A. significant in population specific model + population specific trend is strictly ascending or descending
    ## 3B. AS read t-test p-value is significant + allele specific trend matches with the overall population specific trend
    if (file.exists(RASQUAL_out_sig_loops_donor_genotype_NormCC_file)) {
        sel_idx <- intersect(intersect(def_Sig_Idx, idx_All_Geno_2), union(intersect(pop_Sig_Idx, pop_Trend_idx_2), intersect(ASRead_tTest_Sig_Idx, ASRead_Population_same_Trend_idx_2)))
    } else {
        sel_idx <- intersect(intersect(def_Sig_Idx, idx_All_Geno_2), union(intersect(pop_Sig_Idx, pop_Trend_idx), intersect(ASRead_tTest_Sig_Idx, ASRead_Population_same_Trend_idx)))
    }
    cat(sprintf("\n\n *** Selected final number of loops : %s ", length(sel_idx)))
    if (length(sel_idx) > 0) {
        outfile_rawcc <- paste0(FiltOutDir, "/out_rawCC.txt")        
        write.table(SNP_Loop_Data_RawCC[sel_idx, ], outfile_rawcc, row.names=F, col.names=T, sep="\t", quote=F, append=F)
        if (file.exists(RASQUAL_out_sig_loops_donor_genotype_NormCC_file)) {
            outfile_normcc <- paste0(FiltOutDir, "/out_normCC.txt")        
            write.table(SNP_Loop_Data_NormCC[sel_idx, ], outfile_normcc, row.names=F, col.names=T, sep="\t", quote=F, append=F)
        }            
    }

}   # end file exist condition



