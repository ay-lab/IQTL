#!/usr/bin/env Rscript

##=========================
## IQTL code - sourya - for RASQUAL
## here allele specific reads are re-estimated
## using only reads belonging to the current HiChIP interacting bins (and optionally their neighborhood)
##=========================
library(data.table)
library(parallel)
options(scipen=10)
options(datatable.fread.datatable=FALSE)

# number of processors within the system
ncore <- detectCores()
cat(sprintf("\n Number of cores in the system: %s ", ncore))

## sourya - use 8 threads per RASQUAL command executon
nthread <- 8

##===========================
## functions
##===========================

##======= function to bgzip compress
bgzip_compress_vcf <- function(inpfile) {
    if (file.exists(paste0(inpfile, ".gz"))) {
        system(paste("rm", paste0(inpfile, ".gz")))
    }
    system(paste0("bgzip ", inpfile))
    if (file.exists(paste0(inpfile, ".gz.tbi"))) {
        system(paste("rm", paste0(inpfile, ".gz.tbi")))
    }          
    system(paste0("tabix -p vcf ", inpfile, ".gz"))
}

##======= function to write the interacting bins for the input interactions
Write_Interacting_Intervals <- function(count.dat, start_idx, end_idx, window.size, allintervalfile1, allintervalfile2) {

    for (n.row in seq(start_idx, end_idx)) {
        tmp.count <- count.dat[n.row, ]
        reg.name <- tmp.count$V1
        temp.reg <- sub(":", "_", reg.name)
        temp.reg <- sub("-", "_", temp.reg)
        start1 <- as.numeric( unlist( strsplit( temp.reg, "_"))[2])
        end1 <- (start1 + binsize)
        end2 <- as.numeric( unlist( strsplit( temp.reg, "_"))[3])
        start2 <- (end2 - binsize)
        start1.vcf <- max(0, (start1 - window.size))  # error handle - sourya
        end1.vcf <- (end1 + window.size)
        start2.vcf <- max(0, (start2 - window.size))  # error handle - sourya
        end2.vcf <- (end2 + window.size)
        ## data frame consisting of first interacting intervals
        currdf1 <- data.frame(chr=currchr, start=start1.vcf, end=end1.vcf)
        ## data frame consisting of second interacting intervals
        currdf2 <- data.frame(chr=currchr, start=start2.vcf, end=end2.vcf)
        if (n.row == start_idx) {
            intervaldf1 <- currdf1
            intervaldf2 <- currdf2
        } else {
            intervaldf1 <- rbind.data.frame(intervaldf1, currdf1)
            intervaldf2 <- rbind.data.frame(intervaldf2, currdf2)
        }
    }
    intervaldf1 <- unique(intervaldf1)
    intervaldf2 <- unique(intervaldf2)
    write.table(intervaldf1, allintervalfile1, row.names=F, col.names=F, sep="\t", quote=F, append=F)
    write.table(intervaldf2, allintervalfile2, row.names=F, col.names=F, sep="\t", quote=F, append=F)

}   # end function


##====== function to subsample the BAM reads
##====== only extract the BAM reads whose one end belongs to the intervals in the file "allintervalfile1"
##====== while the other end belongs to the intervals in the file "allintervalfile2"
Subset_BAM_specific_bin_pairs <- function(inphichipalignfile, subsethichipalignfile, allintervalfile1, allintervalfile2) {
    system(paste0("bedtools pairtobed -abam ", inphichipalignfile, " -b ", allintervalfile1, " | bedtools pairtobed -abam stdin -b ", allintervalfile2, " | samtools view -bh -o ", subsethichipalignfile, " - "))
}


##====== function to re-etimate the allele specific reads
##====== with respect to the specified intervals 
ReEst_ASReads <- function(n.row, vcf.z.path, vcf.dat, start1.vcf, end1.vcf, start2.vcf, end2.vcf, curr_SNPset_file, out_currregion_SNPfile, sample.names, out.path.Default, RefGenomeFastaFile) {

    ##===============
    ## step 1: subset VCF by selecting SNPs within these regions
    ##===============
    feat.vcf <- vcf.dat[which(((vcf.dat[, 2] >= start1.vcf) & (vcf.dat[, 2] <= end1.vcf)) | ((vcf.dat[, 2] >= start2.vcf) & (vcf.dat[, 2] <= end2.vcf))), ]
    ## if there is no SNP within the specified intervals, continue
    if (nrow(feat.vcf) == 0) {
        return(0)
    }
    cat(sprintf(" - nrow(feat.vcf) : %s ", nrow(feat.vcf)))

    ##===============
    ## step 2: extract SNPs from VCF file, within this region
    ## Note: this is important since we'll estimate fresh AS reads for these SNPs
    ##===============    
    system(paste0("zcat ", vcf.z.path, " | awk \'{if ((substr($1,1,1)==\"#\") || (($2>=", start1.vcf, ") && ($2<=", end1.vcf, ")) || (($2>=", start2.vcf, ") && ($2<=", end2.vcf, "))) {print $0}}\' - > ", curr_SNPset_file))
    bgzip_compress_vcf(curr_SNPset_file)

    ##===============
    ## step 3: re-estimating AS reads, with respect to the donor specific HiChIP alignments
    ## and for the current HiChIP interaction
    ##===============

    ## parallel processing
    parallel:::mclapply(1:length(sample.names), mc.cores=ncore, function(si) {
    # lapply(1:length(sample.names), function(si) {    

        ## input HiChIP alignment file - already extacted per chromosome, CIS reads 
        ## and only reads within the current set of regions    
        inphichipalignfile <- paste0(out.path.Default, 'out_', sample.names[si], '_reference_all_regions.bam')

        ## output alignment file: to contain the reads for the current interval, and for the current donor
        ## reads whose one end is in the first interval (start1.vcf, end1.vcf), 
        ## and the other end on the second interval (start2.vcf, end2.vcf)
        ## the alignment file needs to be sorted as well (before calling allele specific reads)
        outhichipalignfile <- paste0(out.path.Default, 'out_', n.row, '_', sample.names[si], '.bam')
        outhichipalignfile_sort <- paste0(out.path.Default, 'out_', n.row, '_', sample.names[si], '_sorted.bam')
        system(paste0("samtools view -h ", inphichipalignfile, " | awk \'((substr($1,1,1)==\"@\") || (($4>=", start1.vcf, ") && ($4<=", end1.vcf, ") && ($8>=", start2.vcf, ") && ($8<=", end2.vcf, ")) || (($8>=", start1.vcf, ") && ($8<=", end1.vcf, ") && ($4>=", start2.vcf, ") && ($4<=", end2.vcf, ")))\' - | samtools view -bh -o ", outhichipalignfile, " - "))
        system(paste0("samtools sort -o ", outhichipalignfile_sort, " -@ 8 ", outhichipalignfile))

        ## re-estimate allele specific reads for the current donor 
        ## Note that we have removed duplicate reads before calling allele specific reads (default operation of ASEReadCounter)     
        outASReadFile <- paste0(out.path.Default, 'ASRead_', n.row, '_', sample.names[si], '.csv')

        # system(paste0(GATKExec, " ASEReadCounter -R ", RefGenomeFastaFile, " -O ", outASReadFile, " -I ", outhichipalignfile_sort, " -V ", curr_SNPset_file, ".gz --min-mapping-quality 10 -DF NotDuplicateReadFilter"))
        system(paste0(GATKExec, " ASEReadCounter -R ", RefGenomeFastaFile, " -O ", outASReadFile, " -I ", outhichipalignfile_sort, " -V ", curr_SNPset_file, ".gz --min-mapping-quality 10"))

    } ) # end parallel processing

    ## now accumulate the AS reads for all donors
    for (si in 1:length(sample.names)) {

        outhichipalignfile <- paste0(out.path.Default, 'out_', n.row, '_', sample.names[si], '.bam')
        outhichipalignfile_sort <- paste0(out.path.Default, 'out_', n.row, '_', sample.names[si], '_sorted.bam')

        outASReadFile <- paste0(out.path.Default, 'ASRead_', n.row, '_', sample.names[si], '.csv')
        outASReadData <- data.table::fread(outASReadFile, header=T, sep="\t", stringsAsFactors=F)      

        ## column for the current donor in the reference "feat.vcf" structure
        targetcol <- (9 + si)
        
        ## step 1: first check all SNPs for this donor, in the "feat.vcf" structure
        ## and reset their allele specific counts to 0,0    
        genovec <- sapply(feat.vcf[, targetcol], function(X) unlist(strsplit(X, ":"))[1])
        dummyASvec <- rep(":0,0", nrow(feat.vcf))
        tempdf1 <- data.frame(G=genovec, A=dummyASvec)
        tempdf1$V <- paste(tempdf1$G, tempdf1$A, sep="")
        feat.vcf[, targetcol] <- tempdf1$V

        ## put the allele specific reads for this donor, and for the heterozygous SNPs
        if (nrow(outASReadData) > 0) {

            m <- match(feat.vcf[, 3], outASReadData[, 3])
            idx_s <- which(!is.na(m))
            idx_t <- m[!is.na(m)]         
            v1 <- sapply(feat.vcf[idx_s, targetcol], function(X) unlist(strsplit(X, ":"))[1])

            ## important - sourya
            ## the entries in v1 will have one of the following formats:
            ## 0|0, 1|1, 0|1 and 1|0
            ## outASReadData[, 6] contains count for reference allele; outASReadData[, 7] contains count for alternate allele;
            ## if v1 is 0|0, 1|1, 0|1, use "RefAlleleCount,AltAlleleCount" format
            ## if v1 is 1|0, use "AltAlleleCount,RefAlleleCount" format
            idx_alt_ref <- which(v1 == "1|0")
            idx_ref_alt <- setdiff(seq(1,length(v1)), idx_alt_ref)
            cat(sprintf("\n ==>> adjusting allele specific reads - sample index : %s targetcol : %s donor : %s total allele specific SNPs : %s 1|0 AS : %s 0|1 : %s ", si, targetcol, colnames(feat.vcf)[targetcol], length(v1), length(idx_alt_ref), length(idx_ref_alt)))

            if (length(idx_ref_alt) > 0) {
                ## v1 is 0|0, 1|1, 0|1, use "RefAlleleCount,AltAlleleCount" format
                feat.vcf[idx_s[idx_ref_alt], targetcol] <- paste0(v1[idx_ref_alt], ":", outASReadData[idx_t[idx_ref_alt], 6], ",", outASReadData[idx_t[idx_ref_alt], 7])        
            }
            if (length(idx_alt_ref) > 0) {
                ## v1 is 1|0, use "AltAlleleCount,RefAlleleCount" format
                feat.vcf[idx_s[idx_alt_ref], targetcol] <- paste0(v1[idx_alt_ref], ":", outASReadData[idx_t[idx_alt_ref], 7], ",", outASReadData[idx_t[idx_alt_ref], 6])  
            }

        } # end outASReadData presence condition     
    
        ## delete the temporary files 
        system(paste("rm", outhichipalignfile))
        system(paste("rm", outhichipalignfile_sort))
        system(paste("rm", outASReadFile))
 
    }   # end donor loop - sequential processing

    ## dump the output SNPs with re-estimated allele specific reads
    system(paste0("zcat ", curr_SNPset_file, ".gz | awk \'(substr($1,1,1)==\"#\")\' - > ", out_currregion_SNPfile))
    write.table(feat.vcf, out_currregion_SNPfile, row.names=F, col.names=F, sep="\t", quote=F, append=T)
    bgzip_compress_vcf(out_currregion_SNPfile)

    return(1)
}

##===========================
## main code
##===========================

args <- commandArgs(TRUE)

## Input parameters
currchr <- as.character(args[1])
rasqual.path <- as.character(args[2])
DonorlistFile <- as.character(args[3])
INPDIR <- as.character(args[4])
out.path.Base <- as.character(args[5])
GATKExec <- as.character(args[6])  # GATK executable
RefGenomeFastaFile <- as.character(args[7])  # fasta file for the reference genome
vcf.z.path <- as.character(args[8])
vcf.path <- as.character(args[9])
window.size <- as.numeric(args[10])
binsize <- as.integer(args[11])
chunkID <- as.integer(args[12])
HiChIPAlignmentBaseDir <- as.character(args[13])  ## base directory containing merged HiChIP alignment
MinLoopDist <- as.integer(args[14])

## from the input donor wise list, 
## first column denotes the sample ID
## second column denotes the donor ID
sample.id.col <- 1
donor.id.col <- 2

# Minor allele Frequency used for Rasqual
maf <- 0.05

# Imputation quality used for Rasqual
imputation.quality <- 0.3

y.bin.path <- paste0(INPDIR, '/Y.bin')
y.norm.path <- paste0(INPDIR, '/Y.txt')
k.bin.path <- paste0(INPDIR, '/K.bin')
x.bin.path <- paste0(INPDIR, '/X.bin')
chunkfile <- paste0(INPDIR, '/chunk_info.txt')

## stores default model RASQUAL output (using combined model of population + allele specific)
out.path.Default <- paste0(out.path.Base, '/RASQUAL_Default/', currchr, '/', chunkID, '/')
## stores only population specific RASQUAL output
out.path.Pop <- paste0(out.path.Base, '/RASQUAL_Pop/', currchr, '/', chunkID, '/')
## stores only allele specific RASQUAL output
out.path.AS <- paste0(out.path.Base, '/RASQUAL_AS/', currchr, '/', chunkID, '/')
system(paste("mkdir -p", out.path.Default))
system(paste("mkdir -p", out.path.Pop))
system(paste("mkdir -p", out.path.AS))

## output log file
outlogfile <- paste0(out.path.Default, 'out.log')
sink(outlogfile)

cat("Chromosome:", currchr, "\n")
cat("Rasqual executable:", rasqual.path, "\n")
cat("GATK executable: ", GATKExec, "\n")
cat("Donor list file:", DonorlistFile, "\n")
cat("VCF gzipped path:", vcf.z.path, "\n")
cat("VCF path:", vcf.path, "\n")
cat("RASQUAL input file folder:", INPDIR, "\n")
cat("Output path (default RASQUAL output) :", out.path.Default, "\n")
cat("Output path (population specific RASQUAL output) :", out.path.Pop, "\n")
cat("Output path (allele specific RASQUAL output) :", out.path.AS, "\n")
cat("Window size:", window.size, "\n")
cat("Bin size: ", binsize, "\n")
cat("chunkfile : ", chunkfile, "\n")
cat("chunkID : ", chunkID, "\n")
cat("HiChIPAlignmentBaseDir: ", HiChIPAlignmentBaseDir, "\n")
cat("RefGenomeFastaFile: ", RefGenomeFastaFile, "\n")
cat("MinLoopDist: ", MinLoopDist, "\n")

## Reading data - vcf file (absolute)
vcf.dat <- data.table::fread(vcf.path, header=TRUE)

## Reading data - count file (Y.txt)
count.dat <- data.table::fread(y.norm.path)

## reading data - chunk information file
chunkData <- data.table::fread(chunkfile, header=T)

## Reading data - list of samples
DonorlistData <- data.table::fread(DonorlistFile)
sample.align.dirlist <- as.character(DonorlistData[, sample.id.col])
sample.names <- as.character(DonorlistData[, donor.id.col])
sample.size <- length(sample.names)
cat(sprintf("\n Sample names (donors) : %s ", paste(sample.names, collapse=" ")))

## get the rows to process
start_idx <- as.integer(chunkData[which((chunkData[,1] == currchr) & (chunkData[,2] == chunkID)), 3])
end_idx <- as.integer(chunkData[which((chunkData[,1] == currchr) & (chunkData[,2] == chunkID)), 4])
cat("\n\n rows to be processed : ", start_idx, "to ", end_idx, "\n")

##======= 
## first accumulate all the regions (i.e. interacting bins) 
## belonging to the current set of "count.dat" entries
##======= 
allintervalfile1 <- paste0(out.path.Default, 'complete_set_of_regions_first_interacting_bin.bed')
allintervalfile2 <- paste0(out.path.Default, 'complete_set_of_regions_second_interacting_bin.bed')
Write_Interacting_Intervals(count.dat, start_idx, end_idx, window.size, allintervalfile1, allintervalfile2)

##======= 
## for each donor, get their cis bam files (sorted by read names)
## and extract the reads which overlap with the given set of regions (mentioned in files allintervalfile1 and allintervalfile2)
## use bedtools pairtobed routine, to get reads whose one end overlaps with bins in "allintervalfile1" 
## and the other end overlaps with bins in "allintervalfile2"
## helps to create reference much smaller sized bam files per donor, and reduce the time
##======= 
for (si in 1:length(sample.names)) {
    cat(sprintf("\n ===>> creating bedtools reference for all regions - processing bam : %s ", sample.names[si]))
    ## Note: this bedtools operation requires read alignment file sorted by read name
    inphichipalignfile <- paste0(HiChIPAlignmentBaseDir, '/', sample.align.dirlist[si], '/chrwise/merged_HiChIP_', currchr, '_cis_sorted_readname.bam')
    ## output alignment file
    subsethichipalignfile <- paste0(out.path.Default, 'out_', sample.names[si], '_reference_all_regions.bam')
    if (file.exists(subsethichipalignfile) == FALSE) {
        Subset_BAM_specific_bin_pairs(inphichipalignfile, subsethichipalignfile, allintervalfile1, allintervalfile2)
    }
}

##===== folders to contain the input and output SNPs (with default and re-estimated AS reads)
input_SNP_dir <- paste0(out.path.Default, '/SNPs_default_ASReads')
system(paste("mkdir -p", input_SNP_dir))
output_SNP_dir <- paste0(out.path.Default, '/SNPs_reEst_ASReads')
system(paste("mkdir -p", output_SNP_dir))

##==========================
## Running rasqual
## process individual count entries
for (n.row in seq(start_idx, end_idx)) {

    tmp.count <- count.dat[n.row, ]

    ## sourya - the first field contains chr:start1-end2
    ## for the HiChIP interactions, start1: start coordinate of first interacting bin
    ## for the HiChIP interactions, end2: end coordinate of second interacting bin
    reg.name <- tmp.count$V1
    temp.reg <- sub(":", "_", reg.name)
    temp.reg <- sub("-", "_", temp.reg)
    start1 <- as.numeric( unlist( strsplit( temp.reg, "_"))[2])
    end1 <- (start1 + binsize)
    end2 <- as.numeric( unlist( strsplit( temp.reg, "_"))[3])
    start2 <- (end2 - binsize)

    ## with respect to the given window size (offset from the interacting bin)
    ## define the start and end positions with respect to both interacting bins
    ## this will be useful to define the candidate SNPs
    start1.vcf <- max((start1 - window.size), 0)
    end1.vcf <- (end1 + window.size)
    start2.vcf <- max((start2 - window.size), 0)
    end2.vcf <- (end2 + window.size)

    ## sourya - loop distance condition
    ## check if the distance between end1.vcf and start2.vcf coordinates are < MinLoopDist
    ## in such a case, selectively expand the window so that their distance is >= MinLoopDist   
    if ((start2.vcf - end1.vcf) < MinLoopDist) {
        x <- (start2.vcf - end1.vcf)
        end1.vcf <- as.integer(end1.vcf - ((MinLoopDist - x) / 2))
        start2.vcf <- as.integer(start2.vcf + ((MinLoopDist - x) / 2))
    }

    cat(sprintf("\n ==>> processing row : %s reg.name : %s start1 : %s end1 : %s start2 : %s end2 : %s --- start1.vcf : %s end1.vcf : %s start2.vcf : %s end2.vcf : %s \n", n.row, reg.name, start1, end1, start2, end2, start1.vcf, end1.vcf, start2.vcf, end2.vcf))
    
    ##======================
    ## RASQUAL default command + population specific + allele specific
    ##======================
    ## here use two interacting bins separately as two genes, and employ two RASQUAL commands
    ## for the first interacting bin, employ only the nearby SNPs - first command
    ## for the second interacting bin, employ only the nearby SNPs - second command
    RASQUAL_outfile_1 <- paste0(out.path.Default, "t_", n.row, "_1.txt")
    RASQUAL_outfile_2 <- paste0(out.path.Default, "t_", n.row, "_2.txt")

    RASQUAL_Pop_outfile_1 <- paste0(out.path.Pop, "t_", n.row, "_1.txt")
    RASQUAL_Pop_outfile_2 <- paste0(out.path.Pop, "t_", n.row, "_2.txt")

    RASQUAL_AS_outfile_1 <- paste0(out.path.AS, "t_", n.row, "_1.txt")
    RASQUAL_AS_outfile_2 <- paste0(out.path.AS, "t_", n.row, "_2.txt")

    ## debug - sourya - added this condition - should be removed later
    if (!file.exists(RASQUAL_outfile_1) | !file.exists(RASQUAL_outfile_2) | !file.exists(RASQUAL_Pop_outfile_1) | !file.exists(RASQUAL_Pop_outfile_2) | !file.exists(RASQUAL_AS_outfile_1) | !file.exists(RASQUAL_AS_outfile_2)) {

        ##******************
        ## Re-estimate allele specific reads
        ##******************
        ## the "vcf.dat" structure contains allele specific reads estimated from complete HiChIP coverage
        ## here, we need to re-estimate allele specific reads using only the contacts 
        ## between the regions (start1.vcf, end1.vcf) and (start2.vcf, end2.vcf)
        ## i.e. allowing the window size

        curr_SNPset_file <- paste0(input_SNP_dir, '/input_SNPs_', n.row, '.vcf')
        out_currregion_SNPfile <- paste0(output_SNP_dir, '/out_SNPs_', n.row, '_re_estimated_AS_Reads.vcf')
        status_val <- ReEst_ASReads(n.row, vcf.z.path, vcf.dat, start1.vcf, end1.vcf, start2.vcf, end2.vcf, curr_SNPset_file, out_currregion_SNPfile, sample.names, out.path.Default, RefGenomeFastaFile)
        if (status_val == 0) {
            next
        }        
  
        ##============= process first interacting bin and nearby SNPs
        temp.vcf <- vcf.dat[which((vcf.dat[, 2] >= start1.vcf) & (vcf.dat[, 2] <= end1.vcf)), ]
        if (nrow(temp.vcf) > 0 ) {
            ## Adding snps to avoid error, this parameters are just for memory allocation
            n.snps <- nrow(temp.vcf) + 25
            n.feat.snps <- nrow(temp.vcf) + 5
            cat(sprintf("\n ==>> RASQUAL command - processing row : %s reg.name : %s start1 : %s end1 : %s start2 : %s end2 : %s --- first interacting bin -- start1.vcf : %s end1.vcf : %s nrow(temp.vcf) : %s n.snps : %s n.feat.snps : %s ", n.row, reg.name, start1, end1, start2, end2, start1.vcf, end1.vcf, nrow(temp.vcf), n.snps, n.feat.snps))

            ## print RASQUAL command (default - both population and allele specific operation)
            print( paste0( "tabix ", out_currregion_SNPfile, ".gz ", currchr, ":", start1.vcf, "-", end1.vcf, 
                        " | ", rasqual.path, " -y ",  y.bin.path," -k ", k.bin.path, " -x ", x.bin.path, " -n ", sample.size, 
                        " -j ", n.row," -l ", n.snps," -m ", n.feat.snps ,
                        " -s ", start1, " -e ", end1, " -f ", reg.name, " -z ",
                        " -a ", maf, " --imputation-quality ", imputation.quality ,
                        " --minor-allele-frequency-fsnp ", maf, " --imputation-quality-fsnp ", imputation.quality ,
                        " --n-threads ", nthread, 
                        " > ", RASQUAL_outfile_1))

            ## execute RASQUAL command (default - both population and allele specific operation)
            system( paste0( "tabix ", out_currregion_SNPfile, ".gz ", currchr, ":", start1.vcf, "-", end1.vcf, 
                        " | ", rasqual.path, " -y ",  y.bin.path," -k ", k.bin.path, " -x ", x.bin.path, " -n ", sample.size, 
                        " -j ", n.row," -l ", n.snps," -m ", n.feat.snps ,
                        " -s ", start1, " -e ", end1, " -f ", reg.name, " -z ",
                        " -a ", maf, " --imputation-quality ", imputation.quality ,
                        " --minor-allele-frequency-fsnp ", maf, " --imputation-quality-fsnp ", imputation.quality ,
                        " --n-threads ", nthread, 
                        " > ", RASQUAL_outfile_1))

            ## print RASQUAL command (default - population specific operation)
            print( paste0( "tabix ", out_currregion_SNPfile, ".gz ", currchr, ":", start1.vcf, "-", end1.vcf, 
                        " | ", rasqual.path, " -y ",  y.bin.path," -k ", k.bin.path, " -x ", x.bin.path, " -n ", sample.size, 
                        " -j ", n.row," -l ", n.snps," -m ", n.feat.snps ,
                        " -s ", start1, " -e ", end1, " -f ", reg.name, " -z ",
                        " -a ", maf, " --imputation-quality ", imputation.quality ,
                        " --minor-allele-frequency-fsnp ", maf, " --imputation-quality-fsnp ", imputation.quality ,
                        " --n-threads ", nthread, 
                        " --population-only > ", RASQUAL_Pop_outfile_1))

            ## execute RASQUAL command (default - population specific operation)
            system( paste0( "tabix ", out_currregion_SNPfile, ".gz ", currchr, ":", start1.vcf, "-", end1.vcf, 
                        " | ", rasqual.path, " -y ",  y.bin.path," -k ", k.bin.path, " -x ", x.bin.path, " -n ", sample.size, 
                        " -j ", n.row," -l ", n.snps," -m ", n.feat.snps ,
                        " -s ", start1, " -e ", end1, " -f ", reg.name, " -z ",
                        " -a ", maf, " --imputation-quality ", imputation.quality ,
                        " --minor-allele-frequency-fsnp ", maf, " --imputation-quality-fsnp ", imputation.quality ,
                        " --n-threads ", nthread, 
                        " --population-only > ", RASQUAL_Pop_outfile_1))

            ## print RASQUAL command (default - allele specific operation)
            print( paste0( "tabix ", out_currregion_SNPfile, ".gz ", currchr, ":", start1.vcf, "-", end1.vcf, 
                        " | ", rasqual.path, " -y ",  y.bin.path," -k ", k.bin.path, " -x ", x.bin.path, " -n ", sample.size, 
                        " -j ", n.row," -l ", n.snps," -m ", n.feat.snps ,
                        " -s ", start1, " -e ", end1, " -f ", reg.name, " -z ",
                        " -a ", maf, " --imputation-quality ", imputation.quality ,
                        " --minor-allele-frequency-fsnp ", maf, " --imputation-quality-fsnp ", imputation.quality ,
                        " --n-threads ", nthread, 
                        " --as-only > ", RASQUAL_AS_outfile_1))

            ## execute RASQUAL command (default - allele specific operation)
            system( paste0( "tabix ", out_currregion_SNPfile, ".gz ", currchr, ":", start1.vcf, "-", end1.vcf, 
                        " | ", rasqual.path, " -y ",  y.bin.path," -k ", k.bin.path, " -x ", x.bin.path, " -n ", sample.size, 
                        " -j ", n.row," -l ", n.snps," -m ", n.feat.snps ,
                        " -s ", start1, " -e ", end1, " -f ", reg.name, " -z ",
                        " -a ", maf, " --imputation-quality ", imputation.quality ,
                        " --minor-allele-frequency-fsnp ", maf, " --imputation-quality-fsnp ", imputation.quality ,
                        " --n-threads ", nthread, 
                        " --as-only > ", RASQUAL_AS_outfile_1))    

        } # end if
  
        ##============= process second interacting bin and nearby SNPs
        temp.vcf <- vcf.dat[which((vcf.dat[, 2] >= start2.vcf) & (vcf.dat[, 2] <= end2.vcf)), ]
        if (nrow(temp.vcf) > 0 ) {

            ## Adding snps to avoid error, this parameters are just for memory allocation
            n.snps <- nrow(temp.vcf) + 25
            n.feat.snps <- nrow(temp.vcf) + 5
            cat(sprintf("\n ==>> RASQUAL command - processing row : %s reg.name : %s start1 : %s end1 : %s start2 : %s end2 : %s --- second interacting bin -- start2.vcf : %s end2.vcf : %s nrow(temp.vcf) : %s n.snps : %s n.feat.snps : %s ", n.row, reg.name, start1, end1, start2, end2, start2.vcf, end2.vcf, nrow(temp.vcf), n.snps, n.feat.snps))

            ## print RASQUAL command (default - both population and allele specific operation)
            print( paste0( "tabix ", out_currregion_SNPfile, ".gz ", currchr, ":", start2.vcf, "-", end2.vcf, 
                        " | ", rasqual.path, " -y ",  y.bin.path," -k ", k.bin.path, " -x ", x.bin.path, " -n ", sample.size, 
                        " -j ", n.row," -l ", n.snps," -m ", n.feat.snps ,
                        " -s ", start2, " -e ", end2, " -f ", reg.name, " -z ",
                        " -a ", maf, " --imputation-quality ", imputation.quality ,
                        " --minor-allele-frequency-fsnp ", maf, " --imputation-quality-fsnp ", imputation.quality ,
                        " --n-threads ", nthread, 
                        " > ", RASQUAL_outfile_2))

            ## execute RASQUAL command (default - both population and allele specific operation)
            system( paste0( "tabix ", out_currregion_SNPfile, ".gz ", currchr, ":", start2.vcf, "-", end2.vcf, 
                        " | ", rasqual.path, " -y ",  y.bin.path," -k ", k.bin.path, " -x ", x.bin.path, " -n ", sample.size, 
                        " -j ", n.row," -l ", n.snps," -m ", n.feat.snps ,
                        " -s ", start2, " -e ", end2, " -f ", reg.name, " -z ",
                        " -a ", maf, " --imputation-quality ", imputation.quality ,
                        " --minor-allele-frequency-fsnp ", maf, " --imputation-quality-fsnp ", imputation.quality ,
                        " --n-threads ", nthread, 
                        " > ", RASQUAL_outfile_2))

            ## print RASQUAL command (default - population specific operation)
            print( paste0( "tabix ", out_currregion_SNPfile, ".gz ", currchr, ":", start2.vcf, "-", end2.vcf, 
                        " | ", rasqual.path, " -y ",  y.bin.path," -k ", k.bin.path, " -x ", x.bin.path, " -n ", sample.size, 
                        " -j ", n.row," -l ", n.snps," -m ", n.feat.snps ,
                        " -s ", start2, " -e ", end2, " -f ", reg.name, " -z ",
                        " -a ", maf, " --imputation-quality ", imputation.quality ,
                        " --minor-allele-frequency-fsnp ", maf, " --imputation-quality-fsnp ", imputation.quality ,
                        " --n-threads ", nthread, 
                        " --population-only > ", RASQUAL_Pop_outfile_2))

            ## execute RASQUAL command (default - population specific operation)
            system( paste0( "tabix ", out_currregion_SNPfile, ".gz ", currchr, ":", start2.vcf, "-", end2.vcf, 
                        " | ", rasqual.path, " -y ",  y.bin.path," -k ", k.bin.path, " -x ", x.bin.path, " -n ", sample.size, 
                        " -j ", n.row," -l ", n.snps," -m ", n.feat.snps ,
                        " -s ", start2, " -e ", end2, " -f ", reg.name, " -z ",
                        " -a ", maf, " --imputation-quality ", imputation.quality ,
                        " --minor-allele-frequency-fsnp ", maf, " --imputation-quality-fsnp ", imputation.quality ,
                        " --n-threads ", nthread, 
                        " --population-only > ", RASQUAL_Pop_outfile_2))

            ## print RASQUAL command (default - allele specific operation)
            print( paste0( "tabix ", out_currregion_SNPfile, ".gz ", currchr, ":", start2.vcf, "-", end2.vcf, 
                        " | ", rasqual.path, " -y ",  y.bin.path," -k ", k.bin.path, " -x ", x.bin.path, " -n ", sample.size, 
                        " -j ", n.row," -l ", n.snps," -m ", n.feat.snps ,
                        " -s ", start2, " -e ", end2, " -f ", reg.name, " -z ",
                        " -a ", maf, " --imputation-quality ", imputation.quality ,
                        " --minor-allele-frequency-fsnp ", maf, " --imputation-quality-fsnp ", imputation.quality ,
                        " --n-threads ", nthread, 
                        " --as-only > ", RASQUAL_AS_outfile_2))

            ## execute RASQUAL command (default - allele specific operation)
            system( paste0( "tabix ", out_currregion_SNPfile, ".gz ", currchr, ":", start2.vcf, "-", end2.vcf, 
                        " | ", rasqual.path, " -y ",  y.bin.path," -k ", k.bin.path, " -x ", x.bin.path, " -n ", sample.size, 
                        " -j ", n.row," -l ", n.snps," -m ", n.feat.snps ,
                        " -s ", start2, " -e ", end2, " -f ", reg.name, " -z ",
                        " -a ", maf, " --imputation-quality ", imputation.quality ,
                        " --minor-allele-frequency-fsnp ", maf, " --imputation-quality-fsnp ", imputation.quality ,
                        " --n-threads ", nthread, 
                        " --as-only > ", RASQUAL_AS_outfile_2))                    

        } # end if

    }   # end condition - sourya - remove later        

    gc()

}   # end sequential processing

## delete reference bam files
for (si in 1:length(sample.names)) {
    currdonor <- sample.names[si]
    subsethichipalignfile <- paste0(out.path.Default, 'out_', currdonor, '_reference_all_regions.bam')
    if (file.exists(subsethichipalignfile)) {
        system(paste("rm", subsethichipalignfile))
    }
}

