# IQTL (Interaction QTL)

# Detailed pipeline 
# Quantitative trait loci (QTL) associated with chromatin interaction data.


## ********
## Prior Requirements
## ********

Please check the README of the primary repository link, and install the required packages / libraries.


## ********
## Step-by-step execution
## ********

## ===========
## A. Preprocessing
## ===========


## A.0. Creating the sample list and the metadata for HiChIP datasets

Check the file *Data/DonorList_Annotated.txt* provided along with this repository. 

	- Column 1: Donor ID
	- Column 2: Sample ID. Can be kept same as the donor ID, or user can insert sample specific information. For example, if a donor has multiple HiChIP samples (different sequencing runs) then those individual run information (and corresponding sample ID) can be put here.
	- Column 3: SeqDepth: Sequencing depth of the HiChIP datasets. Initially kept blank. Once HiC-pro is executed, user can fill up this column (see details below).
	- Columns 4 to 7: Sample specific metadata or annotation information, as available for the input set of cohorts. There is no fixed format or header for these fields, can be omitted or inserted as required, depending on the input datasets.


## A.1. Align and pre-process HiChIP datasets using HiC-Pro

Check the folder *Preprocessing/HiCPro* and associated README for details on how to run HiC-pro for one HiChIP sample.

	** Note: We recommend organizing fastq files for individual samples under a common folder (*COMMONDIR*). 

		- For example, the fastq files for a given sample (suppose its name is *sample1*) would be located in the folder *COMMONDIR/sample1/data/rawdata/*.fastq.gz*. 
		- The sample names *sample1* should follow the 2nd column of the sample metadata file *Data/DonorList_Annotated.txt*.

	** Follow the configuration file of HiCPro to run the package. 

	** The output directories for individual samples containing HiCPro output would be *COMMONDIR/sample1/HiCPro/*, *COMMONDIR/sample2/HiCPro/*, where the sample names are *sample1*, *sample2*, etc.


## A.2. Fill the sequencing depth in the sample information metadata file

	** Once HiC-pro is executed for all the samples, we need to compute the sequencing depth for individual samples. In the current study, we only considered CIS valid pairs within the distance range 10 Kb - 3 Mb for calling signfificant loops using FitHiChIP.

	** In view of this, for each sample, we run the following command on the HiC-pro generated valid pairs file to obtain the sequencing depth:

		awk -F'[\t]' -v l="$lowdist" -v h="$highdist" '(($2==$5) && ($2 ~ /^chr([1-9]|2[0-2]|1[0-9])$/ ) && (($6-$3)>=l) && (($6-$3)<=h))' $inpfile

		- where, 
			lowdist = lower distance threshold = 10000 (10 Kb)
			highdist = lower distance threshold = 3000000 (3 Mb)
			We considered the autosomal chromosomes (and corresponding valid pairs) in the current study.
		
	** These sequencing depth values for individual samples need to be put into the sample metadata file *Data/DonorList_Annotated.txt*.


## A.3. Merge ChIP-seq from individual donors and call peaks using MACS2 

Check the file *Preprocessing/ChIP_Peaks/Merge_ChIP_Samples_Infer_Peaks.sh* 

	** Note: this file is a a template script to merge the ChIP-seq alignment (.bam) files from the input donors and call peaks using MACS2. 

	** Although the script shows calling peaks using various p-value or q-value thresholds, our FitHiChIP loop calls and IQTL pipeline (described below) employs peaks with q-value < 0.01 (the file *merged_donors.macs2_peaks.narrowPeak_Q0.01filt*).

	** Note: In our IQTL pipeline, we have reported results (SNPs and loops) from the autosomal chromosomes. If users want to also report the results for autosomal chromosomes only, they can check this generated peak file and discard any entry for the non-autosomal chromosomes.



## A.4. FitHiChIP loops

Check *Preprocessing/FitHiChIP* for details on how to run FitHiChIP for one HiChIP sample.

	** Note: We recommend organizing FitHiChIP loop output directories (named by specific samples / donors) under a base directory *FitHiChIPDir*. This base directory would be useful for the subsequent steps.

		- For example, if the sample names are *sample1*, *sample2*, etc.. then the output directory structure would be *FitHiChIPDir/Sample1*, *FitHiChIPDir/Sample2*, ... where *FitHiChIPDir* is the base directory name.
		
		- The sample names *Sample1*, *Sample2* follow the 2nd column of the sample metadata file *Data/DonorList_Annotated.txt*.

		- The output directory for individual samples can be set in the respective configuration files. We recommend using the following directory for a sample: *FitHiChIPDir/Sample2/FitHiChIPPrefix* where, *FitHiChIPPrefix* is a specific string which should be identical across all the samples.

			- This *FitHiChIPPrefix* value should be provided as an input to the IQTL configuration file (check the IQTL pipeline and the associated configuration file, specifically the parameter *FitHiChIPPrefix*)

	** Check the configuration file and FitHiChIP sample script in the folder *Preprocessing/FitHiChIP* to run FitHiChIP for individual samples. 

		- For details of FitHiChIP, please check its documentation https://github.com/ay-lab/FitHiChIP

	** Note: Check the ChIP-seq peak file used as the input of FitHiChIP. It should be identical across all the input samples. 

		- We used the ChIP-seq peaks derived from the step A.2, specifically peaks with q-value < 0.01 (the file *merged_donors.macs2_peaks.narrowPeak_Q0.01filt*).

	
## A.5. PCA on FitHiChIP loops and for the specific samples

	** Check the folder *Preprocessing/PCA* and the underlying README for details on how to run PCA on the derived FitHiChIP loops, to assess any donor specific variation, or to determine any outlier donors. 

		- In our data, we ran PCA but did not find any specific outliers. 

		- If users find any such outliers, they can remove the corresponding sample, and update the sample metadata file *Data/DonorList_Annotated.txt* as well. 













