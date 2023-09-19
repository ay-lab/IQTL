IQTL (Interaction QTL)
==========================

QTLs (SNPs) or quantitative trait loci associated with chromatin interactions.

*Note*: 

1. Current workflow derives QTLs associated with HiChIP chromatin interactions (details below)

2. The workflow can be adapted to process other types of chromatin interaction data, such as Hi-C, PCHi-C, Micro-C, ChIA-PET, etc.


## Developed by 

Sourya Bhattacharyya

Instructor, Center of Autoimmunity and Inflammation

La Jolla Institute for Immunology, La Jolla, CA 92037, USA


## Background on eQTL

Conventional eQTL (expression quantitative trait loci) studies derive the SNPs associated with the genotype-dependent change of gene expression for a given set of samples, with respect to a specific tissue or cell type.

The bi-allelic SNPs (suppose for a given SNP, we denote its alleles as *X* and *x* where *X* is the reference allele and *x* is the alternate allele) are only considered for eQTL (or any QTL analysis). Samples are grouped according to their genotype: reference homozygous (*X|X*), heterozygous (*X|x*) and alternate homozygous (*x|x*).

The trend of gene expression for different groups of samples according to their genotypes are then assessed for statistical significance, and the significant SNPs are reported as eQTLs. By default, SNPs within 1 Mb from a given gene (TSS or transcription start site) are only assessed for significance.

Below are some of the reference studies in eQTL (and corresponding eQTL databases) that users can look into before diving into the concept of chromatin interaction and IQTLs.

1. GTEx consortium, Science 2020 [Link](https://pubmed.ncbi.nlm.nih.gov/32913098/)

2. DICE eQTL study, Schmiedel et al. Cell 2018 [Link](https://pubmed.ncbi.nlm.nih.gov/30449622/)

3. matrixQTL, a method to derive the significant SNPs associated with genotype-dependent change of gene expression [Link](https://pubmed.ncbi.nlm.nih.gov/22492648/)


## Background on chromatin interactions

Chromatin conformation capture (3C) techniques capture the regulatory interactions between promoters and regulatory elements (like enhancers) which can often be spatially distal (even having more than 1 Mb distance) but comes as 3D-proximal to the promoters due to the loop extrusion principle. Using genome-wide 3C protocols and computational methods, we can decipher the statistically significant regulatory interactions and corresponding interacting Enhancer-Promoter pairs.

Two 3C protocols are widely applied in the reference studies: 
	
1. Hi-C capturing genome-wide all-to-all interactions between every possible fragments. Requires very high sequencing depth. 

2. HiChIP: capturing genome-wide regulatory interactions subject to a particular protein or histone modifications of interest. Requires much lower sequencing depth and has better precision in capturing regularory interactome.

Below are some of the reference studies in chromatin interactions that user can refer for understanding the basics on biological protocols and computational approaches to decode the E-P interactions.

1. *in-situ* HiC protocol: Rao et al. Cell 2014 [Link](https://pubmed.ncbi.nlm.nih.gov/25497547/)

2. HiChIP Protocol: Mumbach et al. Nature Methods 2016 [Link](https://pubmed.ncbi.nlm.nih.gov/27643841/)

3. HiChIP Protocol applied on various immune cells and integrating with GWAS: Mumbach et al. Nature Genetics 2017  [Link](https://pubmed.ncbi.nlm.nih.gov/28945252/)

4. FitHiC - method to identify statistically significant Hi-C interactions [Genome Research, 2014](https://pubmed.ncbi.nlm.nih.gov/24501021/) [Nature Protocols, 2020](https://pubmed.ncbi.nlm.nih.gov/31980751/)

5. FitHiChIP - method to identify statistically significant HiChIP interactions [Nature Communications, 2019](https://pubmed.ncbi.nlm.nih.gov/31530818/)


## Concept of interaction QTL (IQTL)

Interaction QTLs (IQTLs) refer to the SNPs (QTLs) associated with chromatin interaction strength. Both genotype-depedent changes of chromatin contact counts and allele-speific variation of chromatin contacts are taken into account in order to define the interaction QTLs.

In the current study, QTLs associated with HiChIP chromatin interactions are derived. However, the workflow can be adapted to process other types of chromatin interaction data, such as Hi-C, PCHi-C, Micro-C, ChIA-PET, etc.

For example, Fig. 1 shows a SNP rs2305479 as the IQTL for the 40 Kb chromatin loop between the genes *IKZF3* and *GSDMB*. The variation of mean normalized HiChIP contact counts according to different genotypes are indicated by the color scale. Fig. 2 shows the trend of genotype dependent HiChIP contact counts for this SNP and loop.

We employed Naive CD4 H3K27ac HiChIP data of 30 donors (Fig. 3) to derive the IQTLs. To define the IQTLs, we employed RASQUAL using both genotype-dependent variation of chromatin contacts and the allele-specific variation of HiChIP reads (Fig. 4). The outputs of RASQUAL and a separate paired-end t-test of allele-specific reads are combined and filtered to produce the final set of significant SNPs and associated HiChIP loops.

![Figure 1](https://github.com/ay-lab/IQTL/blob/main/images/IQTL_Picture3.png)

*Figure 1: Example IQTL rs2305479 associated with a 40 Kb HiChIP loop between the genes IKZF3 and GSDMB. The color scale indicates the mean normalized contact counts for different genotypes.*

![Figure 2](https://github.com/ay-lab/IQTL/blob/main/images/IQTL_Picture4.png)

*Figure 2: Trend of genotype dependent HiChIP contact counts for rs2305479 and for the 40 Kb HiChIP loop between the genes IKZF3 and GSDMB*

![Figure 3](https://github.com/ay-lab/IQTL/blob/main/images/IQTL_Picture1.png)

*Figure 3: Using Naive CD4 H3K27ac HiChIP data of 30 donors for IQTL derivation*

![Figure 4](https://github.com/ay-lab/IQTL/blob/main/images/IQTL_Picture2.png)

*Figure 4: Schematic of IQTL derivation using genotype-dependent and allele-specific statistics*


## Installation

Following packages (and associated libraries) need to be installed to run IQTL pipeline and other processing scripts:

*Note:* Installation of all these packages should take a maximum of 2 hours.

1. R (version 3.6.1 or higher - we used R 4.1.0). 

	Once R is installed, use R prompt to install the following libraries / packages:

		install.packages(c("optparse", "ggplot2", "data.table", "splines", "fdrtool", "parallel", "tools", "plyr", "dplyr", "R.utils"))

		BiocManager::install("GenomicRanges")

		BiocManager::install("edgeR")		

2. HiC-pro (https://github.com/nservant/HiC-Pro). A package to align the HiChIP datasets with respect to the reference genome. We recommend installing the latest version. Check its associated README to install the package.

3. FitHiChIP (https://github.com/ay-lab/FitHiChIP) and its associated dependencies, to call significant HiChIP loops. Check its documentation for the detailed installation (https://ay-lab.github.io/FitHiChIP/html/usage/installation.html). The significant HiChIP loops will be used as an input to the IQTL derivation pipeline.

4. RASQUAL (https://github.com/natsuhiko/rasqual). Download the source code and install according to the README. The directory containing the downloaded package will be used as an input of the IQTL pipeline.

5. GATK (https://gatk.broadinstitute.org/hc/en-us). Download and install from the mentioned link. The path of the GATK executable needs to be provided as a configuration parameter.

6. bedtools (https://bedtools.readthedocs.io/en/latest/) download and install from the mentioned link.

7. samtools (http://www.htslib.org/) pefereably the latest version (at least >= 1.10). 

	*Note:* Once samtools is installed, run the command "samtools view" in the terminal. If you do not see (-N) option (indicates that a file with the mentioned read names can be provided), you need to upgrade the samtools version.


## Step-by-step execution

### Step A: Creating the sample list and the metadata

Check the file *Data/DonorList_Annotated.txt* provided along with this repository. 

	- Column 1: donor ID

	- Column 2: sample ID

	- Subsequent columns are metadata for individual donors and samples.

*Note*: 

	- There is a difference between donor ID (1st column) and the sample ID (2nd column). 

	- One donor may be associated with multiple samples, like sequences in different runs, resulting Hi-C / HiChIP files for different runs.

	- For example, given a donor ID *D*, there can be multiple samples *D_Run1*, *D_Run2* corresponding to different sequencing runs / batches *Run1* and *Run2*, each corresponding to different HiChIP data.

*Note*: 

	- User should maintain the first 2 columns and enter corresponding information. Duplicate entries are allowed if one donor corresponds to one unique sample.

	- From the column 3, user can put custom columns (and field names) according to input data.


### Step B. Preprocessing HIChIP samples - Running HiC-Pro

Check the section *Preprocessing/HiCPro* and the associated README for details on how to run HiC-pro (to align HiChIP fastq files with respect to the reference genome) on individual HIChIP data.

*Note*: 

	- Once HiC-pro is run for all the samples, user should create a base directory and provide all the HiC-pro output folders (or their links) for all the samples within the base directory. 

	- This base directory will be provided as an input to the IQTL pipeline.


### Step C. Filtering HiChIP datasets (by PCA) (Optional)

If user has any doubt regarding the sample qualities of Hi-C / HiChIP data, or the batch effects, he/she can first run FitHiChIP (https://github.com/ay-lab/FitHiChIP) or FitHiC2 (https://github.com/ay-lab/FitHiC) on the given HiChIP or Hi-C datasets, to identify the significant loops for individual samples. 

User then can check the README of the section *Preprocessing/PCA* to perform principal component analysis (PCA) on the input samples, using the union of (or top-K) significant loops for individual samples, to identify any outliers.

*Note*: 

	- If there is any outlier sample found after PCA, user should manually edit the sample file (*Data/DonorList_Annotated.txt*) or create a separate copy, excluding the outlier samples, and provide it as an input to IQTL pipeline.


### Step D. IQTL pipeline

Check the section *IQTL_Pipeline* for a step-by-step description on how to use the donor / sample specific chromatin interactions and the SNP information to derive the IQTLs.

*Note*: 

	- Due to the constraints of data sharing, we currently could provide a dummy sample of genotype data and corresponding donor information.

	- Once we upload our datasets to suitable repositories, we will share the corresponding links here.


## License

MIT license.


## Reference / Citation

To be provided.


## Release Notes

1. August 24, 2023 - First release, with IQTL pipeline and step-by-step detailed documentation.


## Support / Queries

For any queries, please e-mail: 

Sourya Bhattacharyya <sourya@lji.org>

Ferhat Ay <ferhatay@lji.org>



