# IQTL (Interaction QTL)

## Description

Quantitative trait loci (QTL) associated with chromatin interaction data.

*Note*: Although the current workflow is based on HiChIP data, in practice, the workflow can be adapted to process other types of chromatin interaction data, such as HiC, PCHi-C, Micro-C, ChIA-PET, etc.


## Developed by 

Sourya Bhattacharyya

Instructor, Center of Autoimmunity and Inflammation

La Jolla Institute for Immunology, La Jolla, CA 92037, USA


## ********
## Roadmap
## ********

August 24, 2023 - First release.


## ********
## Installation
## ********

Following packages (and associated libraries) need to be installed to run IQTL pipeline and other processing scripts:

1. R (version 3.6.1 or higher - we used R 4.1.0). 

	Once R is installed, use R prompt to install the following libraries / packages:

		install.packages(c("optparse", "ggplot2", "data.table", "splines", "fdrtool", "parallel", "tools", "plyr", "dplyr", "R.utils"))

		BiocManager::install("GenomicRanges")

		BiocManager::install("edgeR")		

2. HiC-pro (https://github.com/nservant/HiC-Pro). A package to align the HiChIP datasets with respect to the reference genome. We recommend installing the latest version.

3. FitHiChIP (https://github.com/ay-lab/FitHiChIP) and its associated dependencies, to call significant HiChIP loops. Check its documentation for the detailed installation (https://ay-lab.github.io/FitHiChIP/html/usage/installation.html). The significant HiChIP loops will be used as an input to the IQTL derivation pipeline.

4. RASQUAL (https://github.com/natsuhiko/rasqual). Download the source code and install. The directory containing the downloaded package will be used as an input of the IQTL pipeline.

5. GATK (https://gatk.broadinstitute.org/hc/en-us). Download and install. The path of the GATK executable needs to be provided as a configuration parameter.

6. bedtools (https://bedtools.readthedocs.io/en/latest/)

7. samtools (http://www.htslib.org/) pefereably the latest version (at least >= 1.10). 

	** Once samtools is installed, run the command "samtools view" in the terminal. If you do not see (-N) option (indicates that a file with the mentioned read names can be provided), you need to upgrade the samtools version.


## ********
## Step-by-step execution
## ********


## ===========
## A. Creating the sample list and the metadata
## ===========

Check the file *Data/DonorList_Annotated.txt* provided along with this repository. 

The first column is the donor ID, second column is the sample ID, and the subsequent columns are metadata for individual donors and samples.

*Note*: There is a difference between donor ID (1st column) and the sample ID (2nd column). One donor may be associated with multiple samples, like sequences in different runs, resulting HiC / HiChIP files for different runs.

*Note*: For example, given a donor D, there can be multiple samples D_Run1, D_Run2 corresponding to two sequencing runs / batches Run1 and Run2, each generating different HiChIP datasets.

*Note*: User should maintain the above mentioned meta data format of this file, specifically the first two columns. Subsequent columns (and corresponding header information) can be updated / edited / added according to the metadata of the target samples.


## ===========
## B. Preprocessing HIChIP samples and filtering HiChIP datasets and donors (if required)
## ===========

## Running HiC-Pro

Check the section *Preprocessing/HiCPro* and the associated README for details on how to run HiC-pro on individual HIChIP data.

*Note*: Once HiC-pro is run for all the samples, user should create a base directory and provide all the HiC-pro output folders (or their links) for all the samples within the base directory. 

*Note*: This base directory will be provided as an input to the IQTL pipeline.

## Running PCA

If user has any doubt regarding the batch effects of HiC / HiChIP data for individual samples, he/she can first run FitHiChIP (https://github.com/ay-lab/FitHiChIP) or FitHiC2 (https://github.com/ay-lab/FitHiC) on the given HiChIP or HiC datasets for all the samples, to identify the significant loops. 

Once these loops are identified for all the samples, user can check the README of the section *Preprocessing/PCA* to perform principal component analysis on the input samples, using the union (or top-K loops) of significant loops for individual samples, to identify any outlier samples.

*Note*: If there is any outlier sample found after PCA, user should edit the sample file (*Data/DonorList_Annotated.txt*) or create a separate copy, excluding the outlier samples.

*Note*: The modified sample file (after excluding the outlier samples) needs to be provided as an input to the IQTL pipeline.


## ===========
## C. IQTL pipeline
## ===========

Check the section *IQTL_derive* for a step-by-step description on how to use the donor / sample specific chromatin interactions and the SNP information to derive the IQTLs.

*Note*: Due to the constraints of data sharing, we currently could not share the HiChIP loops and also the genotype information. We are currently providing snapshots / examples of data. Once we upload our datasets to suitable repositories, we will share the corresponding links here.


## ********
## License
## ********

MIT license.


## ********
## Reference / Citation
## ********

To be provided.


## ********
## Support / Queries
## ********

For any queries, please e-mail: 

Sourya Bhattacharyya <sourya@lji.org>

Ferhat Ay <ferhatay@lji.org>



