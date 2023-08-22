# IQTL (Interaction QTL)

## Description
Quantitative trait loci (QTL) associated with chromatin interaction data.

	Although the current workflow is based on HiChIP data, in practice, the workflow can be adapted to process other types of chromatin interaction data, such as HiC, PCHi-C, Micro-C, ChIA-PET, etc.


## ********
## Installation
## ********

Following packages (and associated libraries) need to be installed to run IQTL pipeline and other processing scripts:

	1. R (version 3.6.1 or higher) along with the following libraries / packages:

		data.table, edgeR

	2. HiC-pro (https://github.com/nservant/HiC-Pro). A package to align the HiChIP datasets with respect to the reference genome. We recommend installing the latest version.

	3. FitHiChIP (https://github.com/ay-lab/FitHiChIP) and its associated dependencies, to call significant HiChIP loops. Check its documentation for the detailed installation. The significant HiChIP loops will be used as an input to the IQTL derivation pipeline.

	4. RASQUAL (https://github.com/natsuhiko/rasqual). Download the source code and install. The directory containing the downloaded package will be used as an input of the IQTL pipeline.

	5. GATK (https://gatk.broadinstitute.org/hc/en-us). Download and install. The path of the GATK executable needs to be provided as a configuration parameter.

	6. bedtools (https://bedtools.readthedocs.io/en/latest/)


## ********
## Step-by-step execution
## ********

## ===========
## A. Preprocessing
## ===========

## A.0. Creating the sample list and the metadata

Check the file *Data/DonorList_Annotated.txt* provided along with this repository. User should maintain the format of this file, specifying the list of samples and the associated metadata.







## ===========
## B. Running IQTL pipeline (SnakeMake)
## ===========















## A.4. Processing HiChIP alignment files (generated from HiCPro)

Check *Preprocessing/HiChIP_Alignment* for details on how to 


## ===========
## B. Derivation of IQTLs
## ===========

Check the section *IQTL_derive* for a step-by-step description on how to use the donor / sample specific FitHiChIP loops to derive the IQTLs.





## ********
## Support
## ********

For any queries, please e-mail: 

Sourya Bhattacharyya <sourya@lji.org>

Ferhat Ay <ferhatay@lji.org>




