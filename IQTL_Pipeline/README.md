IQTL (Interaction QTL) pipeline
=====================================

# Prior Requirements

1. Please check the README of the primary repository link, and install the required packages / libraries.

2. Execute HiC-pro on individual samples (fastq files), and also organized the output HiCPro directories for individual samples under a common base directory (for details, please check the README associated with the section *Preprocessing/HiCPro*).


# IQTL pipeline

Here we provide a snakemake pipeline, for automated execution of the IQTLs using the input chromatin loops and SNP information.


## Step 1 (Important) - edit the configuration parameters

The file *configfile_IQTL.yaml* contains the configuration parameters. 

The entries are provided as *Key:Value* pairs. Check the corresponding entries and description, and edit accordingly.

*Note:* Missing / incorrect entries in the configuration file will result in errorneous results.


## Step 2:

Edit the file *slurm_submit.sh* (according to the SLURM-based job submission parameters) and edit the corresponding options, if required. 

- Specifically, edit the *code_path* parameter with the directory containing the source code.

- User can edit this file to make it compatible for other job submission environments like *qsub* as well.


## Step 3:

*Snakefile* contains the snakemake pipeline in details.

*Note:* 

- The parameter *configfile* contains the name of the configuration file *configfile_IQTL.yaml* (as provided along with this repository). For any other configiration file, user needs to edit this parameter.

- No other parameters need to be edited.

*Note:* Individual rules / steps of this pipeline are described below for reference:

1. *merge_HiChIP_ValidPairs*: Merge HiChIP valid pairs (CIS) for all samples of a donor. The output is one valid pairs file for a given donor.

2. *calc_seq_depth*: Calculates sequencing depth for individual donors, as the number of CIS valid pairs. Also creates one updated sample file containing the sequencing depth information, which are used later to compute sequencing-depth normalized contact counts.

3. *merge_HiChIP_Alignments*: Merge HiC-pro alignments (bowtie2 generated .bam files) for individual donors and their constituent samples.

4. *subsample_HiChIP*: Subsample the above merged alignment (.bam) files by keeping only those reads which are included in the CIS merged valid pairs file. 

	The alignments are also divided by chromosome, sorted, and indexed. 

	Such distribution helps to quickly compute the allele-specific reads during IQTL derivation.

5. *Call_FitHiChIP*: Call FitHiChIP loops using the merged CIS valid pairs file for individual donors.

6. *MasterSheet*: Generate union of signfificant FitHiChIP loops for all the donors.

7. *RASQUAL_Covariate*: Compute covariates compatible with RASQUAL from the generated master sheet of loops.

8. *run_RASQUAL*: Running RASQUAL for loops from individual chromosomes and chunks (by default, loops from individual chromosomes are distributed in 15 chunks).

9. *summarize_RASQUAL*: Once RASQUAL is run for all the donors, this step summarizes and identifies the significant SNP-loop pairs using the outputs from various models of RASQUAL and also using the paired t-test analysis from the allele-specific reads.

10. *FinalList*: This rule produces the final list of IQTLs, namely the file *Final_IQTLs/Complete_IQTL.txt* within the specified output folder.


## Step 4:

The file *cluster.json* contains the resource allocation (memory, cores, and time) for individual steps. User can look into these rules and modify the resource allotted.

## Step 5

Once the configiration file is edited, run the snakemake pipeline using the command *sbatch slurm_submit.sh* (assuming SLURM job submission).

If everything goes well, the final set of IQTLs would be listed in the file *Final_IQTLs/Complete_IQTL.txt* within the specified output folder.


# Creating genptype file

Due to the constraint of sharing donor-wise genotype information, we cannot share the complete genotype information. 

- So, we provide a sample genotype file in the folder *Data* corresponding to the set of donors *D1* to *D5*, as mentioned in the file *DonorList_Annotated.txt*.

- For each chromosome, we need separate VCF files storing the genotype information. The genotype files should be named as *chr1.vcf.gz* *chr2.vcf.gz* etc. 

- These VCF files should be compressed using bgzip and tabix indexed (using the command *tabix -p vcf chr1.vcf.gz*) as well.

- The format of these VCF files need to be maintained.

*Note:* In future, we'll share the complete chromosome-wise genotype information and VCF files.


## License

MIT license.


## Reference / Citation

To be provided.


## Support / Queries

For any queries, please e-mail: 

Sourya Bhattacharyya <sourya@lji.org>

Ferhat Ay <ferhatay@lji.org>

