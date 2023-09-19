PCA
=========

Performs PCA on the derived FitHiChIP loops, to understand the donor / sample specific variation.

For details on how to run FitHiChIP on individual HiChIP data pre-processed by HiC-pro pipeline, user needs to check the FitHiChIP GitHub repository https://github.com/ay-lab/FitHiChIP and corresponding documentation.


# Step 1

First run the script *Create_MasterSheet_Donor_Specific.sh* by editing the parameters.

By default, top 10K loops from each donor (or sample, depending on the input HiChIP loops) are combined to produce a master sheet of loops.


# Step 2

Then run the script *Perform_PCA.sh* (edit the parameters as required) to perform PCA using the generated master sheet.

Here we use the sample file *DonorList_Annotated.txt* provided within the folder *Data* of this repository.

** Note: In the script *Perform_PCA.sh*, we have mentioned the sample file and the following parameters:

	--colDonor: column containing the donor information.
	
	--colListAnnot: columns (comma or colon separated list) containing the donor / sample specific metadata. These are used to mark the samples and create separate PCA plots.

By default, PCA is computed using three different metrics:

	1. Raw contact count of FitHiChIP contacts. The folder **out_PCA_RawCC** contains corresponding PCA outputs.

	2. FDR (q-value) of FitHiChIP contacts. The folder **out_PCA_QVal** contains corresponding PCA outputs.

	3. Ratio of raw contact count and expected contact count (computed by FitHiChIP). The folder **out_PCA_RawCC_vs_expCC** contains corresponding PCA outputs.

	For individual folders, specifically, check the .pdf files to see the distribution of donors in the reduced feature space and if there is any outlier donor for a specific category.


*Note:* In our data, we ran PCA but did not find any specific outliers. 

*Note:*	If users find any such outliers, they can remove the corresponding sample, and update the sample metadata file *Data/DonorList_Annotated.txt* as well. 



## License

MIT license.


## Reference / Citation

To be provided.


## Support / Queries

For any queries, please e-mail: 

Sourya Bhattacharyya <sourya@lji.org>

Ferhat Ay <ferhatay@lji.org>


