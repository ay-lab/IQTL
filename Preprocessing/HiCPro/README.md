HiCPro
==========

Check the package https://github.com/nservant/HiC-Pro along with the documentation http://nservant.github.io/HiC-Pro/ for details on how to install and run HiCPro on a given HiChIP data.

Here we provided a sample configuration file along with a sample script to run HiC-Pro.

User can edit the configuration file according to the input data and specific parameters.

*Note*: 
	
	We assume for a given sample (suppose its name is *sample1*), the fastq files are located in the folder *sample1/data/rawdata/*.fastq.gz*. 

	This file structure is recommended by HiC-Pro and also used in the given configuration file.

*Note*: 

	We recommend organizing HiC-pro output directories for all the samples under a common directory (*COMMONDIR*). 

	For example, if the sample names are *sample1*, *sample2*, etc.. then the output directory structure would be *COMMONDIR/sample1/HiCPro/*, *COMMONDIR/sample2/HiCPro/*, etc.

	This common directory (*COMMONDIR*) will be provided as an input to the IQTL pipeline.


## License

MIT license.


## Reference / Citation

To be provided.


## Support / Queries

For any queries, please e-mail: 

Sourya Bhattacharyya <sourya@lji.org>

Ferhat Ay <ferhatay@lji.org>


