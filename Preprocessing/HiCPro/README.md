HiCPro
==========

Check the package https://github.com/nservant/HiC-Pro along with the documentation http://nservant.github.io/HiC-Pro/ for details on how to install and run HiCPro on a given HiChIP data.


Here we have provided a sample configuration file along with a sample script to run HiCPro.

User can edit the configuration file according to the input data and specific parameters.

** Note: We assume for a given sample (suppose its name is *sample1*), the fastq files are located in the folder *sample1/data/rawdata/*.fastq.gz*. This file structure is recommended by HiCPro and also used in the given configuration file.

** Note: We recommend organizing HiC-pro output directories for individual samples to be grouped under a common directory (*COMMONDIR*). For example, if the sample names are *sample1*, *sample2*, etc.. then the output directory structure would be *COMMONDIR/sample1/HiCPro/*, *COMMONDIR/sample2/HiCPro/*, etc.

