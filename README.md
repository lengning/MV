# MV
Select genes with high mean high variance

 Wrapper codes to identify genes with high mean and high variance. Example to run the code (from command line):

Rscript MV.R SCexample.csv

or

Rscript MV.R SCexample.csv 100 NULL F

The input values:

The 3rd term indicates the name of the input data set. Currently the program takes csv files or tab delimited file. The input file will be treated as a tab delimited file if the suffix is not '.csv'. Rows are genes and columns are samples. Row names and column names are required in the input file.

The 4th term defines the lower threshold to select genes of interest. Genes with mean > lower threshold will be selected. If set as NULL, then no lower threshold will be defined. Default is 100.  

The 5th term defines the upper threshold to select genes of interest. Genes with mean < upper threshold will be selected. If set as NULL, then no lower threshold will be defined. Default is NULL.  

The 6th term defines whether normalization is needed. If T is specified, median-by-ratio normalization will be performed prior to analysis (default is F).

Outputs:

The high mean genes are defined as above. 
The function will then fit a linear regression on log(variance)~log(mean) on high mean genes. Genes
with variance above this line are considered as the high mean high variance genes.

XX_MV.pdf: mean vs. variance scatter plot. Selected genes are marked in green.

XX_MeanMedVar.txt: mean, median, variance estimates for each gene.

XX_HighMHighV.txt: Selected high mean high variance genes.

XX_expression_HighMHighV.txt: normalized expression of selected high mean high variance genes.
