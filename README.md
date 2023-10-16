# Proteomics_QC
### Running the QC pipeline

"QC_stats.py" contains all the steps to run the pipline and can be used to execute the pipeline on after specifiyin the input files and few parameter. All the back-end methods executed are written in "QC_stats.py" with there details description.
The script includes these steps:
1) density plot to check the distribution for dataset after log2 transformation. 2) mediansubtraction normalization and plotting the distribution. 3) Imputation using down-shift methods
and imputation histogram for evaluating percentage of imputation. 4) Density plot for checking
normal distribution after the imputation. 5) Correlation plots, clustering and principle
component analysis as final step for quality check
