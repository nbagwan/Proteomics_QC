rm(list=ls())

#### this the source where the R package is ####
source("http://www.biostat.jhsph.edu/~kkammers/software/eupa/source.functions.r")

#### these are two packages we installed (LIMMA and Qvalue). See README file
library(limma)
library(qvalue)


#### we are setting up the path where the input file is in setwd command and later reading the file in that location.

#### folder location where the input file is 
setwd("O:/CardiacProteomics/Scripts/NB_Scripts_and Data/Example_data_proteome") 

#### name of the text file in the folder location
x = read.csv("ProteinGroups-filtered-raw-log_norm-imputed.txt",sep="\t",header=TRUE,dec=",")

#### name of the column where the protein/phospho peptide names are.
rownames(x) = x$Leading.protein 

###### next 4 lines converts data into numerical format(common practice not always necessary)
x = x[,-1]
for(i in 1:ncol(x)){
  x[,i] = as.numeric(as.character(x[,i]))
}

### now we are defining the groups we have to compare
### similar to what we did in python
### to do the Ttest we need two groups control and a treatment
### in example below we have 3 replicates for control and 5 for treatment.

#### names of the control samples
controls <- c("Ctrl1","Ctrl2","Ctrl3") 
### name of the treatment samples
treatment <- c("ARVC1", "ARVC2", "ARVC3", "ARVC4","ARVC5") 

### define design according to syntax of limma package
### number of replicate in each condition, in this case we have 3 control that is why  1,1,1. and 5 treatment hence 2,2,2,2,2.
### if we have equal number of replicate in each group, e.g 3 replicates each group the design would look like this:
### design <- model.matrix(~factor(c(1,1,1,2,2,2))) 

### with the 3 vs. 5 layout it looks like this:

design <- model.matrix(~factor(c(1,1,1,2,2,2,2,2))) 
design


### Ttest calculation
### i am also adding the negative log10 of the moderated t-test p-value column to the output file
colnames(design) <- c("Intercept", "Diff")
res.eb <- eb.fit(x[, c(controls,treatment)], design)

###  remember that we are using an moderated t-test, that is why we are using the moderated p-value.
### if you wish to use standard p-value, we will have to change the line below as this: res.eb$log10 <- -log10(res.eb$p.ord)
res.eb$log10 <- -log10(res.eb$p.mod)


### we are writing the output file, of course you can change the file name. 
res.eb.new <- cbind("Leading protein" = rownames(res.eb), res.eb)
write.table(res.eb.new, file = "LIMMA_arvc vs control.txt", sep = "\t", row.names = FALSE) ### name of the output file name.

