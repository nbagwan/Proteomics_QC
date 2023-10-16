__author__ = 'Navratan Bagwan'
###### USER INPUTS ARE GIVEN AT THE END OF THE SCRIPT
##LOADING ALL NECESSARY LIBRARIES
import QC_stats
import pdb
import pandas as pd
import numpy as np
import seaborn as sns

def QC_method(inputFile, outFigs, outFiles,group_names,normlization):

    path_proteingroups_log = inputFile
    norm_method = normlization
    # # set path for output
    folder_output = outFiles
    group_identifier = group_names


    #### load file to pandas data frame
    #### reading input file and saving it as dataframe.
    #### transforming the data into log2
    df_proteingroups = pd.read_csv(path_proteingroups_log, sep="\t", index_col=0)
    df_proteingroups.replace(0, np.nan, inplace=True)
    df_proteingroups_log = df_proteingroups.apply(np.log2)


    #### saving the log2 trasnformed dataframe into a txt file in output folder location
    #### making a density distribution with log2 data
    #### you can always change the output file and figure names if wanted. names between quotes.
    df_proteingroups_log.to_csv(outFiles + "\\" + "ProteinGroups-filtered-raw-log.txt", sep="\t")
    fig_density = QC_stats.plot_density_curves(df_proteingroups_log, title="Intensity distributions before "
                                                                           "" + " normalization")
    QC_stats.save_figure(fig_density, "density " + "before " + "norm", folder_output=outFigs)


    #### normalizing the log2 data by either "median or quantile method" check the comments bellow
    #### saving the normalized data into a txt file in output folder location
    #### making a density distribution with normalized data
    df_proteingroups_log_norm = QC_stats.normalize_prot_int(df_proteingroups_log, method=norm_method)
    df_proteingroups_log_norm.to_csv(outFiles + "\\" + "ProteinGroups-filtered-raw-log_norm.txt", sep="\t")
    fig_density = QC_stats.plot_density_curves(df_proteingroups_log_norm, title="Intensity distributions after "  + "medain normalization")
    QC_stats.save_figure(fig_density, "density " + "after median " + "norm", folder_output=outFigs)


    #### imputing the normalized dataset.
    #### creating the histogram with imputed and non-imputed data for every samaple
    #### saving the imputed data into a txt file in output folder location
    #### making a density distribution with imputed data
    df_proteingroups_log_imputed = QC_stats.impute_gaussian(df_proteingroups_log_norm, path_input=outFiles)
    QC_stats.save_figure(df_proteingroups_log_imputed, "imputation_histograms_gaussian", folder_output=outFigs)
    df_proteingroups_log_imputed.to_csv(outFiles + "\\" + "ProteinGroups-filtered-raw-log_norm-imputed.txt", sep="\t")
    fig_density = QC_stats.plot_density_curves(df_proteingroups_log_imputed, title="Intensity distributions after "  + " normalization_imputed")
    QC_stats.save_figure(fig_density, "density " + "after " + "norm (Imputed)", folder_output=outFigs)

    #### creating Hierarchical Clustering plot using the euclidean metric of all the samples
    fig_clustermap = QC_stats.plot_clustermap(df_proteingroups_log_imputed)
    QC_stats.save_figure(fig_clustermap, "Hierarchical Clustering", folder_output=outFigs)

    #### creating a pearson correlation heatmap
    fig_pearson = QC_stats.plot_pearsonCorr(df_proteingroups_log_imputed)
    QC_stats.save_figure(fig_pearson, "pearson correlation heat-map", folder_output=outFigs)


    #### calculating the principle componant and ploting the PCA based on the groups you have defined.
    #### The groups are defined at the end of the script
    df_proteingroups_log_imputed.reset_index()
    np.where(df_proteingroups_log_imputed.values >= np.finfo(np.float64).max)
    df_PCs, variance_explained, df_feature_importance = QC_stats.calculate_principal_components\
        (df_proteingroups_log_imputed, nr_components=2)
    fig_PCA = QC_stats.plot_principal_components(df_PCs, variance_explained, identifier=group_identifier,
                                                 title="2 Component PCA after " + " normalization", hue=True, annotate=True)
    QC_stats.save_figure(fig_PCA, "PCA_plot_after_median-norm" , folder_output=outFigs)


####### NOTE #######
### parameters bellow needs to be adusted according to your files
## inputfile: path to your inputfile.

## outFigs and outFiles are the path/location where you wanna save your files and figures produced by pipeline.

## groups_names = is the identifier. example bellow is when you have two groups such as control and ARVC group
## assuming you have 3 groups, in that case it would be line this:
## ["Ctrl1|Ctrl2|Ctrl3","ARVC1|ARVC2|ARVC", "HCM1|HCM2|HCM|"],
## make sure the group names are exactly what you have in your input file. check example

## normalization: you have two options 1) median subtraction and 2) quantile.
## median subtraction is often preferred

QC_method(inputFile = r"O:\CardiacProteomics\Scripts\NB_Scripts_and Data\Example_data_proteome\ProteinGroups-filtered-raw.txt",
          outFigs = r"O:\CardiacProteomics\Scripts\NB_Scripts_and Data\Example_data_proteome",
          outFiles = r"O:\CardiacProteomics\Scripts\NB_Scripts_and Data\Example_data_proteome",
          group_names= ["Ctrl1|Ctrl2|Ctrl3","ARVC1|ARVC2|ARVC3|ARVC4|ARVC5"],
          normlization="median subtraction")