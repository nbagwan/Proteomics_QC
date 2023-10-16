__author__ = 'Navratan Bagwan'

# load required libraries
import sys
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde
import pdb
import math

# define required functions
def save_figure(fig, name, folder_output):
    # create new path name for saving the figure as pdf
    path_fig = folder_output + "\\" + name + ".pdf"
    # save figure
    plt.savefig(path_fig, bbox_inches="tight")
    # create new path name for saving the figure as png
    path_fig = folder_output + "\\" + name + ".png"
    # save figure
    plt.savefig(path_fig, bbox_inches="tight", dpi=600)
    return None

## this is the main function which takes the parameter that are defined at the end of the program
def volcanoPlot(inputfile, outputFolder, outFileName,FC_CuttOff, pvalue_CuttOff, Colour_UP, Colour_Down, control_tag, treatment_tag, MapGenes,genesToHighlight):
    path_ttest_results = inputfile

    folder_output = outputFolder
    # # load file to pandas data frame
    df_data = pd.read_table(path_ttest_results)


    # pdb.set_trace()
    # # set protein ids as index
    df_data.set_index("Leading protein", inplace=True)
    df_data.sort_index(axis=0, inplace=True)

    # extract all significant proteins at user defined p-value cut-off calculated value in negative log10 space
    df_signif = df_data[df_data["log10"].astype("float64") >= -math.log10(pvalue_CuttOff)]  # 1.301029996
    # pdb.set_trace()

    # # extract all significant proteins with positive fold change
    df_signif_pos = df_signif[df_signif["logFC"].astype("float64") >= FC_CuttOff]
    # # extract all significant proteins with negative fold change
    df_signif_neg = df_signif[df_signif["logFC"].astype("float64") <= -FC_CuttOff]
    # pdb.set_trace()

    # reduce data frame to columns of interest
    df_data_filtered = df_data[["logFC", "log10"]]
    df_signif_pos = df_signif_pos[["logFC", "log10"]]
    df_signif_neg = df_signif_neg[["logFC", "log10"]]

    # pdb.set_trace()
    ### renaming the coloumns
    df_data_filtered.columns = ["log2(Intensity"+ str(treatment_tag) + " / "+ "Intensity " + str(control_tag), "-log10(p-value)"]
    df_signif_pos.columns = ["log2(Intensity"+ str(treatment_tag) + " / "+ "Intensity " + str(control_tag), "-log10(p-value)"]
    df_signif_neg.columns = ["log2(Intensity"+ str(treatment_tag) + " / "+ "Intensity " + str(control_tag), "-log10(p-value)"]

    # pdb.set_trace()
    #### naming sure that all the numerical data is in float
    df_data_filtered = df_data_filtered.astype("float64")
    df_signif_pos = df_signif_pos.astype("float64")
    df_signif_neg = df_signif_neg.astype("float64")

    ### if we wanna higlight some significant proteins, use True or false as parameter
    ### creating the separate dataframe for the ones we wanna highlight.
    ### if True in parameter, provide a list of genes in text file. (see example data)
    ### in the file containing the IDs to highlight, should be be the one which is in inputput file as first coloumn
    if MapGenes == True:
        geneList = pd.read_table(genesToHighlight, header=None)[0].to_list()
        df_data_filtered_new = df_data_filtered.copy()
        df_data_filtered_new.reset_index(inplace=True)
        # pdb.set_trace()
        df_match = df_data_filtered_new[df_data_filtered_new["Leading protein"].isin(geneList)]
        df_match = df_match.set_index("Leading protein")
        df_match = df_match.astype("float64")
        df_match.to_csv(outputFolder + "/" + "MatchFile.txt", sep="\t")

    num_sign_up = len(df_signif_pos)
    num_sign_down = len(df_signif_neg)
    num_sign = num_sign_up + num_sign_down


    # # create scatter plot with -log10(p-values) on Y axis and log2(difference) on X axis
    # Calculate the point density
    # pdb.set_trace()
    xy = np.vstack([df_data_filtered["log2(Intensity"+ str(treatment_tag) + " / "+ "Intensity " + str(control_tag)],
                    df_data_filtered["-log10(p-value)"]])
    z = gaussian_kde(xy)(xy)

    # pdb.set_trace()
    # get blues colormap and reduce dynamic range (only go from medium dark to dark)
    Blues = cm.get_cmap("Blues", 512)
    Blues_new = ListedColormap(Blues(np.linspace(0.5, 1, 256)))

    # set plot background to plain white
    sns.set_style("white")

    # initiate new figure
    # pdb.set_trace()
    fig, ax = plt.subplots(figsize=(7, 9))
    # pdb.set_trace()
    # create scatter plot of all values colored by density (will be the points seen for the unsignificant proteins later)
    ax.scatter(x="log2(Intensity"+str(treatment_tag) + " / " + "Intensity " + str(control_tag),
               y="-log10(p-value)",
               data=df_data_filtered,
               edgecolor="",
               c=z,
               cmap=Blues_new,
               alpha=0.5)

    # plot significant proteins with positive fold change in specified colour on top of the density scatter plot
    ax.scatter(x="log2(Intensity"+str(treatment_tag) +" / " + "Intensity " + str(control_tag),
               y="-log10(p-value)",
               data=df_signif_pos,
               edgecolor="k",
               color=Colour_UP,
               # label="Significantly higher in Control\nn = " + str(0))
               label="Significantly higher in"+ str(treatment_tag) + "\nn = " + str(num_sign_up))

    # plot significant proteins with negative fold change in specified colour on top of the density scatter plot
    ax.scatter(x="log2(Intensity"+str(treatment_tag) +" / " + "Intensity " + str(control_tag),
               y="-log10(p-value)",
               data=df_signif_neg,
               edgecolor="k",
               color=Colour_Down,
               # label = "Significantly higher in ARVC\nn = " + str(0))
               label="Significantly higher in" + str(control_tag) + "\nn = " + str(num_sign_down))


    ## plot the genes/proteins/unique ids from list to highlight
    if MapGenes == True:
        ax.scatter(x="log2(Intensity"+str(treatment_tag) +" / " + "Intensity " + str(control_tag),
                   y="-log10(p-value)",
                   data=df_match,
                   edgecolor="k",
                   color="#eded02")


    # get legend handels and labels
    handles, labels = ax.get_legend_handles_labels()

    x_limit = max(abs(min(df_data["logFC"])), max(df_data["logFC"])) * 1.05
    plt.xlim(-x_limit, x_limit)

    y_limit_upper = df_data["log10"].max() * 1.05  ### 1.05
    y_limit_lower = -0.05
    plt.ylim(y_limit_lower, y_limit_upper)

    # plt.plot([-x_limit, x_limit], [log_p_cutoff, log_p_cutoff], linestyle="dashed", color="black")
    plt.plot([-x_limit, x_limit], [-math.log10(pvalue_CuttOff), -math.log10(pvalue_CuttOff)], linestyle="dashed", color="black")
    # draw line for fold change cutoff
    plt.plot([-FC_CuttOff, -FC_CuttOff], [-0.05, y_limit_upper], linestyle="dashed", color="black")
    plt.plot([FC_CuttOff, FC_CuttOff], [-0.05, y_limit_upper], linestyle="dashed", color="black")
    # add axis labels
    plt.xlabel("log2(Intensity"+str(treatment_tag) +" / " + "Intensity " + str(control_tag), weight="bold")
    plt.ylabel("-log10(p-value)", weight="bold")

    # add tick marks
    ax.tick_params(axis='both', which='major', direction="out", reset=True, right=False, top=False)

    # add title
    plt.title("Volcano plot using " + str(int(float(pvalue_CuttOff) * 100)) + "% p-value and log2 FC =" + str(FC_CuttOff), weight="bold")

    # create legend stating number of significantly up and down as well as total number of diff. regulated proteins
    legend_elements = [handles[1],
                       handles[2],
                       Line2D([0], [0], marker='o', color='w', markerfacecolor='k', markersize=8,
                              # label="Significant proteins: \n" + str(0))]
                              label="Significant proteins: \n" + str(num_sign))]

    ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    # save figure
    name = outFileName
    save_figure(fig, name, folder_output=outputFolder)
    plt.close()


####NOTE#####
### your inputfile should have these coloumn names. if you are using the output from LIMMA script i made. you will have these names.
## logFC
## log10
## Leading protein

###### change here all the parameters according to your need #####
volcanoPlot(inputfile=r"O:\CardiacProteomics\Scripts\NB_Scripts_and Data\Example_data_proteome\LIMMA_arvc vs control.txt", ## path to the file
            outputFolder=r"O:\CardiacProteomics\Scripts\NB_Scripts_and Data\Example_data_proteome", ## output folder for figure
            outFileName="Test_script", ### figure name
            FC_CuttOff= 0.30, ### remember this foldchange is in log2 scale. (log2fc cutoff of 1 would be foldchange of 2)
            pvalue_CuttOff=0.05,
            Colour_UP= "#68bc69", ## these are colour codes you can find by this pasting '#68bc69" on google
            Colour_Down= "#af67aa",
            control_tag= "Control",
            treatment_tag="ARVC",
            MapGenes= True, ## True or False
            genesToHighlight=r"O:\CardiacProteomics\Scripts\NB_Scripts_and Data\Example_data_proteome\Genes_to_highlight.txt")
            ### if true, please provide the list of proteins in a file see example