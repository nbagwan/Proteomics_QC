__author__ = 'Navratan Bagwan'

import re
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch
import scipy
import sklearn
# from Quantile_Normalize.quantile_norm import quantileNormalize
from scipy.cluster.hierarchy import dendrogram
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from adjustText import adjust_text
from sklearn.cluster import AgglomerativeClustering
from sklearn.impute import KNNImputer
from sklearn.model_selection import ShuffleSplit
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from itertools import combinations
import quantile_norm
import math
from collections import defaultdict
from Bio import SeqIO
import re

def plot_density_curves(df, title="Intensity distribution"):
    """
    Takes a data frame and plots the density curve of each column into one plot

    input: pandas dataframe (all columns are plotted)
    title: str, title of the plot

    returns: figure
    """
    # set plot style to white background
    sns.set_style("white")

    # create new figure
    fig = plt.figure(figsize=(8, 8))

    # iterate through all samples
    for sample in df.columns:

        # extract column of given sample
        df_sample = df[sample]

        # draw density curve
        sns.distplot(df_sample, hist=False, kde=True, kde_kws={'linewidth': 1}, label=sample)

    # Plot formatting
    plt.legend(prop={'size': 8, "weight": "bold"}, title='Sample')
    plt.title(title)
    plt.xlabel('log2(Intensity)', weight="bold")
    plt.ylabel('Density', weight="bold")

    return fig

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

def normalize_prot_int(df, method):
    """
    Normalizes the protein intensities sample wise (coulumns) by using the specified method

    df: pandas dataframe, columns=samples, rows=proteins
    method: quantile, Z or median subtraction

    returns: normalized dataframe
    """
    # choose method
    if method == "quantile":
        df_norm = quantile_norm.quantileNormalize(df)
    elif method == "Z":
        scaler = StandardScaler()
        df_norm = pd.DataFrame(scaler.fit_transform(df), index=df.index, columns=df.columns)
    elif method == "median subtraction":
        df_norm = df.apply(lambda x: x - df.median(axis=0), axis=1)
    else:
        raise ValueError(
            "Unknown normalization method specified. Choose either \"quantile\", \"Z\" or \"median subtraction\"")

    return df_norm

def plot_clustermap(df, title="Hierarchical Clustering"):
    """
    Takes data frame, applies hierarchical clustering to columns and plots the resulting clustermap with dendrogram

    input: pandas dataframe (all columns are plotted)
    title= str, title of the plot

    :return: figure
    """

    # plot clustermap
    fig = sns.clustermap(df, metric="euclidean")

    # remove y axis labels (protein names)
    plt.sca(fig.fig.axes[2])
    plt.yticks([])
    plt.ylabel("")

    # rotate x tick labels
    plt.xticks(rotation=45)

    # set title
    fig.ax_col_dendrogram.set_title(title, weight="bold")

    return fig

def plot_pearsonCorr(df, title="Pearson correlation heatmap"):
    """
        Takes data frame, calculates the pearson correlation and returns a heatmap

        input: pandas dataframe (all columns are plotted)
        title= str, title of the plot

        :return: figure
        """
    corr = df.corr()
    mask = np.triu(np.ones_like(corr, dtype=bool))
    f, ax = plt.subplots(figsize=(11,9))
    colormap = sns.cubehelix_palette(as_cmap=True)
    fig = sns.heatmap(corr, mask=mask, cmap=colormap, square=True, linewidths=0.5)

    return fig

def calculate_principal_components(df, nr_components=2):
    """
    Takes a data frame, performs a principle component analysis on the columns
    and returns PCs, explained variance and feature importance
    """
    # Separating out the features
    x = df.T

    # Separating out the target
    y = x.index.values

    # Standardizing the features
    x = StandardScaler().fit_transform(x)

    # initiate PCA model
    pca = PCA(n_components=nr_components)

    # calculate principal components
    principal_components = pca.fit_transform(x)

    # generate column names
    cols = []
    for i in range(nr_components):
        cols.append("Principal Component " + str(i+1))

    # convert principal components to data frame
    df_principal_components = pd.DataFrame(data=principal_components, columns=cols)

    # add protein names as index to principalDf
    df_principal_components.set_index(y, inplace=True)

    # get variance explained by components
    variance_explained = pca.explained_variance_ratio_

    # get feature importance
    df_feature_importance = pd.DataFrame(data=pca.components_, columns=df.index).T

    # sort features by descending PC1 values
    df_feature_importance.sort_values(by=0, inplace=True, ascending=False)

    return df_principal_components, variance_explained, df_feature_importance


def plot_principal_components(df, variance_explained, identifier=None, title="2 Component PCA", hue=True, annotate=True):
    """
    Takes a data frame of PCs and plots the top two principal components against each other
    """
    if hue is True:
        # add new group column
        # create empty list to fill up with names
        group = pd.Series(np.ones(len(df.index.values)), index=df.index.values)
        # pdb.set_trace()
        # iterate through index and extract
        group_num = 0
        for elm in identifier:
            elm1 = elm.split("|")
            # pdb.set_trace()
            for sample in df.index.values:
                if sample in elm1:
                # if len(re.findall(elm, sample)) != 0:
                #     pdb.set_trace()
                    group[sample] = "Group " + str(group_num)
                else:
                    pass
            group_num += 1

        # create new figure
        fig = plt.figure()

        # plot PCs as scatter plot
        # pdb.set_trace()
        sns.scatterplot(x="Principal Component 1", y="Principal Component 2", hue=group, palette="husl", data=df, s=300)

        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    else:
        # create new figure
        fig = plt.figure()

        # plot PCs as scatter plot
        sns.scatterplot(x="Principal Component 1", y="Principal Component 2", data=df, s=300, palette="husl", legend=False)

    # add title
    plt.title(title, weight="bold")

    # add axis labels
    plt.xlabel("Principal Component 1 (" + str(round(variance_explained[0]*100, 2)) + "%)")
    plt.ylabel("Principal Component 2 (" + str(round(variance_explained[1]*100, 2)) + "%)")

    # add annotations one by one with a loop
    if annotate is True:
        texts = []
        for line in range(len(df)):
            texts.append(plt.text(df["Principal Component 1"][line] + 0.2, df["Principal Component 2"][line],
                         df.index[line], horizontalalignment='left', size='small', color='black'))

        # adjust text to avoid overlapping with arrows
        # adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))

        # adjust text to avoid overlapping with arrows
        adjust_text(texts, expand_points=(1.5, 1.5))

    else:
        pass

    return fig

def impute_KNN(df, neighbors=5):
    """
    Imputes missing values using the kNN imputation algorithm from scikit-learn for each column separately

    df: dataframe (all columns will be imputed)
    neighbors: Number of neighboring samples to use for imputation

    returns: df_imputed, imputed, fig_histograms

    df_imputed: imputed dataframe
    imputed: mask showing which values are imputed (boolean)
    histogram: figure showing the distribution of reals and imputed values
    """
    # keep track of which values are imputed
    imputed = df.isna()

    # initiate KNN imputer
    imputer = KNNImputer(n_neighbors=neighbors, weights="distance")

    # impute missing values
    df_imputed = imputer.fit_transform(df)

    # transpose back again (see above transposing)
    df_imputed = df_imputed.copy()

    # convert to pd dataframe again
    df_imputed = pd.DataFrame(data=df_imputed, columns=df.columns, index=df.index)

    # extract imputed and real values
    df_real_values = df_imputed[~imputed]
    df_imputed_only = df_imputed[imputed]

    # # plot histograms of real and imputed values
    # set background style for the plots
    sns.set_style("white")

    # calculate how many subplots are needed for the amount of samples given (setting 4 columns as fixed)
    num_rows = math.ceil(len(df_imputed.columns)/4)

    # create figure with subplots (layout depends on number of samples)
    fig_histograms, ax = plt.subplots(num_rows, 4, figsize=(10, 6))

    # pdb.set_trace()
    # fill up all the subplots

    for i in range(12):
        # activate next subplot
        row = int(i / 4)
        col = i % 4
        plt.sca(ax[row, col])

        # fill up the first 10 subplots (or change to the number of samples you analyze)
        num_samples = 11
        if i < num_samples:
            # calculate bin edges
            bin_edges = list(np.histogram(a=df_imputed, bins=20)[1])
            # pdb.set_trace()
            # plot histogram of real values
            sns.distplot(df_real_values.iloc[:, i], kde=False, bins=bin_edges, color="#1f77b4",
                         hist_kws=dict(alpha=0.9))

            # plot histogram of imputed values
            sns.distplot(df_imputed_only.iloc[:, i], kde=False, bins=bin_edges, color="#e74c3c",
                         hist_kws=dict(alpha=0.9), ax=ax[row, col])

            # add title
            plt.title(df_real_values.columns[i])

            # set axis limit
            plt.xlim(2.5, 27.5)
            plt.ylim(0, 1100)

            # remove all ticks and labels
            # plt.yticks([], [])
            # plt.xticks([], [])
            plt.ylabel("")
            plt.xlabel("")

        elif i == num_samples:
            legend_elements = [Patch(facecolor="#1f77b4", alpha=0.9, label="Real Values"),
                               Patch(facecolor="#e74c3c", alpha=0.9, label="Imputed Values")]
            ax[row, col].legend(handles=legend_elements, loc="center")

            # remove all ticks and labels
            ax[row, col].axis('off')

        # remove empty subplots
        else:
            plt.axis('off')

    # add title to plot
    plt.suptitle("Distribution of real (blue) and imputed (red) values per sample", weight="bold", y=1.05)

    # make plot look nice
    plt.tight_layout()

    return df_imputed, imputed, fig_histograms



def impute_gaussian(df, path_input, width=0.3, downshift=1.8):
    """

    """
    # keep track of which values are imputed
    imputed = df.isna()

    # copy df
    df_imputed = df.copy()

    # iterate through samples
    for col in df.columns:

        # extract sample
        sample = df[col]

        # calculate mean and standard deviation
        mean = sample.mean()
        std = sample.std()

        # calculate downshifted mean and width of normal distribution
        m = mean - downshift * std
        s = std * width

        # sample from normal distribution if value is 0
        df_imputed[col] = sample.apply(lambda x: np.random.normal(loc=m, scale=s) if str(x) == str(np.nan) else x)

    # extract imputed and real values
    df_real_values = df_imputed[~imputed]
    df_imputed_only = df_imputed[imputed]

    # log transform both
    # df_real_values_log = df_real_values.apply(np.log2)
    # df_imputed_only_log = df_imputed_only.apply(np.log2)

    # plot histograms of real values
    fig = df_real_values.hist(sharex=True, sharey=True, figsize=(16, 10))

    # get axes of figure
    axes = [item for sublist in fig.tolist() for item in sublist]
    empty_subplots = len(df_real_values.columns) - len(axes)

    if empty_subplots is 0:
        # plot histograms of imputed values on same axes
        df_imputed_only.hist(ax=axes)
    else:
        # plot histograms of imputed values on same axes
        df_imputed_only.hist(ax=axes[:empty_subplots])

    # save figure
    # get file name
    fname_input = re.findall("[^\\\\,/]*$", path_input)[0]

    # get input folder
    folder_output = path_input[0:-len(fname_input) - 1]

    # create new path name for saving the figure
    path_fig = folder_output + "\\imputation_histograms_gaussian.pdf"

    # save figure
    plt.savefig(path_fig)

    # create new path name for saving the figure
    path_fig = folder_output + "\\imputation_histograms_gaussian.png"

    # save figure
    plt.savefig(path_fig)

    return df_imputed



def positionList(stringInput):
    outputList = []
    stringInputList = stringInput.split(")")
    for i in stringInputList:
        if "(" in i and len(i) > 1:
            value = i.replace("(", "")
            outputList.append(value)
    return outputList

def singleSTY_cout(stringInput):
    s = ""
    for i in stringInput:
        if (i.isalpha()) == True:
            s+=i
    return len(s)

def find_col_index(string, substring):
    list_index = []
    for i in string.split("\t"):
        if substring in i:
            list_index.append(string.split("\t").index(i))

    return list_index


def peptide_position(substring, position_in_protein):
    start = 0
    diff = 0
    outlist = []
    while start < len(substring) - 1:
        try:
            startInd = substring.index("(", start)
            aa = substring[startInd - 1]
            endInd = substring.index(")", start)
            val = float(substring[substring.index("(", start) + 1:endInd])
            if val >= 0.499:
                # print(aa, val, startInd)
                phospho_position = startInd + position_in_protein
                outlist.append(str(phospho_position))

            substring = substring[0:startInd] + substring[endInd + 1:]
            diff = endInd - startInd
            start = startInd
        except ValueError:
            break

    return outlist

def evidenceTableParse(evidenceFile):
    evedenceDic = {}
    with open(evidenceFile) as eviFile:
        header_names = str(next(eviFile))
        leading_protein_col = find_col_index(string=header_names, substring="Leading protein")  #### coloumn name in expand file, if it is different please change ##
        multiplicity_col = find_col_index(string=header_names, substring="Phospho (STY)")#### coloumn name in expand file, if it is different please change ##
        prob_seq_col = find_col_index(string=header_names, substring="Phospho (STY) Probabilities")#### coloumn name in expand file, if it is different please change ##
        # next(eviFile)
        # pdb.set_trace()
        for line in eviFile:
            if line != "\n":
                splits = line.split("\t")
                if splits[10].strip() == "2" or splits[10].strip() == "3":
                    unique_id = splits[10].strip() + "_" + splits[leading_protein_col[0]].strip().split(";")[0]
                    # pdb.set_trace()
                    if unique_id not in evedenceDic:
                        evedenceDic[unique_id] = [splits[prob_seq_col[0]].strip()]
                    else:
                        evedenceDic[unique_id].append(splits[prob_seq_col[0]].strip())

    return evedenceDic

def parseFasta(fastaQ):
    FastaDic = {}
    for seq_record in SeqIO.parse(fastaQ, "fasta"):
        # print(seq_record.id)
        # print(seq_record.seq)

        if seq_record.id.split("|")[1].strip() not in FastaDic:
            sequence = str(seq_record.seq)
            # pdb.set_trace()
            FastaDic[seq_record.id.split("|")[1].strip()] = str(sequence)

    return FastaDic

def check_evidence_check(id, evi_dict, sequence):
    out = ""
    try:
        list_evidence = evi_dict[id]
        for i in list_evidence:
            # pdb.set_trace()
            if sequence == re.sub("[^a-zA-Z]+", "", i.strip()) or sequence in re.sub("[^a-zA-Z]+", "", i.strip()):
                out+=i.strip()+";"

        # pdb.set_trace()
        if len(out) < 1:
            return "no evidence seq"
        else:
            return out

    except KeyError:
        out = "ID_issue"
        return out
