import sys

# import libraries for plotting
import numpy as np 
from statsmodels.graphics.gofplots import qqplot_2samples
import scipy.stats
import matplotlib.pyplot as plt

print("Hello from script.py Print command line arguments")
print(sys.argv)

# print files in directory 
import os
print(os.listdir())

# define function Load feature Counts
def load_featureCounts(file_name):
    counts = {}
    with open(file_name, 'r') as f:
        # start from 3rd line
        next(f)
        next(f)
        for line in f:
            gene, count = line.strip().split('\t')
            counts[gene] = int(count)
    
    return counts

# define function Load feature Counts HT-seq
def load_featureCounts_htseq(file_name):
    counts = {}
    with open(file_name, 'r') as f:
        for line in f:
            description, gene, count = line.strip().split('\t')
            counts[gene] = int(count)
    
    return counts


def load_geneLength(file_name):
    lengths = {}
    with open(file_name, 'r') as f:
        # start from 3rd line
        next(f)
        next(f)
        for line in f:
            gene, length = line.strip().split('\t')
            lengths[gene] = int(length)

    return lengths

def calculateTPM(counts, lengths):

    # calculate read per kilobase
    RPK = {}
    for gene in counts:
        # Get the counts for this gene
        if gene != '':
            # print("gene " + gene + "counts {}", counts[gene])
            count = counts[gene]
            # Get the length for this gene
            length = lengths[gene]
            if length != 0 :
                RPK[gene] = count / length

    # Calculate “per million” scaling factor.
    scale = sum(RPK.values())

    tpm = {}

    for gene in RPK:
        # Get the counts for this gene
        tpm[gene] = RPK[gene] / scale * 1000000

    return tpm
        

def plot_counts(counts1, countsOld1, title1, title2):
    rez5 = qqplot_2samples(np.array(list(counts1.values())), np.array(list(countsOld1.values())), xlabel="New STAR", ylabel="Olds STAR")

    counts1_array = np.array(list(counts1.values()))
    countsOld1_array = np.array(list(countsOld1.values()))
    pearsonr = scipy.stats.pearsonr(counts1_array, countsOld1_array)[0]

    # add legend with correlations in 3 rows 
    rez5.axes[0].legend(['Pearsonr = ' + str(pearsonr)])


    rez5.show()



def plot_top_10_genes(tpm_dict1, tpm_dict2, title1, title2):
    # sort by value
    sorted_tpm_dict1 = sorted(tpm_dict1.items(), key=lambda kv: kv[1], reverse=True)
    sorted_tpm_dict2 = sorted(tpm_dict2.items(), key=lambda kv: kv[1], reverse=True)

    # get top 10
    top_tpm_dict1 = sorted_tpm_dict1[:10]
    top_tpm_dict2 = sorted_tpm_dict2[:10]

    # print in bar plot where x is gene name and y is tpm

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    axs[0].bar([x[0] for x in top_tpm_dict1], [x[1] for x in top_tpm_dict1], color='g')
    axs[0].set_xticklabels([x[0] for x in top_tpm_dict1], rotation=90)
    axs[0].set_title(title1)

    axs[1].bar([x[0] for x in top_tpm_dict2], [x[1] for x in top_tpm_dict2], color='g')
    axs[1].set_xticklabels([x[0] for x in top_tpm_dict2], rotation=90)
    axs[1].set_title(title2)


    print(top_tpm_dict1)

    plt.show()

# Load featureCounts
counts = load_featureCounts_htseq(sys.argv[1])
print("counts")
print(counts)
