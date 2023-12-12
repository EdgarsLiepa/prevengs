import sys
import argparse

# import libraries for plotting
import numpy as np 
from statsmodels.graphics.gofplots import qqplot_2samples
import scipy.stats
import matplotlib.pyplot as plt

# print files in directory 
import os


# define function Load feature Counts
def load_featureCounts(file_name):
    """
    Load gene count data from a featureCounts output file, starting from the 3rd line.
    
    Parameters:
    file_name (str): The path to the featureCounts file.
    
    Returns:
    dict: A dictionary of gene counts with gene names as keys and counts as values.
    """
    
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
            gene, id, count = line.strip().split('\t')

            if gene != '__no_feature' and gene != '__ambiguous' and gene != '__too_low_aQual' and gene != '__not_aligned' and gene != '__alignment_not_unique':
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
    # save plot to file
    plt.savefig('counts_qqplot.jpg', format='jpg', dpi=1000)   


def calculate_gene_lengths_by_id(gtf_file_path):
    """
    Calculate the lengths of genes from a GTF file using gene_id as the unique identifier.

    :param gtf_file_path: Path to the GTF file.
    :return: A dictionary with gene_ids as keys and their lengths as values.
    """
    gene_lengths = {}

    with open(gtf_file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):  # Skip header lines
                continue
            fields = line.strip().split('\t')
            if fields[2] == 'gene':  # We're only interested in lines describing genes
                gene_info = fields[8]
                gene_id = [info for info in gene_info.split(';') if 'gene_name' in info][0]
                gene_id = gene_id.split('"')[1]

                start_position = int(fields[3])
                end_position = int(fields[4])
                gene_length = end_position - start_position #+ 1  ?

                gene_lengths[gene_id] = gene_length

    return gene_lengths


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
    # save plot to file
    plt.savefig('top_10_genes.jpg', format='jpg', dpi=1000)


def plot_top_10_tpms(tpm,filename):
    """
    Plot the top 10 TPMs

    :param TPM: A dictionary with gene_ids as keys and their TPM as values.
    """
    # Sort the genes by length and get the top 10
    top_10tpms = sorted(tpm.items(), key=lambda x: x[1], reverse=True)[:10]
    gene_ids = [gene[0] for gene in top_10tpms]
    tpm_value = [gene[1] for gene in top_10tpms]

    # Plotting
    plt.figure(figsize=(10, 8))
    bars = plt.bar(gene_ids, tpm_value, color='skyblue')

    # Adding the text on top of each bar
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval + 0.5, f"{yval:.1f}", ha='center', va='bottom')

    # Adding labels and title
    plt.xlabel('Gene ID')
    plt.ylabel('Transcripts per million')
    plt.title('Top 10 TPM')
    plt.xticks(rotation=45, ha='right')  # Rotate the gene ids for better readability
    plt.tight_layout()
    plt.savefig(filename, format='jpg', dpi=1000)

    # Show the plot
    plt.show()
    plt.close()


def main(args):
    """
    Main function to load counts, calculate TPM, and plot results.
    """


    print("args.counts_file")
    print(args.counts_file)
    print("args.gtf_file")
    print(args.gtf_file)

    # Load featureCounts
    try:
        counts = load_featureCounts_htseq(args.counts_file)
        print("Counts loaded successfully.")
    except Exception as e:
        print(f"Error loading counts: {e}")
        sys.exit(1)

    # Calculate gene lengths
    try:
        gene_lengths = calculate_gene_lengths_by_id(args.gtf_file)
        print("Gene lengths calculated successfully.")
    except Exception as e:
        print(f"Error calculating gene lengths: {e}")
        sys.exit(1)


    # get files in directory with .txt extension


    # get sample name from file name
    sample_name = os.path.basename(args.counts_file).split('.')[0]
    
    print(sample_name)
    # Calculate TPM
    tpm = calculateTPM(counts, gene_lengths)

    # Plot the top 10 TPMs
    print("Plotting top 10 TPMs")
    plot_top_10_tpms(tpm, args.out_dir+"/"+sample_name+'_top_10_tpm.jpg')

    # Save TPM to file
    print("Saving TPM to file "+args.out_dir+"/"+sample_name+'_tpm.txt')
    with open(args.out_dir+"/"+sample_name+'_tpm.txt', 'w') as f:
        for gene in tpm:
            f.write(f"{gene}\t{tpm[gene]}\n")
    



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process transcriptome featureCounts.")
    
    parser.add_argument("counts_file", help="The featureCounts file to process.")
    parser.add_argument("gtf_file", help="The GTF file to use for gene length calculation.")
    parser.add_argument("out_dir", help="The directory to save the results in.")

    args = parser.parse_args()
    
    main(args)
    