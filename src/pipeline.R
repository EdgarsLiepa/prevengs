#
# RNA-seq pipeline Cancer gene outlier detection
# from Htseq-counts data.
#
#   INPUT:
#     -  Directory with RNA-seq Htseq-counts results files
#   OUTPUT:
#     - Boxplot of top 5 genes
#     - PCA plot
#     - Sample table
#
#   USAGE:
#     Rscript top5_boxplot.R <input_directory> <output_folder>
#
suppressPackageStartupMessages({
library(dplyr)
library(purrr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(reticulate)
library(factoextra)
})


library(progress)

# Import functions from input_processing.R
source("./src/input_processing.R")
source("./src/dge.R")


#' parse_args function
#'
#' This function parses command line arguments and returns a list
#' containing parameters
#'
#' @return A list with four elements: input_dir, output_folder, reference_gene_annotation_path, and metadata.
#' @param None.
#' parse_args()
parse_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)

    if (length(args) != 4) {
        stop("Usage: Rscript top5_boxplot.R <path to input_directory with ht_seq_files> <path_to_reference_gene_annotation GTF file> <output_folder> <metadata_file>")
    }

    input_dir <- args[1]
    reference_gene_annotation_path <- args[2]
    output_folder <- args[3]
    metadata <- args[4]
    return(list(input_dir=input_dir, output_folder=output_folder, reference_gene_annotation_path=reference_gene_annotation_path, metadata=metadata))
}


main <- function() {
    


    # parse command line arguments and assign to variables
    args <- parse_args()
    output_folder <- args$output_folder

    path <- args$input_dir
    reference_gene_annotation <- args$reference_gene_annotation_path
    metadata <- args$metadata


    # print input arguments
    print(paste0("Input directory: ", path))
    print(paste0("Output folder: ", output_folder))
    print(paste0("Reference gene annotation: ", reference_gene_annotation))
    print(paste0("Metadata: ", metadata))



    # read htseq files to create feature count combined table from all samples.
    htseq_files <- read_htseq_files(path)
    sample_table <- create_sample_table(htseq_files)

    # make output folder if directory does not exist
    if (!dir.exists(output_folder)) {
        dir.create(output_folder)
    }

    # Create plots
    print("Creating box plot...")
    create_box_plot(htseq_files, output_folder)
    print("Creating PCA plot...")
    create_pca_plot(sample_table, output_folder)


    # Save the sample table in the output folder
    print("SAve sample table in output folder /feature_table.csv")
    write.table(sample_table, file = paste0(output_folder, "/feature_table.csv"), sep = "\t", quote = FALSE, row.names = TRUE)


    # run python script
    print("Run Feature Count and TPM analysis script")

    file_list <- list.files(path = path, pattern = "*counts", recursive = TRUE, full.names = TRUE)
    pb <- progress_bar$new(total = length(file_list), format = "[:bar] :percent :eta \n")

    for (i in 1:length(file_list)) {
        system(paste0("python3 src/script.py ", file_list[i], " ", reference_gene_annotation, " ", output_folder))
        pb$tick()
    }

    tpms <- combine_tpms(output_folder)
    write.csv(tpms, file = paste0(output_folder, "/tpms_combined.csv"), row.names = FALSE)


    # # Calculate Z-scores
    # z_score <- calculate_z_scores(tpms)
    # write.csv(z_score, file = paste0(output_folder, "/z-scores.csv"))
    # print(paste0("Z-scores calculated and saved in ",output_folder,"/z-scores.csv"))


    # DGE analysis
    # read feature_table.csv and metadata.tsv
    # run dge.R   
    
    print(paste0("Reading ", metadata," ..."))
    metadata <- tryCatch(
        {
            read.csv(metadata, sep = "\t", header = TRUE, row.names = 1)
        },
        error = function(e) {
            stop("Error reading metadata: ", conditionMessage(e))
        }
    )
    print("Metadata read successfully.")
    
    print("#\n#\n# Startting DGE analysis\n#\n#")
    dge(sample_table, metadata, output_folder)


    #
    # No of 1 Outlier detection
    # 
    # run Outsingle
    #
    #
    
    print("#\n#\n# Startting No of 1 Outlier detection\n#\n#")
        

    print("Running python scripts outsingle/fast_zscore_estimation.py")
    system(paste0("python3 util/outsingle/fast_zscore_estimation.py ", output_folder, "/feature_table.csv"))

    print("Running python scripts outsingle/optht_svd_zs.py")
    system(paste0("python3 util/outsingle/optht_svd_zs.py ", output_folder, "/feature_table-fzse-zs.csv"))

    print("Running python scripts outsingle/combine_outsingle.py")
    system(paste0("python3 src/combine_outsingle.py ", "/feature_table-fzse-zs.csv rez/feature_table-fzse-zs-pv.csv ", output_folder))

    
}

main()
