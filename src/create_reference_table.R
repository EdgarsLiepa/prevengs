source("./src/input_processing.R")

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

    if (length(args) != 2) {
        stop("Usage: Rscript create_reference_table.R <path to input_directory with ht_seq_files> <output folder> \n
              Input directory should contain htseq files for all samples. \n
              Output folder will contain the feature table tsv, bar plot and PCA plot. \n")
    }

    # TODO take files from folder or from pattern
    input_dir <- args[1]
    output_folder <- args[2]

    return(list(input_dir=input_dir, output_folder=output_folder))
}



main <- function()
{
    args <- parse_args()

    ht_seq_dri_path <- args$input_dir
    output_folder <- args$output_folder
    
    # read htseq files to create feature count combined table from all samples.
    htseq_files <- read_htseq_files(ht_seq_dri_path)
    sample_table <- create_sample_table(htseq_files)

    # make output folder if directory does not exist
    if (!dir.exists(output_folder)) {
        dir.create(output_folder)
    }

    # TODO update to nikitas latest version
    # Create plots
    print("Creating bar plot...")
    create_box_plot(htseq_files, output_folder)
    print("Creating PCA plot...")
    create_pca_plot(sample_table, output_folder)

    # Save the sample table in the output folder
    print("SAve sample table in output folder /feature_table.csv")
    write.table(sample_table, file = paste0(output_folder, "/feature_table.csv"), sep = "\t", quote = FALSE, row.names = TRUE)
}

main()

