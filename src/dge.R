suppressPackageStartupMessages({
  library(DESeq2, quietly = TRUE)
})

#' Normalize the counts in a DESeqDataSet object
#'
#' This function takes a DESeqDataSet object and normalizes the counts using the estimateSizeFactors function.
#' It returns the normalized counts.
#'
#' @param dds A DESeqDataSet object containing raw counts.
#' @return A matrix of normalized counts.
#' @export
normalize <- function(dds) {
	dds <- estimateSizeFactors(dds)
	normalized_counts <- counts(dds, normalized=TRUE)
	sizeFactors(dds)
	return(normalized_counts)
}



# Function: dge
# Description: Performs differential gene expression analysis using DESeq2 package.
# Parameters:
#   - feature_table: A matrix of count data representing gene expression levels.
#   - metadata: A data frame containing sample metadata.
# Returns: None
dge <- function(feature_table, metadata, output_folder){

	dds <- DESeqDataSetFromMatrix(countData = feature_table,
								colData = metadata,
								design = ~ sample_type)

	normalize(dds)

	dds <- DESeq(dds)

	res <- results(dds)
	head(results(dds, tidy=TRUE))

	# save results
	write.csv(res, file = paste0(output_folder, "/DseqRes.csv"))
	print(paste("Results saved in output folder:", output_folder, "DseqRes.csv"))

	# plot results
	summary(res)
}



# Function: read_sample_files
# Description: This function reads sample files and returns the data.
# Parameters:
#   - file_paths: A character vector containing the paths of the sample files.
# Returns:
#   - A list containing the data read from the sample files.

read_sample_files <- function(feature_table_path, metadata_path) {
	# Load feature_table.csv
	feature_table <- tryCatch(
		{
			read.csv(feature_table_path, header = TRUE, row.names = 1)
		},
		error = function(e) {
			stop("Error reading feature table: ", conditionMessage(e))
		}
	)
	
	# Load metadata.tsv
	metadata <- tryCatch(
		{
			read.csv(metadata_path, sep = "\t", header = TRUE, row.names = 1)
		},
		error = function(e) {
			stop("Error reading metadata: ", conditionMessage(e))
		}
	)
	
	return(list(feature_table = feature_table, metadata = metadata))
}



#' parse_args function
#'
#' This function parses command line arguments and returns a list
#' containing parameters
#'
#' @return A list with two elements: input_dir and output_folder.
#' @param None.
#' parse_args()
parse_args <- function() {
	args <- commandArgs(trailingOnly = TRUE)

	if (length(args) != 3) {
		stop("Usage: Rscript dge.R <path to input_directory with ht_seq_files> <path_to_reference_gene_annotation GTF file> <output_folder>")
	}

	input_dir <- args[1]
	metadata_path <- args[2]
	output_folder <- args[3]
	reference_gene_annotation_path <- args[4]
	return(list(input_dir=input_dir, metadata_path=metadata_path, output_folder=output_folder, reference_gene_annotation_path=reference_gene_annotation_path))
}

main <- function() {
	
	print("Running DGE analysis dge.R ...")
	args <- parse_args()
	feature_table_path <- args$input_dir
	output_folder <- args$output_folder
	metadata_path <- args$metadata_path
	reference_gene_annotation <- args$reference_gene_annotation_path

	dge()
}
