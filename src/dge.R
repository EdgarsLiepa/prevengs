library(DESeq2)



dge <- function(){

# load feature_table.csv
feature_table <- read.csv("rez/feature_table.csv", header = TRUE, row.names = 1)


# load metadata.tsv
metadata <- read.csv("data/BKUS_SAMPLES/metadata.tsv", sep="\t", header = TRUE, row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = feature_table,
                            colData = metadata,
                            design = ~ condition)

dds <- DESeq(dds)

res <- results(dds)
head(results(dds, tidy=TRUE))

# save results
write.csv(res, file = paste0("rez", "/DseqRes.csv"))

summary(res)


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
        stop("Usage: Rscript top5_boxplot.R <path to input_directory with ht_seq_files> <path_to_reference_gene_annotation GTF file> <output_folder>")
    }

    input_dir <- args[1]
    reference_gene_annotation_path <- args[2]
    output_folder <- args[3]
    return(list(input_dir=input_dir, output_folder=output_folder, reference_gene_annotation_path=reference_gene_annotation_path))
}

main <- function() {
    

    # args <- parse_args()
    # path <- args$input_dir
    # output_folder <- args$output_folder
    # reference_gene_annotation <- args$reference_gene_annotation_path


    dge()
}