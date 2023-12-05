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



library(dplyr)
library(purrr)
library(data.table)
library(ggplot2)
library(factoextra)
library(tidyverse)
library(reticulate)


#' Read HTSeq files from a given path
#'
#' This function reads HTSeq files from a specified path
#' and filter ht-seq summary rows.
#'
#' @param path The path where the HTSeq files are located.
#'
#' @return A list of processed data frames.
#'
#' @examples
#' read_htseq_files("/path/to/files")
#'
read_htseq_files  <- function(path) {
    path %>% 
        list.files(pattern="*counts", recursive = TRUE, full.names = TRUE) %>%
        set_names(tools::file_path_sans_ext(basename(.)))%>%
        lapply(fread) %>%
        map(~ (.x %>% select(-V2)))%>%
        map(~ (.x %>% filter(!V1 %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique"))))
}

process_ht_seq_table <- function(x){
    x[order(x$V3,decreasing=TRUE), ] %>% 
    data.frame(row.names = .[[1]]) %>%
    .[!(row.names(.) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")),] %>% 
    mutate(V1=ifelse(V1 %in% head(V1,5),"Top5",'Not Top5'))%>% 
    mutate(V1=as.character(factor(V1,levels=unique(V1)))) %>%
    group_by(V1) %>% 
    summarize(SUM=sum(V3))
}

create_box_plot <- function(htseq_files) {
    lapply(htseq_files, process_ht_seq_table) %>%
        rbindlist(idcol = "Reduced_name") %>%
        mutate(Reduced_name = gsub(".fastq.genome.htseq_counts", "", Reduced_name)) %>%
        ggplot(aes(fill=V1, y=SUM, x=Reduced_name)) + 
        geom_bar(position="fill", stat="identity")+
        theme_minimal()+
        scale_fill_manual(values = c("Top5" = "#7E22E5", "Not Top5" = "#FFFFFF"))
}

create_sample_table <- function(htseq_files) {
    htseq_files %>%
        rbindlist(idcol = "Sample") %>%
        pivot_wider(names_from = Sample, values_from = "V3") %>%
        `row.names<-`(., NULL) %>% 
        column_to_rownames(var = "V1") %>%
        .[!(row.names(.) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")),]
}

create_pca_plot <- function(transcript_count_table) {
    prcomp(t(transcript_count_table)) %>%
        fviz_pca_ind(repel = TRUE, title = "PCA with all genes")
}

#' Combine TPM files into a single data frame
#'
#' This function takes a path as input and searches for TPM files in the specified directory and its subdirectories.
#' It reads each TPM file into a data frame using the `fread` function from the `data.table` package.
#' The resulting data frames are then combined into a single data frame, where each column represents a TPM file.
#' The column names are derived from the file names without the file extension.
#' The first column of each data frame is excluded from the final combined data frame.
#'
#' @param path The path to the directory containing the TPM files
#' @return A data frame with combined TPM values from all the TPM files
#' @examples
#' combine_tpms("/path/to/tpm/files")
combine_tpms <- function(path) {
    path %>%
        list.files(pattern="*tpm.txt",recursive = T, full.names = TRUE) %>%
        set_names(tools::file_path_sans_ext(basename(.))) %>%
        lapply(fread) %>%
        map(~ (.x %>% select(-V1)))
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
    

    args <- parse_args()
    path <- args$input_dir
    output_folder <- args$output_folder
    reference_gene_annotation <- args$reference_gene_annotation_path
    
    # read htseq files to process
    htseq_files <- read_htseq_files(path)

    # Create data frame with feature counts
    sample_table <- create_sample_table(htseq_files)

    # tpms <- calculate_tpms(htseq_files)


    # Create plots
    box_plot <- create_box_plot(htseq_files)
    pca_plot <- create_pca_plot(sample_table)
    
    # make output folder if directory does not exist
    if (!dir.exists(output_folder)) {
        dir.create(output_folder)
    }
    # Save the sample table in the output folder
    write.csv(sample_table, file = paste0(output_folder, "/feature_table.csv"))
    
    
    # Save the box_plot in output folder
    ggsave(paste0(output_folder, "/box_plot.jpg"), plot = box_plot, width = 10, height = 10, units = "in", dpi = 300)
    # Save the PCA in output folder
    ggsave(paste0(output_folder, "/pca.jpg"), plot = pca_plot, width = 10, height = 10, units = "in", dpi = 300)
    

    # Calculate TPM values


    # Calculate log2FoldChange and p-values


    # Calculate Z-scores


    # run python script

    file_list = list.files(pattern="*counts",recursive = T, full.names = TRUE)
    for (i in 1:length(file_list)) {
        system(paste0("python3 src/script.py ", file_list[i], " ", reference_gene_annotation, " ", output_folder))
    }



}

main()
