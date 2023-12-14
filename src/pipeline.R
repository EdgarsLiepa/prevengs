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
library(tidyverse)
library(reticulate)
library(factoextra)


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
        map(~ if(ncol(.x) == 3) select(.x, -V2) else .x) %>%
        map(~ (.x %>% filter(!V1 %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique"))))
}

process_ht_seq_table <- function(x){
    x[order(x[[ncol(x)]], decreasing=TRUE), ] %>% 
    data.frame(row.names = .[[1]]) %>%
    .[!(row.names(.) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")),] %>% 
    mutate(V1=ifelse(V1 %in% head(V1,5),"Top5",'Not Top5'))%>% 
    mutate(V1=as.character(factor(V1,levels=unique(V1)))) %>%
    group_by(V1) %>% 
    summarize(SUM=sum(x[[ncol(x)]]))
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
        pivot_wider(names_from = Sample, values_from = last_col()) %>%
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
    tpmList <- path %>%
        list.files(pattern="*tpm.txt", recursive = TRUE, full.names = TRUE) %>%
        set_names(tools::file_path_sans_ext(basename(.))) %>%
        lapply(fread) %>%
        map(~ (.x %>% select(1,2)))

    # add first sample to marged list 
    merged_df <- tpmList[[1]] %>% rename(!!names(tpmList[1]) := V2)
    
    # Loop through the remaining data frames and merge
    for (i in 2:length(tpmList)) {
      
      # get list to marge and set collumn name
      df_to_merge <- tpmList[[i]] %>% rename(!!names(tpmList[i]) := V2)
      
      merged_df <- full_join(merged_df, df_to_merge, by = "V1")
    }
    
    # Rename V1 collumn to GeneID
    merged_df <- rename(merged_df, "GeneID" := V1)
    
    return(merged_df) 
}

# input: tpm_df - data frame with combined TPM values
# output: z_score_df - data frame with z-scores
calculate_z_scores <- function(tpm_df) {
  
  # Calculate mean value for collumns in tpm_df
  # Leave for now if per sample mean and deviation is used
  
  # mean_values <- apply(tpm_df[,-1], 2, mean)
  # std <- apply(tpm_df[,-1], 2, sd)
  
  # Calculate mean and standard deviation for data frame with all samples.
  
  tpms_matrix <- as.matrix(tpm_df[,-1])
  stDev <- sd(c(tpms_matrix))
  mean <- mean(c(tpms_matrix))
  
  # print(mean_values)
  # print(std)
  
  print("Calculating z-scores")
  
  # Calculate z-scores
  z_score_df <- (tpm_df[,-1]-mean)/stDev
  
  # # Add gene names back to the z-score data frame
  rownames(z_score_df) <- tpm_df[,1]$GeneID
  
  return(z_score_df)
  
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

    # Create plots
    box_plot <- create_box_plot(htseq_files)
    pca_plot <- create_pca_plot(sample_table)

    # make output folder if directory does not exist
    if (!dir.exists(output_folder)) {
        dir.create(output_folder)
    }
    # Save the sample table in the output folder
    print("SAve sample table in output folder /feature_table.csv")
    write.table(sample_table, file = paste0(output_folder, "/feature_table.csv"), sep = "\t", quote = FALSE, row.names = FALSE)

    # Save the box_plot in output folder
    ggsave(paste0(output_folder, "/box_plot.jpg"), plot = box_plot, width = 10, height = 10, units = "in", dpi = 300)
    # Save the PCA in output folder
    ggsave(paste0(output_folder, "/pca.jpg"), plot = pca_plot, width = 10, height = 10, units = "in", dpi = 300)


    # run python script

    file_list = list.files(pattern="*counts",recursive = T, full.names = TRUE)
    for (i in 1:length(file_list)) {
        system(paste0("python3 src/script.py ", file_list[i], " ", reference_gene_annotation, " ", output_folder))
    }

    tpms <- combine_tpms(output_folder)
    write.csv(tpms, file = paste0(output_folder, "/tpms_combined.csv"), row.names = FALSE)

   # Calculate log2FoldChange and p-values


    # Calculate Z-scores
    z_score <- calculate_z_scores(tpms)
    write.csv(z_score, file = paste0("rez", "/z-scores.csv"))

    # Calculate log2FoldChange and p-values
    # log2FoldChange <- calculate_log2FoldChange(tpms)
    # write.csv(log2FoldChange, file = paste0("rez", "/log2FoldChange.csv"))

    samples_db_path <- paste0(output_folder, "/tpms_combined.csv")
    refSTjude_path <- reference_gene_annotation
    gene_id_column <- "GeneID"
    
    print("Running python scripts outsingle/fast_zscore_estimation.py")
    system(paste0("python3 util/outsingle/fast_zscore_estimation.py ", output_folder, "/feature_table.csv"))
    
    print("Running python scripts outsingle/optht_svd_zs.py")
    system(paste0("python3 util/outsingle/optht_svd_zs.py ", output_folder, "/feature_table-fzse-zs.csv"))
    

}

main()
