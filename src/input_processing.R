suppressPackageStartupMessages({
library(dplyr)
library(purrr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(reticulate)
library(factoextra)
})


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

#' Process HT-seq table
#'
#' This function processes an HT-seq table by ordering it based on the last column in descending order,
#' removing specific rows, categorizing the first column as "Top5" or "Not Top5", and summarizing the
#' last column by category.
#'
#' @param x An HT-seq table
#' @return A processed HT-seq table with summarized values by category
#' @examples
#' process_ht_seq_table(data.frame(V1 = c("A", "B", "C"), V2 = c(10, 20, 30)))
process_ht_seq_table <- function(x){
    x[order(x[[ncol(x)]], decreasing=TRUE), ] %>% 
    data.frame(row.names = .[[1]]) %>%
    .[!(row.names(.) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")),] %>% 
    mutate(V1=ifelse(V1 %in% head(V1,5),"Top5",'Not Top5'))%>% 
    mutate(V1=as.character(factor(V1,levels=unique(V1)))) %>%
    group_by(V1) %>% 
    summarize(SUM=sum(x[[ncol(x)]]))
}

#' Create sample table
#'
#' This function creates a sample table by combining multiple HT-seq files into a single table,
#' pivoting the table to have samples as columns, removing specific rows, and setting the row names.
#'
#' @param htseq_files A list of HT-seq files
#' @return A sample table with samples as columns and specific rows removed
#' @examples
#' create_sample_table(list("file1.txt", "file2.txt"))
create_sample_table <- function(htseq_files) {
    htseq_files %>%
        rbindlist(idcol = "Sample") %>%
        pivot_wider(names_from = Sample, values_from = last_col()) %>%
        `row.names<-`(., NULL) %>% 
        column_to_rownames(var = "V1") %>%
        .[!(row.names(.) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")),]

    
}

# Function: create_box_plot
# Description: Creates a box plot of the summarized expression values for each sample.
# Parameters:
#   - htseq_files: A list of HTSeq count files.
# Returns: A ggplot object representing the box plot.

create_box_plot <- function(htseq_files, output_folder) {
    box_plot <- lapply(htseq_files, process_ht_seq_table) %>%
        rbindlist(idcol = "Reduced_name") %>%
        mutate(Reduced_name = gsub(".fastq.genome.htseq_counts", "", Reduced_name)) %>%
        ggplot(aes(fill=V1, y=SUM, x=Reduced_name)) + 
        geom_bar(position="fill", stat="identity")+
        theme_minimal()+
        scale_fill_manual(values = c("Top5" = "#7E22E5", "Not Top5" = "#FFFFFF"))

    ggsave(paste0(output_folder, "/box_plot.jpg"), plot = box_plot, width = 10, height = 10, units = "in", dpi = 300)
}

# Function: create_pca_plot
# Description: Creates a PCA plot using the transcript count table.
# Parameters:
#   - transcript_count_table: A data frame containing the transcript counts.
# Returns: A ggplot object representing the PCA plot.

create_pca_plot <- function(transcript_count_table, output_folder) {
    pca_plot <- prcomp(t(transcript_count_table)) %>%
        fviz_pca_ind(repel = TRUE, title = "PCA with all genes",
        geom="point")


    ggsave(paste0(output_folder, "/pca_plot.jpg"), plot = pca_plot, width = 10, height = 10, units = "in", dpi = 300)
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
