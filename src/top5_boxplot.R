# This script reads in a list of files containing gene expression data, identifies the top 5 genes with the highest expression levels, and creates a stacked bar chart showing the proportion of each sample's total expression that comes from the top 5 genes.

# Inputs:
# - input_directory: a string specifying the path to the directory containing the input files
# - output_file: a string specifying the path to the output file where the bar chart will be saved

# Outputs:
# - A stacked bar chart saved to the output_file path

# Example usage:
# Rscript top5_boxplot.R /path/to/input/directory /path/to/output/file.pdf

library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(funrar)
# library(tidyverse)
# library(factoextra)

# Define a function to read in a file and return a data frame with a "Top5" column

read_file <- function(file) {
  df <- read.table(file, header=TRUE, sep="\t")
  df[order(df$STDIN,decreasing=TRUE),] %>% 
    mutate(Geneid=ifelse(Geneid %in% head(Geneid,5),"Top5",'Not Top5'))%>% 
    mutate(Geneid=factor(Geneid,levels=unique(Geneid))) %>%
    group_by(Geneid) %>% 
    summarize(SUM=sum(STDIN))
}

# Define a function to process a list of files and return a datagu frame with the results
process_files <- function(files) {
  out <- rbindlist(lapply(files, read_file), fill = TRUE)
  names <- as.data.frame(as.data.frame(files)[rep(seq_len(nrow(as.data.frame(files))), each = 2), ])
  colnames(names) <- "Sample"
  matrix <- as.data.frame(c(out,names))
  return(matrix)
}

# Define a function to create a ggplot2 bar chart and save it to a file
create_bar_chart <- function(matrix, output_file) {
  ggplot(matrix, aes(fill=Geneid, y=SUM, x=Sample)) + 
    geom_bar(position="fill", stat="identity")+
    theme_minimal()+
    scale_fill_manual(values = c("Top5" = "#7E22E5",
                                 "Not Top5" = "#FFFFFF"))
  ggsave(output_file)
}

# Define a function to parse command line arguments
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2) {
    stop("Usage: Rscript top5_boxplot.R <input_directory> <output_file>")
  }
  input_dir <- args[1]
  output_file <- args[2]
  return(list(input_dir=input_dir, output_file=output_file))
}

# Main function
main <- function() {
  # Parse command line arguments
  args <- parse_args()
  input_dir <- args$input_dir
  output_file <- args$output_file
  
  # if in input and output dir not provided quit
  if (input_dir == "" || output_file == "") {
    stop("Usage: Rscript top5_boxplot.R <input_directory> <output_file>")
  }
  

  # print the arguments
  print(paste("Input directory:", input_dir))
  print(paste("Output file:", output_file))


  # Process the files and create the bar chart
  files <- list.files(path=input_dir, pattern="*txt", recursive=TRUE, full.names=TRUE)
  print(paste("Content of current directory:", files))

  matrix <- process_files(files)
  print(paste("Files: ", matrix))
  # create_bar_chart(matrix, output_file)
}

# Call the main function
main()
