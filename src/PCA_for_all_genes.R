library(tidyr)
library(purrr)
library(janitor)
library(tidyverse)
library(factoextra)

# Define a function that takes a path argument
run_pca <- function(path) {
  
  # in one pipeline:
  C_Sycamore  <- path %>% 

    # get csvs full paths. (?i) is for case insentitive
    list.files(pattern="*short",recursive = T, full.names = TRUE) %>% 

    # create a named vector: you need it to assign ids in the next step.
    # and remove file extection to get clean colnames
    set_names(tools::file_path_sans_ext(basename(.))) %>% 

    # read file one by one, bind them in one df and create id column 
    map_dfr(read.table, col.names = c("Geneid", "Sum"), .id = "Sample") %>%

    # pivot to create one column for each .id
    pivot_wider(names_from = Sample, values_from = Sum)%>%
    .[-1, ]%>% 
    column_to_rownames(var="Geneid") %>%
    sapply(., as.numeric)
  
  res.pca <- prcomp(t(C_Sycamore))
  
  pca_fullmatrix <- fviz_pca_ind(res.pca, # Color by the quality of representation
                                 repel = TRUE,     # Avoid text overlapping
                                 title = "PCA with all genes"
  )
  
  #SAVE pca_fullmatrix as file
  write.csv(pca_fullmatrix, file = "pca_fullmatrix.csv")
}

# Define a function to parse command line arguments
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1) {
    stop("Please provide a single path argument")
  }
  if (!file.exists(args[1])) {
    stop("The provided path does not exist")
  }
  return(normalizePath(args[1]))
}

# Call the functions
path <- parse_args()
run_pca(path)

