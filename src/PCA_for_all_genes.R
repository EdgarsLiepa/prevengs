
library(dplyr)
library(tidyverse)
library(data.table)
library(factoextra)

# in one pipeline:

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) != 2) {
    stop("Usage: Rscript PCA_for_all_genes.R <input_directory> <output_file>")
  }

  input_dir <- args[1]
  output_file <- args[2]
  return(list(input_dir=input_dir, output_file=output_file))
}


args <- parse_args()
  path <- args$input_dir
  output_file <- args$output_file

C_Sycamore  <- path %>%

# get csvs full paths. (?i) is for case insentitive
list.files(pattern="*counts",recursive = T, full.names = TRUE) %>%
set_names(tools::file_path_sans_ext(basename(.))) %>%
lapply(fread) %>%
map(~ (.x %>% select(-V2))) %>%
rbindlist(idcol = "Sample") %>%
pivot_wider(names_from = Sample, values_from = "V3") %>%
`row.names<-`(., NULL) %>% 
column_to_rownames(var = "V1") %>%
.[!(row.names(.) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")),]


res.pca <- prcomp(t(C_Sycamore)) %>%
  fviz_pca_ind(repel = TRUE,     # Avoid text overlapping
              title = "PCA with all genes")

# save file at ouptu
pdf(output_file)

res.pca
#SAVE pca_fullmatrix AS PDF
