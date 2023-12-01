library(dplyr)
library(purrr)
library(data.table)
library(ggplot2)
library(factoextra)
library(tidyverse)

parse_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)

    if (length(args) != 2) {
        stop("Usage: Rscript top5_boxplot.R <input_directory> <output_folder>")
    }

    input_dir <- args[1]
    output_folder <- args[2]
    return(list(input_dir=input_dir, output_folder=output_folder))
}

args <- parse_args()
path <- args$input_dir
output_folder <- args$output_folder

testFunc <- function(x){
    x[order(x$V3,decreasing=TRUE), ] %>% 
    data.frame(row.names = .[[1]]) %>%
    .[!(row.names(.) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")),] %>% 
    mutate(V1=ifelse(V1 %in% head(V1,5),"Top5",'Not Top5'))%>% 
    mutate(V1=as.character(factor(V1,levels=unique(V1)))) %>%
    group_by(V1) %>% 
    summarize(SUM=sum(V3))
}

read_htseq_files  <- path %>% 
    list.files(pattern="*counts",recursive = T, full.names = TRUE) %>%
    set_names(tools::file_path_sans_ext(basename(.)))%>%
    lapply(fread) %>%
    map(~ (.x %>% select(-V2)))

box_plot <- lapply(read_htseq_files, testFunc) %>%
    rbindlist(idcol = "Reduced_name") %>%
    mutate(Reduced_name = gsub(".F.fastq.genome.htseq_counts", "", Reduced_name)) %>%
    ggplot(aes(fill=V1, y=SUM, x=Reduced_name)) + 
    geom_bar(position="fill", stat="identity")+
    theme_minimal()+
    scale_fill_manual(values = c("Top5" = "#7E22E5", "Not Top5" = "#FFFFFF"))

C_Sycamore <- path %>%
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
    fviz_pca_ind(repel = TRUE, title = "PCA with all genes")


# Save the box_plot in output folder
ggsave(paste0(output_folder, "/box_plot.jpg"), plot = box_plot, width = 10, height = 10, units = "in", dpi = 300)

# Save the PCA in output folder
ggsave(paste0(output_folder, "/pca.jpg"), plot = res.pca, width = 10, height = 10, units = "in", dpi = 300)


