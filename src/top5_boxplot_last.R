
library(dplyr)
library(purrr)
library(data.table)
library(ggplot2)
library(factoextra)


parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) != 2) {
    stop("Usage: Rscript top5_boxplot.R <input_directory> <output_file>")
  }

  input_dir <- args[1]
  output_file <- args[2]
  return(list(input_dir=input_dir, output_file=output_file))
}


args <- parse_args()
  path <- args$input_dir
  output_file <- args$output_file


# in one pipeline:

#Create a function for the code that we want to run on each file
testFunc<-function(x){
  #Order the rows by the counts column
  x[order(x$V3,decreasing=TRUE), ] %>% 
    #Convert to a data.frame, and use the first column as the row names
    data.frame(row.names = .[[1]]) %>%
    #Remove rows with these names
    .[!(row.names(.) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")),] %>% 
    #Add a column for Top5 or Not Top5
    mutate(V1=ifelse(V1 %in% head(V1,5),"Top5",'Not Top5'))%>% 
    #Convert the V1 column to a factor, with the levels in the same order as they appear in the column
    mutate(V1=as.character(factor(V1,levels=unique(V1)))) %>%
    #Group by the V1 column
    group_by(V1) %>% 
    #Summarize to get the sum of the V3 column
    summarize(SUM=sum(V3))
}

#Get the list of files that we want to run the function on
read_htseq_files  <- path %>% 
  # get csvs full paths. (?i) is for case insentitive
  list.files(pattern="*counts",recursive = T, full.names = TRUE) %>%
  #Set the names of the files to be the file name without the extension
  set_names(tools::file_path_sans_ext(basename(.)))%>%
  #Read in the files as data.tables
  lapply(fread) %>%
  #Run the testFunc function on each data.table
  map(~ (.x %>% select(-V2)))


#Get the results into a data.frame
box_plot <- lapply(read_htseq_files, testFunc) %>%
  #Combine the list of data.frames into one data.frame
  rbindlist(idcol = "Reduced_name") %>%
  #Remove the .F.fastq.genome.htseq_counts from the Reduced_name column
  mutate(Reduced_name = gsub(".F.fastq.genome.htseq_counts", "", Reduced_name)) %>%
  #Create box the plot
  ggplot(aes(fill=V1, y=SUM, x=Reduced_name)) + 
  geom_bar(position="fill", stat="identity")+
  theme_minimal()+
  scale_fill_manual(values = c("Top5" = "#7E22E5", "Not Top5" = "#FFFFFF"))


res.pca <- prcomp(t(read_htseq_files)) %>%
  fviz_pca_ind(repel = TRUE,     # Avoid text overlapping
              title = "PCA with all genes")



pdf(output_file)
#Save as PDF
box_plot
