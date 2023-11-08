a = 1
a
b = "Suh"

print(a)

print(b)

print(a+1)

a_arr=c(1:4)

b_arr=c("1","2","3")

sum(a_arr)


read_expression_counts <- function(folder) {

  file_list <- list.files(path=folder, full.names = TRUE)
  
  file_data_list <- list()
  
  for (file in file_list) {
    print(file)
    file_data <- read.csv(file, sep = "\t", header = FALSE)
    file_name <- basename(file)
    file_data_list[[file_name]] <- file_data  
  }

  return(file_data_list)
  
}

file_data_list_copy <- read_expression_counts("data/")

summary(file_data_list_copy$RNS_FLT3_156.F.fastq.genome.htseq_counts.txt)


z_score <- function(transcripts)
{
  # Calculate mean and standard deviation
  data_mean <- mean(transcripts$V3)
  data_sd <- sd(transcripts$V3)
  
  # Calculate z-scores
  z_scores <- (data - data_mean) / data_sd
  
  # Print z-scores
  print(z_scores)

}

# Here should be data normalization.

z_score(file_data_list_copy$RNS_FLT3_156.F.fastq.genome.htseq_counts.txt)

mean(file_data_list_copy$RNS_FLT3_156.F.fastq.genome.htseq_counts.txt[3])

