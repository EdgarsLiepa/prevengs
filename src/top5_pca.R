library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(factoextra)
library(funrar)
library(tidyverse)


setwd("C:/Users/nikfo/OneDrive/Рабочий стол/IHH/restjudepipeline/")


#Reading all files that matches pattern (e.g. "*short") in multiple directiories
temp = list.files(pattern="*short",recursive = T)


##Show Top5 genes


myfunc <- function(file) {
  df   <- read.table(file, header=TRUE, sep="\t")
  df[order(df$STDIN,decreasing=TRUE),] %>% 
    mutate(Geneid=ifelse(Geneid %in% head(Geneid,5),Geneid,'Not Top5'))%>% 
    mutate(Geneid=factor(Geneid,levels=unique(Geneid))) %>%
    group_by(Geneid) %>% 
    summarize(SUM=sum(STDIN))
}




files = list.files(pattern="*short",recursive = T)

out <- rbindlist(lapply(files, myfunc), fill = TRUE)


names <- as.data.frame(as.data.frame(files)[rep(seq_len(nrow(as.data.frame(files))), each = 6), ])
colnames(names) <- "Sample"

matrix <- as.data.frame(c(out,names))
matrix


matrix1 <- matrix[- grep("Not Top5", matrix$Geneid),]

matrix2 <- dcast(matrix1,Geneid ~ Sample,value.var="SUM",fill=0)
rownames(matrix2) <- matrix2$Geneid
matrix2 <- matrix2[,-1]


# PCA


res.pca <- prcomp(t(matrix2), scale = TRUE)


top5pca <- fviz_pca_ind(res.pca, # Color by the quality of representation
                        repel = TRUE     # Avoid text overlapping
)


# Save PCA