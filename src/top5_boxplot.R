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


##Merge Top5 genes into 1 column


myfunc <- function(file) {
  df   <- read.table(file, header=TRUE, sep="\t")
  df[order(df$STDIN,decreasing=TRUE),] %>% 
    mutate(Geneid=ifelse(Geneid %in% head(Geneid,5),"Top5",'Not Top5'))%>% 
    mutate(Geneid=factor(Geneid,levels=unique(Geneid))) %>%
    group_by(Geneid) %>% 
    summarize(SUM=sum(STDIN))
}

files = list.files(pattern="*short",recursive = T)

out <- rbindlist(lapply(files, myfunc), fill = TRUE)

names <- as.data.frame(as.data.frame(files)[rep(seq_len(nrow(as.data.frame(files))), each = 2), ])
colnames(names) <- "Sample"


matrix <- as.data.frame(c(out,names))


# ggplot2


geom_bar <- ggplot(matrix, aes(fill=Geneid, y=SUM, x=Sample)) + 
  geom_bar(position="fill", stat="identity")+
  theme_minimal()+
  scale_fill_manual(values = c("Top5" = "#7E22E5",
                               "Not Top5" = "#FFFFFF"))


# TO save ggplot2 as file