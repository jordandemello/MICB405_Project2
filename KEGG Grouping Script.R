library(KEGGREST)
library(tidyverse)
library(cowplot)
library(pathview)
library(magrittr)
listDatabases()

raw_data <- read.table("SaanichInlet_100m_ORFs_ko.cleaned.txt", header=TRUE)

grouped_data <- raw_data %>%
  group_by(KO_Number) %>%
  summarize(Count=n()) %>%
  arrange(desc(Count))

ko_numbers = grouped_data[,"KO_Number", drop=FALSE]
ko_names = character(nrow(grouped_data));

grouped_data$ko_name=ko_names

for(i in 1:nrow(grouped_data)){
  ko_names[i] = keggFind("ko", ko_numbers[i,1])
  print(ko_names[i])
}

grouped_data$KO_Name = ko_names
