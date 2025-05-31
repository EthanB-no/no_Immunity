install.packages("ggplot2")
install.packages("tidyverse")
install.packages("immunarch")           # Install the package

library(immunarch)
library(ggplot2)
library(immunarch)

# Load your data from the './Data/' directory
imm_data <- repLoad("./Data/")

repExplore(imm_data$data, "lens") %>% vis()  # Visualise the length distribution of CDR3
repClonality(imm_data$data, "homeo") %>% vis()  # Visualise the relative abundance of clonotypes


library(tidyverse)

imm_data
imm_data$meta

repOverlap(imm_data$data) %>% vis()      # Compute and visualise the most important statistics:
geneUsage(immdata$data[[1]]) %>% vis()  #     public clonotypes, gene usage, sample diversity
repDiversity(immdata$data) %>% vis(.by = "Status", .meta = immdata$meta)      # Group samples