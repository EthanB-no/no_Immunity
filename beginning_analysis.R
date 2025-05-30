library(immunarch)
library(ggplot2)

# Load your data from the './Data/' directory
imm_data <- repLoad("./Data/")

repExplore(imm_data$data, "lens") %>% vis()  # Visualise the length distribution of CDR3
repClonality(imm_data$data, "homeo") %>% vis()  # Visualise the relative abundance of clonotypes



