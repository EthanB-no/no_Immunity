install.packages("ggplot2")
install.packages("tidyverse")
install.packages("immunarch")           # Install the package

library(immunarch)
library(ggplot2)
library(immunarch)

# Load your data from the './Data/' directory
imm_data <- repLoad("./Data2/")

repExplore(imm_data$data, "lens") %>% vis()  # Visualise the length distribution of CDR3
repClonality(imm_data$data, "homeo") %>% vis()  # Visualise the relative abundance of clonotypes


library(tidyverse)

imm_data
imm_data$meta

repOverlap(imm_data$data) %>% vis()      # Compute and visualise the most important statistics:
geneUsage(immdata$data[[1]]) %>% vis()  #     public clonotypes, gene usage, sample diversity
repDiversity(immdata$data) %>% vis(.by = "Status", .meta = immdata$meta)      # Group samples

library(dplyr)

# Assign timepoint group labels to your samples
imm_data$meta <- imm_data$meta %>%
  mutate(Group = case_when(
    grepl("Pre[-]?TX", Sample, ignore.case = TRUE) ~ "Pre-treatment",
    grepl("Post[-]?Rad", Sample, ignore.case = TRUE) ~ "Post-treatment",
    grepl("4weeks", Sample, ignore.case = TRUE) ~ "4 weeks",
    grepl("16week", Sample, ignore.case = TRUE) ~ "16 weeks",
    TRUE ~ "Unknown"
  ))

# Check if the grouping worked
table(imm_data$meta$Group)

repDiversity.methods()

packageVersion("immunarch")

library(immunarch)

div <- repDiversity(imm_data$data, "entropy")

div <- repDiversity(imm_data$data, "norm.entropy")

div <- repDiversity(imm_data$data, "shannon")

?repDiversity

# Hill numbers (effective diversity) across q = 1 to 5
hill_div <- repDiversity(imm_data$data, .method = "hill", .min.q = 1, .max.q = 5)
vis(hill_div, .by = imm_data$meta$Group)

# d50 index (clonality)
d50_div <- repDiversity(imm_data$data, .method = "d50")
vis(d50_div, .by = imm_data$meta$Group)

# Inverse Simpson (evenness)
simp_div <- repDiversity(imm_data$data, .method = "inv.simp")
vis(simp_div, .by = imm_data$meta$Group)

_________________ # Hill test troubleshooting
length(hill_div)
length(imm_data$meta$Group)
names(hill_div)

# First ensure Sample is a column in imm_data$meta
imm_meta <- imm_data$meta
imm_meta$Sample <- rownames(imm_meta)

# Now merge with hill_div by "Sample"
hill_div_merged <- merge(hill_div, imm_meta[, c("Sample", "Group")], by = "Sample")

# Now plot grouped by Group
vis(hill_div_merged, .by = "Group")

library(ggplot2)

ggplot(hill_div_merged, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot() +
  facet_wrap(~ Q, scales = "free_y") +
  labs(
    title = "Hill Diversity by Group and Q",
    y = "Diversity (Hill number)",
    x = "Group"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

str(hill_div_merged)

head(hill_div$Sample)
head(imm_meta$Sample)

imm_meta$Sample <- rownames(imm_meta)
head(imm_meta$Sample)

# Copy sample names from the first column into rownames
rownames(imm_meta) <- imm_meta[[1]]
# Now copy them into a proper column too
imm_meta$Sample <- rownames(imm_meta)
head(imm_meta$Sample)

# Use the names of the data list to set the Sample column
imm_meta <- imm_data$meta
imm_meta$Sample <- names(imm_data$data)

# Now merge
hill_div_merged <- merge(hill_div, imm_meta[, c("Sample", "Group")], by = "Sample")

library(ggplot2)

ggplot(hill_div_merged, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot() +
  facet_wrap(~ Q, scales = "free_y") +
  labs(
    title = "Hill Diversity by Group and Q",
    y = "Diversity (Hill number)",
    x = "Group"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
__________________________________
# Ensure Group is a factor
hill_div_merged$Group <- as.factor(hill_div_merged$Group)

# Run Kruskal-Wallis test for each Q
library(dplyr)

kruskal_results <- hill_div_merged %>%
  group_by(Q) %>%
  summarise(
    p_value = kruskal.test(Value ~ Group)$p.value
  )

print(kruskal_results)

__________________________________

library(vroom)
library(immunarch)
library(dplyr)

# Load immunarch data
immdata <- repLoad("Data2/")

# Add U number to each sample
for (sample_name in names(immdata$data)) {
  
  file_path <- paste0("Data2/", sample_name, ".tsv")
  manual_data <- vroom::vroom(file_path, show_col_types = FALSE)
  
  rep_data <- immdata$data[[sample_name]]
  
  joined_data <- left_join(rep_data,
                           manual_data %>% select(CDR3.aa, V.name, J.name, `U number`),
                           by = c("CDR3.aa", "V.name", "J.name"))
  
  immdata$data[[sample_name]] <- joined_data
}

# ✅ Now inspect outside the loop
head(immdata$data[[1]]$`U number`)  # See if U numbers were added
View(immdata$data[[1]])             # Inspect first sample

colnames(immdata$data[[1]])

sample_name <- names(immdata$data)[1]
file_path <- paste0("Data2/", sample_name, ".tsv")

manual_data <- vroom::vroom(file_path, show_col_types = FALSE)
rep_data <- immdata$data[[sample_name]]

# See top few values for key columns
head(rep_data[, c("CDR3.aa", "V.name", "J.name")])
head(manual_data[, c("CDR3.aa", "V.name", "J.name")])



