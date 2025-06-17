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

__________________________________# Uploading and merging U number table 

# assuming you have meta and patient_map loaded as data.frames
meta_new <- merge(imm_meta, U_numbers, by = "Sample", all.x = TRUE)

colnames(imm_meta)
colnames(U_numbers)

head(meta_new)

imm_data$meta <- meta_new

head(imm_data)

print(imm_data$meta)
________________________________# grouping by patient 

# Group by patient ID (U Number)
immdata_by_patient <- split(imm_data$data, imm_data$meta$`U Number`)

# Each element in immdata_by_patient is now a list of samples for that patient
names(immdata_by_patient)

library(immunarch)

imm_diversity <- repDiversity(imm_data$data, "chao1")

# Compute diversity
imm_diversity <- repDiversity(imm_data$data, "chao1")

# Add metadata
imm_diversity$PatientID <- imm_data$meta$`U Number`
imm_diversity$Group <- imm_data$meta$Group

# Plot diversity by patient over time
library(ggplot2)
ggplot(imm_diversity, aes(x = Group, y = Chao1, color = PatientID, group = PatientID)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Clonal Diversity (Chao1) Over Time by Patient",
       x = "Timepoint", y = "Chao1 Diversity")


____

# Compute diversity
imm_diversity <- repDiversity(imm_data$data, "chao1")

# Convert to data frame
imm_diversity_df <- data.frame(
  Sample = names(imm_diversity),
  Chao1 = as.numeric(imm_diversity)
)

# Merge with metadata
imm_diversity_df <- merge(imm_diversity_df, imm_data$meta, by = "Sample")

# Plot
library(ggplot2)
ggplot(imm_diversity_df, aes(x = Group, y = Chao1, color = `U Number`, group = `U Number`)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Clonal Diversity (Chao1) Over Time by Patient",
       x = "Timepoint", y = "Chao1 Diversity")

str(imm_diversity)


# 1️⃣ Convert matrix to data frame
imm_diversity_df <- as.data.frame(imm_diversity)
imm_diversity_df$Sample <- rownames(imm_diversity_df)

# 2️⃣ Merge with metadata
imm_diversity_df <- merge(imm_diversity_df, imm_data$meta, by = "Sample")

# 3️⃣ Reorder the Group factor
imm_diversity_df$Group <- factor(imm_diversity_df$Group, levels = c(
  "Pre-treatment", "Post-treatment", "4 weeks", "16 weeks"
))

# 4️⃣ Plot
library(ggplot2)
ggplot(imm_diversity_df, aes(x = Group, y = Estimator, color = `U Number`, group = `U Number`)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Clonal Diversity (Chao1) Over Time by Patient",
       x = "Timepoint", y = "Chao1 Diversity")

________#ANOVA test

install.packages("lmerTest")

library(lme4)
# Random intercept for patient
model <- lmer(Estimator ~ Group + (1 | `U Number`), data = imm_diversity_df)
summary(model)

# Check for overall effect of Group

library(lmerTest)
anova(model)

____________________

library(ggplot2)

ggplot(imm_diversity_df, aes(x = Group, y = Estimator, group = `U Number`, color = `U Number`)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "TCR Diversity over Time", y = "Estimator", x = "Timepoint")

______________________ #formatting data for GLIPH2

library(dplyr)
library(stringr)
library(tibble)

# Rename the patient column
meta_data <- meta_new %>%
  rename(Patient = `U Number`)

# Generate GLIPH2 input
gliph_input_list <- lapply(names(imm_data$data), function(s) {
  dat <- imm_data$data[[s]]
  md  <- meta_data %>% filter(Sample == s)
  
  tibble(
    CDR3b               = dat$CDR3.aa,
    Vb                  = str_replace(dat$V.name, "\\*.*", ""),  # drop *01 suffixes
    Jb                  = str_replace(dat$J.name, "\\*.*", ""),
    CDR3a               = NA,  # not used for TCRB, but required
    `Subject:condition` = paste0(md$Patient, ":", md$Group),
    Frequency           = dat$Clones
  ) %>% filter(!is.na(CDR3b) & CDR3b != "")
})

# Combine into one large GLIPH2-compatible table
gliph_input_df <- bind_rows(gliph_input_list) %>%
  select(CDR3b, Vb, Jb, CDR3a, `Subject:condition`, Frequency)

# Write to file
write.table(gliph_input_df, "gliph2_input_merged_full.txt", sep = "\t",
            quote = FALSE, row.names = FALSE)

________________________________ # interpreting GLIPH2_output_data

# Load required libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

# Use your already-loaded data frame
gliph_data <- GLIPH2_output_data

# --- Step 1: Extract patient ID and timepoint ---
gliph_data <- gliph_data %>%
  mutate(
    PatientID = str_extract(Sample, "^.*(?=:)"),  # Everything before the colon
    Timepoint = str_extract(Sample, "(?<=:).*")   # Everything after the colon
  )

# --- Step 2: Filter to one patient ---
patient_id <- "U4564519"  # Change this as needed
patient_data <- gliph_data %>%
  filter(PatientID == patient_id)

# --- Step 3: Identify top clones per timepoint ---
top_clones <- patient_data %>%
  group_by(Timepoint) %>%
  slice_max(order_by = Freq, n = 10, with_ties = FALSE) %>%
  ungroup()

# --- Step 4: Create a unique Clone ID ---
top_clones <- top_clones %>%
  mutate(CloneID = paste(TcRb, V, J, sep = "_"))

# --- Step 5: Plot clone dynamics over time ---
ggplot(top_clones, aes(x = Timepoint, y = Freq, group = CloneID, color = CloneID)) +
  geom_line(linewidth = 1) +  # <-- updated here
  geom_point(size = 2) +
  labs(
    title = paste("Top Clonal Dynamics for Patient", patient_id),
    x = "Timepoint",
    y = "Clone Frequency",
    color = "Clone ID"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))










