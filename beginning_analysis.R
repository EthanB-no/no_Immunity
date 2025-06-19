library(immunarch)
library(ggplot2)
library(dplyr)
library(stringr)

# Load your data from the './Data/' directory
imm_data <- repLoad("./Data/")

repExplore(imm_data$data, "lens") %>% vis()  # Visualise the length distribution of CDR3
repClonality(imm_data$data, "homeo") %>% vis()  # Visualise the relative abundance of clonotypes
repDiversity(imm_data$data) %>% vis()

###################################################################################################
imm_data$meta <- read.csv("./Data/Metadata.csv", row.names = NULL)

meta <- imm_data$meta
# Set row names of meta to the actual sample names
rownames(meta) <- meta$Sample

#########################################################################################
imm_data$meta

imm_data

diversity_results <- repDiversity(imm_data$data, .method = "raref")

#print(diversity_results) # You'll see 'pre-TX', 'post-RAD', '4weeks' rows.

# 2. Visualize these results. Since immdata_obj$meta is populated,
# vis() will automatically use it to group and color the plots if appropriate.
vis(diversity_results, .by = "Time", .meta = meta)

# If you want to explicitly tell vis() which metadata column to use for grouping/coloring:
# vis(diversity_results, .by = "Time") # Use .by if you want specific grouping from meta

# --- Same principle for Clonality ---
clonality_results <- repClonality(imm_data$data, .method = "top", .head = c(10, 100, 1000))
vis(clonality_results, .by = "Time", .meta = meta)

# --- And for Overlap ---
overlap_matrix <- repOverlap(imm_data$data, .method = "overlap", .by = "Time", meta = meta)
vis(overlap_matrix)
# The heatmap will typically label by sample names, which will correspond to your time points.

matched_all <- matched_all %>%
  left_join(meta, by = "Sample")

imm_raref <- repDiversity(imm_data$data, "raref")

p1 <- vis(imm_raref)
p2 <- vis(imm_raref, .by = "Time", .meta = meta)

p1 + p2

# Open a JPEG device
jpeg("raref_plot_1.jpeg", width = 1400, height = 900, quality = 100)

p1 

# Close the device and write the file
dev.off()

# Open a JPEG device
jpeg("raref_plot_2.jpeg", width = 1400, height = 900, quality = 100)

p2

# Close the device and write the file
dev.off()


library(dplyr)
library(purrr)
library(tibble)

# Step 1: Load and clean VDJdb
vdjdb <- read.delim("Data/SearchTable-2025-06-03 19_34_18.85.tsv", sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(
    V_clean = sub("\\*.*", "", V),
    match_key = paste(CDR3, V_clean, sep = "_")
  )

# Step 2: Combine and clean all ImmunoSeq samples with pipe
matched_all <- imm_data$data %>%
  imap_dfr(~ .x %>%
             mutate(
               Sample = .y,
               V_clean = sub("\\*.*", "", V.name),
               match_key = paste(CDR3.aa, V_clean, sep = "_")
             )
  ) %>%
  inner_join(vdjdb, by = "match_key", relationship = "many-to-many") %>%     # Step 3: Match with VDJdb
  left_join(imm_data$meta, by = "Sample", "Time")              # Step 4: Add metadata (e.g., Time)

matched_all %>%
  select(Sample, Time) %>%
  distinct() %>%
  arrange(Sample)


matched_clone_list <- matched_all %>%
  group_by(Sample) %>%
  group_split() %>%
  set_names(map_chr(., ~ .x$Sample[1])) %>%
  map(~ set_names(.x$Clones, .x$match_key))

###############################################################################################

epitope_summary <- matched_all %>%
  group_by(Sample, Epitope.gene, Time) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  )

epitope_summary <- epitope_summary %>%
  left_join(imm_data$meta%>% select(Sample, Time), by = "Time", relationship = "many-to-many")

epitope_time_summary <- epitope_summary %>%
  group_by(Time, Epitope.gene) %>%
  summarise(
    sum_clones = sum(total_clones, na.rm = TRUE),
    .groups = "drop"
  )

head(epitope_time_summary)

epitope_gene_list <- epitope_time_summary$Epitope.gene
head(epitope_gene_list)

print(epitope_gene_list)


ggplot(epitope_time_summary, aes(x = Epitope.gene, y = sum_clones, fill = Epitope.gene)) +
  geom_bar(stat = "identity", position = "dodge") + # stat="identity" uses the actual y value
  facet_wrap(~ Time, scales = "free_x") + # Create separate plots for each 'Time' point
  labs(
    title = "Self Targeted T-Cell Clonotype Counts Over Time",
    x = "Epitope Gene",
    y = "Sum of Clones",
    fill = "Epitope Gene"
  ) +
  theme_minimal() + # A clean theme for the plot
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels for readability
    plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold the title
    legend.position = "none" # Hide the legend if fill is redundant with x-axis labels
  ) +
  scale_fill_viridis_d() # Use a color-blind friendly palette for fills

######################################################################################################################
######################################################################################################################

diversity_results <- repDiversity(imm_data$data, .method = "raref")

#print(diversity_results) # You'll see 'pre-TX', 'post-RAD', '4weeks' rows.

# 2. Visualize these results. Since immdata_obj$meta is populated,
# vis() will automatically use it to group and color the plots if appropriate.
vis(diversity_results, .by = "Patient", .meta = meta)

clonality_results <- repClonality(imm_data$data, .method = "top", .head = c(10, 100, 1000))
vis(clonality_results, .by = "Patient", .meta = meta)

# --- And for Overlap ---
overlap_matrix <- repOverlap(imm_data$data, .method = "overlap", .by = "Patient", meta = meta)
vis(overlap_matrix)


# Step 1: Load and clean VDJdb
vdjdb <- read.delim("Data/SearchTable-2025-06-03 19_34_18.85.tsv", sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(
    V_clean = sub("\\*.*", "", V),
    match_key = paste(CDR3, V_clean, sep = "_")
  )

# Step 2: Combine and clean all ImmunoSeq samples with pipe
matched_all <- imm_data$data %>%
  imap_dfr(~ .x %>%
             mutate(
               Sample = .y,
               V_clean = sub("\\*.*", "", V.name),
               match_key = paste(CDR3.aa, V_clean, sep = "_")
             )
  ) %>%
  inner_join(vdjdb, by = "match_key", relationship = "many-to-many") %>%     # Step 3: Match with VDJdb
  left_join(imm_data$meta, by = "Sample", "Time")              # Step 4: Add metadata (e.g., Time)

matched_all %>%
  select(Sample, Time) %>%
  distinct() %>%
  arrange(Sample)


matched_clone_list <- matched_all %>%
  group_by(Sample) %>%
  group_split() %>%
  set_names(map_chr(., ~ .x$Sample[1])) %>%
  map(~ set_names(.x$Clones, .x$match_key))

###############################################################################################

epitope_summary <- matched_all %>%
  group_by(Sample, Epitope.gene, Time) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  )

epitope_summary <- epitope_summary %>%
  left_join(imm_data$meta%>% select(Sample, Patient), by = "Patient", relationship = "many-to-many")

epitope_time_summary <- epitope_summary %>%
  group_by(Time, Epitope.gene) %>%
  summarise(
    sum_clones = sum(total_clones, na.rm = TRUE),
    .groups = "drop"
  )

head(epitope_time_summary)

epitope_gene_list <- epitope_time_summary$Epitope.gene
head(epitope_gene_list)

print(epitope_gene_list)


ggplot(epitope_time_summary, aes(x = Epitope.gene, y = sum_clones, fill = Epitope.gene)) +
  geom_bar(stat = "identity", position = "dodge") + # stat="identity" uses the actual y value
  facet_wrap(~ Time, scales = "free_x") + # Create separate plots for each 'Time' point
  labs(
    title = "Self Targeted T-Cell Clonotype Counts Over Time",
    x = "Epitope Gene",
    y = "Sum of Clones",
    fill = "Epitope Gene"
  ) +
  theme_minimal() + # A clean theme for the plot
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels for readability
    plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold the title
    legend.position = "none" # Hide the legend if fill is redundant with x-axis labels
  ) +
  scale_fill_viridis_d() # Use a color-blind friendly palette for fills


