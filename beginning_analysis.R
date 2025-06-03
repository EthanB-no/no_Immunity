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


# vis(clonality_results, .by = "Time") # Explicitly group by Time

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


###############################################################################################

library(immunarch)
library(ggplot2)
library(dplyr)
library(stringr)
library(purrr) # For map_chr, imap_dfr
library(tibble) # For tibble::column_to_rownames, if needed


# Step 1: Load and clean VDJdb
vdjdb <- read.delim("Data/SearchTable-2025-06-03 19_34_18.85.tsv", sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(
    V_clean = sub("\\*.*", "", V),
    match_key = paste(CDR3, V_clean, sep = "_")
  )

# Step 2: Combine and clean all ImmunoSeq samples, then match with VDJdb
matched_all <- imm_data$data %>%
  imap_dfr(~ .x %>%
             mutate(
               Sample = .y, # Add a 'Sample' column to each dataframe from its name
               V_clean = sub("\\*.*", "", V.name),
               match_key = paste(CDR3.aa, V_clean, sep = "_")
             )
  ) %>%
  inner_join(vdjdb, by = "match_key", relationship = "many-to-many") %>%
  left_join(imm_data$meta, by = "Sample") # Add metadata (including 'Time')

# --- Sum and Visualize Total Frequency of Matched Clonotypes by Time ---

# Calculate the total frequency of these specific matched clonotypes for each Time point.
# This is where we sum the 'frequency' column *only for the matched clonotypes* within each 'Time' group.
total_frequency_by_time <- matched_all %>%
  group_by(Time) %>%
  summarise(
    Total_Matched_Frequency = sum(frequency), # Sum the frequencies of all matched clonotypes
    .groups = "drop" # Drop the grouping to return a clean dataframe
  ) %>%
  # Ensure the time points are ordered correctly for plotting (e.g., Pre-TX, Post-RAD, 4Weeks)
  # Adjust these levels if your 'Time' values are slightly different
  mutate(Time = factor(Time, levels = c("Pre-TX", "Post-RAD", "4Weeks")))

print(total_frequency_by_time)

# Visualize the total frequency of matched clonotypes per Time Point
p_total_freq_antigens <- ggplot(total_frequency_by_time, aes(x = Time, y = Total_Matched_Frequency, fill = Time)) +
  geom_col(color = "black") + # geom_col creates a bar plot; color adds a border to bars
  # Add labels on top of the bars, formatted to a few decimal places
  geom_text(aes(label = round(Total_Matched_Frequency, 5)), vjust = -0.5, size = 4) +
  labs(
    title = "Total Frequency of Self-Directed/Cancer-Related Clonotypes by Time Point",
    x = "Time Point",
    y = "Total Frequency of Matched Clonotypes"
  ) +
  theme_minimal() + # A clean ggplot theme
  theme(
    legend.position = "none", # No need for a legend if 'fill' is mapped to 'x'
    plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold title
    axis.text.x = element_text(angle = 45, hjust = 1) # Angle x-axis labels if needed
  )

print(p_total_freq_antigens)

# --- Save your plot ---
jpeg("total_matched_antigens_frequency_plot.jpeg", width = 900, height = 700, res = 300) # Increased resolution
print(p_total_freq_antigens)
dev.off()

