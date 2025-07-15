library(readr)
library(ggpicrust2)
library(dplyr)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(glue)

# If you want to analyze KEGG pathway abundance instead of KO within the pathway, turn ko_to_kegg to TRUE.
# KEGG pathways typically have more explainable descriptions.

# Load metadata as a tibble
# data(metadata)
metadata <- read_delim(
  "data/Bombus_2023_metadata_v2.txt",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# Adding a new category of floral species richness:
metadata$category_2 <- ifelse(metadata$flower_abundance <= 27, "low", "high")

# Metadata to compare high and low flower richness from category column
metadata_HL <- metadata %>% 
  filter(category != "Medium")

# Load KEGG pathway abundance
# data(kegg_abundance)
kegg_abundance <- ko2kegg_abundance("data/picrust2_output_0.1/KO_metagenome_out/pred_metagenome_unstrat_clean.tsv")

# Indicate the variable you want to compare
#Variable <- "Time" # Comparing by period early vs. late
Variable <- "category_2" # Comparing by sites and flower species richness

# Perform pathway differential abundance analysis (DAA) using ALDEx2 method
# Please change group to "your_group_column" if you are not using example dataset
#daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata, group = "Time", daa_method = "ALDEx2", select = NULL, reference = NULL)

daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata, group = Variable, daa_method = "ALDEx2", select = NULL, reference = NULL)

# Filter results for ALDEx2_Welch's t test method
# Please check the unique(daa_results_df$method) and choose one
daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "ALDEx2_Wilcoxon rank test", ]

# Annotate pathway results using KO to KEGG conversion
daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df, ko_to_kegg = TRUE)

# Generate pathway error bar plot
# Please change Group to metadata$your_group_column if you are not using example dataset
#p <- pathway_errorbar(abundance = kegg_abundance, daa_results_df = daa_annotated_sub_method_results_df, Group = metadata$Time, p_values_threshold = 0.05, order = "pathway_class", select = NULL, ko_to_kegg = TRUE, p_value_bar = TRUE, colors = NULL, x_lab = "pathway_name")

p <- pathway_errorbar(abundance = kegg_abundance, daa_results_df = daa_annotated_sub_method_results_df, Group = metadata[[Variable]], p_values_threshold = 0.05, order = "pathway_class", select = NULL, ko_to_kegg = TRUE, p_value_bar = TRUE, colors = NULL, x_lab = "pathway_name")

ggsave(glue("results/pathway_errorbar_{Variable}.png"), plot = p, width = 8, height = 6, units = "in", dpi = 300)

# Get significantly different pathways
top_pathways <- daa_annotated_sub_method_results_df %>%
  filter(p_adjust < 0.05) %>%
  arrange(p_adjust) %>%
  slice(1:20) %>%  # select top 20 most significant
  pull(feature)

p_20 <- pathway_errorbar(abundance = kegg_abundance, daa_results_df = daa_annotated_sub_method_results_df, Group = metadata[[Variable]], p_values_threshold = 0.05, order = "pathway_class", select = top_pathways, ko_to_kegg = TRUE, p_value_bar = TRUE, colors = NULL, x_lab = "pathway_name")
ggsave(glue("results/pathway_errorbar_{Variable}_20.png"), plot = p_20, width = 8, height = 6, units = "in", dpi = 300)


# If you want to analyze EC, MetaCyc, and KO without conversions, turn ko_to_kegg to FALSE.

# Load KO abundance as a data.frame
# data(ko_abundance)
ko_abundance <- read_delim("data/picrust2_output_0.1/KO_metagenome_out/pred_metagenome_unstrat_clean.tsv", delim = "\t", col_names = TRUE, trim_ws = TRUE)

### Clean this
#ko_abundance <- "data/picrust2_output_0.1/KO_metagenome_out/pred_metagenome_unstrat_clean.tsv"
#ko_abundance <- read_delim(ko_abundance, delim = "\t", col_names = TRUE, trim_ws = TRUE)

#abundance_file <- "data/picrust2_output_0.1/KO_metagenome_out/pred_metagenome_unstrat_clean.tsv"
#abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)
###

# Perform pathway DAA using ALDEx2 method
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
daa_results_df <- pathway_daa(abundance = ko_abundance %>% column_to_rownames("#function"), metadata = metadata, group = "Time", daa_method = "ALDEx2", select = NULL, reference = NULL)

#########
# Filter results for ALDEx2_Kruskal-Wallace test method
daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "ALDEx2_Wilcoxon rank test", ]

# Annotate pathway results without KO to KEGG conversion
daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df, ko_to_kegg = FALSE)
#########


# Generate pathway error bar plot
# Please change column_to_rownames() to the feature column
# Please change Group to metadata$your_group_column if you are not using example dataset
p2 <- pathway_errorbar(abundance = ko_abundance %>% column_to_rownames("#function"), daa_results_df = daa_annotated_sub_method_results_df, Group = metadata$Time, p_values_threshold = 0.05, order = "group",
                      select = daa_annotated_sub_method_results_df %>% arrange(p_adjust) %>% slice(1:20) %>% dplyr::select(feature) %>% pull(),
                      ko_to_kegg = FALSE,
                      p_value_bar = TRUE,
                      colors = NULL,
                      x_lab = "description")


# Workflow for MetaCyc Pathway and EC

# Load MetaCyc pathway abundance and metadata
#data("metacyc_abundance")
#data("metadata")

# A) MetaCyc Pathway
# Load MetaCyc pathway
metacyc_abundance <- read_delim("data/picrust2_output_0.1/pathways_out/path_abun_unstrat.tsv", delim = "\t", col_names = TRUE, trim_ws = TRUE)

# Perform pathway DAA using LinDA method
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
metacyc_daa_results_df <- pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = Variable, daa_method = "LinDA")

# Annotate MetaCyc pathway results without KO to KEGG conversion
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", daa_results_df = metacyc_daa_results_df, ko_to_kegg = FALSE)

# Generate pathway error bar plot
# Please change column_to_rownames() to the feature column
# Please change Group to metadata$your_group_column if you are not using example dataset
p_metacyc <- pathway_errorbar(abundance = metacyc_abundance %>% column_to_rownames("pathway"), daa_results_df = metacyc_daa_annotated_results_df, Group = metadata[[Variable]], ko_to_kegg = FALSE, p_values_threshold = 0.05, order = "group", select = NULL, p_value_bar = TRUE, colors = NULL, x_lab = "description")

ggsave(glue("results/pathway_errorbar_MEtaCyc_{Variable}.png"), plot = p_metacyc, width = 8, height = 6, units = "in", dpi = 300)

# Get significantly different pathways
top_pathways <- metacyc_daa_annotated_results_df %>%
  filter(p_adjust < 0.05) %>%
  arrange(p_adjust) %>%
  slice(1:20) %>%  # select top 20 most significant
  pull(feature)

p_metacyc_20 <- pathway_errorbar(abundance = metacyc_abundance %>% column_to_rownames("pathway"), daa_results_df = metacyc_daa_annotated_results_df, Group = metadata[[Variable]], ko_to_kegg = FALSE, p_values_threshold = 0.05, order = "group", select = top_pathways, p_value_bar = TRUE, colors = NULL, x_lab = "description")

ggsave(glue("results/pathway_errorbar_MEtaCyc_{Variable}_20.png"), plot = p_metacyc_20, width = 8, height = 6, units = "in", dpi = 300)


# Generate pathway heatmap
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
feature_with_p_0.05 <- metacyc_daa_results_df %>% filter(p_adjust < 0.05)
pathway_heatmap(abundance = metacyc_abundance %>% filter(pathway %in% feature_with_p_0.05$feature) %>% column_to_rownames("pathway"), metadata = metadata, group = Variable)



# Generate pathway PCA plot
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
pathway_pca(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = Variable, colors=c("aquamarine3","plum3"))

# colors=c("green","purple"))
# colors=c("#66C2A5" "#FC8D62"))

# Run pathway DAA for multiple methods
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
methods <- c("ALDEx2", "DESeq2", "edgeR")
daa_results_list <- lapply(methods, function(method) {
  pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = Variable, daa_method = method)
})

# Compare results across different methods
comparison_results <- compare_daa_results(daa_results_list = daa_results_list, method_names = c("ALDEx2_Welch's t test", "ALDEx2_Wilcoxon rank test", "DESeq2", "edgeR"))
# Save comparisons
write.csv(comparison_results, file = glue("results/comparison_results_{Variable}.csv"), row.names = FALSE)

# B) EC
# Load EC
EC_abundance <- read_delim("data/picrust2_output_0.1/EC_metagenome_out/pred_metagenome_unstrat.tsv", delim = "\t", col_names = TRUE, trim_ws = TRUE)
EC_abundance <- EC_abundance %>% rename(EC_function = `function`)

# Perform pathway DAA using LinDA method
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
EC_daa_results_df <- pathway_daa(abundance = EC_abundance %>% column_to_rownames("EC_function"), metadata = metadata, group = Variable, daa_method = "LinDA")

# Annotate EC pathway results without KO to KEGG conversion
EC_daa_annotated_results_df <- pathway_annotation(pathway = "EC", daa_results_df = EC_daa_results_df, ko_to_kegg = FALSE)

## WARNING!! This is not working I get an error in pathway_errorbar function
# Generate pathway error bar plot
# Please change column_to_rownames() to the feature column
# Please change Group to metadata$your_group_column if you are not using example dataset
#pathway_errorbar(abundance = EC_abundance %>% column_to_rownames("EC_function"), daa_results_df = EC_daa_annotated_results_df, Group = metadata[[Variable]], ko_to_kegg = FALSE, p_values_threshold = 0.05, order = "group", select = NULL, p_value_bar = TRUE, colors = NULL, x_lab = "description")

# Generate pathway heatmap
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
feature_with_p_0.05 <- EC_daa_results_df %>% filter(p_adjust < 0.05)
p_path_heatmap_EC <- pathway_heatmap(abundance = EC_abundance %>% filter(EC_function %in% feature_with_p_0.05$feature) %>% column_to_rownames("EC_function"), metadata = metadata, group = Variable)
p_path_heatmap_EC
# Save heatmap
ggsave(glue("results/pathway_heatmap_EC_{Variable}.png"), plot = p_path_heatmap_EC, width = 8, height = 6, units = "in", dpi = 300)


# Generate pathway PCA plot
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
pathway_pca_EC <-pathway_pca(abundance = EC_abundance %>% column_to_rownames("EC_function"), metadata = metadata, group = Variable, colors=c("aquamarine3","plum3"))
pathway_pca_EC

# Run pathway DAA for multiple methods
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
methods <- c("ALDEx2", "DESeq2", "edgeR")
daa_results_list <- lapply(methods, function(method) {
  pathway_daa(abundance = EC_abundance %>% column_to_rownames("EC_function"), metadata = metadata, group = Variable, daa_method = method)
})

# Compare results across different methods
comparison_results_EC <- compare_daa_results(daa_results_list = daa_results_list, method_names = c("ALDEx2_Welch's t test", "ALDEx2_Wilcoxon rank test", "DESeq2", "edgeR"))

#######################

# Split metadata by Time column
metadata_early <- metadata %>% filter(Time == "early")
metadata_late  <- metadata %>% filter(Time == "late")

# Metadata variable: metadata (original); metadata_early or metadata_late
metadata <- "metadata_early"

# Extract sample names from the column sample_name early and late
samples_early <- metadata_early$sample_name
samples_late <- metadata_late$sample_name
# Select columns from ko_abundance: "#function" + those present in metadata_early and metadata_late
ko_abundance_early <- ko_abundance[, c("#function", samples_early)]
ko_abundance_late <- ko_abundance[, c("#function", samples_late)]

# Remove rows with all zero values across samples (excluding the "#function" column)
ko_abundance_early_filtered <- ko_abundance_early[rowSums(ko_abundance_early[,-1]) != 0, ]
ko_abundance_late_filtered <- ko_abundance_late[rowSums(ko_abundance_late[,-1]) != 0, ]

# Select columns kegg_abundance: "#function" + those present in metadata_early and metadata_late
samples_early_in_kegg <- intersect(samples_early, colnames(kegg_abundance))
kegg_abundance_early <- kegg_abundance[, samples_early_in_kegg, drop = FALSE]
samples_late_in_kegg <- intersect(samples_late, colnames(kegg_abundance))
kegg_abundance_late <- kegg_abundance[, samples_late_in_kegg, drop = FALSE]

# Remove rows with all zero values across samples (excluding the "#function" column)
kegg_abundance_early_filtered <- kegg_abundance_early[rowSums(kegg_abundance_early) != 0, ]
kegg_abundance_late_filtered  <- kegg_abundance_late[rowSums(kegg_abundance_late) != 0, ]

# Indicate the variable you want to compare
#Variable <- "Time" # Comparing by period early vs. late
Variable <- "category_2" # Comparing by sites and flower species richness

# Perform pathway differential abundance analysis (DAA) using ALDEx2 method
# Please change group to "your_group_column" if you are not using example dataset
#daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata, group = "Time", daa_method = "ALDEx2", select = NULL, reference = NULL)

daa_results_df <- pathway_daa(abundance = kegg_abundance_early_filtered, metadata = metadata_early, group = Variable, daa_method = "ALDEx2", select = NULL, reference = NULL)

# Filter results for ALDEx2_Welch's t test method
# Please check the unique(daa_results_df$method) and choose one
daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "ALDEx2_Wilcoxon rank test", ]

# Annotate pathway results using KO to KEGG conversion
daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df, ko_to_kegg = TRUE)

# Generate pathway error bar plot
# Please change Group to metadata$your_group_column if you are not using example dataset
#p <- pathway_errorbar(abundance = kegg_abundance, daa_results_df = daa_annotated_sub_method_results_df, Group = metadata$Time, p_values_threshold = 0.05, order = "pathway_class", select = NULL, ko_to_kegg = TRUE, p_value_bar = TRUE, colors = NULL, x_lab = "pathway_name")

p <- pathway_errorbar(abundance = kegg_abundance_early_filtered, daa_results_df = daa_annotated_sub_method_results_df, Group = metadata_early[[Variable]], p_values_threshold = 0.05, order = "pathway_class", select = NULL, ko_to_kegg = TRUE, p_value_bar = TRUE, colors = NULL, x_lab = "pathway_name")

ggsave(glue("results/pathway_errorbar_{Variable}.png"), plot = p, width = 8, height = 6, units = "in", dpi = 300)

# Get significantly different pathways
top_pathways <- daa_annotated_sub_method_results_df %>%
  filter(p_adjust < 0.05) %>%
  arrange(p_adjust) %>%
  slice(1:20) %>%  # select top 20 most significant
  pull(feature)

p_20 <- pathway_errorbar(abundance = kegg_abundance, daa_results_df = daa_annotated_sub_method_results_df, Group = metadata[[Variable]], p_values_threshold = 0.05, order = "pathway_class", select = top_pathways, ko_to_kegg = TRUE, p_value_bar = TRUE, colors = NULL, x_lab = "pathway_name")
ggsave(glue("results/pathway_errorbar_{Variable}_20.png"), plot = p_20, width = 8, height = 6, units = "in", dpi = 300)





#######################  think it's not working

# Compare high and low flower richness from category column. We remove more than half of the data (68 samples). 
metadata_HL 

# Extract sample names from the column sample_name
samples_HL <- metadata_HL$sample_name
# Select columns from ko_abundance: "#function" + those present in metadata_early and metadata_late
ko_abundance_early <- ko_abundance[, c("#function", samples_early)]

# Extract only the columns that match sample_name
kegg_abundance_HL <- kegg_abundance[, colnames(kegg_abundance) %in% samples_HL]
kegg_abundance_HL <- kegg_abundance_HL[rowSums(kegg_abundance_HL != 0) > 0, ]

# Indicate the variable you want to compare
Variable <- "category" 

# Perform pathway differential abundance analysis (DAA) using ALDEx2 method
# Please change group to "your_group_column" if you are not using example dataset
#daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata, group = "Time", daa_method = "ALDEx2", select = NULL, reference = NULL)

daa_results_df <- pathway_daa(abundance = kegg_abundance_HL, metadata = metadata_HL, group = Variable, daa_method = "ALDEx2", select = NULL, reference = NULL)

# Filter results for ALDEx2_Welch's t test method
# Please check the unique(daa_results_df$method) and choose one
daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "ALDEx2_Wilcoxon rank test", ]

# Annotate pathway results using KO to KEGG conversion
daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df, ko_to_kegg = TRUE)

### WARNING: I only recover 1 pathway



# Generate pathway error bar plot
# Please change Group to metadata$your_group_column if you are not using example dataset
#p <- pathway_errorbar(abundance = kegg_abundance, daa_results_df = daa_annotated_sub_method_results_df, Group = metadata$Time, p_values_threshold = 0.05, order = "pathway_class", select = NULL, ko_to_kegg = TRUE, p_value_bar = TRUE, colors = NULL, x_lab = "pathway_name")

p <- pathway_errorbar(abundance = kegg_abundance_HL, daa_results_df = daa_annotated_sub_method_results_df, Group = metadata_HL[[Variable]], p_values_threshold = 0.05, order = "pathway_class", select = NULL, ko_to_kegg = TRUE, p_value_bar = TRUE, colors = NULL, x_lab = "pathway_name")

ggsave(glue("results/pathway_errorbar_{Variable}.png"), plot = p, width = 8, height = 6, units = "in", dpi = 300)

# Get significantly different pathways
top_pathways <- daa_annotated_sub_method_results_df %>%
  filter(p_adjust < 0.05) %>%
  arrange(p_adjust) %>%
  slice(1:20) %>%  # select top 20 most significant
  pull(feature)

p_20 <- pathway_errorbar(abundance = kegg_abundance, daa_results_df = daa_annotated_sub_method_results_df, Group = metadata[[Variable]], p_values_threshold = 0.05, order = "pathway_class", select = top_pathways, ko_to_kegg = TRUE, p_value_bar = TRUE, colors = NULL, x_lab = "pathway_name")
ggsave(glue("results/pathway_errorbar_{Variable}_20.png"), plot = p_20, width = 8, height = 6, units = "in", dpi = 300)
##########


