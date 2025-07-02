library(readr)
library(ggpicrust2)
library(dplyr)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)

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

# Load KEGG pathway abundance
# data(kegg_abundance)
kegg_abundance <- ko2kegg_abundance("data/picrust2_output_0.1/KO_metagenome_out/pred_metagenome_unstrat_clean.tsv")

# Perform pathway differential abundance analysis (DAA) using ALDEx2 method
# Please change group to "your_group_column" if you are not using example dataset
daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata, group = "Time", daa_method = "ALDEx2", select = NULL, reference = NULL)

# Filter results for ALDEx2_Welch's t test method
# Please check the unique(daa_results_df$method) and choose one
daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "ALDEx2_Wilcoxon rank test", ]

# Annotate pathway results using KO to KEGG conversion
daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df, ko_to_kegg = TRUE)

# Generate pathway error bar plot
# Please change Group to metadata$your_group_column if you are not using example dataset
p <- pathway_errorbar(abundance = kegg_abundance, daa_results_df = daa_annotated_sub_method_results_df, Group = metadata$Time, p_values_threshold = 0.05, order = "pathway_class", select = NULL, ko_to_kegg = TRUE, p_value_bar = TRUE, colors = NULL, x_lab = "pathway_name")

# If you want to analyze EC, MetaCyc, and KO without conversions, turn ko_to_kegg to FALSE.

### This was load

# Load metadata as a tibble
# data(metadata)
#metadata <- read_delim("path/to/your/metadata.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

###

# Load KO abundance as a data.frame
# data(ko_abundance)
ko_abundance <- read_delim("data/picrust2_output_0.1/KO_metagenome_out/pred_metagenome_unstrat_clean.tsv", delim = "\t", col_names = TRUE, trim_ws = TRUE)

ko_abundance <- "data/pred_metagenome_unstrat_clean.tsv"
ko_abundance <- read_delim(ko_abundance, delim = "\t", col_names = TRUE, trim_ws = TRUE)

abundance_file <- "data/pred_metagenome_unstrat_clean.tsv"
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)

# Perform pathway DAA using ALDEx2 method
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
daa_results_df <- pathway_daa(abundance = ko_abundance %>% column_to_rownames("#function"), metadata = metadata, group = "Time", daa_method = "ALDEx2", select = NULL, reference = NULL)

# Filter results for ALDEx2_Kruskal-Wallace test method
daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "ALDEx2_Wilcoxon rank test", ]

# Annotate pathway results without KO to KEGG conversion
daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df, ko_to_kegg = FALSE)

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
metacyc_daa_results_df <- pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Time", daa_method = "LinDA")

# Annotate MetaCyc pathway results without KO to KEGG conversion
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", daa_results_df = metacyc_daa_results_df, ko_to_kegg = FALSE)

# Generate pathway error bar plot
# Please change column_to_rownames() to the feature column
# Please change Group to metadata$your_group_column if you are not using example dataset
pathway_errorbar(abundance = metacyc_abundance %>% column_to_rownames("pathway"), daa_results_df = metacyc_daa_annotated_results_df, Group = metadata$Time, ko_to_kegg = FALSE, p_values_threshold = 0.05, order = "group", select = NULL, p_value_bar = TRUE, colors = NULL, x_lab = "description")

# Generate pathway heatmap
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
feature_with_p_0.05 <- metacyc_daa_results_df %>% filter(p_adjust < 0.05)
pathway_heatmap(abundance = metacyc_abundance %>% filter(pathway %in% feature_with_p_0.05$feature) %>% column_to_rownames("pathway"), metadata = metadata, group = "Time")

# Generate pathway PCA plot
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
pathway_pca(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Time", colors=c("aquamarine3","plum3"))

# colors=c("green","purple"))
# colors=c("#66C2A5" "#FC8D62"))

# Run pathway DAA for multiple methods
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
methods <- c("ALDEx2", "DESeq2", "edgeR")
daa_results_list <- lapply(methods, function(method) {
  pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Time", daa_method = method)
})

# Compare results across different methods
comparison_results <- compare_daa_results(daa_results_list = daa_results_list, method_names = c("ALDEx2_Welch's t test", "ALDEx2_Wilcoxon rank test", "DESeq2", "edgeR"))



# B) EC
# Load EC
EC_abundance <- read_delim("data/picrust2_output_0.1/EC_metagenome_out/EC_pred_metagenome_unstrat.tsv", delim = "\t", col_names = TRUE, trim_ws = TRUE)

# Perform pathway DAA using LinDA method
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
EC_daa_results_df <- pathway_daa(abundance = EC_abundance %>% column_to_rownames("EC_function"), metadata = metadata, group = "Time", daa_method = "LinDA")

# Annotate EC pathway results without KO to KEGG conversion
EC_daa_annotated_results_df <- pathway_annotation(pathway = "EC", daa_results_df = EC_daa_results_df, ko_to_kegg = FALSE)

# Generate pathway error bar plot
# Please change column_to_rownames() to the feature column
# Please change Group to metadata$your_group_column if you are not using example dataset
pathway_errorbar(abundance = EC_abundance %>% column_to_rownames("EC_function"), daa_results_df = EC_daa_annotated_results_df, Group = metadata$Time, ko_to_kegg = FALSE, p_values_threshold = 0.05, order = "group", select = NULL, p_value_bar = TRUE, colors = NULL, x_lab = "description")

# Generate pathway heatmap
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
feature_with_p_0.05 <- EC_daa_results_df %>% filter(p_adjust < 0.05)
pathway_heatmap(abundance = EC_abundance %>% filter(EC_function %in% feature_with_p_0.05$feature) %>% column_to_rownames("EC_function"), metadata = metadata, group = "Time")

pathway_heatmap(abundance = metacyc_abundance %>% filter(pathway %in% feature_with_p_0.05$feature) %>% column_to_rownames("pathway"), metadata = metadata, group = "Time")


# Generate pathway PCA plot
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
pathway_pca(abundance = EC_abundance %>% column_to_rownames("EC_function"), metadata = metadata, group = "Time", colors=c("aquamarine3","plum3"))

# colors=c("green","purple"))
# colors=c("#66C2A5" "#FC8D62"))

# Run pathway DAA for multiple methods
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
methods <- c("ALDEx2", "DESeq2", "edgeR")
daa_results_list <- lapply(methods, function(method) {
  pathway_daa(abundance = EC_abundance %>% column_to_rownames("EC_function"), metadata = metadata, group = "Time", daa_method = method)
})

# Compare results across different methods
comparison_results <- compare_daa_results(daa_results_list = daa_results_list, method_names = c("ALDEx2_Welch's t test", "ALDEx2_Wilcoxon rank test", "DESeq2", "edgeR"))
