###
###     Title: Bombus pascuorum microbiome
###     Authors: Chueca Luis J., Poza Jon, Dhami Manpreet P., Donald Marion L., 
###     Rose Jennifer, Salgado-Irazabal Xabier, Hermosilla Brais,
###     Gostout Christian & Magrach Ainhoa.
###     Year: 2025?
###     Related project: GorBEEa
###     Description: Metabarcoding analysis in R
###


# Load libraries
library(decontam); packageVersion("decontam")
library(phyloseq) ; packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(Biostrings); packageVersion("Biostrings")
library(gridExtra); packageVersion("gridExtra")

# Load data
count_tab <- read.table("data/dada2/2023_16S_GorBEEa_prj_ASVs_counts.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")

tax_tab <- as.matrix(read.table("data/dada2/2023_16S_GorBEEa_prj_ASVs_taxonomy.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))

sample_info_tab <- read.delim("data/dada2/sample_info.tsv", header=T, row.names=1, check.names=F, sep="\t")

fasta_seqs <- readDNAStringSet("data/dada2/2023_16S_GorBEEa_prj_ASVs.fa")

# Setting the color column to be of type "character", which helps later
sample_info_tab$color_p <- as.character(sample_info_tab$color_p)
sample_info_tab$color_s <- as.character(sample_info_tab$color_s)

# Create a phyloseq object
otu_table_obj <- otu_table(count_tab, taxa_are_rows = TRUE)
tax_table_obj <- tax_table(as.matrix(tax_tab))
sample_data_obj <- sample_data(sample_info_tab)

# Define the new column names
new_column_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# Assign the new column names to the tax_table
colnames(tax_table_obj) <- new_column_names
# Verify the change
colnames(tax_table_obj)

physeq <- phyloseq(otu_table_obj, tax_table_obj, sample_data_obj)
print(physeq)
head(sample_data(physeq))

#      1. Check for contaminants      ####

## 1.a) Inspect Library Sizes     ####

# Extract sample data into a data frame
df <- as.data.frame(sample_data(physeq))

# Add LibrarySize and Index columns
df$LibrarySize <- sample_sums(physeq)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

# Create a plot for each value in the plate column
plot_list <- lapply(unique(df$plate), function(plate_value) {
  df_subset <- df[df$plate == plate_value, ]
  max_negative <- max(df_subset$quant_reading[df_subset$type == "negative"])
  ggplot(data=df_subset, aes(x=Index, y=LibrarySize, color=type)) + 
    geom_point() +
    geom_hline(yintercept = max_negative, linetype = "dashed", color = "black") +
    ggtitle(paste("Plate:", plate_value)) +
    theme(legend.position = "none") # Remove individual legends
})

# Extract the legend from one of the plots
legend <- get_legend(
  ggplot(data=df[df$plate == unique(df$plate)[1], ], aes(x=Index, y=LibrarySize, color=type)) + 
    geom_point() +
    geom_hline(yintercept = max(df$quant_reading[df$type == "negative"]), linetype = "dashed", color = "black") +
    theme(legend.position = "right")
)

# Combine the plots and the legend using cowplot
combined_plot <- plot_grid(plot_list[[1]], plot_list[[2]], legend, ncol = 3, rel_widths = c(1, 1, 0.3))
combined_plot

ggsave("results/library_size_combined_plot.pdf", plot = combined_plot, device = "pdf", width = 8, height = 6, units = "in")


## 1.b) Identify Contaminants - Prevalence (Separated by Plate) ####

# Create an empty list to store results for each plate
contam_results <- list()

# Define threshold values to test
threshold_values <- c(0.01, 0.05, 0.1, 0.5)

# Loop through each plate separately
for (plate_id in unique(sample_data(physeq)$plate)) {
  
  # Subset phyloseq object for current plate
  physeq_plate <- prune_samples(sample_data(physeq)$plate == plate_id, physeq)
  
  # Identify negative controls only for this plate
  sample_data(physeq_plate)$is.neg <- sample_data(physeq_plate)$type == "negative"
  
  # Function to test multiple thresholds per plate
  test_thresholds <- function(physeq, thresholds, neg_column) {
    results <- list()
    for (threshold in thresholds) {
      contamdf <- isContaminant(physeq, method = "prevalence", neg = neg_column, threshold = threshold)
      results[[as.character(threshold)]] <- list(
        threshold = threshold,
        contaminant_table = table(contamdf$contaminant),
        contaminant_ids = which(contamdf$contaminant),
        dataframe = contamdf
      )
    }
    return(results)
  }
  
  # Run the function for this plate
  results <- test_thresholds(physeq_plate, thresholds = threshold_values, neg_column = "is.neg")
  
  # Store results
  contam_results[[plate_id]] <- results
  
  # Print summarized results for each plate
  cat("\n=== Plate:", plate_id, "===\n")
  for (threshold in names(results)) {
    cat("\nThreshold:", threshold, "\n")
    print(results[[threshold]]$contaminant_table)
  }
}

# Pick a threshold to filter contaminants
Th <- 0.5  # Adjust as needed

# Merge contamination results across plates
contam_combined <- rep(FALSE, ntaxa(physeq))
for (plate_id in names(contam_results)) {
  contam_combined <- contam_combined | contam_results[[plate_id]][[as.character(Th)]]$dataframe$contaminant
}

# Create a dataframe with contamination status
contamdf.prev <- data.frame(contaminant = contam_combined)
rownames(contamdf.prev) <- taxa_names(physeq)

# Check contamination summary
table(contamdf.prev$contaminant)

## 1.c) Remove Contaminants ####

# Remove contaminants from the phyloseq object
ps.nocont <- prune_taxa(!contamdf.prev$contaminant, physeq)

# Create a phyloseq object with only contaminant ASVs
ps.cont <- prune_taxa(contamdf.prev$contaminant, physeq)

## 1.d) Remove contaminants from fasta file     ####
contam_asvs <- taxa_names(ps.cont)  # Get the ASV IDs of contaminants from ps.cont
# Extract the ASV IDs from the headers of the FASTA sequences
fasta_ids <- sapply(names(fasta_seqs), function(x) strsplit(x, " ")[[1]][1])
# Filter out the contaminant ASVs
remaining_ids <- fasta_ids[!fasta_ids %in% contam_asvs]  # IDs to keep
filtered_fasta <- fasta_seqs[which(fasta_ids %in% remaining_ids)]  # Keep only non-contaminant ASVs
# Check the result
length(filtered_fasta)  # Number of remaining sequences after filtering
# Write the filtered sequences to a new FASTA file
fasta_filename <- paste0("data/2023_16S_GorBEEa_prj_ASVs_filt_c.", Th, ".fa")
writeXStringSet(filtered_fasta, file = fasta_filename)

## 1.e) Remove Negative Controls from Final Dataset ####

# Keep only true samples
ps.gbp23 <- prune_samples(sample_data(ps.nocont)$type == "sample", ps.nocont)

#Extract ASV count table and save it
ASV_table <- otu_table(ps.gbp23)
ASV_table <- as.data.frame(ASV_table)
tsv_filename <- paste0("ASV_table_", Th, ".tsv")
write.table(ASV_table, file = tsv_filename, sep = "\t", row.names = FALSE, quote = FALSE)

# Check the number of samples before and after filtering
table(sample_data(ps.nocont)$type)
table(sample_data(ps.gbp23)$type)

# Save the final cleaned phyloseq object
rds_filename <- paste0("data/ps.gbp23_c.", Th, ".RDS")
saveRDS(ps.gbp23, file = rds_filename)