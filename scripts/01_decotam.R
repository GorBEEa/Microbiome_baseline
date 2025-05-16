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

df <- as.data.frame(sample_data(physeq)) # Put sample_data into a ggplot-friendly data.frame
max_negative <- max(df$quant_reading[df$type == "negative"]) # Calculate the highest quant_reading for "negative" samples
df$LibrarySize <- sample_sums(physeq)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
plot <- ggplot(data=df, aes(x=Index, y=LibrarySize, color=type)) + 
  geom_point() +
  geom_hline(yintercept = max_negative, linetype = "dashed", color = "black")
ggsave("results/library_size_plot.pdf", plot = plot, device = "pdf", width = 8, height = 6, units = "in")
plot

###########################################################################################################################
###     
###                         WARNING


# Identify Contaminants - Frequency 
# LJ --> I think is not doing correctly. not with read frequency, is has to be done with DNA concentration

#contamdf.freq <- isContaminant(physeq, method="frequency", conc="quant_reading")
#head(contamdf.freq)
#table(contamdf.freq$contaminant)
#head(which(contamdf.freq$contaminant))

#plot_frequency(physeq, taxa_names(physeq)[c(1,6)], conc="quant_reading") + 
#xlab("DNA Concentration (PicoGreen fluorescent intensity)")
##########################################################################################################################


## 1.b) Identify Contaminants - Prevalence      ####

# Check different threshold values

# set negative or control samples
sample_data(physeq)$is.neg <- sample_data(physeq)$type == "negative"
# Function to test multiple thresholds
test_thresholds <- function(physeq, thresholds, neg_column) {
  # List to store the results
  results <- list()
  
  # Iterate over each threshold value
  for (threshold in thresholds) {
    # Apply isContaminant with the current threshold
    contamdf <- isContaminant(physeq, method = "prevalence", neg = neg_column, threshold = threshold)
    
    # Save the results in the list
    results[[as.character(threshold)]] <- list(
      threshold = threshold,
      contaminant_table = table(contamdf$contaminant),
      contaminant_ids = which(contamdf$contaminant),
      dataframe = contamdf
    )
  }
  
  # Return all results
  return(results)
}

# Threshold values to test
threshold_values <- c(0.01, 0.05, 0.1, 0.5)

# Use the function
results <- test_thresholds(physeq, thresholds = threshold_values, neg_column = "is.neg")

# Explore results:
# Ptint all summarized results:
for (threshold in names(results)) {
  cat("\nThreshold:", threshold, "\n")
  print(results[[threshold]]$contaminant_table)
}

# Check each threshold individually:
Th <- 0.5 # Replace by your value
results[[as.character(Th)]]$contaminant_table
results[[as.character(Th)]]$contaminant_ids
head(results[[as.character(Th)]]$dataframe, n = 20)

# Run with your selected threshold value
contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev$contaminant)
head(contamdf.prev, n=20)
head(which(contamdf.prev$contaminant))

# Make phyloseq object of presence-absence in negative controls and true samples
physeq.pa <- transform_sample_counts(physeq, function(abund) 1*(abund>0))
physeq.pa.neg <- prune_samples(sample_data(physeq.pa)$type == "negative", physeq.pa)
physeq.pa.pos <- prune_samples(sample_data(physeq.pa)$type == "sample", physeq.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(physeq.pa.pos), pa.neg=taxa_sums(physeq.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")


## 1.c) Remove contaminants     ####

# create phyloseq object with contaminant ASVs removed  
ps.nocont <- prune_taxa(!contamdf.prev$contaminant, physeq)
# create a phyloseq object with only contaminant ASVs
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
writeXStringSet(filtered_fasta, file = "data/2023_16S_GorBEEa_prj_ASVs_filt_0.5.fa")

## 1.e) Remove negative controls from phyloseq object     ####

# Filter out "negative" samples based on the "type" variable in sample_data
ps.gbp23 <- prune_samples(sample_data(ps.nocont)$type == "sample", ps.nocont)
# Verify the number of samples and variables in the new object
ps.gbp23
# Check the original sample types
table(sample_data(ps.nocont)$type)
# Check the filtered sample types
table(sample_data(ps.gbp23)$type)

# Save the ps.gbp23 object to an RDS file
saveRDS(ps.gbp23, file = "data/ps.gbp23.0.5.RDS")

