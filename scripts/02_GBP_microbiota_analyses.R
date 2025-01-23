###################################################################
###
###     Title: Bombus pascuorum microbiome
###     Authors: Chueca Luis J., Poza Jon, Dhami Manpreet P., Donald Marion L., 
###     Rose Jennifer, Salgado-Irazabal Xabier, Hermosilla Brais,
###     Gostout Christian & Magrach Ainhoa.
###     Year: 2025?
###     Related project: GorBEEa
###     Description: Metabarcoding analysis in R
###
##############################################################


# Load libraries
library(Biostrings) ; packageVersion("Biostrings")
library(tidyverse) ; packageVersion("tidyverse")
library(phyloseq) ; packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(ggstatsplot) ; packageVersion("ggstatsplot")
library(vegan) ; packageVersion("vegan") 
library(DESeq2) ; packageVersion("DESeq2") 
library(dendextend) ; packageVersion("dendextend") 
library(viridis) ; packageVersion("viridis")
library(patchwork) ; packageVersion("patchwork")
library (dplyr) ; packageVersion("dplyr")

# Load the phyloseq object
ps.gbp23 <- readRDS("data/ps.gbp23.RDS")
print(ps.gbp23)

####################################################
###                                              ###
###       2. Comparison between periods          ###
###                                              ###
####################################################

# From the phyloseq contaminant-free file
# Extract the ASV count table, sample_data and tax_table
count_tab.cl <- otu_table(ps.gbp23)
count_tab.cl <- as.data.frame(count_tab.cl)
sample_data_tab.cl <- sample_data(ps.gbp23)
tax_table.cl <- tax_table(ps.gbp23)

deseq_counts <- DESeqDataSetFromMatrix(count_tab.cl, colData = sample_data_tab.cl, design = ~period) 
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)
euc_dist <- dist(t(vst_trans_count_tab))

# Hierarchical clustering for 
euc_clust <- hclust(euc_dist, method="ward.D2")

# hclust objects like this can be plotted with the generic plot() function
plot(euc_clust) 
# but i like to change them to dendrograms for two reasons:
# 1) it's easier to color the dendrogram plot by groups
# 2) if wanted you can rotate clusters with the rotate() 
#    function of the dendextend package

euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- as.character(sample_data_tab.cl$color_p[order.dendrogram(euc_dend)])
labels_colors(euc_dend) <- dend_cols

plot(euc_dend, ylab="VST Euc. dist.")


# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_data_tab.cl)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

# generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

plot_ordination(vst_physeq, vst_pcoa, color="period") + 
  geom_point(size=1) + labs(col="period") + 
  geom_text(aes(label=rownames(sample_data_tab.cl), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values=unique(sample_data_tab.cl$color_p[order(sample_data_tab.cl$period)])) + 
  theme_bw() + theme(legend.position="none")

####################################################
###                                              ###
###       3. Comparison between sites            ###
###                                              ###
####################################################

deseq_counts <- DESeqDataSetFromMatrix(count_tab.cl, colData = sample_data_tab.cl, design = ~site) 
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)
euc_dist <- dist(t(vst_trans_count_tab))

# Hierarchical clustering for 
euc_clust <- hclust(euc_dist, method="ward.D2")

# hclust objects like this can be plotted with the generic plot() function
plot(euc_clust) 
# but i like to change them to dendrograms for two reasons:
# 1) it's easier to color the dendrogram plot by groups
# 2) if wanted you can rotate clusters with the rotate() 
#    function of the dendextend package

euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- as.character(sample_data_tab.cl$color_s[order.dendrogram(euc_dend)])
labels_colors(euc_dend) <- dend_cols

plot(euc_dend, ylab="VST Euc. dist.")


# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_data_tab.cl)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

# generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

plot_ordination(vst_physeq, vst_pcoa, color="site") + 
  geom_point(size=1) + labs(col="site") + 
  geom_text(aes(label=rownames(sample_data_tab.cl), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values=unique(sample_data_tab.cl$color_s[order(sample_data_tab.cl$site)])) + 
  theme_bw() + theme(legend.position="none")


####################################################
###                                              ###
###       4. Alpha diversity                     ###
###                                              ###
####################################################

# Rarefaction curves 
rarecurve(t(count_tab.cl), step=100, col=sample_data_tab.cl$color_p, lwd=2, ylab="ASVs", label=F)
# and adding a vertical line at the fewest seqs in any sample
abline(v=(min(rowSums(t(count_tab)))))


###     4.a) Richness and diversity estimates
###

################################## Perhaps not neccessary
# first we need to create a phyloseq object using our un-transformed count table
count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)

ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)




ps.gbp23
sample_data_tab.cl
##################################


plot_richness(ps.gbp23, color="period", measures=c("Chao1", "Shannon")) + 
  scale_color_manual(values=unique(sample_data_tab.cl$color_p[order(sample_data_tab.cl$period)])) +
  theme_bw() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

period_richness <- plot_richness(ps.gbp23, x="period", color="period", measures=c("Chao1", "Shannon")) + 
  scale_color_manual(values=unique(sample_data_tab.cl$color_p[order(sample_data_tab.cl$period)])) +
  theme_bw() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_violin() + geom_jitter(height = 0, width = 0.1)

period_richness + scale_fill_manual(values=c("#440154FF","#414487FF","#2A788EFF","#22A884FF","#7AD151FF","#FDE725FF","red"))

site_richness <- plot_richness(ps.gbp23, x="site", color="site", measures=c("Chao1", "Shannon")) + 
  scale_color_manual(values=unique(sample_data_tab.cl$color_s[order(sample_data_tab.cl$site)])) +
  theme_bw() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

site_richness + geom_violin() + geom_jitter(height = 0, width = 0.1)


# Statistical comparisons

# Estimate richness (Chao1 and Shannon index) by using phyloseq
richness_data <- estimate_richness(ps.gbp23, measures = c("Chao1", "Shannon"))
# Add period and site metadata
richness_data$period <- sample_data(ps.gbp23)$period
richness_data$site <- sample_data(ps.gbp23)$site

# Shapiro-Wilk normality test
shapiro.test(richness_data$Shannon)
shapiro.test(richness_data$Chao1)

# Define our custom color palette
viridis_palette <- c("#440154FF", "#414487FF", "#2A788EFF", "#22A884FF", "#7AD151FF", "#FDE725FF")


plt1 <- ggbetweenstats(data = richness_data, x = period, y = Chao1, type = "nonparametric"
) +
  scale_color_manual(values = viridis_palette) +
  scale_fill_manual(values = viridis_palette)

plt2 <-ggbetweenstats(data = richness_data, x = period, y = Shannon, type = "nonparametric"
) +
  scale_color_manual(values = viridis_palette) +
  scale_fill_manual(values = viridis_palette)

plt3 <- ggbetweenstats(data = richness_data, x = site, y = Chao1, type = "nonparametric",
                       package = "colorBlindness",
                       palette = "Blue2DarkOrange18Steps",
)

plt4 <- ggbetweenstats(data = richness_data, x = site, y = Shannon, type = "nonparametric",
                       package = "colorBlindness",
                       palette = "Blue2DarkOrange18Steps",
)

combined_plot <- (plt1 + plt2) / (plt3 + plt4)
combined_plot

# Not significant difference between periods and sites.

####
####      4.b) TAXONOMIC SUMMARIES
####


ps.gbp23
sample_data_tab.cl
tax_table.cl
count_tab.cl

# using phyloseq to make a count table that has summed all ASVs
# that were in the same phylum
phyla_counts_tab <- otu_table(tax_glom(ps.gbp23, taxrank="Phylum")) 

# making a vector of phyla names to set as row names
phyla_tax_vec <- as.vector(tax_table(tax_glom(ps.gbp23, taxrank="Phylum"))[,"Phylum"]) 
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

# we also have to account for sequences that weren't assigned any
# taxonomy even at the phylum level 
# these came into R as 'NAs' in the taxonomy table, but their counts are
# still in the count table
# so we can get that value for each sample by subtracting the column sums
# of this new table (that has everything that had a phylum assigned to it)
# from the column sums of the starting count table (that has all
# representative sequences)
unclassified_tax_counts <- colSums(count_tab.cl) - colSums(phyla_counts_tab)
# and we'll add this row to our phylum count table:
phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unclassified"=unclassified_tax_counts)

# now we'll remove the Proteobacteria, so we can next add them back in
# broken down by class
temp_major_taxa_counts_tab <- phyla_and_unidentified_counts_tab[!row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ]

# making count table broken down by class (contains classes beyond the
# Proteobacteria too at this point)
class_counts_tab <- otu_table(tax_glom(ps.gbp23, taxrank="Class")) 

# making a table that holds the phylum and class level info
class_tax_phy_tab <- tax_table(tax_glom(ps.gbp23, taxrank="Class")) 

phy_tmp_vec <- class_tax_phy_tab[,2]
class_tmp_vec <- class_tax_phy_tab[,3]
rows_tmp <- row.names(class_tax_phy_tab)
class_tax_tab <- data.frame("Phylum"=phy_tmp_vec, "Class"=class_tmp_vec, row.names = rows_tmp)

# making a vector of just the Proteobacteria classes
proteo_classes_vec <- as.vector(class_tax_tab[class_tax_tab$Phylum == "Proteobacteria", "Class"])

# changing the row names like above so that they correspond to the taxonomy,
# rather than an ASV identifier
rownames(class_counts_tab) <- as.vector(class_tax_tab$Class) 

# making a table of the counts of the Proteobacterial classes
proteo_class_counts_tab <- class_counts_tab[row.names(class_counts_tab) %in% proteo_classes_vec, ] 

# there are also possibly some some sequences that were resolved to the level
# of Proteobacteria, but not any further, and therefore would be missing from
# our class table
# we can find the sum of them by subtracting the proteo class count table
# from just the Proteobacteria row from the original phylum-level count table
proteo_no_class_annotated_counts <- phyla_and_unidentified_counts_tab[row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ] - colSums(proteo_class_counts_tab)

# now combining the tables:
major_taxa_counts_tab <- rbind(temp_major_taxa_counts_tab, proteo_class_counts_tab, "Unresolved_Proteobacteria"=proteo_no_class_annotated_counts)

# and to check we didn't miss any other sequences, we can compare the column
# sums to see if they are the same
# if "TRUE", we know nothing fell through the cracks
identical(colSums(major_taxa_counts_tab), colSums(count_tab)) 

# now we'll generate a proportions table for summarizing:
major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)

# if we check the dimensions of this table at this point
dim(major_taxa_proportions_tab)
# we see there are currently 42 rows, which might be a little busy for a
# summary figure
# many of these taxa make up a very small percentage, so we're going to
# filter some out
# this is a completely arbitrary decision solely to ease visualization and
# intepretation, entirely up to your data and you
# here, we'll only keep rows (taxa) that make up greater than 5% in any
# individual sample
temp_filt_major_taxa_proportions_tab <- data.frame(major_taxa_proportions_tab[apply(major_taxa_proportions_tab, 1, max) > 5, ])
# checking how many we have that were above this threshold
dim(temp_filt_major_taxa_proportions_tab) 
# now we have 12, much more manageable for an overview figure

# though each of the filtered taxa made up less than 5% alone, together they
# may add up and should still be included in the overall summary
# so we're going to add a row called "Other" that keeps track of how much we
# filtered out (which will also keep our totals at 100%)
filtered_proportions <- colSums(major_taxa_proportions_tab) - colSums(temp_filt_major_taxa_proportions_tab)
filt_major_taxa_proportions_tab <- rbind(temp_filt_major_taxa_proportions_tab, "Other"=filtered_proportions)

## don't worry if the numbers or taxonomy vary a little, this might happen due to different versions being used 
## from when this was initially put together







# first let's make a copy of our table that's safe for manipulating
filt_major_taxa_proportions_tab_for_plot <- filt_major_taxa_proportions_tab

# and add a column of the taxa names so that it is within the table, rather
# than just as row names (this makes working with ggplot easier)
filt_major_taxa_proportions_tab_for_plot$Major_Taxa <- row.names(filt_major_taxa_proportions_tab_for_plot)

# now we'll transform the table into narrow, or long, format (also makes
# plotting easier)
filt_major_taxa_proportions_tab_for_plot.g <- pivot_longer(filt_major_taxa_proportions_tab_for_plot, !Major_Taxa, names_to = "Sample", values_to = "Proportion") %>% data.frame()

# take a look at the new table and compare it with the old one
head(filt_major_taxa_proportions_tab_for_plot.g)
head(filt_major_taxa_proportions_tab_for_plot)
# manipulating tables like this is something you may need to do frequently in R

# now we want a table with "color" and "characteristics" of each sample to
# merge into our plotting table so we can use that more easily in our plotting
# function
# here we're making a new table by pulling what we want from the sample
# information table
sample_info_for_merge<-data.frame("Sample"=row.names(sample_info_tab), "period"=sample_info_tab$period, "color"=sample_info_tab$color_p, stringsAsFactors=F)

# and here we are merging this table with the plotting table we just made
# (this is an awesome function!)
filt_major_taxa_proportions_tab_for_plot.g2 <- merge(filt_major_taxa_proportions_tab_for_plot.g, sample_info_for_merge)

# and now we're ready to make some summary figures with our wonderfully
# constructed table

## a good color scheme can be hard to find, i included the viridis package
## here because it's color-blind friendly and sometimes it's been really
## helpful for me, though this is not demonstrated in all of the following :/ 

# one common way to look at this is with stacked bar charts for each taxon per sample:
ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="All samples")

ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(Major_Taxa, Proportion)) +
  geom_jitter(aes(color=factor(period), shape=factor(period)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_tab_for_plot.g2$color[order(filt_major_taxa_proportions_tab_for_plot.g2$period)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.title=element_blank()) +
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="All samples")


# Function to generate a ggplot for a given period
generate_major_taxa_plot <- function(period, data, sample_info_tab) {
  # Get sample IDs for the specified period
  sample_IDs <- row.names(sample_info_tab)[sample_info_tab$period == period]
  
  # Filter the data for the given period
  filt_data <- data[data$Sample %in% sample_IDs, ]
  
  # Generate the plot
  plot <- ggplot(filt_data, aes(Major_Taxa, Proportion)) +
    scale_y_continuous(limits = c(0, 100)) +
    geom_jitter(aes(color = factor(period)), size = 2, width = 0.15, height = 0) +
    scale_color_manual(values = unique(filt_data$color[order(filt_data$period)])) +
    geom_boxplot(fill = NA, outlier.color = NA) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(
      x = "Major Taxa", 
      y = "% of 16S rRNA gene copies recovered", 
      title = paste("Period", gsub("p", "", period), "samples")
    )
  
  return(plot)
}

# List of periods to analyze
periods <- c("p01", "p02", "p03", "p04", "p05", "p06")

# Store all generated plots in a list
plots <- list()
for (period in periods) {
  plots[[period]] <- generate_major_taxa_plot(period, filt_major_taxa_proportions_tab_for_plot.g2, sample_info_tab)
}


# Combine all plots into a grid with patchwork package
major_taxa_combined <- (plots[["p01"]] + plots[["p02"]] + plots[["p03"]]) / 
  (plots[["p04"]] + plots[["p05"]] + plots[["p06"]])
major_taxa_combined

# Save the plots. !LJC Should I include this?
# Save the combined plot
#ggsave("major_taxa_combined_plot.png", major_taxa_combined, width = 10, height = 8)

# Save individual plots
#for (period in periods) {
#  ggsave(paste0("major_taxa_", period, ".png"), plots[[period]], width = 5, height = 4)
#}


# LCJ # This correspond to the previous previous loop. Perhaps better remove
###########################################################

# let's set some helpful variables first:
p01_sample_IDs <- row.names(sample_info_tab)[sample_info_tab$period == "p01"]
p02_sample_IDs <- row.names(sample_info_tab)[sample_info_tab$period == "p02"]
p03_sample_IDs <- row.names(sample_info_tab)[sample_info_tab$period == "p03"]
p04_sample_IDs <- row.names(sample_info_tab)[sample_info_tab$period == "p04"]
p05_sample_IDs <- row.names(sample_info_tab)[sample_info_tab$period == "p05"]
p06_sample_IDs <- row.names(sample_info_tab)[sample_info_tab$period == "p06"]



# first we need to subset our plotting table to include just the rock samples to plot
filt_major_taxa_proportions_p02_only_tab_for_plot.g <- filt_major_taxa_proportions_tab_for_plot.g2[filt_major_taxa_proportions_tab_for_plot.g2$Sample %in% p02_sample_IDs, ]
# and then just the water samples
filt_major_taxa_proportions_p01_only_tab_for_plot.g <- filt_major_taxa_proportions_tab_for_plot.g2[filt_major_taxa_proportions_tab_for_plot.g2$Sample %in% p01_sample_IDs, ]

filt_major_taxa_proportions_p03_only_tab_for_plot.g <- filt_major_taxa_proportions_tab_for_plot.g2[filt_major_taxa_proportions_tab_for_plot.g2$Sample %in% p03_sample_IDs, ]

filt_major_taxa_proportions_p04_only_tab_for_plot.g <- filt_major_taxa_proportions_tab_for_plot.g2[filt_major_taxa_proportions_tab_for_plot.g2$Sample %in% p04_sample_IDs, ]

filt_major_taxa_proportions_p05_only_tab_for_plot.g <- filt_major_taxa_proportions_tab_for_plot.g2[filt_major_taxa_proportions_tab_for_plot.g2$Sample %in% p05_sample_IDs, ]

filt_major_taxa_proportions_p06_only_tab_for_plot.g <- filt_major_taxa_proportions_tab_for_plot.g2[filt_major_taxa_proportions_tab_for_plot.g2$Sample %in% p06_sample_IDs, ]


# and now we can use the same code as above just with whatever minor alterations we want

# Period p01 samples
major_taxa_p01 <- ggplot(filt_major_taxa_proportions_p01_only_tab_for_plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,100)) + # adding a setting for the y axis range so the rock and water plots are on the same scale
  geom_jitter(aes(color=factor(period)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_p01_only_tab_for_plot.g$color[order(filt_major_taxa_proportions_p01_only_tab_for_plot.g$period)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none") +
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="Period 01 samples")

# Period p02 samples
major_taxa_p02 <- ggplot(filt_major_taxa_proportions_p02_only_tab_for_plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,100)) + # adding a setting for the y axis range so the rock and water plots are on the same scale
  geom_jitter(aes(color=factor(period)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_p02_only_tab_for_plot.g$color[order(filt_major_taxa_proportions_p02_only_tab_for_plot.g$period)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none") +
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="Period 02 samples")

# Period p03 samples
major_taxa_p03 <- ggplot(filt_major_taxa_proportions_p03_only_tab_for_plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,100)) + # adding a setting for the y axis range so the rock and water plots are on the same scale
  geom_jitter(aes(color=factor(period)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_p03_only_tab_for_plot.g$color[order(filt_major_taxa_proportions_p03_only_tab_for_plot.g$period)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none") +
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="Period 03 samples")

# Period p04 samples
major_taxa_p04 <- ggplot(filt_major_taxa_proportions_p04_only_tab_for_plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,100)) + # adding a setting for the y axis range so the rock and water plots are on the same scale
  geom_jitter(aes(color=factor(period)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_p04_only_tab_for_plot.g$color[order(filt_major_taxa_proportions_p04_only_tab_for_plot.g$period)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none") +
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="Period 04 samples")

# Period p05 samples
major_taxa_p05 <- ggplot(filt_major_taxa_proportions_p05_only_tab_for_plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,100)) + # adding a setting for the y axis range so the rock and water plots are on the same scale
  geom_jitter(aes(color=factor(period)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_p05_only_tab_for_plot.g$color[order(filt_major_taxa_proportions_p05_only_tab_for_plot.g$period)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none") +
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="Period 05 samples")

# Period p06 samples
major_taxa_p06 <- ggplot(filt_major_taxa_proportions_p06_only_tab_for_plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,100)) + # adding a setting for the y axis range so the rock and water plots are on the same scale
  geom_jitter(aes(color=factor(period)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_p06_only_tab_for_plot.g$color[order(filt_major_taxa_proportions_p06_only_tab_for_plot.g$period)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none") +
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="Period 06 samples")
#################################################################




################################################################
################################################################

# https://benjjneb.github.io/dada2/tutorial.html

# Transform data to proportions as appropriate for Bray-Curtis distances
ASV_physeq.prop <- transform_sample_counts(ps.gbp23, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ASV_physeq.prop, method="NMDS", distance="bray")

# Plot with explicit colors
plot_ordination(ASV_physeq.prop, ord.nmds.bray, color="period", title="Bray NMDS") +
  scale_color_manual(values = viridis_palette)

# If we want to focusing without outliers
plot_ordination(ASV_physeq.prop, ord.nmds.bray, color="period", title="Bray NMDS") +
  scale_color_manual(values = viridis_palette) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))


# Bar plot
top20 <- names(sort(taxa_sums(ps.gbp23), decreasing=TRUE))[1:20]
ASV_physeq.top20 <- transform_sample_counts(ps.gbp23, function(OTU) OTU/sum(OTU))
ASV_physeq.top20 <- prune_taxa(top20, ASV_physeq.top20)
plot_bar(ASV_physeq.top20, x="period", fill="Genus") #+ facet_wrap(~When, scales="free_x")

plot_bar(ASV_physeq.top20, x="Sample", fill="Genus") #+ facet_wrap(~When, scales="free_x")


################################################################
# filter out samples we don’t want to include in our analysis such as the extraction and pcr blanks

# LJC   # We removed this negative controls before
# No Whites (NW)
#ASV_physeq.NW <- ASV_physeq %>%
#  subset_samples(type == "sample") %>%
#  prune_taxa(taxa_sums(.) > 0, .)
################################################################

# Now we will filter out Eukaryotes, Archaea, chloroplasts and mitochondria, because we only intended to amplify bacterial sequences
ASV_physeq.NW.Bact <- ps.gbp23 %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Order   != "Chloroplast"
  )

ASV_physeq.NW.Bact

# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(ASV_physeq.NW.Bact))

ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

# mean, max and min of sample read counts
smin <- min(sample_sums(ASV_physeq.NW.Bact))
smean <- mean(sample_sums(ASV_physeq.NW.Bact))
smax <- max(sample_sums(ASV_physeq.NW.Bact))


# melt to long format (for ggploting) 
# prune out phyla below 2% in each sample

GorBEEa_2023_phylum <- ASV_physeq.NW.Bact %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum



ggplot(GorBEEa_2023_phylum, aes(x = period, y = Abundance, fill = Phylum)) + 
  #facet_grid(site~.) +
  geom_bar(stat = "identity", position="fill") +
  scale_fill_viridis_d(option = "viridis") +
  # scale_fill_manual(values =sample_info_tab$color_p) +
  #scale_x_discrete(
  # breaks = c("7/8", "8/4", "9/2", "10/6"),
  # labels = c("Jul", "Aug", "Sep", "Oct"), 
  #drop = FALSE
  # ) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition of GorBEEa \n Bacterial Communities by Period") 



# melt to long format (for ggploting) 
# prune out Genera below 5% in each sample

# Set colors for plotting
genera_colors <- c("#D32F2F", # Deep Red
                   "#2C7975", # Teal Green
                   "#FFC107", # Amber Yellow
                   "#9967CE", # Lavender Purple
                   "#CD9BCD", # Lilac
                   "#43978D", # Deep Aqua
                   "#E91E63", # Pink
                   "#522157", # Deep Purple
                   "#AF4474", # Rose Wine
                   "#1F2F98", # Royal Blue
                   "#F06292", # Light Pink
                   "#FD8F52", # Coral Orange
                   "#7BE495", # Mint Green
                   "#0191B4", # Cyan
                   "#A5CAD2", # Soft Blue
                   "#D3DD18", # Lime Green
                   "#4378A2"  # Steel Blue
)

GorBEEa_2023_genus <- ASV_physeq.NW.Bact %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.05) %>%                         # Filter out low abundance taxa
  arrange(Genus)                                      # Sort data frame alphabetically by phylum

genus_composition_indv <- ggplot(GorBEEa_2023_genus, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity") +
  #scale_fill_viridis_d(option = "viridis") + # I do not recommend this color palette for many genera
  scale_fill_manual(values =genera_colors) + # Color palette customized for this dataset and 0.05 value
  # Remove x axis title
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)) +
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genus > 5%) \n") +
  ggtitle("Genus Composition of GorBEEa \n Bacterial Communities by specimen")

plot(genus_composition_indv)

jpeg(filename = "genus_composition_indv.jpeg", width = 800, height = 600, quality = 100)
print(genus_composition_indv)
dev.off()

ggsave("/Users/luisja/Downloads/genus_composition_indv.png", plot = genus_composition_period, width = 8, height = 6, units = "in", dpi = 300)

###############################
###
###     Group by periods
###
##############################

sample_counts <- GorBEEa_2023_genus %>%
  group_by(period) %>%
  summarise(n_samples = n_distinct(Sample))

genus_composition_period <- ggplot(GorBEEa_2023_genus, aes(x = period, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position="fill") +
  #scale_fill_viridis_d(option = "viridis") + # I do not recommend this color palette for many genera
  scale_fill_manual(values =genera_colors) + # Color palette customized for this dataset and 0.05 value
  geom_text(data = sample_counts, aes(x = period, y = 0, label = paste("n=", n_samples)), 
            inherit.aes = FALSE, vjust = 1.5, size = 3, color = "black") +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genus > 5%) \n") +
  ggtitle("Genus Composition of GorBEEa \n Bacterial Communities by Period")

plot(genus_composition_period)

jpeg(filename = "genus_composition_period.jpeg", width = 800, height = 600, quality = 100)
print(genus_composition_period)
dev.off()

ggsave("/Users/luisja/Downloads/genus_composition_period.pdf", plot = genus_composition_period, width = 8, height = 6, units = "in", dpi = 300)

###
###   Relative Abundance by site
###

# Calcular el número de muestras por 'site'

sample_counts <- GorBEEa_2023_genus %>%
  group_by(site) %>%
  summarise(n_samples = n_distinct(Sample))

genus_composition_site <- ggplot(GorBEEa_2023_genus, aes(x = site, y = Abundance, fill = Genus)) + 
  #facet_grid(site~.) +
  geom_bar(stat = "identity", position="fill") +
  geom_text(data = sample_counts, aes(x = site, y = 0, label = paste("n=", n_samples)), 
            inherit.aes = FALSE, vjust = 1.5, size = 3, color = "black") +
  #scale_fill_viridis_d(option = "viridis") +
  scale_fill_manual(values =genera_colors) +
  #scale_x_discrete(
  # breaks = c("7/8", "8/4", "9/2", "10/6"),
  # labels = c("Jul", "Aug", "Sep", "Oct"), 
  #drop = FALSE
  # ) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genus > 5%) \n") +
  ggtitle("Genus Composition of GorBEEa \n Bacterial Communities by Site") 

plot(genus_composition_site)

jpeg(filename = "genus_composition_site.jpeg", width = 800, height = 600, quality = 100)
print(genus_composition_site)
dev.off()

ggsave("genus_composition_site.png", plot = genus_composition_site, width = 8, height = 6, units = "in", dpi = 300)


### Group by site and period

# Contar muestras por sitio y periodo

sample_counts <- GorBEEa_2023_genus %>%
  group_by(site, period) %>%
  summarise(n_samples = n_distinct(Sample))

genus_composition_site_period <- ggplot(GorBEEa_2023_genus, aes(x = period, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position="fill") +
  facet_wrap(site~.) +
  #geom_bar(stat = "identity") +
  #scale_fill_viridis_d(option = "viridis") +
  scale_fill_manual(values = genera_colors) +
  geom_text(data = sample_counts, aes(x = period, y = -0.05, label = paste("n=", n_samples)), inherit.aes = FALSE, vjust = 0.5, size = 3) +
  #scale_x_discrete(
  # breaks = c("7/8", "8/4", "9/2", "10/6"),
  # labels = c("Jul", "Aug", "Sep", "Oct"), 
  #drop = FALSE
  # ) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genus > 5%) \n") +
  ggtitle("Genus Composition of GorBEEa \n Bacterial Communities by Site & Period") 

plot(genus_composition_site_period)

jpeg(filename = "genus_composition_site_period.jpeg", width = 800, height = 600, quality = 100)
print(genus_composition_site_period)
dev.off()

ggsave("genus_composition_site_period.png", plot = genus_composition_site_period, width = 8, height = 6, units = "in", dpi = 300)