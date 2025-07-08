data / ASV_table_0.1.tsv
data / ps.gbp23_c.0.1.RDS

# Load necessary library
library(dplyr)
library(iNEXT)
library(ggplot2)

# Read the TSV file
df_asv <- read.delim("data/ASV_table_0.1.tsv", header = TRUE, sep = "\t")

# Take a quick look at the structure of the data
glimpse(df_asv)


# Assuming df_asv is already loaded, and ASVs are rows, samples are columns.

# 1. Set ASVs as rownames and remove ASV column (you did this already)
rownames(df_asv) <- df_asv$ASV
df_asv <- df_asv[, -1]

# 2. Transpose so samples are rows, ASVs are columns
df_asv_t <- t(df_asv)

# 3. Convert each sample row into a numeric vector (not dataframe!)
asv_list <- apply(df_asv_t, 1, function(x) as.numeric(x))

# 4. `apply` returns a matrix for numeric input, so convert to list
asv_list <- split(asv_list, seq(nrow(df_asv_t)))

# But `split` won't work properly here; better do:
asv_list <- lapply(1:nrow(df_asv_t), function(i) as.numeric(df_asv_t[i, ]))
names(asv_list) <- rownames(df_asv_t)

# Now run iNEXT
asv_inext <- iNEXT(asv_list, q = 0, datatype = "abundance", size = NULL)

# Plot results
ggiNEXT(asv_inext) + theme(legend.position = "none")


# Plot the rarefaction curves
ggiNEXT(asv_inext, type = 1) +
  ggtitle("Rarefaction and Extrapolation Curves") +
  theme_minimal() +
  theme(legend.position = "none")


library(ggplot2)

p <- ggiNEXT(asv_inext, type = 1) +
  ylab("ASV diversity") +
  theme_minimal() +
  theme(legend.position = "none")

# Modify line size in the plot object
p$layers <- lapply(p$layers, function(layer) {
  if ("GeomLine" %in% class(layer$geom)) {
    layer$aes_params$linewidth <- 1
  }
  layer
})

p # now plot with thinner lines


# Plot the sampling completeness curves
ggiNEXT(asv_inext, type = 2) +
  ylab("ASV diversity") +
  theme_minimal() +
  theme(legend.position = "none")

asv_inext$iNextEst
str(asv_inext)
colnames(asv_inext$iNextEst)

#This shows the 10 samples with highest completeness.

asv_inext$DataInfo %>%
  select(Assemblage, n, S.obs, SC) %>%
  arrange(desc(SC)) %>%
  head(10)


summary(asv_inext$iNextEst$coverage_based)

mean(asv_inext$DataInfo$SC)
summary(asv_inext$DataInfo$SC)


asv_inext$AsyEst$completeness <- asv_inext$AsyEst$Observed /
  asv_inext$AsyEst$Estimator

mean(asv_inext$AsyEst$completeness)
range(asv_inext$AsyEst$completeness)


## change through time in bacterial diversity
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install("phyloseq")

library(phyloseq)
packageVersion("phyloseq")

ps.gbp23 <- readRDS("data/ps.gbp23_c.0.1.RDS")
print(ps.gbp23)
sample_data_tab.cl <- as(sample_data(ps.gbp23), "data.frame")
glimpse(sample_data_tab.cl)

rownames(sample_data_tab.cl)
# Convert row names in df2 to a column
sample_data_tab.cl$Assemblage <- rownames(sample_data_tab.cl)

#join both datasets

asv_inext$AsyEst

# Merge by species_id
merged_df <- merge(asv_inext$AsyEst, sample_data_tab.cl, by = "Assemblage")
glimpse(merged_df)

str(merged_df)

merged_df$Diversity <- as.factor(merged_df$Diversity)
merged_df$period <- as.factor(merged_df$period)
merged_df$site <- as.factor(merged_df$site)
library(ggplot2)

str(merged_df)

ggplot(merged_df, aes(x = period, y = Observed)) +
  geom_violin() +
  facet_wrap(~Diversity, scales = "free_y") +
  theme_minimal() +
  labs(x = "Period", y = "Diversity")

ggplot(merged_df, aes(x = site, y = Observed)) +
  geom_violin() +
  facet_wrap(~Diversity, scales = "free_y") +
  theme_minimal() +
  labs(x = "Period", y = "Diversity")


# library(RColorBrewer)
# pal<- brewer.pal(n = 6, name = "Set2")
pal <- c("#3A9AB2", "#91BAB6", "#BDC881", "#E3B710", "#EC7A05", "#F11B00")
p1 <- ggplot(merged_df, aes(x = period, y = Observed, fill = period)) +
  geom_violin(alpha = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = pal) +
  theme_minimal() +
  labs(x = "Period", y = "Diversity") +
  facet_wrap(~Diversity, scales = "free_y")
p1

# Save the plot
ggsave(
  "results/violin_plot_diversity_period.png",
  plot = p1,
  width = 8,
  height = 5,
  dpi = 400
)

ggplot(merged_df, aes(x = site, y = Observed, fill = site)) +
  geom_violin(alpha = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
  theme_minimal() +
  labs(x = "Site", y = "Diversity") +
  facet_wrap(~Diversity, scales = "free_y")


#change to numeric to quantitatively assess changes through time

merged_df$period_num <- as.numeric(sub(
  "^p",
  "",
  as.character(merged_df$period)
))


library(dplyr)

model_results <- merged_df %>%
  group_by(Diversity) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(Observed ~ period_num, data = .x)),
    tidy_model = map(model, tidy)
  ) %>%
  unnest(tidy_model)


model_results %>%
  filter(term == "period_num") %>%
  select(Diversity, estimate, std.error, statistic, p.value)


ggplot(merged_df, aes(x = period_num, y = Observed)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~Diversity, scales = "free_y") +
  theme_minimal() +
  labs(x = "Period", y = "Diversity")


# merge with actual names of ASVs

df_asv

ps.gbp23 <- readRDS("data/ps.gbp23_c.0.1.RDS")
count_tab.cl <- otu_table(ps.gbp23)
count_tab.cl <- as.data.frame(count_tab.cl)
sample_data_tab.cl <- as(sample_data(ps.gbp23), "data.frame")
tax_table.cl <- as.data.frame(tax_table(ps.gbp23))


rownames(tax_table.cl)

tax_table.cl$ASV_code <- rownames(tax_table.cl)
df_asv$ASV_code <- rownames(df_asv)

# Merge by species_id

merged_tax <- merge(df_asv, tax_table.cl, by = "ASV_code")


glimpse(merged_df)
glimpse(merged_tax)
head(merged_tax)

library(dplyr)
library(tidyr)

# Filtrar órdenes no nulos
merged_tax_filtered <- merged_tax %>%
  filter(!is.na(Order))

# Seleccionar columnas de muestras (que comienzan con 'GBP')
sample_columns <- grep("^GBP", names(merged_tax_filtered), value = TRUE)

# Calcular el número total de muestras
n_samples <- length(sample_columns)


# Convertir los rownames de sample_data_tab.cl a una columna para hacer el join
sample_data_tab.cl <- sample_data_tab.cl %>%
  tibble::rownames_to_column(var = "Sample")

long_df <- merged_tax_filtered %>%
  pivot_longer(
    cols = all_of(sample_columns),
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  filter(Abundance > 0) %>%
  select(Order, Sample) %>%
  left_join(
    sample_data_tab.cl %>% select(Sample, period, site),
    by = c("Sample" = "Sample")
  ) %>%
  distinct()


# Calcular prevalencia por orden
order_prevalence <- long_df %>%
  group_by(Order, period) %>%
  summarise(n_samples_present = n_distinct(Sample), .groups = "drop") %>%
  mutate(prevalence = (n_samples_present / n_samples) * 100) %>%
  arrange(desc(prevalence))

top_orders <- order_prevalence %>%
  group_by(Order) %>%
  summarise(max_prevalence = max(prevalence)) %>%
  slice_max(max_prevalence, n = 15) %>%
  pull(Order)

order_prevalence_top <- order_prevalence %>%
  filter(Order %in% top_orders)

library(ggplot2)
pal <- c("#3A9AB2", "#91BAB6", "#BDC881", "#E3B710", "#EC7A05", "#F11B00")


ggplot(
  order_prevalence_top,
  aes(x = reorder(Order, prevalence), y = prevalence, fill = period)
) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_text(
    aes(label = paste0(round(prevalence, 1), "%")),
    position = position_dodge(width = 0.9),
    hjust = -0.1,
    size = 3
  ) +
  coord_flip() +
  scale_fill_manual(values = pal) +
  labs(
    x = "Order",
    y = "Prevalence (%)",
    fill = "Period"
  ) +
  theme_minimal(base_size = 14)


# Filtrar familias no nulos
merged_tax_filtered <- merged_tax %>%
  filter(!is.na(Family))

# Seleccionar columnas de muestras (que comienzan con 'GBP')
sample_columns <- grep("^GBP", names(merged_tax_filtered), value = TRUE)

# Calcular el número total de muestras
n_samples <- length(sample_columns)


# # Convertir los rownames de sample_data_tab.cl a una columna para hacer el join
# sample_data_tab.cl <- sample_data_tab.cl %>%
#   tibble::rownames_to_column(var = "Sample_ID")

long_df <- merged_tax_filtered %>%
  pivot_longer(
    cols = all_of(sample_columns),
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  filter(Abundance > 0) %>%
  select(Family, Sample) %>%
  left_join(
    sample_data_tab.cl %>% select(Sample, period, site),
    by = c("Sample" = "Sample")
  ) %>%
  distinct()


# Calcular prevalencia por orden
Family_prevalence <- long_df %>%
  group_by(Family, period) %>%
  summarise(n_samples_present = n_distinct(Sample), .groups = "drop") %>%
  mutate(prevalence = (n_samples_present / n_samples) * 100) %>%
  arrange(desc(prevalence))

top_Familys <- Family_prevalence %>%
  group_by(Family) %>%
  summarise(max_prevalence = max(prevalence)) %>%
  slice_max(max_prevalence, n = 15) %>%
  pull(Family)

Family_prevalence_top <- Family_prevalence %>%
  filter(Family %in% top_Familys)

library(ggplot2)
pal <- c("#3A9AB2", "#91BAB6", "#BDC881", "#E3B710", "#EC7A05", "#F11B00")


ggplot(
  Family_prevalence_top,
  aes(x = reorder(Family, prevalence), y = prevalence, fill = period)
) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_text(
    aes(label = paste0(round(prevalence, 1), "%")),
    position = position_dodge(width = 0.9),
    hjust = -0.1,
    size = 3
  ) +
  coord_flip() +
  scale_fill_manual(values = pal) +
  labs(
    x = "Family",
    y = "Prevalence (%)",
    fill = "Period"
  ) +
  theme_minimal(base_size = 14)


# Filtrar class no nulos
merged_tax_filtered <- merged_tax %>%
  filter(!is.na(Class))

# Seleccionar columnas de muestras (que comienzan con 'GBP')
sample_columns <- grep("^GBP", names(merged_tax_filtered), value = TRUE)

# Calcular el número total de muestras
n_samples <- length(sample_columns)


# Convertir los rownames de sample_data_tab.cl a una columna para hacer el join
# sample_data_tab.cl <- sample_data_tab.cl %>%
#   tibble::rownames_to_column(var = "Sample_ID")

long_df <- merged_tax_filtered %>%
  pivot_longer(
    cols = all_of(sample_columns),
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  filter(Abundance > 0) %>%
  select(Class, Sample) %>%
  left_join(
    sample_data_tab.cl %>% select(Sample, period),
    by = c("Sample" = "Sample")
  ) %>%
  distinct()


# Calcular prevalencia por orden
Class_prevalence <- long_df %>%
  group_by(Class, period) %>%
  summarise(n_samples_present = n_distinct(Sample), .groups = "drop") %>%
  mutate(prevalence = (n_samples_present / n_samples) * 100) %>%
  arrange(desc(prevalence))

top_Classs <- Class_prevalence %>%
  group_by(Class) %>%
  summarise(max_prevalence = max(prevalence)) %>%
  slice_max(max_prevalence, n = 15) %>%
  pull(Class)

Class_prevalence_top <- Class_prevalence %>%
  filter(Class %in% top_Classs)

library(ggplot2)
pal <- c("#3A9AB2", "#91BAB6", "#BDC881", "#E3B710", "#EC7A05", "#F11B00")


ggplot(
  Class_prevalence_top,
  aes(x = reorder(Class, prevalence), y = prevalence, fill = period)
) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_text(
    aes(label = paste0(round(prevalence, 1), "%")),
    position = position_dodge(width = 0.9),
    hjust = -0.1,
    size = 3
  ) +
  coord_flip() +
  scale_fill_manual(values = pal) +
  labs(
    x = "Class",
    y = "Prevalence (%)",
    fill = "Period"
  ) +
  theme_minimal(base_size = 14)


##NMDS to see changes in composition through time

long_df_comp <- merged_tax_filtered %>%
  pivot_longer(
    cols = all_of(sample_columns),
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  filter(Abundance > 0) %>%
  left_join(
    sample_data_tab.cl %>% select(Sample, period, site),
    by = c("Sample" = "Sample")
  ) %>%
  distinct()


library(vegan)
library(tidyr)
library(dplyr)
library(tibble)


# Step 1: Sum abundances per sample-ASV combo (if there are duplicates)
long_df_comp_sum <- long_df_comp %>%
  group_by(Sample, ASV_code) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Step 2: Pivot to wide format
community_matrix <- long_df_comp_sum %>%
  pivot_wider(
    names_from = ASV_code,
    values_from = Abundance,
    values_fill = 0
  )

# Step 3: Set rownames to Sample
library(tibble)
community_matrix <- community_matrix %>%
  tibble::column_to_rownames("Sample")


library(vegan)

# Step 4: NMDS
nmds_result <- metaMDS(community_matrix, distance = "bray", k = 2, trymax = 100)

# Step 5: Add metadata (e.g., period)
metadata <- sample_data_tab.cl %>%
  filter(Sample %in% rownames(community_matrix)) %>%
  tibble::column_to_rownames("Sample")

# Step 6: PERMANOVA

adonis_result <- adonis2(
  community_matrix ~ period + site,
  data = metadata,
  method = "bray",
  permutations = 999,
  by = "margin"
)

adonis_result

# Get NMDS scores (site/samples coordinates)
nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))

# Add sample names
nmds_scores$Sample <- rownames(nmds_scores)

# Add metadata (like 'period')
metadata <- sample_data_tab.cl %>%
  filter(Sample %in% rownames(community_matrix)) # ensure same samples

# Merge NMDS coords with metadata
nmds_df <- left_join(nmds_scores, metadata, by = "Sample")


library(ggplot2)

# Your custom palette
pal <- c("#3A9AB2", "#91BAB6", "#BDC881", "#E3B710", "#EC7A05", "#F11B00")

p2 <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = period)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = pal) +
  labs(title = "NMDS of ASV Composition", color = "Period") +
  theme_minimal(base_size = 14)

p2
# Save the plot
ggsave("results/nmds.png", plot = p2, width = 8, height = 5, dpi = 400)

#remove outliers
nmds_df_sub <- subset(nmds_df, nmds_df$NMDS1 > -1.5)

p3 <- ggplot(nmds_df_sub, aes(x = NMDS1, y = NMDS2, color = period)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = pal) +
  labs(title = "NMDS of ASV Composition", color = "Period") +
  theme_minimal(base_size = 14)
p3

# Save the plot
ggsave("results/nmds_subset.png", plot = p3, width = 8, height = 5, dpi = 400)


ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = site)) +
  geom_point(size = 3, alpha = 0.8) +
  #scale_color_manual(values = pal) +
  labs(title = "NMDS of ASV Composition", color = "Period") +
  theme_minimal(base_size = 14)

#add information on plant species diversity per site

#use data from flower transects
pacman::p_load(rgbif, assertr, dplyr, tidyverse, ggplot2, lme4, car) #instala y lee paquetes a la vez; sustituye a Install.packages() y library()

df_fl <- read.csv("data/Floral_resources_2023.csv")
print(df_fl)
## Check taxonomy
df_fl2 <- as.data.frame(
  rgbif::name_backbone_checklist(df_fl$Planta) |>
    assert(in_set("EXACT", "FUZZY", "HIGHERRANK", "NONE"), matchType)
)

none <- subset(df_fl2, df_fl2$matchType == "NONE")

#clean those that can be cleaned and re-run
df_fl$Planta <- dplyr::recode(
  df_fl$Planta,
  "Erica sp." = "Erica",
  "Hutchinisia alpina" = "Hutchinsia alpina",
  "Ericea sp." = "Erica",
  "Trifolium platino" = "Trifolium pratense",
  "Globulana sp." = "Globularia vulgaris",
  "Cragaetus sp." = "Crataegus monogyna",
  "Thymus sp." = "Thymus praecox",
  "Cerasium sp." = "Cerastium fontanum"
)


df_fl3 <- cbind(df_fl, df_fl2)


df_fl4 <- subset(
  df_fl3,
  df_fl3$matchType == "EXACT" |
    df_fl3$matchType == "FUZZY" |
    df_fl3$matchType == "HIGHERRANK"
)
head(df_fl4)


#get plant richness per site and period
df_fl5 <- df_fl4 %>%
  group_by(Periodo, Sitio) %>%
  summarise(n_plants = n_distinct(species))
head(df_fl5)


df_fl5$j <- paste(df_fl5$Periodo, df_fl5$Sitio)


merged_df$site_num <- as.numeric(sub(
  "^s",
  "",
  as.character(merged_df$site)
))


merged_df$j <- paste(merged_df$period_num, merged_df$site_num)


merged_df_fl <- left_join(merged_df, df_fl5, by = "j")


pal <- c("#3A9AB2", "#91BAB6", "#BDC881", "#E3B710", "#EC7A05", "#F11B00")

p4 <- ggplot(merged_df_fl, aes(x = n_plants, y = Observed)) +
  geom_point(aes(color = period)) +
  geom_smooth(method = "lm", color = "black") + # single line, not per period
  facet_wrap(~Diversity, scales = "free_y") +
  scale_color_manual(values = pal) +
  xlab("Plant species richness") +
  ylab("Observed ASV diversity/richness") +
  theme_minimal()

p4
ggsave(
  "results/flower_richness.png",
  plot = p4,
  width = 8,
  height = 5,
  dpi = 400
)


ggplot(merged_df_fl, aes(x = n_plants, y = Observed)) +
  geom_point(aes(color = site)) +
  geom_smooth(method = "lm", color = "black") + # single line, not per period
  facet_wrap(~Diversity, scales = "free_y") +
  xlab("Plant species richness") +
  ylab("Observed ASV diversity/richness") +
  # scale_color_manual(values = pal) +
  theme_minimal()

ggsave(
  "results/flower_richness.png",
  plot = p4,
  width = 8,
  height = 5,
  dpi = 400
)

library(dplyr)

model_results_fl <- merged_df_fl %>%
  group_by(Diversity) %>%
  nest() %>%
  mutate(
    model = map(
      data,
      ~ lm(Observed ~ period_num + site_num + n_plants, data = .x)
    ),
    tidy_model = map(model, tidy)
  ) %>%
  unnest(tidy_model)


model_results_fl %>%
  filter(term %in% c("period_num", "site_num", "n_plants")) %>%
  select(Diversity, estimate, std.error, statistic, p.value)

#
#
# base_data$j <- paste(base_data$period, base_data$site)
#
#
# base_data$plant_sps_flor_res <- df_fl5$n_plants[match(base_data$j, df_fl5$j)]
#
#
# ggplot(
#   base_data,
#   aes(
#     x = plant_sps_flor_res,
#     y = presence_count,
#     color = factor(site),
#     group = site
#   )
# ) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, aes(group = site))
# # Add a line of fit without confidence interval shading
#
# ggplot(base_data, aes(x = plant_sps_flor_res, y = presence_count)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE) +
#   facet_grid(~site)
# # Add a line of fit without confidence interval shading
#
# m1 <- lmer(presence_count ~ plant_sps_flor_res + (1 | site), data = base_data)
# summary(m1)
# car::Anova(m1)
#
#
# #keep only sites with at least measures for 3 periods of time
#
# # Filter the dataset
# filtered_data <- base_data %>%
#   group_by(site) %>% # Group by site
#   filter(n_distinct(period) >= 2) %>% # Keep sites with at least 3 unique periods
#   ungroup() # Ungroup the data
#
# # View the filtered dataset
# filtered_data
#
# ggplot(filtered_data, aes(x = plant_sps_flor_res, y = presence_count)) +
#   geom_point() +
#   geom_smooth(method = "loess", se = FALSE)
# # Add a line of fit without confidence interval shading
#
# m1 <- lmer(
#   presence_count ~ plant_sps_flor_res + period + (1 | site),
#   data = filtered_data
# )
# summary(m1)
# car::Anova(m1)
#
#
# #stack data to calculate richness per order
#
# # Step 1: Keep columns 1 to 10
# df_micro_fixed <- df_micro[, 1:10]
#
# # Stack columns from the 11th onward
# df_micro_stacked <- df_micro[, 11:ncol(df_micro)]
#
# # Combine the fixed columns with the stacked data
# df_micro_long <- df_micro_stacked %>%
#   pivot_longer(
#     cols = everything(), # Pivot all remaining columns
#     names_to = "X", # New column for the original column names
#     values_to = "value" # New column for the values
#   ) %>%
#   bind_cols(df_micro_fixed[
#     rep(1:nrow(df_micro_fixed), each = ncol(df_micro_stacked)),
#   ]) # Repeat the first 10 columns for each stacked row
#
# head(df_micro_long)
#
# #combine with second dataset that has taxomonic information for each OTU
#
# df_micro_tax <- read.csv("data/microbiome/GPB23_tax_table_v2.csv")
# head(df_micro_tax)
# str(df_micro_tax)
#
# merged_df_micro <- merge(df_micro_long, df_micro_tax, by = "X", all.x = TRUE)
#
# head(merged_df_micro)
# summary(merged_df_micro)
# str(merged_df_micro)
#
#
# #remove Order Chloroplast might be plant material
#
# merged_df_micro2 <- subset(
#   merged_df_micro,
#   merged_df_micro$Order != "Chloroplast"
# )
# str(merged_df_micro2)
#
# # Filter rows where quant_reading > 0 and count occurrences by Order
# summary_by_order <- merged_df_micro2 %>%
#   filter(value > 0) %>% # Keep only rows where value > 0
#   group_by(Row.names, year, period, site, Order) %>% # Group by Row.names, year, period, site, and Order
#   summarise(order_count = n(), .groups = "drop") %>% # Count occurrences of each Order
#   ungroup() %>% # Remove grouping to prepare for stacking
#   pivot_longer(
#     cols = order_count, # Pivot the order_count column
#     names_to = "order", # Keep a column for the order names
#     values_to = "value"
#   ) # Keep a column for the counts
#
#
# # View the summarized dataset
# print(summary_by_order)
# str(summary_by_order)
# summary_by_order2 <- as.data.frame(summary_by_order)
# str(summary_by_order2)
#
# # Update 'period' and 'site' columns to keep only numbers
# summary_by_order2$period <- as.numeric(gsub(
#   "\\D",
#   "",
#   summary_by_order2$period
# )) # Remove non-digits
# summary_by_order2$site <- as.numeric(gsub("\\D", "", summary_by_order2$site)) # Remove non-digits
#
#
# df_fl5$j <- paste(df_fl5$Periodo, df_fl5$Sitio)
# summary_by_order2$j <- paste(summary_by_order2$period, summary_by_order2$site)
#
#
# summary_by_order2$plant_sps_flor_res <- df_fl5$n_plants[match(
#   summary_by_order2$j,
#   df_fl5$j
# )]
# summary_by_order2$plant_sps <- df5$n_plants[match(summary_by_order2$j, df5$j)]
#
#
# str(summary_by_order2)
#
# filtered_orders <- summary_by_order2 %>%
#   group_by(Order) %>%
#   summarise(n_specimens = n_distinct(Row.names), .groups = "drop") %>%
#   filter(n_specimens >= 60) # Keep orders with at least 10 distinct Row.names
#
# # Now, filter your original dataset to keep only those orders
# summary_by_order2_filtered <- summary_by_order2 %>%
#   filter(Order %in% filtered_orders$Order)
#
# # View the result
# head(summary_by_order2_filtered)
#
#
# ggplot(
#   summary_by_order2_filtered,
#   aes(x = plant_sps_flor_res, y = value, color = factor(Order), group = Order)
# ) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, aes(group = Order)) +
#   theme(legend.position = "none") # This removes the legend
#
#
# ggplot(
#   summary_by_order2_filtered,
#   aes(x = plant_sps, y = value, color = factor(Order), group = Order)
# ) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, aes(group = Order))
#
# library(ggplot2)
#
# # Manually define 6 color values
# colors <- c("#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20", "#29211F")
#
#
# ggplot(
#   summary_by_order2_filtered,
#   aes(x = plant_sps_flor_res, y = value, color = factor(Order), group = Order)
# ) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, aes(group = Order)) +
#   scale_color_manual(values = colors) + # Manually defined 6 colors
#   labs(
#     x = "Plant species richness", # Change x-axis label
#     y = "Bacteria species richness", # Change y-axis label
#     color = "Bacteria Order" # Change legend title
#   ) +
#   theme_classic() +
#   theme(
#     axis.title = element_text(size = 14), # Increase axis title font size
#     axis.text = element_text(size = 12), # Increase axis text font size
#     legend.title = element_text(size = 14), # Increase legend title font size
#     legend.text = element_text(size = 12), # Increase legend text font siz
#   ) +
#   facet_grid(~period, scales = "free")
#
#
# str(summary_by_order2_filtered)
#
# m1 <- lmer(
#   value ~ plant_sps * Order + (1 | site),
#   data = summary_by_order2_filtered
# )
# summary(m1)
# car::Anova(m1)
#
# m1 <- lm(value ~ plant_sps + period, data = summary_by_order2_filtered)
# summary(m1)
#
#
# m1 <- lmer(
#   value ~ plant_sps_flor_res + (1 | Order),
#   data = summary_by_order2_filtered
# )
# summary(m1)
# car::Anova(m1)
#
# m1 <- lm(value ~ plant_sps_flor_res, data = summary_by_order2_filtered)
# summary(m1)
#
# ##look at rickettsia sps
#
# rick <- subset(merged_df_micro2, merged_df_micro2$Order == "Rickettsiales")
#
# #keep only sites with at least measures for 3 periods of time
#
# # Filter the dataset
# filtered_data <- summary_by_order2 %>%
#   group_by(site) %>% # Group by site
#   filter(n_distinct(period) >= 2) %>% # Keep sites with at least 3 unique periods
#   ungroup() # Ungroup the data
#
# # View the filtered dataset
# filtered_data <- droplevels(filtered_data)
#
#
# ggplot(
#   filtered_data,
#   aes(x = plant_sps, y = value, color = factor(Order), group = Order)
# ) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, aes(group = Order)) +
#   theme(legend.position = "none") # This removes the legend
#
# m1 <- lmer(
#   value ~ plant_sps_flor_res + period + (1 | site),
#   data = filtered_data
# )
# summary(m1)
# car::Anova(m1)
#
#
# ##same but using Family
#
# # Filter rows where value > 0 and count occurrences by Order
# summary_by_fam <- merged_df_micro2 %>%
#   filter(value > 0) %>% # Keep only rows where value > 0
#   group_by(Row.names, year, period, site, Family) %>% # Group by Row.names, year, period, site, and Order
#   summarise(family_count = n(), .groups = "drop") %>% # Count occurrences of each Order
#   ungroup() %>% # Remove grouping to prepare for stacking
#   pivot_longer(
#     cols = family_count, # Pivot the order_count column
#     names_to = "family", # Keep a column for the order names
#     values_to = "value"
#   ) # Keep a column for the counts
#
#
# # View the summarized dataset
# print(summary_by_fam)
#
# summary_by_fam2 <- as.data.frame(summary_by_fam)
#
#
# # Update 'period' and 'site' columns to keep only numbers
# summary_by_fam2$period <- as.numeric(gsub("\\D", "", summary_by_fam2$period)) # Remove non-digits
# summary_by_fam2$site <- as.numeric(gsub("\\D", "", summary_by_fam2$site)) # Remove non-digits
#
#
# df_fl5$j <- paste(df_fl5$Periodo, df_fl5$Sitio)
# summary_by_fam2$j <- paste(summary_by_fam2$period, summary_by_fam2$site)
#
#
# summary_by_fam2$plant_sps_flor_res <- df_fl5$n_plants[match(
#   summary_by_fam2$j,
#   df_fl5$j
# )]
# summary_by_fam2$plant_sps <- df5$n_plants[match(summary_by_fam2$j, df5$j)]
#
#
# str(summary_by_fam2)
#
# summary_by_fam3 <- summary_by_fam2 %>%
#   drop_na()
#
#
# ggplot(
#   summary_by_fam3,
#   aes(x = plant_sps_flor_res, y = value, color = factor(Family), group = Family)
# ) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, aes(group = Family)) +
#   theme(legend.position = "none") # This removes the legend
#
#
# ggplot(
#   summary_by_fam3,
#   aes(x = plant_sps, y = value, color = factor(Family), group = Family)
# ) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, aes(group = Family)) +
#   theme(legend.position = "none") # This removes the legend
#
#
# m1 <- lmer(value ~ plant_sps_flor_res + (1 | site), data = summary_by_fam3)
# summary(m1)
# car::Anova(m1)
#
#
# filtered_fams <- summary_by_fam3 %>%
#   group_by(Family) %>%
#   summarise(n_specimens = n_distinct(Row.names), .groups = "drop") %>%
#   filter(n_specimens >= 60) # Keep fams with at least 10 distinct Row.names
#
# # Now, filter your original dataset to keep only those fams
# summary_by_fam2_filtered <- summary_by_fam2 %>%
#   filter(Family %in% filtered_fams$Family)
#
# # View the result
# head(summary_by_fam2_filtered)
#
#
# ggplot(
#   summary_by_fam2_filtered,
#   aes(x = plant_sps_flor_res, y = value, color = factor(Family), group = Family)
# ) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, aes(group = Family)) +
#   theme(legend.position = "none") # This removes the legend
#
#
# ggplot(
#   summary_by_fam2_filtered,
#   aes(x = plant_sps, y = value, color = factor(Family), group = Family)
# ) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, aes(group = Family)) +
#   theme(legend.position = "none") # This removes the legend
#
#
# m1 <- lm(value ~ plant_sps_flor_res, data = summary_by_fam2_filtered)
# summary(m1)
# car::Anova(m1)
#
# m1 <- lm(value ~ plant_sps, data = summary_by_fam2_filtered)
# summary(m1)
# car::Anova(m1)
#
#
# #keep only sites with at least measures for 3 periods of time
#
# # Filter the dataset
# filtered_data <- summary_by_fam3 %>%
#   group_by(site) %>% # Group by site
#   filter(n_distinct(period) >= 2) %>% # Keep sites with at least 3 unique periods
#   ungroup() # Ungroup the data
#
# # View the filtered dataset
# filtered_data
#
#
# ggplot(
#   filtered_data,
#   aes(x = plant_sps, y = value, color = factor(Family), group = Family)
# ) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, aes(group = Family)) +
#   theme(legend.position = "none") # This removes the legend
#
# ggplot(
#   filtered_data,
#   aes(x = plant_sps_flor_res, y = value, color = factor(Family), group = Family)
# ) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, aes(group = Family)) +
#   theme(legend.position = "none") # This removes the legend
#
# m1 <- lmer(
#   value ~ plant_sps_flor_res + period + (1 | site),
#   data = filtered_data
# )
# summary(m1)
# car::Anova(m1)
#
#
# ##same but using genus
#
# # Filter rows where value > 0 and count occurrences by Order
# summary_by_gen <- merged_df_micro2 %>%
#   filter(value > 0) %>% # Keep only rows where value > 0
#   group_by(Row.names, year, period, site, Genus) %>% # Group by Row.names, year, period, site, and Order
#   summarise(genus_count = n(), .groups = "drop") %>% # Count occurrences of each Order
#   ungroup() %>% # Remove grouping to prepare for stacking
#   pivot_longer(
#     cols = genus_count, # Pivot the order_count column
#     names_to = "genus", # Keep a column for the order names
#     values_to = "value"
#   ) # Keep a column for the counts
#
#
# # View the summarized dataset
# print(summary_by_gen)
#
# summary_by_gen2 <- as.data.frame(summary_by_gen)
#
#
# # Update 'period' and 'site' columns to keep only numbers
# summary_by_gen2$period <- as.numeric(gsub("\\D", "", summary_by_gen2$period)) # Remove non-digits
# summary_by_gen2$site <- as.numeric(gsub("\\D", "", summary_by_gen2$site)) # Remove non-digits
#
#
# df_fl5$j <- paste(df_fl5$Periodo, df_fl5$Sitio)
# summary_by_gen2$j <- paste(summary_by_gen2$period, summary_by_gen2$site)
#
#
# summary_by_gen2$plant_sps_flor_res <- df_fl5$n_plants[match(
#   summary_by_gen2$j,
#   df_fl5$j
# )]
# summary_by_gen2$plant_sps <- df5$n_plants[match(summary_by_gen2$j, df5$j)]
#
#
# str(summary_by_gen2)
#
# summary_by_gen3 <- summary_by_gen2 %>%
#   drop_na()
#
# summary_by_gen3 <- droplevels(summary_by_gen3)
#
# ggplot(
#   summary_by_gen3,
#   aes(x = plant_sps_flor_res, y = value, color = factor(Genus), group = Genus)
# ) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, aes(group = Genus)) +
#   theme(legend.position = "none") # This removes the legend
#
#
# ggplot(
#   summary_by_gen3,
#   aes(x = plant_sps, y = value, color = factor(Genus), group = Genus)
# ) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, aes(group = Genus)) +
#   theme(legend.position = "none") # This removes the legend
#
#
# m1 <- lmer(value ~ plant_sps_flor_res + (1 | site), data = summary_by_gen3)
# summary(m1)
# car::Anova(m1)
#
#
# filtered_gens <- summary_by_gen3 %>%
#   group_by(Genus) %>%
#   summarise(n_specimens = n_distinct(Row.names), .groups = "drop") %>%
#   filter(n_specimens >= 30) # Keep fams with at least 10 distinct Row.names
#
# # Now, filter your original dataset to keep only those fams
# summary_by_gen2_filtered <- summary_by_gen3 %>%
#   filter(Genus %in% filtered_gens$Genus)
#
# # View the result
# head(summary_by_gen2_filtered)
#
#
# ggplot(
#   summary_by_gen2_filtered,
#   aes(x = plant_sps_flor_res, y = value, color = factor(Genus), group = Genus)
# ) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, aes(group = Genus)) +
#   theme(legend.position = "none") + # This removes the legend
#   facet_grid(~period)
#
#
# ggplot(
#   summary_by_gen2_filtered,
#   aes(x = plant_sps, y = value, color = factor(Genus), group = Genus)
# ) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, aes(group = Genus)) +
#   theme(legend.position = "none") # This removes the legend
#
#
# m1 <- lm(value ~ plant_sps_flor_res, data = summary_by_gen2_filtered)
# summary(m1)
# car::Anova(m1)
#
# m1 <- lm(value ~ plant_sps, data = summary_by_gen2_filtered)
# summary(m1)
# car::Anova(m1)
#
#
# #keep only sites with at least measures for 3 periods of time
#
# # Filter the dataset
# filtered_data <- summary_by_gen3 %>%
#   group_by(site) %>% # Group by site
#   filter(n_distinct(period) >= 3) %>% # Keep sites with at least 3 unique periods
#   ungroup() # Ungroup the data
#
# # View the filtered dataset
# filtered_data
#
#
# ggplot(
#   filtered_data,
#   aes(x = plant_sps, y = value, color = factor(Genus), group = Genus)
# ) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, aes(group = Genus)) +
#   theme(legend.position = "none") # This removes the legend
#
# ggplot(
#   filtered_data,
#   aes(x = plant_sps_flor_res, y = value, color = factor(Genus), group = Genus)
# ) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, aes(group = Genus)) +
#   theme(legend.position = "none") # This removes the legend
#
# m1 <- lmer(
#   value ~ plant_sps_flor_res + period + (1 | site),
#   data = filtered_data
# )
# summary(m1)
# car::Anova(m1)
#
#
# m1 <- lm(value ~ plant_sps_flor_res, data = filtered_data)
# summary(m1)
# car::Anova(m1)
