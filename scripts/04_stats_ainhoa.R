# Load necessary libraries
library(dplyr)
library(iNEXT)
library(ggplot2)

# Read the TSV file
df_asv <- read.delim("data/ASV_table_0.1.tsv", header = TRUE, sep = "\t")

# Take a quick look at the structure of the data
glimpse(df_asv)


###############################
####SAMPLING COMPLETENESS######
##############################

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


ggsave(
  "results/sampling_completeness.png",
  plot = p,
  width = 4,
  height = 5,
  dpi = 400
)


save(p, file = "./results/RData/sampling_complete.RData")


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


############################################################
###ASSESS CHANGE THROUGH TIME AND WITH FLORAL DIVERSITY######
###########################################################

# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# BiocManager::install("phyloseq")

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


##keep only shannon diversity

merged_df_sub <- subset(merged_df, merged_df$Diversity == "Shannon diversity")


# First, convert period to numeric for the trend line
merged_df_sub$period_num <- as.numeric(factor(merged_df_sub$period)) # preserves order


p1_sh <- ggplot(merged_df_sub, aes(x = period, y = Observed, fill = period)) +
  geom_violin(alpha = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
  geom_smooth(
    data = merged_df_sub,
    aes(x = period_num, y = Observed),
    method = "lm",
    se = FALSE,
    color = "darkgrey",
    linewidth = 0.5,
    linetype = 2,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = pal) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Period", y = "ASV Diversity")
p1_sh

# Save the plot
ggsave(
  "results/violin_plot_SHANNON_diversity_period.png",
  plot = p1_sh,
  width = 4,
  height = 5,
  dpi = 400
)


save(p1_sh, file = "./results/RData/violin_plots_SHANNON_ASV_DIv.RData")


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
library(broom)

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

library(lme4)

m1 <- lmer(Observed ~ period_num + (1 | site), data = merged_df_sub)
summary(m1)
library(car)
car::Anova(m1)


m1 <- lm(Observed ~ period_num + site, data = merged_df_sub)
summary(m1)
library(car)
car::Anova(m1)


ggplot(merged_df, aes(x = period_num, y = Observed)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~Diversity, scales = "free_y") +
  theme_minimal() +
  labs(x = "Period", y = "Diversity")


#add information on plant species diversity per site

#use data from flower transects
pacman::p_load(rgbif, assertr, dplyr, tidyverse, ggplot2, lme4, car) #instala y lee paquetes a la vez; sustituye a Install.packages() y library()

df_fl <- read.csv("data/Floral_resources_2023.csv")
print(df_fl)
## Check taxonomy
# df_fl2 <- as.data.frame(
#   rgbif::name_backbone_checklist(df_fl$Planta) |>
#     assert(in_set("EXACT", "FUZZY", "HIGHERRANK", "NONE"), matchType)
# )

library(rgbif)
library(purrr)
library(dplyr)
library(stringi)

names_vec <- df_fl$Planta |> as.character() |> stri_trim_both() |> enc2utf8()

chunks <- split(names_vec, ceiling(seq_along(names_vec) / 100))

safe_fetch <- safely(
  function(v)
    name_backbone_checklist(
      v,
      curlopts = list(timeout = 120, connecttimeout = 60)
    ),
  otherwise = tibble()
)

res_list <- map(
  chunks,
  ~ {
    out <- safe_fetch(.x)
    if (!is.null(out$error)) {
      Sys.sleep(2)
      out <- safe_fetch(.x)
    }
    out$result
  }
)

df_fl2 <- bind_rows(res_list)

# sanity check: should match original length
stopifnot(nrow(df_fl2) == length(names_vec))


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

# write.csv(merged_df_fl, "data/merged_ASV_floral_diversity.csv", row.names = FALSE)

merged_df_fl <- read.csv("data/merged_ASV_floral_diversity.csv")

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


##only shannon div

merged_df_fl_sub <- subset(
  merged_df_fl,
  merged_df_fl$Diversity == "Shannon diversity"
)
p4_sh <- ggplot(merged_df_fl_sub, aes(x = n_plants, y = Observed)) +
  geom_point(aes(color = period)) +
  geom_smooth(method = "lm", color = "black") + # single line, not per period
  scale_color_manual(values = pal) +
  xlab("Plant species richness") +
  ylab("Observed ASV diversity") +
  theme_minimal()

p4_sh


save(p4_sh, file = "./results/RData/scatter_plant_ASV_Div.RData")


ggsave(
  "results/flower_richness_SHANNON.png",
  plot = p4_sh,
  width = 4,
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


##We include time and floral diversity in the same mixed model,
# separate within-site temporal effects from between-site differences,
# handling their correlation, and checking whether the “time effect” is
#actually mediated by floral diversity.

# We use a linear mixed model for Shannon diversity (it’s typically fine with Gaussian
#errors), with site as a random effect and random slopes for time. Add an AR(1)
#correlation because we have repeated measures through time.
#Cluster-mean centering (within–between decomposition) for floral diversity to
#reduce confounding with time and to interpret effects cleanly:
#This lets us answer: “When floral diversity at a site is higher than its own
#long-term average, is bacterial Shannon higher?” while also asking “Do sites
#with overall higher floral diversity also have higher bacterial Shannon?”

library(dplyr)

glimpse(merged_df_fl)

merged_df_fl_sub <- subset(
  merged_df_fl,
  merged_df_fl$Diversity == "Shannon diversity"
)

glimpse(merged_df_fl_sub)

merged_df_fl_sub2 <- merged_df_fl_sub %>%
  group_by(site_num) %>%
  mutate(
    floral_between = mean(n_plants, na.rm = TRUE),
    floral_within = n_plants - floral_between
  ) %>%
  ungroup()


#check for collinearity

m_in <- lm(
  Observed ~ Periodo + floral_within + floral_between,
  data = merged_df_fl_sub2
)


car::vif(m_in) #no collinearity


library(lme4)
m1 <- lmer(
  Observed ~
    Periodo + floral_within + floral_between + (1 + Periodo | site_num),
  data = merged_df_fl_sub2
)

summary(m1)
car::Anova(m1)

# If temporal autocorrelation matters or you want AR(1):
library(dplyr)

merged_df_fl_sub2 <- merged_df_fl_sub2 %>%
  arrange(site_num, Periodo) %>%
  mutate(
    # ensure Periodo is numeric (e.g., date -> as.numeric(date))
    Periodo_num = as.numeric(Periodo),
    # standardize predictors (helps optimizer & identifiability)
    Periodo_sc = as.numeric(scale(Periodo_num)),
    floral_within_sc = as.numeric(scale(floral_within)),
    floral_between_sc = as.numeric(scale(floral_between))
  )


library(glmmTMB)

m0 <- glmmTMB(
  Observed ~ Periodo_sc + floral_within_sc + floral_between_sc + (1 | site_num),
  data = merged_df_fl_sub2,
  family = gaussian(),
  ziformula = ~0,
  dispformula = ~1,
  REML = TRUE
)
summary(m0) # Should give SEs and no singular warning
m1 <- glmmTMB(
  Observed ~
    Periodo_sc +
      floral_within_sc +
      floral_between_sc +
      (1 + Periodo_sc || site_num), # note the double pipes
  data = merged_df_fl_sub2,
  family = gaussian(),
  ziformula = ~0,
  dispformula = ~1,
  REML = TRUE
)
summary(m1)
m2 <- glmmTMB(
  Observed ~
    Periodo_sc +
      floral_within_sc +
      floral_between_sc +
      (1 + Periodo_sc | site_num),
  data = merged_df_fl_sub2,
  family = gaussian(),
  ziformula = ~0,
  dispformula = ~1,
  REML = TRUE
)
summary(m2) # If this reintroduces singularity, stick with m1.


# Standardized model already fit as m1
confint(m2) # 95% CI
MuMIn::r.squaredGLMM(m2) # R2m, R2c

# Partial effect plot
library(ggeffects)
plot(ggpredict(m2, terms = "floral_within_sc [all]")) +
  ggplot2::labs(
    x = "Floral richness (within-site, z)",
    y = "Predicted Shannon H"
  )

###add depth of sequencing
merged_df_fl_sub2 <- merged_df_fl_sub2 |>
  dplyr::mutate(
    log_reads = scale(log(quant_reading)) # or whatever your depth column is
  )

m2_depth <- glmmTMB(
  Observed ~
    Periodo_sc +
      floral_within_sc +
      floral_between_sc +
      log_reads +
      (1 + Periodo_sc || site_num),
  data = merged_df_fl_sub2,
  family = gaussian(),
  ziformula = ~0,
  dispformula = ~1,
  REML = TRUE
)
summary(m2_depth)


# Standardized model already fit as m1
confint(m2_depth) # 95% CI
MuMIn::r.squaredGLMM(m2_depth) # R2m, R2c

# Partial effect plot
library(ggeffects)
plot(ggpredict(m2_depth, terms = "floral_within_sc [all]")) +
  ggplot2::labs(
    x = "Floral richness (within-site, z)",
    y = "Predicted Shannon H"
  )


p4_sh2 <- ggplot(merged_df_fl_sub2, aes(x = floral_within_sc, y = Observed)) +
  geom_point(aes(color = period)) +
  geom_smooth(method = "lm", color = "black") + # single line, not per period
  scale_color_manual(values = pal) +
  xlab("Plant species richness") +
  ylab("Observed ASV diversity") +
  theme_minimal()

p4_sh2


ggsave(
  "results/flower_richness_temporal_autocorrelation.png",
  plot = p4,
  width = 8,
  height = 5,
  dpi = 400
)

save(
  p4_sh2,
  file = "./results/RData/flower_richness_temporal_autocorrelation.RData"
)


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


plot_order <- ggplot(
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


ggsave(
  "results/barchart_order.png",
  plot = plot_order,
  width = 8,
  height = 5,
  dpi = 400
)

save(
  plot_order,
  file = "./results/RData/barchart_order.RData"
)


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


#####NMDS########
#################################################
## 0) Packages
#################################################
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)
library(permute)
library(ggplot2)

## Optional: create output dirs
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("results/RData", showWarnings = FALSE, recursive = TRUE)

#################################################
## 1) Build community matrix + aligned metadata
#################################################

# long table with ASV abundances and metadata
long_df_comp <- merged_tax_filtered %>%
  pivot_longer(
    cols = all_of(sample_columns), # <- your sample columns
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  filter(Abundance > 0) %>%
  left_join(
    sample_data_tab.cl %>% select(Sample, period, site),
    by = "Sample"
  ) %>%
  distinct()

# sum duplicates (Sample × ASV_code)
long_df_comp_sum <- long_df_comp %>%
  group_by(Sample, ASV_code) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# wide community matrix (samples × ASVs)
community_matrix <- long_df_comp_sum %>%
  pivot_wider(
    names_from = ASV_code,
    values_from = Abundance,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("Sample")

# aligned metadata for the same samples
metadata <- sample_data_tab.cl %>%
  filter(Sample %in% rownames(community_matrix)) %>%
  mutate(
    period = factor(trimws(as.character(period))),
    site = factor(trimws(as.character(site)))
  ) %>%
  tibble::column_to_rownames("Sample")

# final alignment check
stopifnot(all(rownames(community_matrix) == rownames(metadata)))

#################################################
## 2) NMDS (as before)
#################################################
set.seed(1)

# square-root transform abundances
community_matrix_sqrt <- sqrt(community_matrix)

# then use this transformed matrix
nmds_result <- metaMDS(
  community_matrix_sqrt,
  distance = "bray",
  k = 2,
  trymax = 100
)


# nmds_result <- metaMDS(community_matrix, distance = "bray", k = 2, trymax = 100)

# NMDS scores + metadata for plotting
nmds_df <- as.data.frame(scores(nmds_result, display = "sites")) %>%
  rownames_to_column("Sample") %>%
  left_join(
    sample_data_tab.cl %>%
      filter(Sample %in% rownames(community_matrix)) %>%
      select(Sample, period, site),
    by = "Sample"
  )

# palette (your colors)
pal <- c("#3A9AB2", "#91BAB6", "#BDC881", "#E3B710", "#EC7A05", "#F11B00")

p2 <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = period)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = pal) +
  labs(color = "Period") +
  theme_minimal(base_size = 14)

# convex hulls per period
hulls <- nmds_df %>%
  group_by(period) %>%
  slice(chull(NMDS1, NMDS2))

p2_hull <- p2 +
  geom_polygon(
    data = hulls,
    aes(x = NMDS1, y = NMDS2, fill = period, group = period),
    alpha = 0.2,
    color = NA
  ) +
  scale_fill_manual(values = pal) +
  guides(fill = "none")

ggsave("results/nmds.png", plot = p2_hull, width = 8, height = 5, dpi = 400)
save(p2_hull, file = "results/RData/nmds.RData")


####graphs per site

p2_site <- ggplot(nmds_df, aes(NMDS1, NMDS2, color = site)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_viridis_d() + # needs library(viridis)
  labs(color = "Site") +
  theme_minimal(base_size = 14)

# # convex hulls per period
# hulls_site <- nmds_df %>%
#   group_by(site) %>%
#   slice(chull(NMDS1, NMDS2))
#
# p2_site_hull <- p2_site +
#   geom_polygon(
#     data = hulls_site,
#     aes(x = NMDS1, y = NMDS2, fill = site, group = site),
#     alpha = 0.2,
#     color = NA
#   ) +
#   scale_color_viridis_d() +
#   guides(fill = "none")

ggsave(
  "results/nmds_site.png",
  plot = p2_site,
  width = 8,
  height = 5,
  dpi = 400
)
save(p2_site, file = "results/RData/nmds_site.RData")


#################################################
## 3) Global PERMANOVA (period + site)
#################################################
adonis_result <- adonis2(
  community_matrix_sqrt ~ period + site,
  data = metadata,
  method = "bray",
  permutations = 999,
  by = "margin" # marginal tests (like your original)
)
print(adonis_result)

#################################################
## 4) Dispersion check for 'period' (assumption)
#################################################
dist_bc <- vegdist(community_matrix_sqrt, method = "bray")
bd <- betadisper(dist_bc, group = metadata$period)
disp_test <- permutest(bd)
print(disp_test)
# If significant, note as a caveat in the text.
# plot(bd) # optional visualization

#################################################
## 5) Pairwise PERMANOVA between periods (vegan-only)
##    - permutations blocked within site when possible
#################################################

pairwise_permanova_period <- function(
  community_matrix_sqrt,
  metadata,
  period_col = "period",
  site_col = "site",
  method = "bray",
  nperm = 999,
  p_adjust = "BH",
  min_per_group = 2
) {
  stopifnot(all(rownames(community_matrix) == rownames(metadata)))

  meta <- metadata
  meta[[period_col]] <- droplevels(factor(meta[[period_col]]))
  meta[[site_col]] <- droplevels(factor(meta[[site_col]]))

  levs <- levels(meta[[period_col]])
  if (length(levs) < 2) {
    return(data.frame(
      pair = NA,
      n1 = NA,
      n2 = NA,
      F = NA,
      R2 = NA,
      p = NA,
      p_adj = NA,
      note = "Need at least 2 levels in period after alignment"
    ))
  }

  # Precompute Bray–Curtis once
  Dm <- as.matrix(vegdist(community_matrix, method = method))

  pairs <- t(combn(levs, 2))
  out <- vector("list", nrow(pairs))

  for (i in seq_len(nrow(pairs))) {
    lev <- pairs[i, ]
    idx <- meta[[period_col]] %in% lev
    sub_meta_all <- droplevels(meta[idx, , drop = FALSE])

    # need >= min_per_group per period
    tabP <- table(sub_meta_all[[period_col]])
    n1 <- as.integer(tabP[1])
    n2 <- as.integer(tabP[2])
    if (any(tabP < min_per_group)) {
      out[[i]] <- data.frame(
        pair = paste(lev, collapse = " vs "),
        n1 = n1,
        n2 = n2,
        F = NA,
        R2 = NA,
        p = NA,
        p_adj = NA,
        note = "Skipped: < min_per_group samples in a period"
      )
      next
    }

    # try blocked permutations within site: drop 1-sample blocks (cannot permute within them)
    tabS_all <- table(sub_meta_all[[site_col]])
    keep_sites <- names(tabS_all)[tabS_all >= 2]
    sub_meta <- droplevels(sub_meta_all[
      sub_meta_all[[site_col]] %in% keep_sites,
      ,
      drop = FALSE
    ])

    use_blocks <- TRUE
    if (nrow(sub_meta) < 2 || length(unique(sub_meta[[period_col]])) < 2) {
      sub_meta <- sub_meta_all
      use_blocks <- FALSE
    } else {
      tabP2 <- table(sub_meta[[period_col]])
      if (any(tabP2 < min_per_group)) {
        sub_meta <- sub_meta_all
        use_blocks <- FALSE
      }
    }

    # subset distance matrix and convert back to 'dist'
    ids <- rownames(sub_meta)
    Dij <- as.dist(Dm[ids, ids, drop = FALSE])

    # permutation control: permute within site blocks when possible
    ctrl <- if (use_blocks)
      permute::how(nperm = nperm, blocks = sub_meta[[site_col]]) else
      permute::how(nperm = nperm)

    # run adonis2 with distance response; formula references a local column
    a2 <- vegan::adonis2(
      Dij ~ sub_period,
      data = data.frame(
        sub_period = sub_meta[[period_col]],
        sub_site = sub_meta[[site_col]],
        row.names = ids
      ),
      permutations = ctrl,
      by = "terms"
    )

    Fval <- a2$F[1]
    R2 <- a2$R2[1]
    pval <- a2$`Pr(>F)`[1]
    note <- if (use_blocks) "Permutations blocked within site" else
      "Unconstrained permutations"

    out[[i]] <- data.frame(
      pair = paste(lev, collapse = " vs "),
      n1 = as.integer(table(sub_meta[[period_col]])[1]),
      n2 = as.integer(table(sub_meta[[period_col]])[2]),
      F = Fval,
      R2 = R2,
      p = pval,
      note = note
    )
  }

  res <- bind_rows(out)
  if (nrow(res) > 0 && any(!is.na(res$p))) {
    res$p_adj <- p.adjust(res$p, method = p_adjust)
  } else {
    res$p_adj <- NA
  }
  res[order(res$p_adj, res$p), , drop = FALSE]
}

# Run pairwise PERMANOVA
pw <- pairwise_permanova_period(
  community_matrix = community_matrix_sqrt,
  metadata = metadata,
  period_col = "period",
  site_col = "site",
  method = "bray",
  nperm = 999,
  p_adjust = "BH"
)
print(pw)

write.csv(pw, "results/pairwise_permanova_period.csv", row.names = FALSE)

#################################################
## 6) (Optional) quick heatmap of pairwise p_adj
#################################################
if (nrow(pw) > 0 && !all(is.na(pw$p_adj))) {
  # expand pairs for plotting
  split_pairs <- do.call(rbind, strsplit(pw$pair, " vs "))
  hp <- data.frame(
    A = split_pairs[, 1],
    B = split_pairs[, 2],
    minuslog10padj = -log10(pw$p_adj)
  )
  gg <- ggplot(hp, aes(A, B, fill = minuslog10padj)) +
    geom_tile() +
    scale_fill_viridis_c(name = "-log10(FDR p)") +
    theme_minimal(base_size = 12) +
    labs(x = "Period", y = "Period")
  ggsave(
    "results/pairwise_permanova_heatmap.png",
    gg,
    width = 6,
    height = 5,
    dpi = 300
  )
}


# install.packages("patchwork") # if not yet installed
library(patchwork)


# add panel tags A) and B)
combined <- p2_hull + gg + plot_annotation(tag_levels = 'A')

# save combined figure
ggsave(
  "results/nmds_pairwise_combined.png",
  combined,
  width = 12,
  height = 5,
  dpi = 300
)


save(combined, file = "results/RData/combined_nmds.RData")

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
