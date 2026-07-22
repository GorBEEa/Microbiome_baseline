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
  ylab("Bacterial ASV diversity") +
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
merged_df_rich <- subset(merged_df, merged_df$Diversity == "Species richness")

#use richness o shannon?
summary_by_period <- merged_df_sub %>%
  group_by(period) %>%
  summarise(
    n = n(),
    min = min(Observed, na.rm = TRUE),
    max = max(Observed, na.rm = TRUE),
    range = max - min,
    mean = mean(Observed, na.rm = TRUE),
    median = median(Observed, na.rm = TRUE),
    sd = sd(Observed, na.rm = TRUE)
  )
print(summary_by_period)


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
  ylab("Observed bacterial diversity") +
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

p4_sh_site <- ggplot(
  merged_df_fl_sub,
  aes(x = n_plants, y = Observed, color = site)
) +
  # geom_point(aes(color = site)) +
  geom_point() +
  #geom_smooth(method = "lm", color = "black") + # single line, not per period
  geom_smooth(method = "lm", se = FALSE) + # single line, not per period
  #scale_color_manual(values = pal) +
  xlab("Plant species richness") +
  ylab("Observed bacterial diversity") +
  theme_minimal()

p4_sh_site

ggsave(
  "results/flower_richness_bacterial_site_SHANNON.png",
  plot = p4_sh_site,
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


###betadiversity

########################
##########FUNGI#########
########################

# Read the TSV file
df_fung <- read.delim("data/ASV_table_fungi.tsv", header = TRUE, sep = "\t")

# Take a quick look at the structure of the data
glimpse(df_fung)


###############################
####SAMPLING COMPLETENESS######
##############################

# Assuming df_asv is already loaded, and ASVs are rows, samples are columns.

# 1. Set ASVs as rownames and remove ASV column (you did this already)
rownames(df_fung) <- df_fung$ASV
df_fung <- df_fung[, -1]

# 2. Transpose so samples are rows, ASVs are columns
df_fung_t <- t(df_fung)

# 3. Convert each sample row into a numeric vector (not dataframe!)
fung_list <- apply(df_fung_t, 1, function(x) as.numeric(x))

# 4. `apply` returns a matrix for numeric input, so convert to list
fung_list <- split(fung_list, seq(nrow(df_fung_t)))

# But `split` won't work properly here; better do:
fung_list <- lapply(1:nrow(df_fung_t), function(i) as.numeric(df_fung_t[i, ]))
names(fung_list) <- rownames(df_fung_t)

# Now run iNEXT
fung_inext <- iNEXT(fung_list, q = 0, datatype = "abundance", size = NULL)

# Plot results
ggiNEXT(fung_inext) + theme(legend.position = "none")


# Plot the rarefaction curves
ggiNEXT(fung_inext, type = 1) +
  ggtitle("Rarefaction and Extrapolation Curves") +
  theme_minimal() +
  theme(legend.position = "none")


library(ggplot2)

p_fung <- ggiNEXT(fung_inext, type = 1) +
  ylab("Fungal ASV diversity") +
  theme_minimal() +
  theme(legend.position = "none")

# Modify line size in the plot object
p_fung$layers <- lapply(p_fung$layers, function(layer) {
  if ("GeomLine" %in% class(layer$geom)) {
    layer$aes_params$linewidth <- 1
  }
  layer
})

p_fung # now plot with thinner lines


ggsave(
  "results/sampling_completeness_fungi.png",
  plot = p_fung,
  width = 6,
  height = 4,
  dpi = 400
)


save(p_fung, file = "./results/RData/sampling_complete_fungi.RData")


# Plot the sampling completeness curves
ggiNEXT(fung_inext, type = 2) +
  ylab("ASV diversity") +
  theme_minimal() +
  theme(legend.position = "none")

fung_inext$iNextEst
str(fung_inext)
colnames(fung_inext$iNextEst)

#This shows the 10 samples with highest completeness.

fung_inext$DataInfo %>%
  select(Assemblage, n, S.obs, SC) %>%
  arrange(desc(SC)) %>%
  head(10)


summary(fung_inext$iNextEst$coverage_based)

mean(fung_inext$DataInfo$SC)
summary(fung_inext$DataInfo$SC)


fung_inext$AsyEst$completeness <- fung_inext$AsyEst$Observed /
  fung_inext$AsyEst$Estimator

mean(fung_inext$AsyEst$completeness)
range(fung_inext$AsyEst$completeness)
sort(fung_inext$AsyEst$completeness)

############################################################
###ASSESS CHANGE THROUGH TIME AND WITH FLORAL DIVERSITY######
###########################################################

# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# BiocManager::install("phyloseq")

library(phyloseq)
packageVersion("phyloseq")

ps.gbp23_fun <- readRDS("data/ps.gbp23.f.RDS")

print(ps.gbp23_fun)
sample_data_tab.cl_fun <- as(sample_data(ps.gbp23_fun), "data.frame")
glimpse(sample_data_tab.cl_fun)

rownames(sample_data_tab.cl_fun)
# Convert row names in df2 to a column
sample_data_tab.cl_fun$Assemblage <- rownames(sample_data_tab.cl_fun)

#join both datasets

fung_inext$AsyEst

# Merge by species_id
merged_df_fung <- merge(
  fung_inext$AsyEst,
  sample_data_tab.cl_fun,
  by = "Assemblage"
)
glimpse(merged_df_fung)

str(merged_df_fung)

merged_df_fung$Diversity <- as.factor(merged_df_fung$Diversity)
merged_df_fung$period <- as.factor(merged_df_fung$period)
merged_df_fung$site <- as.factor(merged_df_fung$site)
library(ggplot2)

str(merged_df_fung)

ggplot(merged_df_fung, aes(x = period, y = Observed)) +
  geom_violin() +
  facet_wrap(~Diversity, scales = "free_y") +
  theme_minimal() +
  labs(x = "Period", y = "Diversity")

ggplot(merged_df_fung, aes(x = site, y = Observed)) +
  geom_violin() +
  facet_wrap(~Diversity, scales = "free_y") +
  theme_minimal() +
  labs(x = "Period", y = "Diversity")


# library(RColorBrewer)
# pal<- brewer.pal(n = 6, name = "Set2")
pal <- c("#3A9AB2", "#91BAB6", "#BDC881", "#E3B710", "#EC7A05", "#F11B00")
p1_fung <- ggplot(
  merged_df_fung,
  aes(x = period, y = Observed, fill = period)
) +
  geom_violin(alpha = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = pal) +
  theme_minimal() +
  labs(x = "Period", y = "Fungal diversity") +
  facet_wrap(~Diversity, scales = "free_y")
p1_fung

# Save the plot
ggsave(
  "results/violin_plot_diversity_fungi_period.png",
  plot = p1_fung,
  width = 8,
  height = 5,
  dpi = 400
)


##keep only shannon diversity

merged_df_fung_sub <- subset(
  merged_df_fung,
  merged_df_fung$Diversity == "Shannon diversity"
)


# First, convert period to numeric for the trend line
merged_df_fung_sub$period_num <- as.numeric(factor(merged_df_fung_sub$period)) # preserves order


p1_sh_fung <- ggplot(
  merged_df_fung_sub,
  aes(x = period, y = Observed, fill = period)
) +
  geom_violin(alpha = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
  geom_smooth(
    data = merged_df_fung_sub,
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
  labs(x = "Period", y = "Fungi Diversity")
p1_sh_fung

# Save the plot
ggsave(
  "results/violin_plot_SHANNON_diversity_fungi_period.png",
  plot = p1_sh_fung,
  width = 4,
  height = 5,
  dpi = 400
)


save(p1_sh_fung, file = "./results/RData/violin_plots_SHANNON_fung_DIv.RData")


ggplot(merged_df_fung, aes(x = site, y = Observed, fill = site)) +
  geom_violin(alpha = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
  theme_minimal() +
  labs(x = "Site", y = "Fungi diversity") +
  facet_wrap(~Diversity, scales = "free_y")


#change to numeric to quantitatively assess changes through time

merged_df_fung$period_num <- as.numeric(sub(
  "^p",
  "",
  as.character(merged_df_fung$period)
))

library(dplyr)
library(broom)

model_results_fung <- merged_df_fung %>%
  group_by(Diversity) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(Observed ~ period_num, data = .x)),
    tidy_model = map(model, tidy)
  ) %>%
  unnest(tidy_model)


model_results_fung %>%
  filter(term == "period_num") %>%
  select(Diversity, estimate, std.error, statistic, p.value)

library(lme4)

m1_fung <- lmer(Observed ~ period_num + (1 | site), data = merged_df_fung_sub)
summary(m1_fung)
library(car)
car::Anova(m1_fung)


m1_fung <- lm(Observed ~ period_num + site, data = merged_df_fung_sub)
summary(m1_fung)
library(car)
car::Anova(m1_fung)


ggplot(merged_df_fung, aes(x = period_num, y = Observed)) +
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


merged_df_fung$site_num <- as.numeric(sub(
  "^s",
  "",
  as.character(merged_df_fung$site)
))


merged_df_fung$j <- paste(merged_df_fung$period_num, merged_df_fung$site_num)


merged_df_fung_fl <- left_join(merged_df_fung, df_fl5, by = "j")

# write.csv(merged_df_fung_fl, "data/merged_fungi_floral_diversity.csv", row.names = FALSE)

merged_df_fung_fl <- read.csv("data/merged_fungi_floral_diversity.csv")

pal <- c("#3A9AB2", "#91BAB6", "#BDC881", "#E3B710", "#EC7A05", "#F11B00")

p4_fung <- ggplot(merged_df_fung_fl, aes(x = n_plants, y = Observed)) +
  geom_point(aes(color = period)) +
  geom_smooth(method = "lm", color = "black") + # single line, not per period
  facet_wrap(~Diversity, scales = "free_y") +
  scale_color_manual(values = pal) +
  xlab("Plant species richness") +
  ylab("Observed fungi diversity/richness") +
  theme_minimal()

p4_fung


ggsave(
  "results/flower_richness_fungi.png",
  plot = p4_fung,
  width = 8,
  height = 5,
  dpi = 400
)


##only shannon div

merged_df_fung_fl_sub <- subset(
  merged_df_fung_fl,
  merged_df_fung_fl$Diversity == "Shannon diversity"
)
p4_sh_fung <- ggplot(merged_df_fung_fl_sub, aes(x = n_plants, y = Observed)) +
  geom_point(aes(color = period)) +
  geom_smooth(method = "lm", color = "black") + # single line, not per period
  scale_color_manual(values = pal) +
  xlab("Plant species richness") +
  ylab("Observed fungi diversity") +
  theme_minimal()

p4_sh_fung


save(p4_sh_fung, file = "./results/RData/scatter_plant_fung_Div.RData")


ggsave(
  "results/flower_richness_SHANNON_fung.png",
  plot = p4_sh_fung,
  width = 4,
  height = 5,
  dpi = 400
)


ggplot(merged_df_fung_fl, aes(x = n_plants, y = Observed)) +
  geom_point(aes(color = site)) +
  geom_smooth(method = "lm", color = "black") + # single line, not per period
  facet_wrap(~Diversity, scales = "free_y") +
  xlab("Plant species richness") +
  ylab("Observed fungi diversity/richness") +
  # scale_color_manual(values = pal) +
  theme_minimal()

ggsave(
  "results/flower_richness_fung.png",
  plot = p4_fung,
  width = 8,
  height = 5,
  dpi = 400
)

library(dplyr)

model_results_fl_fung <- merged_df_fung_fl %>%
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


model_results_fl_fung %>%
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

glimpse(merged_df_fung_fl)

merged_df_fung_fl_sub <- subset(
  merged_df_fung_fl,
  merged_df_fung_fl$Diversity == "Shannon diversity"
)

glimpse(merged_df_fung_fl_sub)

merged_df_fung_fl_sub2 <- merged_df_fung_fl_sub %>%
  group_by(site_num) %>%
  mutate(
    floral_between = mean(n_plants, na.rm = TRUE),
    floral_within = n_plants - floral_between
  ) %>%
  ungroup()


#check for collinearity

m_in_fun <- lm(
  Observed ~ Periodo + floral_within + floral_between,
  data = merged_df_fung_fl_sub2
)


car::vif(m_in_fun) #no collinearity


library(lme4)
m1_fung <- lmer(
  Observed ~
    Periodo + floral_within + floral_between + (1 + Periodo | site_num),
  data = merged_df_fung_fl_sub2
)

summary(m1_fung)
car::Anova(m1_fung)

# If temporal autocorrelation matters or you want AR(1):
library(dplyr)

merged_df_fung_fl_sub2 <- merged_df_fung_fl_sub2 %>%
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

m0_fung <- glmmTMB(
  Observed ~ Periodo_sc + floral_within_sc + floral_between_sc + (1 | site_num),
  data = merged_df_fung_fl_sub2,
  family = gaussian(),
  ziformula = ~0,
  dispformula = ~1,
  REML = TRUE
)
summary(m0_fung) # Should give SEs and no singular warning
m1_fung <- glmmTMB(
  Observed ~
    Periodo_sc +
      floral_within_sc +
      floral_between_sc +
      (1 + Periodo_sc || site_num), # note the double pipes
  data = merged_df_fung_fl_sub2,
  family = gaussian(),
  ziformula = ~0,
  dispformula = ~1,
  REML = TRUE
)
summary(m1_fung)
m2_fung <- glmmTMB(
  Observed ~
    Periodo_sc +
      floral_within_sc +
      floral_between_sc +
      (1 + Periodo_sc | site_num),
  data = merged_df_fung_fl_sub2,
  family = gaussian(),
  ziformula = ~0,
  dispformula = ~1,
  REML = TRUE
)
summary(m2_fung) # If this reintroduces singularity, stick with m1.


# Standardized model already fit as m1
confint(m2_fung) # 95% CI
MuMIn::r.squaredGLMM(m2_fung) # R2m, R2c

# Partial effect plot
library(ggeffects)
plot(ggpredict(m2_fung, terms = "floral_within_sc [all]")) +
  ggplot2::labs(
    x = "Floral richness (within-site, z)",
    y = "Predicted Shannon H"
  )

###add depth of sequencing
merged_df_fung_fl_sub2 <- merged_df_fung_fl_sub2 |>
  dplyr::mutate(
    log_reads = scale(log(quant_reading)) # or whatever your depth column is
  )

m2_depth_fung <- glmmTMB(
  Observed ~
    Periodo_sc +
      floral_within_sc +
      floral_between_sc +
      log_reads +
      (1 + Periodo_sc || site_num),
  data = merged_df_fung_fl_sub2,
  family = gaussian(),
  ziformula = ~0,
  dispformula = ~1,
  REML = TRUE
)
summary(m2_depth_fung)


# Standardized model already fit as m1
confint(m2_depth_fung) # 95% CI
MuMIn::r.squaredGLMM(m2_depth_fung) # R2m, R2c

# Partial effect plot
library(ggeffects)
plot(ggpredict(m2_depth_fung, terms = "floral_within_sc [all]")) +
  ggplot2::labs(
    x = "Floral richness (within-site, z)",
    y = "Predicted Shannon H"
  )


p4_sh_fung2 <- ggplot(
  merged_df_fung_fl_sub2,
  aes(x = floral_within_sc, y = Observed)
) +
  geom_point(aes(color = period)) +
  geom_smooth(method = "lm", color = "black") + # single line, not per period
  scale_color_manual(values = pal) +
  xlab("Plant species richness") +
  ylab("Observed Fungi diversity") +
  theme_minimal()

p4_sh_fung2


ggsave(
  "results/flower_richness_temporal_autocorrelation_fungi.png",
  plot = p4_sh_fung2,
  width = 8,
  height = 5,
  dpi = 400
)

save(
  p4_sh_fung2,
  file = "./results/RData/flower_richness_temporal_autocorrelation_fungi.RData"
)


# merge with actual names of ASVs

df_fung

# ps.gbp23 <- readRDS("data/ps.gbp23_c.0.1.RDS")
ps.gbp23_fung <- readRDS("data/ps.gbp23.f.RDS")

count_tab.cl_fung <- otu_table(ps.gbp23_fung)
count_tab.cl_fung <- as.data.frame(count_tab.cl_fung)
sample_data_tab.cl_fung <- as(sample_data(ps.gbp23_fung), "data.frame")
tax_table.cl_fung <- as.data.frame(tax_table(ps.gbp23_fung))


rownames(tax_table.cl_fung)

tax_table.cl_fung$ASV_code <- rownames(tax_table.cl_fung)
df_fung$ASV_code <- rownames(df_fung)

# Merge by species_id

merged_tax_fung <- merge(df_fung, tax_table.cl_fung, by = "ASV_code")


glimpse(merged_df_fung)
glimpse(merged_tax_fung)
head(merged_tax_fung)

library(dplyr)
library(tidyr)

# Filtrar órdenes no nulos
merged_tax_filtered_fung <- merged_tax_fung %>%
  filter(!is.na(Order))

# Seleccionar columnas de muestras (que comienzan con 'GBP')
sample_columns_fung <- grep(
  "^GBP",
  names(merged_tax_filtered_fung),
  value = TRUE
)

# Calcular el número total de muestras
n_samples_fung <- length(sample_columns_fung)


# Convertir los rownames de sample_data_tab.cl a una columna para hacer el join
sample_data_tab.cl_fung <- sample_data_tab.cl_fung %>%
  tibble::rownames_to_column(var = "Sample")


long_df_fung <- merged_tax_filtered_fung %>%
  pivot_longer(
    cols = all_of(sample_columns_fung),
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  filter(Abundance > 0) %>%
  select(Order, Sample) %>%
  left_join(
    sample_data_tab.cl_fung %>% select(Sample, period, site),
    by = c("Sample" = "Sample")
  ) %>%
  distinct()


# Calcular prevalencia por orden
order_prevalence_fung <- long_df_fung %>%
  group_by(Order, period) %>%
  summarise(n_samples_present = n_distinct(Sample), .groups = "drop") %>%
  mutate(prevalence = (n_samples_present / n_samples) * 100) %>%
  arrange(desc(prevalence))

top_orders_fung <- order_prevalence_fung %>%
  group_by(Order) %>%
  summarise(max_prevalence = max(prevalence)) %>%
  slice_max(max_prevalence, n = 15) %>%
  pull(Order)

order_prevalence_top_fung <- order_prevalence_fung %>%
  filter(Order %in% top_orders_fung)

library(ggplot2)
pal <- c("#3A9AB2", "#91BAB6", "#BDC881", "#E3B710", "#EC7A05", "#F11B00")


plot_order_fung <- ggplot(
  order_prevalence_top_fung,
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
  "results/barchart_order_fungi.png",
  plot = plot_order_fung,
  width = 8,
  height = 5,
  dpi = 400
)

save(
  plot_order_fung,
  file = "./results/RData/barchart_order_fung.RData"
)


# Filtrar familias no nulos
merged_tax_filtered_fung <- merged_tax_fung %>%
  filter(!is.na(Family))

# Seleccionar columnas de muestras (que comienzan con 'GBP')
sample_columns_fung <- grep(
  "^GBP",
  names(merged_tax_filtered_fung),
  value = TRUE
)

# Calcular el número total de muestras
n_samples_fung <- length(sample_columns_fung)


# # Convertir los rownames de sample_data_tab.cl a una columna para hacer el join
# sample_data_tab.cl <- sample_data_tab.cl %>%
#   tibble::rownames_to_column(var = "Sample_ID")

long_df_fung <- merged_tax_filtered_fung %>%
  pivot_longer(
    cols = all_of(sample_columns_fung),
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  filter(Abundance > 0) %>%
  select(Family, Sample) %>%
  left_join(
    sample_data_tab.cl_fung %>% select(Sample, period, site),
    by = c("Sample" = "Sample")
  ) %>%
  distinct()


# Calcular prevalencia por orden
Family_prevalence_fung <- long_df_fung %>%
  group_by(Family, period) %>%
  summarise(n_samples_present = n_distinct(Sample), .groups = "drop") %>%
  mutate(prevalence = (n_samples_present / n_samples) * 100) %>%
  arrange(desc(prevalence))

top_Familys_fung <- Family_prevalence_fung %>%
  group_by(Family) %>%
  summarise(max_prevalence = max(prevalence)) %>%
  slice_max(max_prevalence, n = 15) %>%
  pull(Family)

Family_prevalence_top_fung <- Family_prevalence_fung %>%
  filter(Family %in% top_Familys_fung)

library(ggplot2)
pal <- c("#3A9AB2", "#91BAB6", "#BDC881", "#E3B710", "#EC7A05", "#F11B00")


ggplot(
  Family_prevalence_top_fung,
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
merged_tax_filtered_fung <- merged_tax_fung %>%
  filter(!is.na(Class))

# Seleccionar columnas de muestras (que comienzan con 'GBP')
sample_columns_fung <- grep(
  "^GBP",
  names(merged_tax_filtered_fung),
  value = TRUE
)

# Calcular el número total de muestras
n_samples_fung <- length(sample_columns_fung)


# Convertir los rownames de sample_data_tab.cl a una columna para hacer el join
# sample_data_tab.cl <- sample_data_tab.cl %>%
#   tibble::rownames_to_column(var = "Sample_ID")

long_df_fung <- merged_tax_filtered_fung %>%
  pivot_longer(
    cols = all_of(sample_columns_fung),
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  filter(Abundance > 0) %>%
  select(Class, Sample) %>%
  left_join(
    sample_data_tab.cl_fung %>% select(Sample, period),
    by = c("Sample" = "Sample")
  ) %>%
  distinct()


# Calcular prevalencia por orden
Class_prevalence_fung <- long_df_fung %>%
  group_by(Class, period) %>%
  summarise(n_samples_present = n_distinct(Sample), .groups = "drop") %>%
  mutate(prevalence = (n_samples_present / n_samples) * 100) %>%
  arrange(desc(prevalence))

top_Classs_fung <- Class_prevalence_fung %>%
  group_by(Class) %>%
  summarise(max_prevalence = max(prevalence)) %>%
  slice_max(max_prevalence, n = 15) %>%
  pull(Class)

Class_prevalence_top_fung <- Class_prevalence_fung %>%
  filter(Class %in% top_Classs_fung)

library(ggplot2)
pal <- c("#3A9AB2", "#91BAB6", "#BDC881", "#E3B710", "#EC7A05", "#F11B00")


ggplot(
  Class_prevalence_top_fung,
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
long_df_comp_fung <- merged_tax_filtered_fung %>%
  pivot_longer(
    cols = all_of(sample_columns_fung), # <- your sample columns
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  filter(Abundance > 0) %>%
  left_join(
    sample_data_tab.cl_fung %>% select(Sample, period, site),
    by = "Sample"
  ) %>%
  distinct()

# sum duplicates (Sample × ASV_code)
long_df_comp_sum_fung <- long_df_comp_fung %>%
  group_by(Sample, ASV_code) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# wide community matrix (samples × ASVs)
community_matrix_fung <- long_df_comp_sum_fung %>%
  pivot_wider(
    names_from = ASV_code,
    values_from = Abundance,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("Sample")

# aligned metadata for the same samples
metadata_fung <- sample_data_tab.cl_fung %>%
  filter(Sample %in% rownames(community_matrix_fung)) %>%
  mutate(
    period = factor(trimws(as.character(period))),
    site = factor(trimws(as.character(site)))
  ) %>%
  tibble::column_to_rownames("Sample")

# final alignment check
stopifnot(all(rownames(community_matrix_fung) == rownames(metadata_fung)))

#################################################
## 2) NMDS (as before)
#################################################
set.seed(1)

# square-root transform abundances
community_matrix_sqrt_fung <- sqrt(community_matrix_fung)

# then use this transformed matrix
nmds_result_fung <- metaMDS(
  community_matrix_sqrt_fung,
  distance = "bray",
  k = 2,
  trymax = 100
)


# nmds_result <- metaMDS(community_matrix, distance = "bray", k = 2, trymax = 100)

# NMDS scores + metadata for plotting
nmds_df_fung <- as.data.frame(scores(nmds_result_fung, display = "sites")) %>%
  rownames_to_column("Sample") %>%
  left_join(
    sample_data_tab.cl_fung %>%
      filter(Sample %in% rownames(community_matrix_fung)) %>%
      select(Sample, period, site),
    by = "Sample"
  )


df <- nmds_df_fung

# --- Robust outlier removal (per-axis MAD threshold; tweak k if needed) ---
k <- 6
med1 <- median(df$NMDS1, na.rm = TRUE)
mad1 <- mad(df$NMDS1, na.rm = TRUE)
med2 <- median(df$NMDS2, na.rm = TRUE)
mad2 <- mad(df$NMDS2, na.rm = TRUE)

clean <- df %>%
  mutate(z1 = abs((NMDS1 - med1) / mad1), z2 = abs((NMDS2 - med2) / mad2)) %>%
  filter(z1 <= k, z2 <= k) %>%
  select(-z1, -z2)

# palette (your colors)
pal <- c("#3A9AB2", "#91BAB6", "#BDC881", "#E3B710", "#EC7A05", "#F11B00")

p2_fung <- ggplot(clean, aes(x = NMDS1, y = NMDS2, color = period)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = pal) +
  labs(color = "Period") +
  theme_minimal(base_size = 14)


p2_fung

# convex hulls per period
hulls_fung <- clean %>%
  group_by(period) %>%
  slice(chull(NMDS1, NMDS2))

p2_hull_fung <- p2_fung +
  geom_polygon(
    data = hulls_fung,
    aes(x = NMDS1, y = NMDS2, fill = period, group = period),
    alpha = 0.2,
    color = NA
  ) +
  scale_fill_manual(values = pal) +
  guides(fill = "none")


p2_hull_fung


ggsave(
  "results/nmds_fung.png",
  plot = p2_fung,
  width = 8,
  height = 5,
  dpi = 400
)
save(p2_fung, file = "results/RData/nmds_fung.RData")


####graphs per site

p2_site_fung <- ggplot(clean, aes(NMDS1, NMDS2, color = site)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_viridis_d() + # needs library(viridis)
  labs(color = "Site") +
  theme_minimal(base_size = 14)

# # convex hulls per period
hulls_site <- nmds_df %>%
  group_by(site) %>%
  slice(chull(NMDS1, NMDS2))
#
p2_site_hull <- p2_site +
  geom_polygon(
    data = hulls_site,
    aes(x = NMDS1, y = NMDS2, fill = site, group = site),
    alpha = 0.2,
    color = NA
  ) +
  scale_color_viridis_d() +
  guides(fill = "none")

ggsave(
  "results/nmds_site_fung.png",
  plot = p2_site_fung,
  width = 8,
  height = 5,
  dpi = 400
)
save(p2_site_fung, file = "results/RData/nmds_site_fung.RData")


#################################################
## 3) Global PERMANOVA (period + site)
#################################################
adonis_result_fung <- adonis2(
  community_matrix_sqrt_fung ~ period + site,
  data = metadata_fung,
  method = "bray",
  permutations = 999,
  by = "margin" # marginal tests (like your original)
)
print(adonis_result_fung)

#################################################
## 4) Dispersion check for 'period' (assumption)
#################################################
dist_bc_fung <- vegdist(community_matrix_sqrt_fung, method = "bray")
bd_fung <- betadisper(dist_bc_fung, group = metadata_fung$period)
disp_test_fung <- permutest(bd_fung)
print(disp_test_fung)
# If significant, note as a caveat in the text.
# plot(bd) # optional visualization

#################################################
## 5) Pairwise PERMANOVA between periods (vegan-only)
##    - permutations blocked within site when possible
#################################################

pairwise_permanova_period_fung <- function(
  community_matrix_sqrt_fung,
  metadata_fung,
  period_col = "period",
  site_col = "site",
  method = "bray",
  nperm = 999,
  p_adjust = "BH",
  min_per_group = 2
) {
  stopifnot(all(rownames(community_matrix_fung) == rownames(metadata_fung)))

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
pw_fung <- pairwise_permanova_period_fung(
  community_matrix = community_matrix_sqrt_fung,
  metadata = metadata_fung,
  period_col = "period",
  site_col = "site",
  method = "bray",
  nperm = 999,
  p_adjust = "BH"
)
print(pw_fung)

write.csv(
  pw_fung,
  "results/pairwise_permanova_period_fung.csv",
  row.names = FALSE
)

#################################################
## 6) (Optional) quick heatmap of pairwise p_adj
#################################################
if (nrow(pw_fung) > 0 && !all(is.na(pw_fung$p_adj))) {
  # expand pairs for plotting
  split_pairs_fung <- do.call(rbind, strsplit(pw_fung$pair, " vs "))
  hp_fung <- data.frame(
    A = split_pairs[, 1],
    B = split_pairs[, 2],
    minuslog10padj = -log10(pw_fung$p_adj)
  )
  gg_fung <- ggplot(hp_fung, aes(A, B, fill = minuslog10padj)) +
    geom_tile() +
    scale_fill_viridis_c(name = "-log10(FDR p)") +
    theme_minimal(base_size = 12) +
    labs(x = "Period", y = "Period")

  ggsave(
    "results/pairwise_permanova_heatmap_fung.png",
    gg,
    width = 6,
    height = 5,
    dpi = 300
  )
}


# install.packages("patchwork") # if not yet installed
library(patchwork)


# add panel tags A) and B)
combinedgg_fung <- p2_fung + gg_fung + plot_annotation(tag_levels = 'A')

# save combined figure
ggsave(
  "results/nmds_pairwise_combined_fung.png",
  combinedgg_fung,
  width = 12,
  height = 5,
  dpi = 300
)


save(combinedgg_fung, file = "results/RData/combined_nmds_fung.RData")

######################################
########Plant gut content##############
#######################################

##add data on plant species found in gut contents and relate the diversity
##of plant gut contents with that of bacteria and fungi

# Read the TSV file
df_asv_pl <- read.delim(
  "data/2023_plant_GorBEEa_ASVs_counts.tsv",
  header = TRUE,
  sep = "\t"
)

# Take a quick look at the structure of the data
glimpse(df_asv_pl)


###############################
####SAMPLING COMPLETENESS######
##############################

# Assuming df_asv is already loaded, and ASVs are rows, samples are columns.

# 1. Set ASVs as rownames and remove ASV column (you did this already)
# df_asv_pl$ASV <- df_asv_pl$X
colnames(df_asv_pl)
df_asv_pl <- df_asv_pl[, -(128:135)]
rownames(df_asv_pl) <- df_asv_pl$X
df_asv_pl <- df_asv_pl[, -1]

# 2. Transpose so samples are rows, ASVs are columns
df_asv_pl_t <- t(df_asv_pl)

# 3. Convert each sample row into a numeric vector (not dataframe!)
asv_list_pl <- apply(df_asv_pl_t, 1, function(x) as.numeric(x))

# 4. `apply` returns a matrix for numeric input, so convert to list
asv_list_pl <- split(asv_list_pl, seq(nrow(df_asv_pl_t)))

# But `split` won't work properly here; better do:
asv_list_pl <- lapply(
  1:nrow(df_asv_pl_t),
  function(i) as.numeric(df_asv_pl_t[i, ])
)
names(asv_list_pl) <- rownames(df_asv_pl_t)


# 0) Make sure it's a numeric matrix with ASVs in rows, samples in columns
mat <- as.matrix(df_asv_pl[, -1]) # drop ASV column
rownames(mat) <- df_asv_pl$ASV
storage.mode(mat) <- "numeric"

# 1) Sanitize values
mat[!is.finite(mat)] <- 0 # replace NA/NaN/Inf with 0
mat[mat < 0] <- 0 # no negatives
# If your data are relative abundances (sum≈1), iNEXT expects counts, not proportions.

# 2) Drop empty samples (sum == 0)
keep <- colSums(mat) > 0
if (!all(keep)) {
  message(
    "Dropping ",
    sum(!keep),
    " empty sample(s): ",
    paste(colnames(mat)[!keep], collapse = ", ")
  )
}
mat <- mat[, keep, drop = FALSE]

# 3) Build the list that iNEXT expects (one numeric vector per sample)
asv_list_pl <- lapply(seq_len(ncol(mat)), function(j) as.numeric(mat[, j]))
names(asv_list_pl) <- colnames(mat)

# 4) Run iNEXT
library(iNEXT)
asv_pl_inext <- iNEXT(asv_list_pl, q = 0, datatype = "abundance")

# 5) Plot
library(ggplot2)
ggiNEXT(asv_pl_inext) + theme(legend.position = "none")


# Plot the rarefaction curves
ggiNEXT(asv_pl_inext, type = 1) +
  ggtitle("Rarefaction and Extrapolation Curves") +
  theme_minimal() +
  theme(legend.position = "none")


library(ggplot2)

p_pl <- ggiNEXT(asv_pl_inext, type = 1) +
  ylab("Plant ASV diversity") +
  theme_minimal() +
  theme(legend.position = "none")

# Modify line size in the plot object
p_pl$layers <- lapply(p_pl$layers, function(layer) {
  if ("GeomLine" %in% class(layer$geom)) {
    layer$aes_params$linewidth <- 1
  }
  layer
})

p_pl # now plot with thinner lines


ggsave(
  "results/sampling_completeness_plants.png",
  plot = p_pl,
  width = 4,
  height = 5,
  dpi = 400
)


save(p_pl, file = "./results/RData/sampling_complete_plants.RData")


# Plot the sampling completeness curves
ggiNEXT(asv_pl_inext, type = 2) +
  ylab("ASV diversity") +
  theme_minimal() +
  theme(legend.position = "none")

asv_pl_inext$iNextEst
str(asv_pl_inext)
colnames(asv_pl_inext$iNextEst)

#This shows the 10 samples with highest completeness.

asv_pl_inext$DataInfo %>%
  select(Assemblage, n, S.obs, SC) %>%
  arrange(desc(SC)) %>%
  head(10)


summary(asv_pl_inext$iNextEst$coverage_based)

mean(asv_pl_inext$DataInfo$SC)
summary(asv_pl_inext$DataInfo$SC)


asv_pl_inext$AsyEst$completeness <- asv_pl_inext$AsyEst$Observed /
  asv_pl_inext$AsyEst$Estimator

mean(asv_pl_inext$AsyEst$completeness)
range(asv_pl_inext$AsyEst$completeness)

############################################################
###ASSESS RELATIONSHIP BETWEEN GUT PLANT AND MICROBE DIVERSITY######
###########################################################

# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# BiocManager::install("phyloseq")

library(phyloseq)
packageVersion("phyloseq")

ps.gbp23_pl <- readRDS("data/gbp23.plant.decontam.0.5.RDS")
print(ps.gbp23_pl)
sample_data_tab.cl_pl <- as(sample_data(ps.gbp23_pl), "data.frame")
glimpse(sample_data_tab.cl_pl)

rownames(sample_data_tab.cl_pl)
# Convert row names in df2 to a column
sample_data_tab.cl_pl$Assemblage <- rownames(sample_data_tab.cl_pl)

#join both datasets

asv_pl_inext$AsyEst

# Merge by species_id
merged_df_pl <- merge(
  asv_pl_inext$AsyEst,
  sample_data_tab.cl_pl,
  by = "Assemblage"
)
glimpse(merged_df_pl)

str(merged_df_pl)

merged_df_pl$Diversity <- as.factor(merged_df_pl$Diversity)
merged_df_pl$period <- as.factor(merged_df_pl$period)
merged_df_pl$site <- as.factor(merged_df_pl$site)
library(ggplot2)

str(merged_df_pl)

ggplot(merged_df_pl, aes(x = period, y = Observed)) +
  geom_violin() +
  facet_wrap(~Diversity, scales = "free_y") +
  theme_minimal() +
  labs(x = "Period", y = "Diversity")

ggplot(merged_df_pl, aes(x = site, y = Observed)) +
  geom_violin() +
  facet_wrap(~Diversity, scales = "free_y") +
  theme_minimal() +
  labs(x = "Period", y = "Diversity")


# library(RColorBrewer)
# pal<- brewer.pal(n = 6, name = "Set2")
pal <- c("#3A9AB2", "#91BAB6", "#BDC881", "#E3B710", "#EC7A05", "#F11B00")
p1_pl <- ggplot(merged_df_pl, aes(x = period, y = Observed, fill = period)) +
  geom_violin(alpha = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = pal) +
  theme_minimal() +
  labs(x = "Period", y = "Diversity") +
  facet_wrap(~Diversity, scales = "free_y")
p1_pl

# Save the plot
ggsave(
  "results/violin_plot_diversity_period_plant.png",
  plot = p1_pl,
  width = 8,
  height = 5,
  dpi = 400
)


##keep only shannon diversity

merged_df_sub_pl <- subset(
  merged_df_pl,
  merged_df_pl$Diversity == "Shannon diversity"
)


# First, convert period to numeric for the trend line
merged_df_sub_pl$period_num <- as.numeric(factor(merged_df_sub_pl$period)) # preserves order


p1_sh_pl <- ggplot(
  merged_df_sub_pl,
  aes(x = period, y = Observed, fill = period)
) +
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
  labs(x = "Period", y = "Gut content plant Diversity")
p1_sh_pl

# Save the plot
ggsave(
  "results/violin_plot_SHANNON_diversity_period_plant.png",
  plot = p1_sh_pl,
  width = 4,
  height = 5,
  dpi = 400
)


save(p1_sh_pl, file = "./results/RData/violin_plots_SHANNON_ASV_DIv_pl.RData")


ggplot(merged_df_pl, aes(x = site, y = Observed, fill = site)) +
  geom_violin(alpha = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
  theme_minimal() +
  labs(x = "Site", y = "Diversity") +
  facet_wrap(~Diversity, scales = "free_y")


#change to numeric to quantitatively assess changes through time

merged_df_pl$period_num <- as.numeric(sub(
  "^p",
  "",
  as.character(merged_df_pl$period)
))

library(dplyr)
library(broom)

model_results_pl <- merged_df_pl %>%
  group_by(Diversity) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(Observed ~ period_num, data = .x)),
    tidy_model = map(model, tidy)
  ) %>%
  unnest(tidy_model)


model_results_pl %>%
  filter(term == "period_num") %>%
  select(Diversity, estimate, std.error, statistic, p.value)

library(lme4)

m1_pl <- lmer(Observed ~ period_num + (1 | site), data = merged_df_sub_pl)
summary(m1_pl)
library(car)
car::Anova(m1_pl)


m1_pl <- lm(Observed ~ period_num + site, data = merged_df_sub_pl)
summary(m1_pl)
library(car)
car::Anova(m1_pl)


ggplot(merged_df_pl, aes(x = period_num, y = Observed)) +
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


merged_df_pl$site_num <- as.numeric(sub(
  "^s",
  "",
  as.character(merged_df_pl$site)
))


merged_df_pl$j <- paste(merged_df_pl$period_num, merged_df_pl$site_num)


merged_df_fl_pl <- left_join(merged_df_pl, df_fl5, by = "j")

#write.csv(merged_df_fl_pl, "data/merged_ASV_floral_diversity_gut_plants.csv", row.names = FALSE)

merged_df_fl_pl <- read.csv("data/merged_ASV_floral_diversity_gut_plants.csv")

pal <- c("#3A9AB2", "#91BAB6", "#BDC881", "#E3B710", "#EC7A05", "#F11B00")

p4_pl <- ggplot(merged_df_fl_pl, aes(x = n_plants, y = Observed)) +
  geom_point(aes(color = period)) +
  geom_smooth(method = "lm", color = "black") + # single line, not per period
  facet_wrap(~Diversity, scales = "free_y") +
  scale_color_manual(values = pal) +
  xlab("Plant species richness") +
  ylab("Observed gut plant richness") +
  theme_minimal()

p4_pl


ggsave(
  "results/flower_richness_gut_plants.png",
  plot = p4_pl,
  width = 8,
  height = 5,
  dpi = 400
)


##only shannon div

merged_df_fl_sub_pl <- subset(
  merged_df_fl_pl,
  merged_df_fl_pl$Diversity == "Shannon diversity"
)
p4_sh_pl <- ggplot(merged_df_fl_sub_pl, aes(x = n_plants, y = Observed)) +
  geom_point(aes(color = period)) +
  geom_smooth(method = "lm", color = "black") + # single line, not per period
  scale_color_manual(values = pal) +
  xlab("Plant species richness") +
  ylab("Observed gut plant diversity") +
  theme_minimal()

p4_sh_pl


save(p4_sh_pl, file = "./results/RData/scatter_plant_ASV_Div_gut_plants.RData")


ggsave(
  "results/flower_richness_SHANNON_gut_plants.png",
  plot = p4_sh_pl,
  width = 4,
  height = 5,
  dpi = 400
)


###hasta aqui!!

ggplot(merged_df_fl_pl, aes(x = n_plants, y = Observed)) +
  geom_point(aes(color = site)) +
  geom_smooth(method = "lm", color = "black") + # single line, not per period
  facet_wrap(~Diversity, scales = "free_y") +
  xlab("Plant species richness") +
  ylab("Observed gut plant species richness") +
  # scale_color_manual(values = pal) +
  theme_minimal()


library(dplyr)

model_results_fl_pl <- merged_df_fl_pl %>%
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


model_results_fl_pl %>%
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

glimpse(merged_df_fl_pl)

merged_df_fl_sub_pl <- subset(
  merged_df_fl_pl,
  merged_df_fl_pl$Diversity == "Shannon diversity"
)

glimpse(merged_df_fl_sub_pl)

merged_df_fl_sub_pl2 <- merged_df_fl_sub_pl %>%
  group_by(site_num) %>%
  mutate(
    floral_between = mean(n_plants, na.rm = TRUE),
    floral_within = n_plants - floral_between
  ) %>%
  ungroup()


#check for collinearity

m_in_pl <- lm(
  Observed ~ Periodo + floral_within + floral_between,
  data = merged_df_fl_sub_pl2
)


car::vif(m_in_pl) #no collinearity


library(lme4)
m1_pl <- lmer(
  Observed ~
    Periodo + floral_within + floral_between + (1 + Periodo | site_num),
  data = merged_df_fl_sub_pl2
)

summary(m1_pl)
car::Anova(m1_pl)

# If temporal autocorrelation matters or you want AR(1):
library(dplyr)

merged_df_fl_sub_pl2 <- merged_df_fl_sub_pl2 %>%
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

m0_pl <- glmmTMB(
  Observed ~ Periodo_sc + floral_within_sc + floral_between_sc + (1 | site_num),
  data = merged_df_fl_sub_pl2,
  family = gaussian(),
  ziformula = ~0,
  dispformula = ~1,
  REML = TRUE
)
summary(m0_pl) # Should give SEs and no singular warning
m1_pl <- glmmTMB(
  Observed ~
    Periodo_sc +
      floral_within_sc +
      floral_between_sc +
      (1 + Periodo_sc || site_num), # note the double pipes
  data = merged_df_fl_sub_pl2,
  family = gaussian(),
  ziformula = ~0,
  dispformula = ~1,
  REML = TRUE
)
summary(m1_pl)
m2_pl <- glmmTMB(
  Observed ~
    Periodo_sc +
      floral_within_sc +
      floral_between_sc +
      (1 + Periodo_sc | site_num),
  data = merged_df_fl_sub_pl2,
  family = gaussian(),
  ziformula = ~0,
  dispformula = ~1,
  REML = TRUE
)
summary(m2_pl) # If this reintroduces singularity, stick with m1.


# Standardized model already fit as m1
confint(m2_pl) # 95% CI
MuMIn::r.squaredGLMM(m2_pl) # R2m, R2c

# Partial effect plot
library(ggeffects)
plot(ggpredict(m2_pl, terms = "floral_within_sc [all]")) +
  ggplot2::labs(
    x = "Floral richness (within-site, z)",
    y = "Predicted Shannon H"
  )

###add depth of sequencing
merged_df_fl_sub_pl2 <- merged_df_fl_sub_pl2 |>
  dplyr::mutate(
    log_reads = scale(log(quant_reading)) # or whatever your depth column is
  )

m2_depth_pl <- glmmTMB(
  Observed ~
    Periodo_sc +
      floral_within_sc +
      floral_between_sc +
      log_reads +
      (1 + Periodo_sc || site_num),
  data = merged_df_fl_sub_pl2,
  family = gaussian(),
  ziformula = ~0,
  dispformula = ~1,
  REML = TRUE
)
summary(m2_depth_pl) #this is the model we keep


# Standardized model already fit as m1
confint(m2_depth_pl) # 95% CI
MuMIn::r.squaredGLMM(m2_depth_pl) # R2m, R2c

# Partial effect plot
library(ggeffects)
plot(ggpredict(m2_depth_pl, terms = "floral_within_sc [all]")) +
  ggplot2::labs(
    x = "Floral richness (within-site, z)",
    y = "Predicted Shannon H"
  )


p4_sh_pl2 <- ggplot(
  merged_df_fl_sub_pl2,
  aes(x = floral_within_sc, y = Observed)
) +
  geom_point(aes(color = period)) +
  geom_smooth(method = "lm", color = "black") + # single line, not per period
  scale_color_manual(values = pal) +
  xlab("Plant species richness") +
  ylab("Observed gut content plant species richness") +
  theme_minimal()

p4_sh_pl2


ggsave(
  "results/flower_richness_temporal_autocorrelation_gut plants.png",
  plot = p4_sh_pl2,
  width = 8,
  height = 5,
  dpi = 400
)

save(
  p4_sh_pl2,
  file = "./results/RData/flower_richness_temporal_autocorrelation_gut plants.RData"
)


# merge with actual names of ASVs

df_asv_pl

ps.gbp23_pl <- readRDS("data/gbp23.plant.decontam.0.5.RDS")
count_tab.cl_pl <- otu_table(ps.gbp23_pl)
count_tab.cl_pl <- as.data.frame(count_tab.cl_pl)
sample_data_tab.cl_pl <- as(sample_data(ps.gbp23_pl), "data.frame")
tax_pln <- as.data.frame(tax_table(ps.gbp23_pl))
tax_pln$ASV_code <- row.names(tax_pln)


rownames(tax_pln)

tax_pln$ASV_code <- rownames(tax_pln)
df_asv_pl$ASV_code <- rownames(df_asv_pl)

# Merge by species_id

merged_tax_pl <- merge(df_asv_pl, tax_pln, by = "ASV_code")


glimpse(merged_tax_pl)
head(merged_tax_pl)

library(dplyr)
library(tidyr)

# Filtrar órdenes no nulos
merged_tax_pl_filtered <- merged_tax_pl %>%
  filter(!is.na(Order))

# Seleccionar columnas de muestras (que comienzan con 'GBP')
sample_columns_pl <- grep("^GBP", names(merged_tax_pl_filtered), value = TRUE)


#####AQUI#####

# Calcular el número total de muestras
n_samples_pl <- length(sample_columns_pl)


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


#############################
####ASSESS RELATIONSHIP BETWEEN GUT PLANT DIVERSITY AND BACTERIAL
####AND FUNGAL DIVERSITY

merged_df_fl_sub_pl
library(dplyr)

core <- function(x) sub("_.*$", "", as.character(x))

# Right table: collapse to one row per Core and name the value as gut_plants
lk <- merged_df_fl_sub_pl2 %>%
  mutate(Core = core(Assemblage)) %>%
  group_by(Core) %>%
  summarise(
    gut_plants = dplyr::first(Observed[!is.na(Observed)]),
    .groups = "drop"
  )

# Left join and you’re done
merged_df_fl_sub <- merged_df_fl_sub %>%
  mutate(Core = core(Assemblage)) %>%
  left_join(lk, by = "Core")


ggplot(
  merged_df_fl_sub,
  aes(x = gut_plants.y, y = Observed)
) +
  geom_point(aes(color = period)) +
  geom_smooth(method = "lm", color = "black") + # single line, not per period
  scale_color_manual(values = pal) +
  xlab("Gut Plant species richness") +
  ylab("Observed bacterial ASV richness") +
  theme_minimal()


m2_gut_pl <- glmmTMB(
  Observed ~
    period +
      gut_plants.y +
      (1 + period || site),
  data = merged_df_fl_sub,
  family = gaussian(),
  ziformula = ~0,
  dispformula = ~1,
  REML = TRUE
)
summary(m2_gut_pl) #this is the model we keep


###fungi

merged_df_fung_fl_sub
library(dplyr)

core <- function(x) sub("_.*$", "", as.character(x))

# Right table: collapse to one row per Core and name the value as gut_plants
lk <- merged_df_fl_sub_pl2 %>%
  mutate(Core = core(Assemblage)) %>%
  group_by(Core) %>%
  summarise(
    gut_plants = dplyr::first(Observed[!is.na(Observed)]),
    .groups = "drop"
  )

# Left join and you’re done
merged_df_fung_fl_sub <- merged_df_fung_fl_sub %>%
  mutate(Core = core(Assemblage)) %>%
  left_join(lk, by = "Core")


head(merged_df_fung_fl_sub)

ggplot(
  merged_df_fung_fl_sub,
  aes(x = gut_plants, y = Observed)
) +
  geom_point(aes(color = period)) +
  geom_smooth(method = "lm", color = "black", se = FALSE) + # single line, not per period
  scale_color_manual(values = pal) +
  xlab("Gut Plant species richness") +
  ylab("Observed fungal diversity") +
  theme_minimal()


m2_gut_pl_fun <- glmmTMB(
  Observed ~
    period +
      gut_plants +
      (1 + period || site),
  data = merged_df_fung_fl_sub,
  family = gaussian(),
  ziformula = ~0,
  dispformula = ~1,
  REML = TRUE
)
summary(m2_gut_pl_fun) #this is the model we keep


##SAVE DATA FOR BETADIVERSITY CODE
save(
  community_matrix,
  metadata,
  community_matrix_fung,
  metadata_fung,
  file = "./results/RData/data_betadiv.RData"
)


##SAVE DATA FOR CO-OCCURRENCE CODE
save(
  df_asv_t,
  df_fung_t,
  df_asv_pl_t,
  tax_table.cl,
  tax_table.cl_fung,
  file = "./results/RData/data_co-occurrence.RData"
)
