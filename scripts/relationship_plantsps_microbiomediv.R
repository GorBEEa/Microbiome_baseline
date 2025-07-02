###relate plant species richness to microbiome richness in 2023


#packages
pacman::p_load(rgbif, assertr, dplyr, tidyverse, ggplot2, lme4, car) #instala y lee paquetes a la vez; sustituye a Install.packages() y library()


df<-read.csv("data/Datos_transecto_2023.csv")

## Check taxonomy
df2 <- as.data.frame(rgbif::name_backbone_checklist(df$Planta) |>
  assert(in_set("EXACT", "FUZZY", "HIGHERRANK", "NONE"), matchType))

none<-subset(df2, df2$matchType=="NONE")

#clean those that can be cleaned and re-run
df$Planta <- dplyr::recode(df$Planta,
                                "Hutchinisia alpina" = "Hutchinsia alpina",
                           "Hutchsinsia alpina" =  "Hutchinsia alpina",
                                 "Hippocrepis nomosa" = "Hippocrepis comosa",
                           "Anthillys vulneraria" = "Anthyllis vulneraria",
                           "Ranuculus sp." = "Ranunculus sp.",
                           "Thymus sp." = "Thymus praecox",
                           "Thymus sp. " = "Thymus praecox",
                           "Thymus sp. 2"= "Thymus praecox",
                           "Myositis sp." = "Myosotis sp.",
                           "Thesium sp." = "Thesium pyrenaicum",
                           "Globularia sp." = "Globularia vulgaris",
                           "Carduus sp. cf"= "Cardus",
                           "Carduus"= "Cardus",
                           "Carduus sp."= "Cardus",
                           "Campanula sp." = "Campanula",
                           "Falso Taraxacum sp."= "Taraxacum",
                           "Cordo sp." ="Cardus",
                           "Fake manzanita" = "Chamaemelum nobile",
                           "Manzanilla falso" = "Chamaemelum nobile",
                           "Manzanita falso" = "Chamaemelum nobile",
                           "Manzanilla sp."= "Chamaemelum nobile",
                           "Manzania sp."= "Chamaemelum nobile",
                           "Manzanita sp."= "Chamaemelum nobile",
                           "Fake manzania" = "Chamaemelum nobile",
                           "Manzanilla" = "Chamaemelum nobile",
                           "Manzalnilla" = "Chamaemelum nobile")


df3<-cbind(df, df2)


df4<-subset(df3, df3$matchType=="EXACT"| df3$matchType=="FUZZY"| df3$matchType=="HIGHERRANK")
head(df4)


#get plant richness per site and period
df5 <- df4 %>%
  group_by(Periodo, Sitio) %>%
  summarise(n_plants = n_distinct(species))
head(df5)

#clean microbiome data

df_micro<-read.csv("data/microbiome/GBP23_ASVs.csv")
head(df_micro)
str(df_micro)
# Keep the first 10 columns
base_data <- df_micro[, 1:10]

# Count the number of columns > 0 for each row in the remaining columns
presence_count <- rowSums(df_micro[, 11:ncol(df_micro)] > 0)

# Add the count as a new column
base_data$presence_count <- presence_count

# View the transformed data
head(base_data)

# Update 'period' and 'site' columns to keep only numbers
base_data$period <- as.numeric(gsub("\\D", "", base_data$period)) # Remove non-digits
base_data$site <- as.numeric(gsub("\\D", "", base_data$site))     # Remove non-digits

print(base_data)


df5$j<-paste(df5$Periodo, df5$Sitio)
base_data$j<-paste(base_data$period, base_data$site)


base_data$plant_sps<-df5$n_plants[match(base_data$j, df5$j)]


ggplot(base_data, aes(x = plant_sps, y = presence_count, color = factor(site), group = site)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(group = site)) 
# Add a line of fit without confidence interval shading

ggplot(base_data, aes(x = plant_sps, y = presence_count)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + facet_grid(~period, scales="free")
# Add a line of fit without confidence interval shading


m1<-lmer(presence_count~plant_sps + (1|site), data=base_data)
summary(m1)
car::Anova(m1)


#keep only sites with at least measures for 3 periods of time

# Filter the dataset
filtered_data <- base_data %>%
  group_by(site) %>%                      # Group by site
  filter(n_distinct(period) >= 3) %>%     # Keep sites with at least 3 unique periods
  ungroup()                               # Ungroup the data

# View the filtered dataset
filtered_data

ggplot(filtered_data, aes(x = plant_sps, y = presence_count)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + facet_grid(~period, scales="free") +
  labs(
    x = "Plant species richness",  # Change x-axis label
    y = "Bacteria species richness") + theme_classic() +
  theme(
    axis.title = element_text(size = 14),  # Increase axis title font size
    axis.text = element_text(size = 12),   # Increase axis text font size
    legend.title = element_text(size = 14),  # Increase legend title font size
    legend.text = element_text(size = 12),   # Increase legend text font siz
  ) + facet_grid(~period, scales= "free")




m1<-lm(presence_count~plant_sps , data=filtered_data)
summary(m1)








#use data from flower transects

df_fl<-read.csv("data/Floral_resources_2023.csv")
print(df_fl)
## Check taxonomy
df_fl2 <- as.data.frame(rgbif::name_backbone_checklist(df_fl$Planta) |>
                          assert(in_set("EXACT", "FUZZY", "HIGHERRANK", "NONE"), matchType))

none<-subset(df_fl2, df_fl2$matchType=="NONE")

#clean those that can be cleaned and re-run
df_fl$Planta <- dplyr::recode(df_fl$Planta,
                              "Erica sp." = "Erica",
                              "Hutchinisia alpina" = "Hutchinsia alpina",
                              "Ericea sp." = "Erica",
                              "Trifolium platino"="Trifolium pratense",
                              "Globulana sp." = "Globularia vulgaris",
                              "Cragaetus sp."= "Crataegus monogyna",
                              "Thymus sp." = "Thymus praecox",
                              "Cerasium sp." = "Cerastium fontanum")


df_fl3<-cbind(df_fl, df_fl2)


df_fl4<-subset(df_fl3, df_fl3$matchType=="EXACT"| df_fl3$matchType=="FUZZY"| df_fl3$matchType=="HIGHERRANK")
head(df_fl4)


#get plant richness per site and period
df_fl5 <- df_fl4 %>%
  group_by(Periodo, Sitio) %>%
  summarise(n_plants = n_distinct(species))
head(df_fl5)



df_fl5$j<-paste(df_fl5$Periodo, df_fl5$Sitio)
base_data$j<-paste(base_data$period, base_data$site)


base_data$plant_sps_flor_res<-df_fl5$n_plants[match(base_data$j, df_fl5$j)]


ggplot(base_data, aes(x = plant_sps_flor_res, y = presence_count, color = factor(site), group = site)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(group = site)) 
# Add a line of fit without confidence interval shading

ggplot(base_data, aes(x = plant_sps_flor_res, y = presence_count)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + facet_grid(~site)
# Add a line of fit without confidence interval shading


m1<-lmer(presence_count~plant_sps_flor_res + (1|site), data=base_data)
summary(m1)
car::Anova(m1)


#keep only sites with at least measures for 3 periods of time

# Filter the dataset
filtered_data <- base_data %>%
  group_by(site) %>%                      # Group by site
  filter(n_distinct(period) >= 2) %>%     # Keep sites with at least 3 unique periods
  ungroup()                               # Ungroup the data

# View the filtered dataset
filtered_data

ggplot(filtered_data, aes(x = plant_sps_flor_res, y = presence_count)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) 
# Add a line of fit without confidence interval shading


m1<-lmer(presence_count~plant_sps_flor_res + period + (1|site), data=filtered_data)
summary(m1)
car::Anova(m1)





#stack data to calculate richness per order

# Step 1: Keep columns 1 to 10
df_micro_fixed <- df_micro[, 1:10]

# Stack columns from the 11th onward
df_micro_stacked <- df_micro[, 11:ncol(df_micro)]

# Combine the fixed columns with the stacked data
df_micro_long <- df_micro_stacked %>%
  pivot_longer(
    cols = everything(), # Pivot all remaining columns
    names_to = "X",  # New column for the original column names
    values_to = "value"  # New column for the values
  ) %>%
  bind_cols(df_micro_fixed[rep(1:nrow(df_micro_fixed), each = ncol(df_micro_stacked)), ]) # Repeat the first 10 columns for each stacked row

head(df_micro_long)

#combine with second dataset that has taxomonic information for each OTU

df_micro_tax<-read.csv("data/microbiome/GPB23_tax_table_v2.csv")
head(df_micro_tax)
str(df_micro_tax)

merged_df_micro <- merge(df_micro_long, df_micro_tax, by = "X", all.x = TRUE)

head(merged_df_micro)
summary(merged_df_micro)
str(merged_df_micro)


#remove Order Chloroplast might be plant material

merged_df_micro2<-subset(merged_df_micro, merged_df_micro$Order!="Chloroplast")
str(merged_df_micro2)

# Filter rows where quant_reading > 0 and count occurrences by Order
summary_by_order <- merged_df_micro2 %>%
  filter(value > 0) %>%  # Keep only rows where value > 0
  group_by(Row.names, year, period, site, Order) %>%  # Group by Row.names, year, period, site, and Order
  summarise(order_count = n(), .groups = "drop") %>%  # Count occurrences of each Order
  ungroup() %>%  # Remove grouping to prepare for stacking
  pivot_longer(cols = order_count,  # Pivot the order_count column
               names_to = "order",   # Keep a column for the order names
               values_to = "value")  # Keep a column for the counts




# View the summarized dataset
print(summary_by_order)
str(summary_by_order)
summary_by_order2<-as.data.frame(summary_by_order)
str(summary_by_order2)

# Update 'period' and 'site' columns to keep only numbers
summary_by_order2$period <- as.numeric(gsub("\\D", "", summary_by_order2$period)) # Remove non-digits
summary_by_order2$site <- as.numeric(gsub("\\D", "", summary_by_order2$site))     # Remove non-digits




df_fl5$j<-paste(df_fl5$Periodo, df_fl5$Sitio)
summary_by_order2$j<-paste(summary_by_order2$period, summary_by_order2$site)


summary_by_order2$plant_sps_flor_res<-df_fl5$n_plants[match(summary_by_order2$j, df_fl5$j)]
summary_by_order2$plant_sps<-df5$n_plants[match(summary_by_order2$j, df5$j)]


str(summary_by_order2)

filtered_orders <- summary_by_order2 %>%
  group_by(Order) %>%
  summarise(n_specimens = n_distinct(Row.names), .groups = "drop") %>%
  filter(n_specimens >= 60)  # Keep orders with at least 10 distinct Row.names

# Now, filter your original dataset to keep only those orders
summary_by_order2_filtered <- summary_by_order2 %>%
  filter(Order %in% filtered_orders$Order)

# View the result
head(summary_by_order2_filtered)






ggplot(summary_by_order2_filtered, aes(x = plant_sps_flor_res, y = value, color = factor(Order), group = Order)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(group = Order)) +
  theme(legend.position = "none")  # This removes the legend




ggplot(summary_by_order2_filtered, aes(x = plant_sps, y = value, color = factor(Order), group = Order)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(group = Order)) 

library(ggplot2)

# Manually define 6 color values
colors <- c("#DD8D29", "#E2D200" ,"#46ACC8" ,"#E58601", "#B40F20","#29211F")


ggplot(summary_by_order2_filtered, aes(x = plant_sps_flor_res, y = value, color = factor(Order), group = Order)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(group = Order)) +
  scale_color_manual(values = colors) + # Manually defined 6 colors
  labs(
    x = "Plant species richness",  # Change x-axis label
    y = "Bacteria species richness",               # Change y-axis label
    color = "Bacteria Order"       # Change legend title
  ) + theme_classic() +
  theme(
    axis.title = element_text(size = 14),  # Increase axis title font size
    axis.text = element_text(size = 12),   # Increase axis text font size
    legend.title = element_text(size = 14),  # Increase legend title font size
    legend.text = element_text(size = 12),   # Increase legend text font siz
  ) + facet_grid(~period, scales= "free")


str(summary_by_order2_filtered)

m1<-lmer(value~plant_sps *  Order+ (1|site) , data=summary_by_order2_filtered)
summary(m1)
car::Anova(m1)

m1<-lm(value~plant_sps + period, data=summary_by_order2_filtered)
summary(m1)


m1<-lmer(value~plant_sps_flor_res + (1|Order), data=summary_by_order2_filtered)
summary(m1)
car::Anova(m1)

m1<-lm(value~plant_sps_flor_res , data=summary_by_order2_filtered)
summary(m1)

##look at rickettsia sps

rick<-subset(merged_df_micro2, merged_df_micro2$Order=="Rickettsiales")

#keep only sites with at least measures for 3 periods of time

# Filter the dataset
filtered_data <- summary_by_order2 %>%
  group_by(site) %>%                      # Group by site
  filter(n_distinct(period) >= 2) %>%     # Keep sites with at least 3 unique periods
  ungroup()                               # Ungroup the data

# View the filtered dataset
filtered_data<-droplevels(filtered_data)


ggplot(filtered_data, aes(x = plant_sps, y = value, color = factor(Order), group = Order)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(group = Order)) +
  theme(legend.position = "none")  # This removes the legend

m1<-lmer(value~plant_sps_flor_res + period + (1|site), data=filtered_data)
summary(m1)
car::Anova(m1)



##same but using Family

# Filter rows where value > 0 and count occurrences by Order
summary_by_fam <- merged_df_micro2 %>%
  filter(value > 0) %>%  # Keep only rows where value > 0
  group_by(Row.names, year, period, site, Family) %>%  # Group by Row.names, year, period, site, and Order
  summarise(family_count = n(), .groups = "drop") %>%  # Count occurrences of each Order
  ungroup() %>%  # Remove grouping to prepare for stacking
  pivot_longer(cols = family_count,  # Pivot the order_count column
               names_to = "family",   # Keep a column for the order names
               values_to = "value")  # Keep a column for the counts




# View the summarized dataset
print(summary_by_fam)

summary_by_fam2<-as.data.frame(summary_by_fam)


# Update 'period' and 'site' columns to keep only numbers
summary_by_fam2$period <- as.numeric(gsub("\\D", "", summary_by_fam2$period)) # Remove non-digits
summary_by_fam2$site <- as.numeric(gsub("\\D", "", summary_by_fam2$site))     # Remove non-digits




df_fl5$j<-paste(df_fl5$Periodo, df_fl5$Sitio)
summary_by_fam2$j<-paste(summary_by_fam2$period, summary_by_fam2$site)


summary_by_fam2$plant_sps_flor_res<-df_fl5$n_plants[match(summary_by_fam2$j, df_fl5$j)]
summary_by_fam2$plant_sps<-df5$n_plants[match(summary_by_fam2$j, df5$j)]


str(summary_by_fam2)

summary_by_fam3 <- summary_by_fam2 %>%
  drop_na()



ggplot(summary_by_fam3, aes(x = plant_sps_flor_res, y = value, color = factor(Family), group = Family)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(group = Family)) +
  theme(legend.position = "none")  # This removes the legend




ggplot(summary_by_fam3, aes(x = plant_sps, y = value, color = factor(Family), group = Family)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(group = Family)) +
  theme(legend.position = "none")  # This removes the legend




m1<-lmer(value~plant_sps_flor_res + (1|site), data=summary_by_fam3)
summary(m1)
car::Anova(m1)


filtered_fams <- summary_by_fam3 %>%
  group_by(Family) %>%
  summarise(n_specimens = n_distinct(Row.names), .groups = "drop") %>%
  filter(n_specimens >= 60)  # Keep fams with at least 10 distinct Row.names

# Now, filter your original dataset to keep only those fams
summary_by_fam2_filtered <- summary_by_fam2 %>%
  filter(Family %in% filtered_fams$Family)

# View the result
head(summary_by_fam2_filtered)




ggplot(summary_by_fam2_filtered, aes(x = plant_sps_flor_res, y = value, color = factor(Family), group = Family)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(group = Family)) +
  theme(legend.position = "none")  # This removes the legend




ggplot(summary_by_fam2_filtered, aes(x = plant_sps, y = value, color = factor(Family), group = Family)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(group = Family)) +
  theme(legend.position = "none")  # This removes the legend




m1<-lm(value~plant_sps_flor_res , data=summary_by_fam2_filtered)
summary(m1)
car::Anova(m1)

m1<-lm(value~plant_sps , data=summary_by_fam2_filtered)
summary(m1)
car::Anova(m1)




#keep only sites with at least measures for 3 periods of time

# Filter the dataset
filtered_data <- summary_by_fam3 %>%
  group_by(site) %>%                      # Group by site
  filter(n_distinct(period) >= 2) %>%     # Keep sites with at least 3 unique periods
  ungroup()                               # Ungroup the data

# View the filtered dataset
filtered_data


ggplot(filtered_data, aes(x = plant_sps, y = value, color = factor(Family), group = Family)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(group = Family)) +
  theme(legend.position = "none")  # This removes the legend

ggplot(filtered_data, aes(x = plant_sps_flor_res, y = value, color = factor(Family), group = Family)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(group = Family)) +
  theme(legend.position = "none")  # This removes the legend

m1<-lmer(value~plant_sps_flor_res + period + (1|site), data=filtered_data)
summary(m1)
car::Anova(m1)





##same but using genus

# Filter rows where value > 0 and count occurrences by Order
summary_by_gen <- merged_df_micro2 %>%
  filter(value > 0) %>%  # Keep only rows where value > 0
  group_by(Row.names, year, period, site, Genus) %>%  # Group by Row.names, year, period, site, and Order
  summarise(genus_count = n(), .groups = "drop") %>%  # Count occurrences of each Order
  ungroup() %>%  # Remove grouping to prepare for stacking
  pivot_longer(cols = genus_count,  # Pivot the order_count column
               names_to = "genus",   # Keep a column for the order names
               values_to = "value")  # Keep a column for the counts




# View the summarized dataset
print(summary_by_gen)

summary_by_gen2<-as.data.frame(summary_by_gen)


# Update 'period' and 'site' columns to keep only numbers
summary_by_gen2$period <- as.numeric(gsub("\\D", "", summary_by_gen2$period)) # Remove non-digits
summary_by_gen2$site <- as.numeric(gsub("\\D", "", summary_by_gen2$site))     # Remove non-digits




df_fl5$j<-paste(df_fl5$Periodo, df_fl5$Sitio)
summary_by_gen2$j<-paste(summary_by_gen2$period, summary_by_gen2$site)


summary_by_gen2$plant_sps_flor_res<-df_fl5$n_plants[match(summary_by_gen2$j, df_fl5$j)]
summary_by_gen2$plant_sps<-df5$n_plants[match(summary_by_gen2$j, df5$j)]


str(summary_by_gen2)

summary_by_gen3 <- summary_by_gen2 %>%
  drop_na()

summary_by_gen3<-droplevels(summary_by_gen3)

ggplot(summary_by_gen3, aes(x = plant_sps_flor_res, y = value, color = factor(Genus), group = Genus)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(group = Genus)) +
  theme(legend.position = "none")  # This removes the legend




ggplot(summary_by_gen3, aes(x = plant_sps, y = value, color = factor(Genus), group = Genus)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(group = Genus)) +
  theme(legend.position = "none")  # This removes the legend




m1<-lmer(value~plant_sps_flor_res + (1|site), data=summary_by_gen3)
summary(m1)
car::Anova(m1)


filtered_gens <- summary_by_gen3 %>%
  group_by(Genus) %>%
  summarise(n_specimens = n_distinct(Row.names), .groups = "drop") %>%
  filter(n_specimens >= 30)  # Keep fams with at least 10 distinct Row.names

# Now, filter your original dataset to keep only those fams
summary_by_gen2_filtered <- summary_by_gen3 %>%
  filter(Genus %in% filtered_gens$Genus)

# View the result
head(summary_by_gen2_filtered)




ggplot(summary_by_gen2_filtered, aes(x = plant_sps_flor_res, y = value, color = factor(Genus), group = Genus)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(group = Genus)) +
  theme(legend.position = "none") +  # This removes the legend 
facet_grid(~period)



ggplot(summary_by_gen2_filtered, aes(x = plant_sps, y = value, color = factor(Genus), group = Genus)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(group = Genus)) +
  theme(legend.position = "none")  # This removes the legend




m1<-lm(value~plant_sps_flor_res , data=summary_by_gen2_filtered)
summary(m1)
car::Anova(m1)

m1<-lm(value~plant_sps , data=summary_by_gen2_filtered)
summary(m1)
car::Anova(m1)




#keep only sites with at least measures for 3 periods of time

# Filter the dataset
filtered_data <- summary_by_gen3 %>%
  group_by(site) %>%                      # Group by site
  filter(n_distinct(period) >= 3) %>%     # Keep sites with at least 3 unique periods
  ungroup()                               # Ungroup the data

# View the filtered dataset
filtered_data


ggplot(filtered_data, aes(x = plant_sps, y = value, color = factor(Genus), group = Genus)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(group = Genus)) +
  theme(legend.position = "none")  # This removes the legend

ggplot(filtered_data, aes(x = plant_sps_flor_res, y = value, color = factor(Genus), group = Genus)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(group = Genus)) +
  theme(legend.position = "none")  # This removes the legend

m1<-lmer(value~plant_sps_flor_res + period + (1|site), data=filtered_data)
summary(m1)
car::Anova(m1)


m1<-lm(value~plant_sps_flor_res , data=filtered_data)
summary(m1)
car::Anova(m1)



