load("./results/RData/data_betadiv.RData")

# Install once if needed
# install.packages(c("betapart", "vegan", "dplyr", "tidyr", "readr", "ggplot2", "geodist"))

library(betapart)
library(vegan)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(geodist)

# ==== 1) Load your data ====
# Replace with your file paths / data frames
comm <- community_matrix # rows = samples, cols = ASVs, first col = sample_id
comm$sample_id <- rownames(comm)
meta <- metadata # columns: sample_id, lat, lon, [habitat, date, ...]
meta$sample_id <- rownames(meta)
##add lat and lon

coor <- read.csv("data/coor.csv")
head(coor)
head(meta)

# make sure the column exists
stopifnot("site" %in% names(coor))

# normalize to lowercase, strip "site_", pad to 2 digits, and prefix with "s"
num <- sub("^site_", "", tolower(coor$site))
num <- as.integer(num) # turns "01" or "1" into 1
coor$site <- sprintf("s%02d", num) # -> "s01", "s02", ...


meta$lat <- coor$Latitude[match(meta$site, coor$site)]
meta$lon <- coor$Longitude[match(meta$site, coor$site)]


# Ensure sample_id is the key
stopifnot("sample_id" %in% names(comm), "sample_id" %in% names(meta))

# Move sample_id to rownames and keep only ASV columns
# Start from your 'comm' tibble/data.frame that has a 'sample_id' column
comm_mat <- comm %>%
  as.data.frame(stringsAsFactors = FALSE) %>% # ensure plain data.frame
  tibble::remove_rownames() %>% # <-- clear any existing rownames
  {
    # ensure the key exists and is usable
    stopifnot("sample_id" %in% names(.))
    .$sample_id <- as.character(.$sample_id)
    if (anyNA(.$sample_id)) stop("sample_id contains NA values.")
    if (anyDuplicated(.$sample_id)) {
      warning(
        "Duplicate sample_id values found; making them unique with make.unique()."
      )
      .$sample_id <- make.unique(.$sample_id, sep = "_dup")
    }
    .
  } %>%
  tibble::column_to_rownames("sample_id") # now it will work


# Optional: pool replicates by site (if you have a 'site' column in metadata)
# meta must include a 'site' column for this step; comment this block out if not pooling
if ("site" %in% names(meta)) {
  # align rows
  meta_aligned <- meta %>% filter(sample_id %in% rownames(comm_mat))
  comm_mat <- comm_mat[meta_aligned$sample_id, , drop = FALSE]
  # pool by site (sum counts across replicates)
  comm_mat <- meta_aligned %>%
    select(sample_id, site) %>%
    cbind(comm_mat) %>%
    group_by(site) %>%
    summarise(across(-sample_id, sum), .groups = "drop") %>%
    as.data.frame() %>%
    tibble::column_to_rownames("site")
  # collapse metadata to unique site rows (mean lat/lon if slight jitter)
  meta <- meta %>%
    group_by(site) %>%
    summarise(
      lat = mean(lat, na.rm = TRUE),
      lon = mean(lon, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(sample_id = site)
}

# Keep only samples present in both tables, same order
meta <- meta %>% filter(sample_id %in% rownames(comm_mat))
comm_mat <- comm_mat[meta$sample_id, , drop = FALSE]

# ==== 2) Choose presence–absence or abundance path ====
# Presence–absence (Jaccard/Sørensen; turnover + nestedness)
comm_pa <- (comm_mat > 0) * 1

# Abundance (Bray-Curtis; splits into balanced variation vs abundance gradient)
# Keep comm_mat as is for abundance partitioning

# ==== 3) betapart: presence–absence partition (turnover & nestedness) ====
core_pa <- betapart.core(comm_pa)
# Jaccard family (index.family="jac") is common; Sørensen ("sor") also ok.
pair_pa <- beta.pair(core_pa, index.family = "jac")
# pair_pa is a list of dist objects:
#   beta.jtu = turnover (Jaccard turnover)
#   beta.jne = nestedness-resultant
#   beta.jac = total Jaccard dissimilarity

# ==== 4) betapart: abundance partition (Bray-Curtis components) ====
pair_ab <- beta.pair.abund(comm_mat, index.family = "bray")
# pair_ab contains:
#   beta.bray.bal = balanced variation in abundance (analogous to turnover)
#   beta.bray.gra = abundance gradient (loss/gain)
#   beta.bray     = total Bray-Curtis dissimilarity

# ==== 5) Build a tidy pairwise table and add geographic distances ====
# Geographic distance matrix (meters) using haversine
coords <- meta %>%
  select(sample_id, lon, lat) %>%
  arrange(match(sample_id, rownames(comm_mat)))

gdist <- geodist(coords[, c("lon", "lat")], measure = "haversine") # meters
rownames(gdist) <- coords$sample_id
colnames(gdist) <- coords$sample_id

# Helper to melt a 'dist' object to a 3-col data.frame
dist_to_df <- function(d, name) {
  m <- as.matrix(d)
  tibble::tibble(
    sample_i = rep(rownames(m), times = ncol(m)),
    sample_j = rep(colnames(m), each = nrow(m)),
    value = as.vector(m)
  ) %>%
    filter(sample_i < sample_j) %>%
    rename(!!name := value)
}

df <- dist_to_df(pair_pa$beta.jtu, "beta_jacc_turnover") %>%
  left_join(
    dist_to_df(pair_pa$beta.jne, "beta_jacc_nestedness"),
    by = c("sample_i", "sample_j")
  ) %>%
  left_join(
    dist_to_df(pair_pa$beta.jac, "beta_jacc_total"),
    by = c("sample_i", "sample_j")
  ) %>%
  left_join(
    dist_to_df(pair_ab$beta.bray.bal, "beta_bray_balanced"),
    by = c("sample_i", "sample_j")
  ) %>%
  left_join(
    dist_to_df(pair_ab$beta.bray.gra, "beta_bray_gradient"),
    by = c("sample_i", "sample_j")
  ) %>%
  left_join(
    dist_to_df(pair_ab$beta.bray, "beta_bray_total"),
    by = c("sample_i", "sample_j")
  )

# Add geographic distances (km)
gdist_df <- gdist %>%
  as.matrix() %>%
  {
    tibble::tibble(
      sample_i = rep(rownames(.), times = ncol(.)),
      sample_j = rep(colnames(.), each = nrow(.)),
      dist_km = as.vector(.) / 1000
    )
  } %>%
  filter(sample_i < sample_j)

df <- df %>% left_join(gdist_df, by = c("sample_i", "sample_j"))

# ==== 6) Distance–decay visualization (turnover vs distance) ====
ggplot(df, aes(dist_km, beta_jacc_turnover)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Geographic distance (km)",
    y = "Jaccard turnover (β_jtu)",
    title = "Distance–decay of bacterial community turnover"
  )

# (Optional) same for Bray–Curtis balanced change
ggplot(df, aes(dist_km, beta_bray_balanced)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Geographic distance (km)",
    y = "Bray–Curtis balanced component",
    title = "Distance–decay of abundance-based turnover"
  )

# ==== 7) Simple tests of distance–decay ====
# Mantel test (note: interpret with care; spatial autocorrelation caveats)
mantel_turnover <- mantel(
  pair_pa$beta.jtu,
  dist(as.matrix(gdist)),
  method = "pearson",
  permutations = 999
)
mantel_turnover

# Linear model on pairwise data (pseudo-replication caveats)
summary(lm(beta_jacc_turnover ~ dist_km, data = df))

# ==== 8) Multi-site β (summary across all samples) ====
multi_pa <- beta.multi(core_pa, index.family = "jac")
multi_pa # yields overall turnover, nestedness, total dissimilarity

###through time

# =========================
# TEMPORAL β-diversity (sample-level only)
# =========================

# Packages
library(dplyr)
library(tibble)
library(stringr)
library(lubridate)
library(ggplot2)
library(betapart)
library(vegan)

# -------------------------
# 0) Inputs and guards
# -------------------------

# Replace with your file paths / data frames
comm <- community_matrix # rows = samples, cols = ASVs, first col = sample_id
comm$sample_id <- rownames(comm)
meta <- metadata # columns: sample_id, lat, lon, [habitat, date, ...]
meta$sample_id <- rownames(meta)
##add lat and lon

stopifnot(exists("metadata"), exists("comm_mat"))


# ensure sample_id column in metadata
meta <- meta_sample %>% as.data.frame(stringsAsFactors = FALSE)
if (!"sample_id" %in% names(meta)) {
  if (!is.null(rownames(meta))) {
    meta$sample_id <- rownames(meta)
  } else {
    stop("metadata must have rownames or a 'sample_id' column.")
  }
}

# -------------------------
# 1) Standardize period/site & join period–site dates
# -------------------------
fmt_period <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("^period[o]?_?", "", x)
  x <- gsub("^p_?", "p", x)
  x <- gsub("^p$", "", x)
  x <- ifelse(grepl("^p", x), sub("^p", "", x), x)
  sprintf("p%02d", as.integer(x))
}
fmt_site <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("^site_?", "", x)
  x <- gsub("^sitio_?", "", x)
  x <- gsub("^s_?", "s", x)
  x <- gsub("^s$", "", x)
  x <- ifelse(grepl("^s", x), sub("^s", "", x), x)
  sprintf("s%02d", as.integer(x))
}

# guess period/site columns in meta (English/Spanish)
period_col <- names(meta)[str_detect(tolower(names(meta)), "period")][1]
site_col <- names(meta)[str_detect(tolower(names(meta)), "site|sitio")][1]
stopifnot(!is.na(period_col), !is.na(site_col))

meta <- meta %>%
  mutate(
    period_std = fmt_period(.data[[period_col]]),
    site_std = fmt_site(.data[[site_col]]),
    key = paste0(period_std, "_", site_std)
  )

# read the small period–site date table and standardize it too
sample_dates <- read.csv("data/sampling_dates.csv", stringsAsFactors = FALSE)

pd_col <- names(sample_dates)[str_detect(
  tolower(names(sample_dates)),
  "^period"
)][1]
st_col <- names(sample_dates)[str_detect(
  tolower(names(sample_dates)),
  "^sitio|^site"
)][1]
dt_col <- names(sample_dates)[str_detect(tolower(names(sample_dates)), "date")][
  1
]
stopifnot(!is.na(pd_col), !is.na(st_col), !is.na(dt_col))

sample_dates <- sample_dates %>%
  mutate(
    period_std = fmt_period(.data[[pd_col]]),
    site_std = fmt_site(.data[[st_col]]),
    key = paste0(period_std, "_", site_std),
    Representative_Date = as.Date(.data[[dt_col]])
  ) %>%
  select(key, Representative_Date)

# join date into meta (prefer existing meta$date if present+parseable)
if (!"date" %in% names(meta)) meta$date <- NA
meta <- meta %>%
  left_join(sample_dates, by = "key") %>%
  mutate(date = coalesce(as.Date(date), Representative_Date))

# fallback: try auto-parse any Fecha/Date column if still NA
if (all(is.na(meta$date))) {
  cand <- names(meta)[str_detect(
    tolower(names(meta)),
    "fecha|date|sampling|collection"
  )]
  if (length(cand) > 0) {
    parsed <- parse_date_time(
      meta[[cand[1]]],
      orders = c(
        "ymd",
        "dmy",
        "mdy",
        "Y-m-d",
        "d/m/Y",
        "m/d/Y",
        "ymd HMS",
        "dmy HMS",
        "mdy HMS"
      ),
      tz = "UTC"
    )
    meta$date <- as.Date(parsed)
  }
}
stopifnot(any(!is.na(meta$date)))

# -------------------------
# 2) Align with comm_mat_sample and keep valid-dated samples
# -------------------------
# comm_mat_sample must have rownames = sample_id
if (is.null(rownames(comm))) stop("comm_mat must have rownames = sample_id.")

meta_t <- meta %>% filter(sample_id %in% rownames(comm), !is.na(date))
comm_mat_t <- comm[meta_t$sample_id, , drop = FALSE]
if (nrow(comm_mat_t) < 2) {
  stop(
    "Fewer than 2 samples with valid dates after alignment. Check IDs and date joins."
  )
}

# -------------------------
# 3) Compute β (presence–absence; Jaccard family)
# -------------------------
comm_pa_t <- (comm_mat_t > 0) * 1
core_pa_t <- betapart.core(comm_pa_t)
pair_pa_t <- beta.pair(core_pa_t, index.family = "jac")

# -------------------------
# 4) Helper: dist -> tidy df (upper triangle) with safe labels
# -------------------------
dist_to_df <- function(d, name) {
  m <- as.matrix(d)
  labs <- rownames(m)
  if (is.null(labs) || any(labs == "")) labs <- attr(d, "Labels")
  if (is.null(labs) || anyNA(labs)) labs <- as.character(seq_len(nrow(m)))
  rownames(m) <- labs
  colnames(m) <- labs
  idx <- which(upper.tri(m), arr.ind = TRUE)
  tibble(
    sample_i = labs[idx[, 1]],
    sample_j = labs[idx[, 2]],
    !!name := m[idx]
  )
}

# pairwise β table
df_time <- dist_to_df(pair_pa_t$beta.jtu, "beta_jacc_turnover") %>%
  left_join(
    dist_to_df(pair_pa_t$beta.jne, "beta_jacc_nestedness"),
    by = c("sample_i", "sample_j")
  ) %>%
  left_join(
    dist_to_df(pair_pa_t$beta.jac, "beta_jacc_total"),
    by = c("sample_i", "sample_j")
  )

# -------------------------
# 5) Add dates and within-site temporal differences
# -------------------------
# choose a site column to carry forward (use standardized if needed)
site_to_use <- if ("site" %in% names(meta_t)) "site" else "site_std"
if (!site_to_use %in% names(meta_t)) {
  meta_t$site <- meta_t$site_std
  site_to_use <- "site"
}

meta_min <- meta_t %>% select(sample_id, date, site = all_of(site_to_use))

df_time <- df_time %>%
  left_join(
    meta_min %>% rename(date_i = date, site_i = site),
    by = c("sample_i" = "sample_id")
  ) %>%
  left_join(
    meta_min %>% rename(date_j = date, site_j = site),
    by = c("sample_j" = "sample_id")
  ) %>%
  mutate(
    time_diff_days = abs(as.numeric(difftime(date_j, date_i, units = "days")))
  )

# all pairs with valid time differences
df_time_all <- df_time %>% filter(!is.na(time_diff_days))
# within-site only (recommended to avoid spatial confounding)
df_time_within <- df_time_all %>% filter(site_i == site_j)

message(
  "Pairs: all = ",
  nrow(df_time_all),
  "; within-site = ",
  nrow(df_time_within)
)

# -------------------------
# 6) Plots & stats (within-site recommended)
# -------------------------
ggplot(df_time_within, aes(time_diff_days, beta_jacc_turnover)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Temporal distance (days)",
    y = "Jaccard turnover (β_JTU)",
    title = "Temporal distance–decay of bacterial community turnover (within-site)"
  )

# Mantel test vs temporal distance (matrix form; all samples)
dates_vec <- meta_t$date
names(dates_vec) <- meta_t$sample_id
tdiff_mat <- outer(dates_vec, dates_vec, function(a, b) abs(as.numeric(a - b)))
tdiff_dist <- as.dist(tdiff_mat[rownames(comm_mat_t), rownames(comm_mat_t)])

mantel_time <- mantel(
  pair_pa_t$beta.jtu,
  tdiff_dist,
  method = "pearson",
  permutations = 999
)
print(mantel_time)

# Pairwise LM (note non-independence)
print(summary(lm(beta_jacc_turnover ~ time_diff_days, data = df_time_within)))

# -------------------------
# 7) Adjacent-time turnover within sites (clear interpretation)
# -------------------------
adjacent_turnover <- meta_t %>%
  arrange(.data[[site_to_use]], date) %>%
  group_by(.data[[site_to_use]]) %>%
  mutate(
    next_sample = lead(sample_id),
    next_diff_days = as.numeric(difftime(lead(date), date, units = "days"))
  ) %>%
  ungroup() %>%
  select(
    site = all_of(site_to_use),
    sample_i = sample_id,
    sample_j = next_sample,
    time_diff_days = next_diff_days
  ) %>%
  filter(!is.na(sample_j), !is.na(time_diff_days)) %>%
  left_join(
    df_time %>% select(sample_i, sample_j, beta_jacc_turnover),
    by = c("sample_i", "sample_j")
  )

ggplot(adjacent_turnover, aes(time_diff_days, beta_jacc_turnover)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Days between successive samples (within-site)",
    y = "Jaccard turnover (β_JTU)",
    title = "Within-site temporal turnover (adjacent timepoints)"
  )

# Mantel test (β_JTU vs temporal distance)
# Uses the temporal distance matrix tdiff_dist already built from meta_t$date
mantel_turnover_time <- mantel(
  pair_pa_t$beta.jtu,
  tdiff_dist,
  method = "pearson",
  permutations = 999
)
mantel_turnover_time

# Linear model on pairwise data (use within-site to avoid spatial confounding)
summary(lm(beta_jacc_turnover ~ time_diff_days, data = df_time_within))


# Presence–absence (Jaccard family)
multi_pa_time <- beta.multi(core_pa_t, index.family = "jac")
multi_pa_time # $beta.JTU (turnover), $beta.JNE (nestedness), $beta.JAC (total)


####FUNGIIIII

# ==== 1) Load your data ====
# Replace with your file paths / data frames
comm <- community_matrix_fung # rows = samples, cols = ASVs, first col = sample_id
comm$sample_id <- rownames(comm)
meta <- metadata_fung # columns: sample_id, lat, lon, [habitat, date, ...]
meta$sample_id <- rownames(meta)
##add lat and lon

coor <- read.csv("data/coor.csv")
head(coor)
head(meta)

# make sure the column exists
stopifnot("site" %in% names(coor))

# normalize to lowercase, strip "site_", pad to 2 digits, and prefix with "s"
num <- sub("^site_", "", tolower(coor$site))
num <- as.integer(num) # turns "01" or "1" into 1
coor$site <- sprintf("s%02d", num) # -> "s01", "s02", ...


meta$lat <- coor$Latitude[match(meta$site, coor$site)]
meta$lon <- coor$Longitude[match(meta$site, coor$site)]


# Ensure sample_id is the key
stopifnot("sample_id" %in% names(comm), "sample_id" %in% names(meta))

# Move sample_id to rownames and keep only ASV columns
# Start from your 'comm' tibble/data.frame that has a 'sample_id' column
comm_mat <- comm %>%
  as.data.frame(stringsAsFactors = FALSE) %>% # ensure plain data.frame
  tibble::remove_rownames() %>% # <-- clear any existing rownames
  {
    # ensure the key exists and is usable
    stopifnot("sample_id" %in% names(.))
    .$sample_id <- as.character(.$sample_id)
    if (anyNA(.$sample_id)) stop("sample_id contains NA values.")
    if (anyDuplicated(.$sample_id)) {
      warning(
        "Duplicate sample_id values found; making them unique with make.unique()."
      )
      .$sample_id <- make.unique(.$sample_id, sep = "_dup")
    }
    .
  } %>%
  tibble::column_to_rownames("sample_id") # now it will work


# Optional: pool replicates by site (if you have a 'site' column in metadata)
# meta must include a 'site' column for this step; comment this block out if not pooling
if ("site" %in% names(meta)) {
  # align rows
  meta_aligned <- meta %>% filter(sample_id %in% rownames(comm_mat))
  comm_mat <- comm_mat[meta_aligned$sample_id, , drop = FALSE]
  # pool by site (sum counts across replicates)
  comm_mat <- meta_aligned %>%
    select(sample_id, site) %>%
    cbind(comm_mat) %>%
    group_by(site) %>%
    summarise(across(-sample_id, sum), .groups = "drop") %>%
    as.data.frame() %>%
    tibble::column_to_rownames("site")
  # collapse metadata to unique site rows (mean lat/lon if slight jitter)
  meta <- meta %>%
    group_by(site) %>%
    summarise(
      lat = mean(lat, na.rm = TRUE),
      lon = mean(lon, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(sample_id = site)
}

# Keep only samples present in both tables, same order
meta <- meta %>% filter(sample_id %in% rownames(comm_mat))
comm_mat <- comm_mat[meta$sample_id, , drop = FALSE]

# ==== 2) Choose presence–absence or abundance path ====
# Presence–absence (Jaccard/Sørensen; turnover + nestedness)
comm_pa <- (comm_mat > 0) * 1

# Abundance (Bray-Curtis; splits into balanced variation vs abundance gradient)
# Keep comm_mat as is for abundance partitioning

# ==== 3) betapart: presence–absence partition (turnover & nestedness) ====
core_pa <- betapart.core(comm_pa)
# Jaccard family (index.family="jac") is common; Sørensen ("sor") also ok.
pair_pa <- beta.pair(core_pa, index.family = "jac")
# pair_pa is a list of dist objects:
#   beta.jtu = turnover (Jaccard turnover)
#   beta.jne = nestedness-resultant
#   beta.jac = total Jaccard dissimilarity

# ==== 4) betapart: abundance partition (Bray-Curtis components) ====
pair_ab <- beta.pair.abund(comm_mat, index.family = "bray")
# pair_ab contains:
#   beta.bray.bal = balanced variation in abundance (analogous to turnover)
#   beta.bray.gra = abundance gradient (loss/gain)
#   beta.bray     = total Bray-Curtis dissimilarity

# ==== 5) Build a tidy pairwise table and add geographic distances ====
# Geographic distance matrix (meters) using haversine
coords <- meta %>%
  select(sample_id, lon, lat) %>%
  arrange(match(sample_id, rownames(comm_mat)))

gdist <- geodist(coords[, c("lon", "lat")], measure = "haversine") # meters
rownames(gdist) <- coords$sample_id
colnames(gdist) <- coords$sample_id

# Helper to melt a 'dist' object to a 3-col data.frame
dist_to_df <- function(d, name) {
  m <- as.matrix(d)
  tibble::tibble(
    sample_i = rep(rownames(m), times = ncol(m)),
    sample_j = rep(colnames(m), each = nrow(m)),
    value = as.vector(m)
  ) %>%
    filter(sample_i < sample_j) %>%
    rename(!!name := value)
}

df <- dist_to_df(pair_pa$beta.jtu, "beta_jacc_turnover") %>%
  left_join(
    dist_to_df(pair_pa$beta.jne, "beta_jacc_nestedness"),
    by = c("sample_i", "sample_j")
  ) %>%
  left_join(
    dist_to_df(pair_pa$beta.jac, "beta_jacc_total"),
    by = c("sample_i", "sample_j")
  ) %>%
  left_join(
    dist_to_df(pair_ab$beta.bray.bal, "beta_bray_balanced"),
    by = c("sample_i", "sample_j")
  ) %>%
  left_join(
    dist_to_df(pair_ab$beta.bray.gra, "beta_bray_gradient"),
    by = c("sample_i", "sample_j")
  ) %>%
  left_join(
    dist_to_df(pair_ab$beta.bray, "beta_bray_total"),
    by = c("sample_i", "sample_j")
  )

# Add geographic distances (km)
gdist_df <- gdist %>%
  as.matrix() %>%
  {
    tibble::tibble(
      sample_i = rep(rownames(.), times = ncol(.)),
      sample_j = rep(colnames(.), each = nrow(.)),
      dist_km = as.vector(.) / 1000
    )
  } %>%
  filter(sample_i < sample_j)

df <- df %>% left_join(gdist_df, by = c("sample_i", "sample_j"))

# ==== 6) Distance–decay visualization (turnover vs distance) ====
ggplot(df, aes(dist_km, beta_jacc_turnover)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Geographic distance (km)",
    y = "Jaccard turnover (β_jtu)",
    title = "Distance–decay of bacterial community turnover"
  )

# (Optional) same for Bray–Curtis balanced change
ggplot(df, aes(dist_km, beta_bray_balanced)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Geographic distance (km)",
    y = "Bray–Curtis balanced component",
    title = "Distance–decay of abundance-based turnover"
  )

# ==== 7) Simple tests of distance–decay ====
# Mantel test (note: interpret with care; spatial autocorrelation caveats)
mantel_turnover <- mantel(
  pair_pa$beta.jtu,
  dist(as.matrix(gdist)),
  method = "pearson",
  permutations = 999
)
mantel_turnover

# Linear model on pairwise data (pseudo-replication caveats)
summary(lm(beta_jacc_turnover ~ dist_km, data = df))

# ==== 8) Multi-site β (summary across all samples) ====
multi_pa <- beta.multi(core_pa, index.family = "jac")
multi_pa # yields overall turnover, nestedness, total dissimilarity

###through time

# =========================
# TEMPORAL β-diversity (sample-level only)
# =========================

# Packages
library(dplyr)
library(tibble)
library(stringr)
library(lubridate)
library(ggplot2)
library(betapart)
library(vegan)

# -------------------------
# 0) Inputs and guards
# -------------------------

# Replace with your file paths / data frames
comm <- community_matrix_fung # rows = samples, cols = ASVs, first col = sample_id
comm$sample_id <- rownames(comm)
meta <- metadata_fung # columns: sample_id, lat, lon, [habitat, date, ...]
meta$sample_id <- rownames(meta)
##add lat and lon

# =======================================================
# Fungi — Temporal β-diversity (sample-level, robust)
# Inputs:
#   - comm: community matrix (rows = sample_id-like strings, cols = ASVs/OTUs)
#   - metadata: data.frame with period/site columns; sample_id as column or rownames
#   - data/sampling_dates.csv with Periodo/Sitio + Representative_Date (YYYY-MM-DD)
# =======================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
  library(lubridate)
  library(betapart)
  library(vegan)
  library(ggplot2)
})

# ---------------- 0) Inputs & guards ----------------
stopifnot(exists("comm"), exists("metadata"))
if (is.null(rownames(comm)))
  stop("comm must have rownames (sample-level identifiers).")
if (all(grepl("^s\\d{2}$", rownames(comm)))) {
  stop(
    "comm rownames look like site IDs (e.g., s01). Provide the *sample-level* fungi matrix, not pooled-by-site."
  )
}

# Metadata with sample_id
meta <- metadata %>% as.data.frame(stringsAsFactors = FALSE)
if (!"sample_id" %in% names(meta)) {
  if (!is.null(rownames(meta))) meta$sample_id <- rownames(meta) else
    stop("metadata must have rownames or a 'sample_id' column.")
}

# If 'sample_id' column is inside comm, move it to rownames
if ("sample_id" %in% colnames(comm)) {
  rn <- as.character(comm$sample_id)
  comm <- as.data.frame(comm, stringsAsFactors = FALSE)
  rownames(comm) <- rn
  comm$sample_id <- NULL
}

# --------------- 1) Standardize period/site & join dates ---------------
fmt_period <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("^period[o]?_?", "", x)
  x <- gsub("^p_?", "p", x)
  x <- gsub("^p$", "", x)
  x <- ifelse(grepl("^p", x), sub("^p", "", x), x)
  sprintf("p%02d", as.integer(x))
}
fmt_site <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("^site_?", "", x)
  x <- gsub("^sitio_?", "", x)
  x <- gsub("^s_?", "s", x)
  x <- gsub("^s$", "", x)
  x <- ifelse(grepl("^s", x), sub("^s", "", x), x)
  sprintf("s%02d", as.integer(x))
}

# Detect period/site columns (English/Spanish)
pcol <- names(meta)[grepl("period", tolower(names(meta)))]
scol <- names(meta)[grepl("site|sitio", tolower(names(meta)))]
if (length(pcol) == 0 || length(scol) == 0)
  stop(
    "metadata must contain period and site columns (e.g., 'period'/'Periodo' and 'site'/'Sitio')."
  )
pcol <- pcol[1]
scol <- scol[1]

meta <- meta %>%
  mutate(
    period_std = fmt_period(.data[[pcol]]),
    site_std = fmt_site(.data[[scol]]),
    key = paste0(period_std, "_", site_std)
  )

# Join representative dates if date missing
if (!"date" %in% names(meta) || all(is.na(meta$date))) {
  if (file.exists("data/sampling_dates.csv")) {
    sd <- read.csv("data/sampling_dates.csv", stringsAsFactors = FALSE)
    pd <- names(sd)[grepl("^period", tolower(names(sd)))][1]
    st <- names(sd)[grepl("^sitio|^site", tolower(names(sd)))][1]
    dt <- names(sd)[grepl("date", tolower(names(sd)))][1]
    stopifnot(!is.na(pd), !is.na(st), !is.na(dt))
    sd <- sd %>%
      mutate(
        period_std = fmt_period(.data[[pd]]),
        site_std = fmt_site(.data[[st]]),
        key = paste0(period_std, "_", site_std),
        Representative_Date = as.Date(.data[[dt]])
      ) %>%
      select(key, Representative_Date)
    meta <- left_join(meta, sd, by = "key")
  }
}

# Normalize/parse meta$date and coalesce with Representative_Date
if (!"date" %in% names(meta)) meta$date <- as.Date(NA)
if (is.list(meta$date) || is.data.frame(meta$date))
  meta$date <- unlist(meta$date, use.names = FALSE)
if (is.numeric(meta$date)) {
  meta$date <- as.Date(meta$date, origin = "1899-12-30") # Excel serials
} else {
  meta$date <- as.character(meta$date)
  meta$date <- suppressWarnings(as.Date(parse_date_time(
    meta$date,
    orders = c(
      "ymd",
      "dmy",
      "mdy",
      "Y-m-d",
      "d/m/Y",
      "m/d/Y",
      "ymd HMS",
      "dmy HMS",
      "mdy HMS"
    ),
    tz = "UTC"
  )))
}
if ("Representative_Date" %in% names(meta)) {
  if (!inherits(meta$Representative_Date, "Date"))
    meta$Representative_Date <- as.Date(meta$Representative_Date)
  meta$date <- coalesce(meta$date, meta$Representative_Date)
  meta$Representative_Date <- NULL
}
if (all(is.na(meta$date))) {
  cand <- names(meta)[grepl(
    "fecha|date|sampling|collection",
    tolower(names(meta))
  )]
  if (length(cand) > 0) {
    parsed <- parse_date_time(
      meta[[cand[1]]],
      orders = c(
        "ymd",
        "dmy",
        "mdy",
        "Y-m-d",
        "d/m/Y",
        "m/d/Y",
        "ymd HMS",
        "dmy HMS",
        "mdy HMS"
      ),
      tz = "UTC"
    )
    meta$date <- as.Date(parsed)
  }
}
if (all(is.na(meta$date)))
  stop(
    "No parseable dates in metadata. Supply a date or fix data/sampling_dates.csv."
  )

# --------------- 2) CORE-ID RECONCILIATION & ALIGNMENT ---------------
# Many pipelines differ only by suffix (e.g., _ITS_F96 vs _16S_B96); match on the core before "_"
get_core_id <- function(x) sub("_.*$", "", as.character(x))

# Ensure site column for later within-site analysis
if (!"site" %in% names(meta)) meta$site <- meta$site_std

# Build core IDs
comm_core <- tibble(
  rn_comm = rownames(comm),
  core_id = get_core_id(rownames(comm))
)
meta <- meta %>%
  mutate(sample_id = as.character(sample_id), core_id = get_core_id(sample_id))

# Uniqueness checks
dups_comm <- comm_core %>% count(core_id) %>% filter(n > 1)
dups_meta <- meta %>% count(core_id) %>% filter(n > 1)
if (nrow(dups_comm) > 0)
  stop(
    "Non-unique core_id in comm. Examples: ",
    paste(head(dups_comm$core_id, 5), collapse = ", ")
  )
if (nrow(dups_meta) > 0)
  stop(
    "Non-unique core_id in metadata. Examples: ",
    paste(head(dups_meta$core_id, 5), collapse = ", "),
    "\nDisambiguate metadata (e.g., keep the fungi/ITS row for each core)."
  )

# Align by core_id
map3 <- left_join(
  comm_core,
  meta %>% select(core_id, sample_id, date, site, site_std),
  by = "core_id"
)
matched3 <- sum(!is.na(map3$sample_id))
if (matched3 == 0) {
  stop(
    "Could not align comm to metadata via core_id.\n",
    "Examples comm IDs: ",
    paste(head(rownames(comm), 5), collapse = ", "),
    "\n",
    "Examples meta IDs: ",
    paste(head(meta$sample_id, 5), collapse = ", "),
    "\n",
    "Check that the part before '_' matches across data."
  )
}

# Build aligned meta_t and comm_mat_t (index comm by core_id; label rows with full sample_id)
rownames(comm) <- comm_core$core_id
meta_t <- map3 %>%
  filter(!is.na(sample_id), !is.na(date)) %>%
  transmute(
    core_id,
    sample_id,
    date = as.Date(date),
    site = coalesce(site, site_std)
  )
stopifnot(all(meta_t$core_id %in% rownames(comm)))
comm_mat_t <- comm[meta_t$core_id, , drop = FALSE]
rownames(comm_mat_t) <- meta_t$sample_id # for readability

if (nrow(comm_mat_t) < 2)
  stop(
    "Fewer than 2 fungi samples with valid dates after alignment by core_id."
  )

# --------------- 3) Temporal β (presence–absence; Jaccard) ---------------
comm_pa_t <- (comm_mat_t > 0) * 1
core_pa_t <- betapart::betapart.core(comm_pa_t)
pair_pa_t <- betapart::beta.pair(core_pa_t, index.family = "jac")

# Helper: dist/matrix -> tidy upper-triangle
dist_to_df <- function(d, name) {
  m <- as.matrix(d)
  labs <- rownames(m)
  if (is.null(labs) || any(labs == "")) labs <- attr(d, "Labels")
  if (is.null(labs) || anyNA(labs)) labs <- as.character(seq_len(nrow(m)))
  rownames(m) <- labs
  colnames(m) <- labs
  idx <- which(upper.tri(m), arr.ind = TRUE)
  tibble(sample_i = labs[idx[, 1]], sample_j = labs[idx[, 2]], !!name := m[idx])
}

df_time <- dist_to_df(pair_pa_t$beta.jtu, "beta_jacc_turnover") %>%
  left_join(
    dist_to_df(pair_pa_t$beta.jne, "beta_jacc_nestedness"),
    by = c("sample_i", "sample_j")
  ) %>%
  left_join(
    dist_to_df(pair_pa_t$beta.jac, "beta_jacc_total"),
    by = c("sample_i", "sample_j")
  )

# Add dates & within-site temporal differences
meta_min <- meta_t %>% select(sample_id, date, site)
df_time <- df_time %>%
  left_join(
    meta_min %>% rename(date_i = date, site_i = site),
    by = c("sample_i" = "sample_id")
  ) %>%
  left_join(
    meta_min %>% rename(date_j = date, site_j = site),
    by = c("sample_j" = "sample_id")
  ) %>%
  mutate(
    time_diff_days = abs(as.numeric(difftime(date_j, date_i, units = "days")))
  )
df_time_within <- df_time %>% filter(site_i == site_j, !is.na(time_diff_days))

# --------------- 4) Tests, plots, and summaries (with fallback) ---------------
# Temporal distance matrix (comm_mat_t order)
dates_vec <- meta_t$date
names(dates_vec) <- meta_t$sample_id
tdiff_mat <- outer(dates_vec, dates_vec, function(a, b) abs(as.numeric(a - b)))
tdiff_dist <- as.dist(tdiff_mat[rownames(comm_mat_t), rownames(comm_mat_t)])

use_all_pairs <- FALSE
if (nrow(df_time_within) == 0) {
  warning(
    "No within-site temporal pairs (site_i == site_j). Falling back to all pairs."
  )
  df_time_within <- df_time %>% filter(!is.na(time_diff_days))
  use_all_pairs <- TRUE
}

# Mantel (β_JTU vs temporal distance)
mantel_time_fungi <- vegan::mantel(
  pair_pa_t$beta.jtu,
  tdiff_dist,
  method = "pearson",
  permutations = 999
)
print(mantel_time_fungi)

# Pairwise LM
lm_time_fungi <- lm(beta_jacc_turnover ~ time_diff_days, data = df_time_within)
print(summary(lm_time_fungi))

# Plot
ttl <- if (use_all_pairs) "Fungi: temporal distance–decay (all pairs)" else
  "Fungi: temporal distance–decay (within-site)"
ggplot(df_time_within, aes(time_diff_days, beta_jacc_turnover)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Temporal distance (days)",
    y = "Jaccard turnover (β_JTU)",
    title = ttl
  )

# Multi-site β (presence–absence; Jaccard)
multi_pa_time_fungi <- betapart::beta.multi(core_pa_t, index.family = "jac")
print(multi_pa_time_fungi) # $beta.JTU, $beta.JNE, $beta.JAC

# (Optional) Adjacent-time turnover (within-site)
adjacent_turnover_fungi <- meta_t %>%
  arrange(site, date) %>%
  group_by(site) %>%
  mutate(
    next_sample = lead(sample_id),
    next_diff_days = as.numeric(difftime(lead(date), date, units = "days"))
  ) %>%
  ungroup() %>%
  select(
    site,
    sample_i = sample_id,
    sample_j = next_sample,
    time_diff_days = next_diff_days
  ) %>%
  filter(!is.na(sample_j), !is.na(time_diff_days)) %>%
  left_join(
    df_time %>% select(sample_i, sample_j, beta_jacc_turnover),
    by = c("sample_i", "sample_j")
  )

if (nrow(adjacent_turnover_fungi) > 0) {
  print(head(adjacent_turnover_fungi))
  ggplot(adjacent_turnover_fungi, aes(time_diff_days, beta_jacc_turnover)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(
      x = "Days between successive samples (within-site)",
      y = "Jaccard turnover (β_JTU)",
      title = "Fungi: within-site temporal turnover (adjacent timepoints)"
    )
} else {
  message("No sites with ≥2 temporal samples — skipping adjacent-time plot.")
}
