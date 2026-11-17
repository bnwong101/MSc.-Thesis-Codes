#Descriptive Statistics - Thesis Brandon Wong

library(dplyr)
library ("readxl")
#Clean dataset

# 1 ── File inspection & set-up -------------------------------------------
setwd("/Users/Brandon W/Documents/Western/THESIS/Analysis/FINAL ANALYSIS")  
csv <- "stroke_pm25_temp_imputed.csv"

# 2 ── Load data -----------------------------------------------------------
df <- read.csv(csv, stringsAsFactors = FALSE)

# 3 ── Inspect data --------------------------------------------------------
colnames(df)

head(df, 10)

names(df) [names(df) == "stoke_type"] <- "stroke_type"

# 4 ── Filter data ----------------------------------------------------------
thesis_data <- df[,c("case", "match", "fis_yr", "visit_d", "fsa", "age_yr", 
                     "male", "stroke_type", "pm25_lag0", "pm25_lag1", 
                     "pm25_lag2", "pm25_lag3", "pm25_lag4", "pm25_lag5", 
                     "pm25_lag6", "pm25_lag7", "temp_avg_lag0", 
                     "temp_avg_lag1", "temp_avg_lag2", "temp_avg_lag3", 
                     "temp_avg_lag4", "temp_avg_lag5", "temp_avg_lag6", 
                     "temp_avg_lag7", "temp_max_lag0", "temp_max_lag1", 
                     "temp_max_lag2", "temp_max_lag3", "temp_max_lag4", 
                     "temp_max_lag5", "temp_max_lag6", "temp_max_lag7", 
                     "humidity_lag0", "humidity_lag1", "humidity_lag2", 
                     "humidity_lag3", "humidity_lag4", "humidity_lag5", 
                     "humidity_lag6", "humidity_lag7")]
colnames(thesis_data)
head(thesis_data, 10)

#Categorical age variable
thesis_data$age_cat <- cut(
  thesis_data$age_yr,
  breaks = c(0, 65, 75, Inf),
  labels = FALSE,
  right = TRUE
) - 1

# Dataset includes only those >=18
# 1 = 18–65
# 2 = 65–74
# 3 = ≥75

# Define binary region indicator
thesis_data <- thesis_data %>%
  mutate(
    region_north = ifelse(substr(fsa, 1, 1) == "P", 1, 0)
  )

thesis_data$region_north <- factor(
  thesis_data$region_north,
  levels = c(0, 1),
  labels = c("Southern Ontario", "Northern Ontario")
)

## Attach health-region codes
fsa_hr <- read_xlsx("Canada_FSA_to_HR18_Correspondence.xlsx") %>% 
  rename(fsa = FSA)
thesis_data  <- merge(thesis_data, fsa_hr[,c("fsa","HR_UID")], by = "fsa")

#Descriptive Statistics for characteristics of cases
df1 <- thesis_data %>% filter(case==1)
colnames(df1)
head(df1, 10)

#HR
cases_per_HR <- df1 %>% 
  count(HR_UID, name = "n_cases")
n_HR <- nrow(cases_per_HR)
print(n_HR) # Number of HR
HR_summary <- cases_per_HR %>% 
  summarise(
    mean_cases = mean(n_cases),
    sd_cases   = sd(n_cases),
    min_cases  = min(n_cases),
    max_cases  = max(n_cases)
  )
print(HR_summary)


#Age
age_table <- table(df1$age_cat)
age_percentage_table <- prop.table(age_table)*100

print(age_table)
print(age_percentage_table)

mean(df1$age_yr)
min(df1$age_yr)
max(df1$age_yr)
sd(df1$age_yr)

#Sex
sex_table <- table(df1$male)
sex_percentage_table <- prop.table(sex_table)*100

print(sex_table)
print(sex_percentage_table)

n_cases <- nrow(df1)                      # total number of case==1 rows
n_missing_sex <- sum(is.na(df1$male))     # how many of those have male==NA
n_cases; n_missing_sex

tablex<- table(df1$male, useNA="ifany")
tablexp<-prop.table(sex_table)*100
print(tablexp)


#Stroke types
stroke_table <- table(df1$stroke_type)
stroke_percentage_table <- prop.table(stroke_table)*100

print(stroke_table)
print(stroke_percentage_table)

#FSAs
cases_per_fsa <- df1 %>% 
  count(fsa, name = "n_cases")
n_fsa <- nrow(cases_per_fsa)
print(n_fsa) # Number of FSAs
fsa_summary <- cases_per_fsa %>% 
  summarise(
    mean_cases = mean(n_cases),
    sd_cases   = sd(n_cases),
    min_cases  = min(n_cases),
    max_cases  = max(n_cases)
  )
print(fsa_summary)

#Region
region_table      <- table(df1$region_north)
region_pct        <- prop.table(region_table) * 100
print(region_table)
print(region_pct)

#Descriptive Statistics for environmental variables


# Unweighted statistics
# Create FSA-Year variable
thesis_data <- thesis_data %>%
  mutate(fsa_year = paste(fsa, fis_yr, sep = "_"))

# Calculate FSA-year level means for each environmental variable
fsa_year_means_all <- thesis_data %>%
  group_by(fsa, fis_yr) %>%
  summarise(
    temp_avg_lag0 = mean(temp_avg_lag0, na.rm = TRUE),
    humidity_lag0 = mean(humidity_lag0, na.rm = TRUE),
    pm25_lag0 = mean(pm25_lag0, na.rm = TRUE),
    .groups = "drop"
  )

# Descriptive stats for each environmental variable
summary_stats <- function(x) {
  c(
    mean   = mean(x, na.rm = TRUE),
    sd     = sd(x, na.rm = TRUE),
    min    = min(x, na.rm = TRUE),
    p25    = quantile(x, 0.25, na.rm = TRUE),
    median = quantile(x, 0.50, na.rm = TRUE),
    p75    = quantile(x, 0.75, na.rm = TRUE),
    max    = max(x, na.rm = TRUE)
  )
}

# Apply function to each environmental variable
env_stats <- bind_rows(
  temp_avg    = summary_stats(fsa_year_means_all$temp_avg_lag0),
  humidity    = summary_stats(fsa_year_means_all$humidity_lag0),
  pm25        = summary_stats(fsa_year_means_all$pm25_lag0),
  .id = "variable"
)

# Show results
print(env_stats)

# Weighted statistics (by case)
# Summarize environmental means and case counts at the FSA-year level
fsa_year_summary <- df1 %>%
  group_by(fsa, fis_yr) %>%
  summarise(
    n_cases         = n(),
    temp_avg_lag0   = mean(temp_avg_lag0, na.rm = TRUE),
    humidity_lag0   = mean(humidity_lag0, na.rm = TRUE),
    pm25_lag0       = mean(pm25_lag0, na.rm = TRUE),
    .groups = "drop"
  )

# Define weighted stat functions manually
weighted_mean <- function(x, w) {
  sum(w * x, na.rm = TRUE) / sum(w[!is.na(x)])
}

weighted_var <- function(x, w) {
  m <- weighted_mean(x, w)
  sum(w * (x - m)^2, na.rm = TRUE) / sum(w[!is.na(x)])
}

weighted_sd <- function(x, w) {
  sqrt(weighted_var(x, w))
}

weighted_quantile <- function(x, w, probs) {
  # Sort x and corresponding weights
  ord <- order(x)
  x_sorted <- x[ord]
  w_sorted <- w[ord]
  cum_weights <- cumsum(w_sorted) / sum(w_sorted)
  return(approx(cum_weights, x_sorted, xout = probs)$y)
}

# Apply these to each variable
get_weighted_stats <- function(var, weights) {
  c(
    weighted_mean   = weighted_mean(var, weights),
    weighted_sd     = weighted_sd(var, weights),
    weighted_min    = min(var, na.rm = TRUE),
    weighted_p25    = weighted_quantile(var, weights, 0.25),
    weighted_median = weighted_quantile(var, weights, 0.50),
    weighted_p75    = weighted_quantile(var, weights, 0.75),
    weighted_max    = max(var, na.rm = TRUE)
  )
}

# Compile results
weighted_env_stats <- bind_rows(
  temp_avg = get_weighted_stats(fsa_year_summary$temp_avg_lag0, fsa_year_summary$n_cases),
  humidity = get_weighted_stats(fsa_year_summary$humidity_lag0, fsa_year_summary$n_cases),
  pm25     = get_weighted_stats(fsa_year_summary$pm25_lag0, fsa_year_summary$n_cases),
  .id = "variable"
)

# View results
print(weighted_env_stats)


# 0) Install/load required packages
# install.packages(c("sf","dplyr","ggplot2","viridis"))

library(sf)
library(dplyr)
library(ggplot2)
library(viridis)

# 1) Read in FSA boundary shapefile (downloaded from Statistics Canada)
setwd("/Users/Brandon W/Documents/Western/THESIS/Analysis/FINAL ANALYSIS/FSA") 
fsa_shapes <- st_read("lfsa000b16a_e.shp")

# 2) Keep just CFSAUID and rename it to match your data
fsa_shapes <- fsa_shapes %>%
  select(CFSAUID) %>%
  rename(fsa = CFSAUID)

# 2.5) Filter to Ontario FSAs only (prefix K, L, M, N, P)  ### <— added
ont_prefixes <- c("K","L","M","N","P")
fsa_shapes <- fsa_shapes %>%
  filter(substr(fsa, 1, 1) %in% ont_prefixes)

# 3) Compute stroke counts per FSA from your case data (df1)
stroke_by_fsa <- df1 %>%
  count(fsa, name = "n_strokes")


# 1) Read in the three census files
pop11 <- read.csv("2011 data.CSV", stringsAsFactors = FALSE)
pop16 <- read.csv("2016 data.CSV", stringsAsFactors = FALSE)
pop21 <- read.csv("2021 data.csv", stringsAsFactors = FALSE)

colnames(pop11)
colnames(pop16)
colnames(pop21)


# 2) Select & rename the key columns (adjust “CFSAUID” & “Population” 
#    to whatever your files actually call them)
pop11 <- pop11 %>% 
  select(fsa = Geographic.name, pop2011 = Population..2011)
pop16 <- pop16 %>% 
  select(fsa = Geographic.code, pop2016 = Population..2016)
pop21 <- pop21 %>% 
  select(fsa = Geographic.name, pop2021 = Population..2021)

# 1) Clean each pop column: remove all non-digits, then convert
pop11 <- pop11 %>% 
  mutate(pop2011 = as.numeric(gsub("[^0-9]", "", pop2011)))

pop16 <- pop16 %>% 
  mutate(pop2016 = as.numeric(gsub("[^0-9]", "", pop2016)))

pop21 <- pop21 %>% 
  mutate(pop2021 = as.numeric(gsub("[^0-9]", "", pop2021)))

# 2) (Optional) check how many NAs remain
pop11 %>% summarise(n_na = sum(is.na(pop2011)))
pop16 %>% summarise(n_na = sum(is.na(pop2016)))
pop21 %>% summarise(n_na = sum(is.na(pop2021)))

# 3) Re-join and compute average population
pop_all <- pop11 %>%
  full_join(pop16, by = "fsa") %>%
  full_join(pop21, by = "fsa") %>%
  mutate(pop_avg = rowMeans(across(starts_with("pop")), na.rm = TRUE))

# 4) Combine with your stroke counts and compute rate per 1,000
stroke_rate_by_fsa <- stroke_by_fsa %>%
  left_join(pop_all, by = "fsa") %>%
  mutate(
    rate_per_1000 = n_strokes / pop_avg * 1000
  )

# 5) Join back to the Ontario‐only shapes
map_rate <- fsa_shapes %>%
  left_join(stroke_rate_by_fsa, by = "fsa")

# 6) Plot the choropleth of stroke rates
ggplot(map_rate) +
  geom_sf(aes(fill = rate_per_1000), color = NA) +
  scale_fill_viridis(
    na.value   = "grey90",
    name       = "Stroke rate\n(per 1,000)"
  ) +
  theme_minimal() +
  coord_sf(expand = FALSE)

# 4) Join counts back to the Ontario‐only spatial object
map_data <- fsa_shapes %>%
  left_join(stroke_by_fsa, by = "fsa")

# 5) Plot choropleth of stroke counts (Ontario only)
ggplot(map_data) +
  geom_sf(aes(fill = n_strokes), color = NA) +
  scale_fill_viridis(na.value = "grey90", name = "Stroke\ncount") +
  theme_minimal() +
  coord_sf(expand = FALSE)  # ensures tight bounding box

# 6) Likewise, for temperature distribution:
temp_by_fsa <- thesis_data %>%
  group_by(fsa) %>%
  summarise(mean_temp = mean(temp_avg_lag0, na.rm = TRUE), .groups="drop")

map_temp <- fsa_shapes %>%
  left_join(temp_by_fsa, by = "fsa")

ggplot(map_temp) +
  geom_sf(aes(fill = mean_temp), color = NA) +
  scale_fill_viridis(na.value = "grey90", name = "Mean\nTemp (°C)") +
  theme_minimal() +
  coord_sf(expand = FALSE)  # ensure it’s zoomed to Ontario