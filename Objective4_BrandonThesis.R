# Brandon Wong Thesis Objective 4
# Exposure of interest: temperature
# Outcome of interest: stroke-related ED visits
# Potential confounder: humidity

# 1 Data preparation

# Upload packages
library ("dlnm")
library ("mixmeta")
library ('splines')
library ('survival')
library ("boot")
library ("lme4")
library ("dplyr")
library ("readxl")
library ("timeDate")
library ("lubridate")

# Set working directory
setwd("/Users/Brandon W/Documents/Western/THESIS/Analysis/FINAL ANALYSIS")

# Load dataset
df <- read.csv("stroke_pm25_temp_imputed.csv", stringsAsFactors = FALSE)

df  <- df %>% rename(stroke_type = stoke_type)                 # fix typo
keep_vars <- c(
  "case","match","fis_yr","visit_d","fsa","age_yr","male","stroke_type",
  paste0("temp_avg_lag",0:7), paste0("temp_max_lag",0:7),
  paste0("humidity_lag",0:7)
)
thesis_data <- df[, keep_vars]

## Attach health-region codes
fsa_hr <- read_xlsx("Canada_FSA_to_HR18_Correspondence.xlsx") %>% 
  rename(fsa = FSA)
colnames(fsa_hr)
final  <- merge(thesis_data, fsa_hr[,c("fsa","ENGNAME")], by = "fsa")

## Add visit date and holiday flag
final$visit_date <- dmy(final$visit_d)
holiday_dates    <- as.Date(holidayTSX(unique(year(final$visit_date))))
final$stat_hol   <- as.integer(final$visit_date %in% holiday_dates)

colnames(final)

# ─────────────────────────────────────────────────────────────────────────────
# 4) Heatwave indicator + simple regression (Objective 4)
# ─────────────────────────────────────────────────────────────────────────────

library(data.table)  # for rleid()


# a) Compute heatwave indicator
setDT(final)  # convert to data.table in place
final[, thr_95 := quantile(temp_avg_lag0, probs = 0.95, na.rm = TRUE),
      by = ENGNAME]
final[, is_hot := temp_avg_lag0 > thr_95]
# identify consecutive runs of same is_hot value
final[, run_id := rleid(is_hot)]
# compute run length within each run
final[, run_len := .N, by = .(ENGNAME, run_id)]
# flag heatwave days: part of a hot run of length >= 2
final[, heatwave := as.integer(is_hot & run_len >= 2)]
# back to data.frame if you prefer
setDF(final)

# ─────────────────────────────────────────────────────────────────────────────
# 4b) Region‐specific conditional logistic models
# ─────────────────────────────────────────────────────────────────────────────

library(broom)    # for tidy()
library(mixmeta)  # for mixmeta()

# 1) Get the list of regions
regions <- unique(final$ENGNAME)

# 2) Loop over regions, fit clogit, extract heatwave log‐OR & SE
region_fit <- lapply(regions, function(hr) {
  dat <- subset(final, ENGNAME == hr)
  # only fit if there’s enough data
  if (nrow(dat) < 10) return(NULL)
  
  m <- clogit(case ~ heatwave + humidity_lag0 + strata(match),
              data = dat)
  td <- tidy(m, exponentiate = FALSE, conf.int = FALSE) 
  hw <- td[td$term == "heatwave", ]
  
  data.frame(
    ENGNAME = hr,
    beta   = hw$estimate,
    se     = hw$std.error
  )
})
region_df <- do.call(rbind, region_fit)
region_df$var <- region_df$se^2

# ─────────────────────────────────────────────────────────────────────────────
# 6) Pooled random‐effects meta‐analysis + 95% CI (fix for missing pool_mod)
# ─────────────────────────────────────────────────────────────────────────────

# build a list of 1×1 variance–covariance matrices, one per region
S_list <- lapply(region_df$var, function(v) matrix(v, 1, 1))

# fit random‐effects meta‐analysis of the heatwave log‐ORs
pool_mod <- mixmeta(
  formula = beta ~ 1,
  S       = S_list,
  data    = region_df,
  method  = "reml"
)

# extract fixed‐effect variance (on log‐OR scale) and between‐region τ²
var_fixed <- vcov(pool_mod)[1,1]
tau2      <- pool_mod$Psi[1,1]

# total variance on the log scale
total_var <- var_fixed + tau2
pooled_se <- sqrt(total_var)

# pooled log‐OR → OR and 95% CI
pooled_logOR <- coef(pool_mod)
pooled_OR    <- exp(pooled_logOR)
ci_low       <- exp(pooled_logOR - 1.96 * pooled_se)
ci_high      <- exp(pooled_logOR + 1.96 * pooled_se)

cat(sprintf(
  "Pooled OR for heatwave: %.3f (95%% CI: %.3f – %.3f)\n",
  pooled_OR, ci_low, ci_high
))


# ─────────────────────────────────────────────────────────────────────────────
# 7) Forest plot with ggplot2 (fixed)
# ─────────────────────────────────────────────────────────────────────────────

library(ggplot2)
library(dplyr)

plot_df <- region_df %>%
  # turn ENGNAME into a character region name
  transmute(
    region = as.character(ENGNAME),
    beta,
    se
  ) %>%
  # bind on the pooled estimate
  bind_rows(
    tibble(
      region = "Pooled",
      beta   = pooled_logOR,
      se     = pooled_se
    )
  ) %>%
  mutate(
    OR    = exp(beta),
    lower = exp(beta - 1.96 * se),
    upper = exp(beta + 1.96 * se),
    # define factor levels so Pooled comes last
    region = factor(region, levels = c(sort(unique(region)[region != "Pooled"]), "Pooled"))
  )

ggplot(plot_df, aes(x = OR, y = region)) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(
    x     = "Odds Ratio (95% CI)",
    y     = "Health Region",
    title = "Region-specific and Pooled Heatwave Effect"
  ) +
  theme_minimal(base_size = 14)
