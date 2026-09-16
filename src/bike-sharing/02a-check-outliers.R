rm(list = ls())

# PACKAGES AND SOURCE  ----
suppressPackageStartupMessages(library(tidyverse))
source("src/01a-utils.R")
source("src/01b-gibbs.R")

# LOAD DATA AND FILTER ----
## Load data ----
data = readRDS("data/bike-sharing/milan_bikemi_data_for_regression.rds") |> 
  mutate(Weekend = as.numeric(weekdays(commit_at) %in% c("Saturday", "Sunday")),
         Rain = as.numeric(Precipitation > 0))
## Filter 2025 and no bikes ----
nobikes2025 = data |> filter(year(commit_at) == 2025 & unavail_type == "no_bikes")

# Exclude stations close to Pero: 334, 500 and 501. Actually 500 does not have any observation in 2025.
D = nobikes2025 |> filter(ID != 334, ID != 501) |> 
  select(-commit_at, -Precipitation, -Temperature, -Rain)


covariates = D |> 
  select(station, population_density, starts_with("count"), starts_with("density"), starts_with("closest")) |> 
  distinct() |> 
  pivot_longer(cols = -station) |> 
  mutate(log_value = log1p(value)) |> 
  pivot_longer(cols = c(value, log_value), names_to = "scale", values_to = "value") |>
  group_by(name, scale) |> 
  mutate(low_out = quantile(value, 0.25) - 5 * IQR(value),
         high_out = quantile(value, 0.75) + 5 * IQR(value)) |>
  ungroup() |> 
  mutate(outlier = ifelse(value < low_out | value > high_out, "yes", "no")) |>
  separate(name, into = c("measure", "type"), sep = "_") |> 
  mutate(type = if_else(measure == "population", "pop_dens", type))


covariates  |>
  filter(measure != "closest") |> 
  ggplot() +
  geom_boxplot(aes(y = value, col = scale)) +
  geom_text(aes(y = value, label = station, alpha = outlier), x = 0.2,
            size = 2) +
  scale_alpha_manual(values = c("no" = 0, "yes" = 1), guide = "none") +
  ggh4x::facet_nested(measure+scale~type, scales = "free")

ggsave("application/bike-sharing/pics/data/check_outliers.pdf",
       width = 11, height = 7)
