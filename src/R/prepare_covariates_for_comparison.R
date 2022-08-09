library(tidyverse)

comparators_df <- tribble(
  ~epiweek, ~relative_self_mobility,
  ~number_introduction, ~number_introduction_lag1,
  ~number_introduction_lag2, ~proportion_cum_1st_dose_vaccinated,
  ~proportion_cum_2nd_dose_vaccinated,
  0, 0, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 0, 0, 0,
  0, 1, 0, 0, 0, 0, 0,
  0, 0, 1, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 0,
  0, 0, 0, 0, 1, 0, 0,
  0, 0, 0, 0, 0, 1, 0,
  0, 0, 0, 0, 0, 0, 1,
  1, 1, 0, 0, 0, 0, 0,
  1, 0, 1, 0, 0, 0, 0,
  1, 0, 0, 1, 0, 0, 0,
  1, 0, 0, 0, 1, 0, 0,
  1, 0, 0, 0, 0, 1, 0,
  1, 0, 0, 0, 0, 0, 1,
  1, 1, 0, 0, 0, 0, 1
  )  %>% 
  mutate(id=seq_along(epiweek))

saveRDS(comparators_df, "data/processed/covariates_to_include.rds")  
