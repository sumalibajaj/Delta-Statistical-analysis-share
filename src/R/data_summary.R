library(tidyverse)

df <- readRDS("data/processed/England_weekly_processed.rds")

n_data_points <- df %>% 
  group_by(location) %>% 
  summarise(n=sum(!is.na(total_count))) %>% 
  ungroup() %>% 
  summarise(lower=min(n),
            middle=mean(n),
            upper=max(n)) %>% 
  mutate(name="# data points")

delta_positive <- df %>% 
  filter(!is.na(total_count)) %>% 
  mutate(delta=delta_count/total_count) %>% 
  group_by(location) %>% 
  summarise(n=mean(delta) * 100) %>% 
  ungroup() %>% 
  summarise(lower=min(n),
            middle=mean(n),
            upper=max(n)) %>% 
  mutate(name="average % Delta positive (of non-missing samples)")

n_samples <- df %>% 
  filter(!is.na(total_count)) %>% 
  group_by(location) %>% 
  summarise(n=mean(total_count)) %>% 
  ungroup() %>% 
  summarise(lower=min(n),
            middle=mean(n),
            upper=max(n)) %>% 
  mutate(name="average total # samples per data point")

covars <- df %>% 
  mutate(proportion_cum_1st_dose_vaccinated=proportion_cum_1st_dose_vaccinated * 100,
         proportion_cum_2nd_dose_vaccinated=proportion_cum_2nd_dose_vaccinated * 100) %>% 
  select(epiweek, location,
         relative_self_mobility,
         proportion_cum_1st_dose_vaccinated,
         proportion_cum_2nd_dose_vaccinated,
         number_introduction) %>% 
  pivot_longer(-c(epiweek, location)) %>% 
  group_by(location, name) %>% 
  summarise(n=mean(value)) %>% 
  ungroup() %>% 
  group_by(name) %>% 
  summarise(lower=min(n),
            middle=mean(n),
            upper=max(n))

stacked_df <- n_data_points %>% 
  bind_rows(delta_positive) %>% 
  bind_rows(n_samples) %>% 
  bind_rows(covars) %>% 
  relocate(name) %>% 
  rename(Min=lower,
         Mean=middle,
         Max=upper) %>% 
  mutate(Min=round(Min, 2),
         Mean=round(Mean, 2),
         Max=round(Max, 2)) %>% 
  mutate(name=gsub("relative_self_mobility",
                   "average within UTLA mobility", name)) %>% 
  mutate(name=gsub("proportion_cum_2nd_dose_vaccinated",
                   "average cumulative % 2nd dose vacc.", name)) %>%
  mutate(name=gsub("proportion_cum_1st_dose_vaccinated",
                   "average cumulative % 1st dose vacc.", name)) %>% 
  mutate(name=gsub("number_introduction",
                   "average # introductions", name)) %>% 
  rename(Variable=name)
write.csv(stacked_df, "data/processed/data_summary.csv",
          row.names = FALSE)
