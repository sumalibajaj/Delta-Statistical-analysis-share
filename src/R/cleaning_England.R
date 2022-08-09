library(dplyr)
library(lubridate) # to use epiweek function

################################################################################
# UK NUTS1 level data and covariates
# First keep only UTLA level UK data
# Keep data from week 66 (starting date = 28 Mar 2021), and subtract 53 from all epiweek to only 2021 data
# Calculate the delta percentage for every week and each UTLA
# removing weeks following the first time 95% delta samples of total is reached
# Keep only those UTLAs which have data on delta for (doesn't matter if it's zero samples) more than three weeks
################################################################################
rm(list = ls())

dat <- read.csv(file = "data/raw/updated_delta.csv")

# To find earliest week with a delta sample
dat_min <-  dat %>%
  filter(level == "UK_utla") %>%
  group_by(code) %>%
  arrange(code, epiweek) %>%
  mutate(delta_yn = ifelse(delta_count >0 , 1, 0),
         delta_yn_cu = cumsum(delta_yn)) %>%
  filter(delta_yn_cu == 1 & delta_yn == 1) 

# Fill in missing NUTS1 / level
lookup <- dat %>% 
  select(code, NUTS1, level, population) %>% 
  unique() %>% 
  na.omit()

dat <- dat %>% 
  select(-NUTS1) %>%
  select(-level) %>% 
  select(-population) %>% 
  left_join(lookup)

# introductions: if NA, this means there were no estimated introductions
dat <- dat %>% 
  mutate(number_introduction=ifelse(is.na(number_introduction), 0, number_introduction))

dat_short <-  dat %>%
  filter(level == "UK_utla") %>%
  filter(!(NUTS1 %in% c("Northern_Ireland", "Wales", "Scotland"))) %>%
  mutate(prop_delta = delta_count/total_count,
         yn95 = ifelse(prop_delta >= 0.95, 1, 0)) %>%
  dplyr::rename(region = code) %>%
  group_by(region) %>%
  arrange(region, epiweek) %>%
  mutate(cu_yn95 = cumsum(yn95),
         prop_cu_1vac_utla_week = cumsum(proportion_population_new_1st_dose),
         prop_cu_2vac_utla_week = cumsum(proportion_population_new_2nd_dose)) %>%
  select(cu_yn95, epiweek, region, population, delta_count, total_count, location, 
         new_cases_reported, relative_self_mobility, number_introduction,
         prop_cu_1vac_utla_week, prop_cu_2vac_utla_week, 
         new_1st_dose_vaccine, new_2nd_dose_vaccine,
         NUTS1) %>% 
  mutate(number_introduction_lag1 = lag(number_introduction, default = NA),
         number_introduction_lag2 = lag(number_introduction, n=2, default = NA)) 

dat_final_short <- dat_short %>%
  arrange(region, epiweek) %>%
  filter(epiweek >= 66 & epiweek <= max(epiweek, na.rm = TRUE)) %>%
  mutate(epiweek = epiweek-53) %>%
  filter(epiweek <= 23) %>%# after week 23 >95% of samples are Delta
  group_by(region) %>%
  mutate(delta_utla_no = sum(!is.na(total_count))) %>%
  filter(delta_utla_no>=9) %>% 
  droplevels() %>%
  ungroup() %>% 
  mutate(nuts=as.numeric(as.factor(NUTS1)))

saveRDS(dat_final_short, "data/processed/England_weekly_cleaned.rds")

