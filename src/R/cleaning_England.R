library(dplyr)
library(lubridate) # to use epiweek function
library(ggplot2)

################################################################################
# UK NUTS1 level data and covariates
# First keep only UTLA level UK data
# Keep data from week 66 (starting date = 28 Mar 2021), and subtract 53 from all epiweek to only 2021 data
# Calculate the delta percentage for every week and each UTLA
# removing weeks following the first time 95% delta samples of total is reached
# Keep only those UTLAs which have data on delta for (doesn't matter if it's zero samples) more than three weeks
################################################################################
rm(list = ls())

# dat <- read.csv(file = "data/raw/master_with_number_introduction_and_lag.csv")
dat <- read.csv(file = "data/raw/master_total_updated.csv")

# To find earliest week with a delta sample
dat_min <-  dat %>%
  filter(level == "UK_utla") %>%
  group_by(code) %>%
  arrange(code, epiweek) %>%
  mutate(delta_yn = ifelse(delta_count >0 , 1, 0),
         delta_yn_cu = cumsum(delta_yn)) %>%
  filter(delta_yn_cu == 1 & delta_yn == 1) 

dat_vac <-  dat %>%
  filter(total_count > 0) %>%
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
         new_cases_reported,
         prop_cu_1vac_utla_week, prop_cu_2vac_utla_week, 
         new_1st_dose_vaccine, new_2nd_dose_vaccine,
         NUTS1) %>%
  na.omit()

dat_mob <-  dat %>%
  filter(total_count > 0) %>%
  filter(level == "UK_utla") %>%
  filter(!(NUTS1 %in% c("Northern_Ireland", "Wales", "Scotland"))) %>%
  mutate(prop_delta = delta_count/total_count,
         yn95 = ifelse(prop_delta >= 0.95, 1, 0)) %>%
  dplyr::rename(region = code) %>%
  group_by(region) %>%
  arrange(region, epiweek) %>%
  mutate(cu_yn95 = cumsum(yn95)) %>%
  select(cu_yn95, epiweek, region, population, delta_count, total_count, location, 
         new_cases_reported, relative_self_mobility,
         NUTS1) %>%
  na.omit()

dat_intro <-  dat %>%
  filter(total_count > 0) %>%
  filter(level == "UK_utla") %>%
  filter(!(NUTS1 %in% c("Northern_Ireland", "Wales", "Scotland"))) %>%
  mutate(prop_delta = delta_count/total_count,
         yn95 = ifelse(prop_delta >= 0.95, 1, 0)) %>%
  dplyr::rename(region = code) %>%
  arrange(region, epiweek) %>%
  group_by(region) %>%
  mutate(cu_yn95 = cumsum(yn95)) %>%
  select(cu_yn95, epiweek, region, population, delta_count, total_count, location, 
         new_cases_reported, number_introduction, relative_self_mobility,
         NUTS1) %>%
  na.omit()

dat_intro <- dat_intro %>% 
              mutate(number_introduction_lag1 = lag(number_introduction, default = NA),
                     number_introduction_lag2 = lag(number_introduction, n=2, default = NA)) 

dat_final_vac <- dat_vac %>%
  arrange(region, epiweek) %>%
  filter(cu_yn95 <1) %>% # first time 95% of samples become delta only
  filter(epiweek >= 66 & epiweek <= max(epiweek, na.rm = TRUE)) %>%
  mutate(epiweek = epiweek-53) %>%
  filter(epiweek <= 23) %>% # after week 23 most are missing values for delta
 group_by(region) %>%
  mutate(delta_utla_no = n()) %>%
  filter(delta_utla_no>=3) %>% 
  droplevels()

dat_final_vac$nuts <- as.numeric(as.factor(dat_final_vac$NUTS1))


dat_final_mob <- dat_mob %>%
  arrange(region, epiweek) %>%
  filter(cu_yn95 <1) %>% # first time 95% of samples become delta only
  filter(epiweek >= 66 & epiweek <= max(epiweek, na.rm = TRUE)) %>%
  mutate(epiweek = epiweek-53) %>%
  filter(epiweek <= 23) %>%# after week 23 most are missing values for delta
  group_by(region) %>%
  mutate(delta_utla_no = n()) %>%
  filter(delta_utla_no>=3) %>% 
  droplevels()

dat_final_mob$nuts <- as.numeric(as.factor(dat_final_mob$NUTS1))


dat_final_intro <- dat_intro %>%
  arrange(region, epiweek) %>%
  filter(cu_yn95 <1) %>% # first time 95% of samples become delta only
  filter(epiweek >= 66 & epiweek <= max(epiweek, na.rm = TRUE)) %>%
  mutate(epiweek = epiweek-53) %>%
  filter(epiweek <= 23) %>% # after week 23 most are missing values for delta
 group_by(region) %>%
  mutate(delta_utla_no = n()) %>%
  filter(delta_utla_no>=3) %>% 
  droplevels()

dat_final_intro$nuts <- as.numeric(as.factor(dat_final_intro$NUTS1))

# 574 observations (06 August 2021)
# 549 observations (17 August 2021)
# 528 observations (17 August 2021) # adding new_lincount to the data
# 588 observations (06 September 2021)


# 707 vac observations (14 September 2021)
# 707 mob observations (14 September 2021)
# 313 intro observations (14 September 2021)

# 599 vac observations (14 September 2021)
# 599 mob observations (14 September 2021)
# 313 intro observations (14 September 2021)

# 590 vac observations (12 October 2021)
# 590 mob observations (12 October 2021)
# 332 intro observations (12 October 2021)

# 590 vac observations (22 October 2021)
# 590 mob observations (22 October 2021)
# 299 intro observations (22 October 2021)

saveRDS(dat_final_vac, "data/processed/England_weekly_vac.rds")
saveRDS(dat_final_mob, "data/processed/England_weekly_mob.rds")
saveRDS(dat_final_intro, "data/processed/England_weekly_intro.rds")

