library(dplyr)
library(lubridate) # to use epiweek function
library(ggplot2)

################################################################################
# USA state level data and covariates
# First keep only UTLA level UK data
# Calculate weekly and regional population, vaccinations (1st and 2nd prop)
# Calculate the delta percentage for every week and each UTLA
# Keep data from week 64 (starting date = 14 Mar 2021), and subtract 53 from all epiweek to only 2021 data
# removing weeks following the first time 95% delta samples of total is reached
# Keep only those UTLAs which have data on delta for (doesn't matter if it's zero samples) more than three weeks
################################################################################
############
# DONT USE TOTAL IMPOTATION INTENSITY IN V8 OF THE DATA
############

rm(list = ls())

# dat <- read.csv(file = "data/raw/master_with_number_introduction_and_lag.csv")
dat <- read.csv(file = "data/raw/master_total_updated.csv")

# To find earliest week with a delta sample
dat_min <-  dat %>%
            filter(level == "US_state") %>%
            group_by(code) %>%
            arrange(code, epiweek) %>%
            mutate(delta_yn = ifelse(delta_count >0 , 1, 0),
                   delta_yn_cu = cumsum(delta_yn)) %>%
            filter(delta_yn_cu == 1 & delta_yn == 1) 


dat_new <-  dat %>%
    filter(total_count > 0) %>%
    filter(level == "US_state") %>%
    mutate(prop_delta = delta_count/total_count,
           yn95 = ifelse(prop_delta >= 0.95, 1, 0)) %>%
    rename(region = code) %>%
    group_by(region) %>%
    arrange(region, epiweek) %>%
    mutate(cu_yn95 = cumsum(yn95),
           prop_cu_1vac_utla_week = cumsum(proportion_population_new_1st_dose),
           prop_cu_2vac_utla_week = cumsum(proportion_population_new_2nd_dose),
           prop_cu_2vac_utla_week_lag1 = lag(prop_cu_2vac_utla_week),
           prop_cu_2vac_utla_week_lag2 = lag(prop_cu_2vac_utla_week, n=2)) %>%
    select(prop_delta, yn95, cu_yn95, epiweek, region, population, delta_count, total_count, 
           new_cases_reported, relative_self_mobility,
           proportion_population_new_1st_dose, proportion_population_new_2nd_dose,
           prop_cu_1vac_utla_week, prop_cu_2vac_utla_week, 
           new_1st_dose_vaccine, new_2nd_dose_vaccine,
           NUTS1) %>%
    na.omit()

# cumsum(ifelse(is.na(x), 0, x)) + x*0
dat_blimmunity <- dat %>%
  filter(level == "US_state", epiweek == 63) %>%
  select(bl_immunity = immunity_level, region = code)


dat_blsusc <- dat %>%
  filter(level == "US_state") %>%
  filter(!is.na(susceptibility)) %>%
  group_by(code) %>%
  filter(epiweek == max(epiweek)) %>%
  mutate(immunity1 = 1 - susceptibility) %>%
  select(code, bl_immunity1 = immunity1, bl_susc = susceptibility, region = code)


dat_final <- dat_new %>%
  arrange(region, epiweek) %>%
  filter(cu_yn95 <1) %>% # first time 95% of samples become delta only
  filter(epiweek >= 64 & epiweek <= max(epiweek, na.rm = TRUE)) %>%
  mutate(epiweek = epiweek-53) %>%
  filter(epiweek <= 25) %>% # after week 25 most are missing values for delta
  group_by(region) %>%
  mutate(delta_utla_no = n()) %>%
  filter(delta_utla_no>3) %>% 
  droplevels()
  

dat_final$nuts <- as.numeric(as.factor(dat_final$NUTS1))
dat_final <- merge(dat_final, dat_blimmunity, by = "region")
# dat_final <- merge(dat_final, dat_blvac, by = "region")
dat_final <- merge(dat_final, dat_blsusc, by = "region")

dat_final <- dat_final %>% 
  group_by(region) %>%
  arrange(region, epiweek) %>%
  mutate(prop_cu_2vac_utla_week_lag1 = lag(prop_cu_2vac_utla_week, default = NA),
         prop_cu_2vac_utla_week_lag2 = lag(prop_cu_2vac_utla_week, n=2, default = NA),
         prop_cu_2vac_utla_week_lag4 = lag(prop_cu_2vac_utla_week, n=4, default = NA)) 


# 739 observations (09 August 2021)
# 735observations (10 August 2021)
# 732 observations (12 August 2021)
# 735 observations (16 August 2021)
# 735 observations (19 August 2021)
# 735 observations (25 August 2021 v4) 
# 765 observation (01 September 2021 v8)
# 735 observation (02 September 2021 v9)
# 735 observation (07 September 2021 v9)

# 735 observation (12 October 2021 master)
# 765 observation (20 October 2021 master updated)
# 735 observation (22 October 2021 master updated and total count corrected)


# saveRDS(dat_final, "data/processed/US_weekly_old.rds")
saveRDS(dat_final, "data/processed/US_weekly.rds")

################################################################################
################################################################################
# DATA WITH IMPORTATION (REDUCES TO ALMOST HALF)
################################################################################
################################################################################

rm(list = ls())

# dat <- read.csv(file = "data/raw/master_with_number_introduction_and_lag.csv")
dat <- read.csv(file = "data/raw/master_total_updated.csv")

# To find earliest week with a delta sample
dat_min <-  dat %>%
  filter(level == "US_state") %>%
  group_by(code) %>%
  arrange(code, epiweek) %>%
  mutate(delta_yn = ifelse(delta_count >0 , 1, 0),
         delta_yn_cu = cumsum(delta_yn)) %>%
  filter(delta_yn_cu == 1 & delta_yn == 1) 


dat_new <-  dat %>%
  filter(total_count > 0) %>%
  filter(level == "US_state") %>%
  mutate(prop_delta = delta_count/total_count,
         yn95 = ifelse(prop_delta >= 0.95, 1, 0)) %>%
  rename(region = code) %>%
  group_by(region) %>%
  arrange(region, epiweek) %>%
  mutate(cu_yn95 = cumsum(yn95),
         prop_cu_1vac_utla_week = cumsum(proportion_population_new_1st_dose),
         prop_cu_2vac_utla_week = cumsum(proportion_population_new_2nd_dose),
         prop_cu_2vac_utla_week_lag1 = lag(prop_cu_2vac_utla_week),
         prop_cu_2vac_utla_week_lag2 = lag(prop_cu_2vac_utla_week, n=2)) %>%
  select(cu_yn95, epiweek, region, population, delta_count, total_count, 
         new_cases_reported, relative_self_mobility,
         proportion_population_new_1st_dose, proportion_population_new_2nd_dose,
         prop_cu_1vac_utla_week, prop_cu_2vac_utla_week, 
         new_1st_dose_vaccine, new_2nd_dose_vaccine,
         total_importation_intensity,
         NUTS1) %>%
  na.omit()


dat_blimmunity <- dat %>%
  filter(level == "US_state", epiweek == 63) %>%
  select(bl_immunity = immunity_level, region = code)


dat_blsusc <- dat %>%
  filter(level == "US_state") %>%
  filter(!is.na(susceptibility)) %>%
  group_by(code) %>%
  filter(epiweek == max(epiweek)) %>%
  mutate(immunity1 = 1 - susceptibility) %>%
  select(code, bl_immunity1 = immunity1, bl_susc = susceptibility, region = code)


dat_final <- dat_new %>%
  arrange(region, epiweek) %>%
  filter(cu_yn95 <1) %>% # first time 95% of samples become delta only
  filter(epiweek >= 64 & epiweek <= max(epiweek, na.rm = TRUE)) %>%
  mutate(epiweek = epiweek-53) %>%
  filter(epiweek <= 25) %>% # after week 25 most are missing values for delta
  group_by(region) %>%
  mutate(delta_utla_no = n()) %>%
  filter(delta_utla_no>3) %>% 
  droplevels()


dat_final$nuts <- as.numeric(as.factor(dat_final$NUTS1))
dat_final <- merge(dat_final, dat_blimmunity, by = "region")
# dat_final <- merge(dat_final, dat_blvac, by = "region")
dat_final <- merge(dat_final, dat_blsusc, by = "region")

dat_final <- dat_final %>% 
  group_by(region) %>%
  arrange(region, epiweek) %>%
  mutate(prop_cu_2vac_utla_week_lag1 = lag(prop_cu_2vac_utla_week, default = NA),
         prop_cu_2vac_utla_week_lag2 = lag(prop_cu_2vac_utla_week, n=2, default = NA),
         prop_cu_2vac_utla_week_lag4 = lag(prop_cu_2vac_utla_week, n=4, default = NA)) 

# 389 observation (05 October 2021 master)
# 387 observation (12 October 2021 master)
# 396 observation (12 October 2021 master updated)
# 387 observation (22 October 2021 master updated and corrected)

saveRDS(dat_final, "data/processed/US_weekly_intro.rds")
  
