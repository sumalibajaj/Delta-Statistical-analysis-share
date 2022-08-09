library(tidyverse)
library(lubridate)

dat_final_short <- readRDS("data/processed/England_weekly_cleaned.rds")

# look at cumulative vaccinations: from
# https://api.coronavirus.data.gov.uk/v2/data?areaType=utla&metric=cumPeopleVaccinatedSecondDoseByVaccinationDate&format=csv
vaccs <- read.csv("data/raw/utla_2022-04-13.csv") %>% 
  rename(location=areaName) %>% 
  rename(region=areaCode) %>% 
  rename(cum_2nd_dose_vaccinated=cumPeopleVaccinatedSecondDoseByVaccinationDate) %>% 
  rename(cum_1st_dose_vaccinated=cumPeopleVaccinatedFirstDoseByVaccinationDate) %>% 
  mutate(date=as_date(date)) %>% 
  filter(date >= "2021-01-03") %>% 
  filter(date <= "2021-12-31") %>% 
  mutate(epiweek=epiweek(date)) %>%
  arrange(region, date) %>% 
  group_by(location, region, epiweek) %>% 
  summarise(cum_2nd_dose_vaccinated=mean(cum_2nd_dose_vaccinated, rm.na=TRUE),
            cum_1st_dose_vaccinated=mean(cum_1st_dose_vaccinated, rm.na=TRUE))

# group UTLAs to match genome regions
greater_london <- c('E09000007', 'E09000011','E09000012',
                    'E09000013','E09000019','E09000020',
                    'E09000022','E09000023', 'E09000028',
                    'E09000030', 'E09000032', 'E09000033',
                    'E09000001', 'E09000002', 'E09000003',
                    'E09000004', 'E09000005', 'E09000006',
                    'E09000008', 'E09000009', 'E09000010',
                    'E09000014', 'E09000015', 'E09000016',
                    'E09000017', 'E09000018',  'E09000021',
                    'E09000027', 'E09000024', 'E09000025',
                    'E09000026', 'E09000029', 'E09000031')
greater_london <- tibble(region=greater_london,
                         new_region="E13000001|E13000002",
                         location="Greater London")
merseyside <- c('E08000012', 'E08000014', 'E08000011', 'E08000013', 'E08000015')
merseyside <- tibble(region=merseyside,
                     new_region="E11000002",
                     location="Merseyside")
tyneandweir <- c('E08000037','E08000021', 'E08000022', 'E08000023', 'E08000024')
tyneandweir <- tibble(region=tyneandweir,
                      new_region="E11000007",
                      location="Tyne and Wear")
greatermanchester <- c('E08000003', 'E08000007', 'E08000008',
                       'E08000004', 'E08000005', 'E08000002',
                       'E08000001', 'E08000010', 'E08000006',
                       'E08000009')
greatermanchester <- tibble(region=greatermanchester,
                            new_region="E11000001",
                            location="Greater Manchester")
westyorkshire <- c('E08000035', 'E08000036', 'E08000034', 'E08000033', 'E08000032')
westyorkshire <- tibble(region=westyorkshire,
                        new_region="E11000006",
                        location="WEST YORKSHIRE")
southyorkshire <- c('E08000019', 'E08000018', 'E08000017', 'E08000016')
southyorkshire <- tibble(region=southyorkshire,
                         new_region="E11000003",
                         location="SOUTH YORKSHIRE")
westmidlands <- c('E08000026', 'E08000029', 'E08000025',
                  'E08000028', 'E08000030', 'E08000027',
                  'E08000031')
westmidlands <- tibble(region=westmidlands,
                       new_region="E11000005",
                       location="WEST MIDLANDS")
new_categorisations <- greater_london %>% 
  bind_rows(merseyside) %>% 
  bind_rows(tyneandweir) %>% 
  bind_rows(greatermanchester) %>% 
  bind_rows(westyorkshire) %>% 
  bind_rows(southyorkshire) %>% 
  bind_rows(westmidlands) %>% 
  rename(new_location=location)

vaccs <- vaccs %>% 
  left_join(new_categorisations) %>% 
  mutate(region=if_else(is.na(new_region), region, new_region)) %>% 
  mutate(location=if_else(is.na(new_location), location, new_location))

# clean area codes as there are some differences in region names in vaccines data
area_lookup <- dat_final_short %>% 
  select(region, location) %>% 
  unique() %>% 
  arrange(location)
cornwall <- vaccs %>%
  filter(str_detect(location, "Cornwall")) %>%
  mutate(location="Cornwall")
vaccs <- vaccs %>%
  filter(!str_detect(location, "Cornwall")) %>%
  bind_rows(cornwall)
vaccs <- vaccs %>% 
  ungroup() %>% 
  select(-region) %>% 
  left_join(area_lookup)

dat_with_vaccs <- dat_final_short %>% 
  full_join(vaccs) %>% 
  select(-one_of("prop_cu_1vac_utla_week", "prop_cu_2vac_utla_week")) %>% 
  arrange(region, epiweek) %>% 
  filter(region %in% dat_final_short$region)

# look at unmatched vaccination places
test <- dat_with_vaccs[is.na(dat_with_vaccs$cum_1st_dose_vaccinated), ]
if(nrow(test) > 0)
  stop("There are remaining unmatched regions.")

dat_with_vaccs <- dat_with_vaccs %>% 
  select(-new_region) %>% 
  select(-new_location)

# check that those regions that were aggregated
# had zero sd in relative self mobility (since this
# variable was created at this grouping)
test <- dat_with_vaccs %>% 
  group_by(epiweek, region, location) %>% 
  summarise(relative_self_mobility_sd=var(relative_self_mobility))
if(max(test$relative_self_mobility_sd, na.rm = T) > 0)
  stop("Relative self mobility has variance when it shouldn't.")

# for those regions that are grouped together,
# aggregate
aggd <- dat_with_vaccs %>% 
  group_by(epiweek, region, location, nuts, NUTS1) %>% 
  summarise(
    population=sum(population),
    delta_count=sum(delta_count),
    total_count=sum(total_count),
    new_cases_reported=sum(new_cases_reported),
    relative_self_mobility=mean(relative_self_mobility),
    number_introduction=sum(number_introduction),
    number_introduction_lag1=sum(number_introduction_lag1),
    number_introduction_lag2=sum(number_introduction_lag2),
    cum_1st_dose_vaccinated=sum(cum_1st_dose_vaccinated),
    cum_2nd_dose_vaccinated=sum(cum_2nd_dose_vaccinated)
  ) %>% 
  ungroup() %>% 
  mutate(proportion_cum_1st_dose_vaccinated=cum_1st_dose_vaccinated / population,
         proportion_cum_2nd_dose_vaccinated=cum_2nd_dose_vaccinated / population) %>% 
  filter(epiweek >= 13) %>% 
  filter(epiweek <= 23) %>% 
  arrange(region, epiweek)

# check vaccine coverages
a_max <- max(aggd$proportion_cum_1st_dose_vaccinated)
a_min <- min(aggd$proportion_cum_1st_dose_vaccinated)
if(a_max > 1 || a_min < 0)
  stop("Vaccine coverages not correct.")
a_max <- max(aggd$proportion_cum_2nd_dose_vaccinated)
a_min <- min(aggd$proportion_cum_2nd_dose_vaccinated)
if(a_max > 1 || a_min < 0)
  stop("Vaccine coverages not correct.")

saveRDS(aggd, "data/processed/England_weekly_processed.rds")
