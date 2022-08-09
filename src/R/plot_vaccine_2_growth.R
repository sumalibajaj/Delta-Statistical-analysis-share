library(tidyverse)

dat_og <- readRDS("data/processed/England_weekly_processed.rds") %>%
  group_by(region) %>% 
  mutate(proportion_cum_2nd_dose_vaccinated=mean(proportion_cum_2nd_dose_vaccinated))

dat_og <- dat_og %>% 
  mutate(fraction_delta=delta_count / total_count)

dat_growth <- dat_og %>% 
  group_by(location) %>% 
  arrange(epiweek) %>% 
  mutate(growth=c(NA, diff(fraction_delta)))

lookup <- tibble(NUTS1=unique(dat_growth$NUTS1)) %>% 
  mutate(new_name=gsub("_", " ", NUTS1))

dat_growth <- dat_growth %>% 
  left_join(lookup)

g <- dat_growth %>% 
  ggplot(aes(x=proportion_cum_2nd_dose_vaccinated, y=growth)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~new_name, scales="free") +
  xlab("Proportion vaccinated (averaged)") +
  ylab("Delta growth, weekly change in proportion")
ggsave("outputs/vaccines_2_growth.pdf", g, width = 12, height = 8)
