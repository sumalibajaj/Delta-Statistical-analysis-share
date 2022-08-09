library(tidyverse)

dat_og <- readRDS("data/processed/England_weekly_processed.rds")

dat_og <- dat_og %>% 
  mutate(fraction_delta=delta_count / total_count)

averaged <- dat_og %>% 
  group_by(location) %>% 
  summarise(proportion_cum_2nd_dose_vaccinated=mean(proportion_cum_2nd_dose_vaccinated))

totals <- dat_og %>% 
  group_by(location) %>% 
  summarise(total_count=sum(total_count))

locations_with_high <- totals %>% 
  filter(total_count >= 500) %>% 
  pull(location)

overall_averaged <- dat_og %>% 
  filter(location %in% locations_with_high) %>% 
  ungroup() %>% 
  summarise(proportion_cum_2nd_dose_vaccinated=mean(proportion_cum_2nd_dose_vaccinated)) %>% 
  pull(proportion_cum_2nd_dose_vaccinated)

averaged <- averaged %>% 
  filter(location %in% locations_with_high) %>% 
  arrange(desc(proportion_cum_2nd_dose_vaccinated)) %>% 
  mutate(top_index=seq_along(location)) %>% 
  arrange(proportion_cum_2nd_dose_vaccinated) %>% 
  mutate(bottom_index=seq_along(location))

n_extreme <- 1:5
categorised_df <-  dat_og %>%
  select(-proportion_cum_2nd_dose_vaccinated) %>% 
  filter(location %in% locations_with_high) %>% 
  left_join(averaged) %>% 
  ungroup() %>% 
  mutate(top_5=if_else(top_index %in% n_extreme, 1, 0)) %>% 
  mutate(bottom_5=if_else(bottom_index %in% n_extreme, 1, 0)) %>% 
  mutate(other=if_else(top_5==1, "high", if_else(bottom_5==1, "low", "medium"))) %>% 
  mutate(other=as.factor(other)) %>% 
  mutate(other=fct_relevel(other, "high", "medium", "low")) %>% 
  select(epiweek, fraction_delta, other, location)

g <- categorised_df %>% 
  ggplot(aes(x=epiweek, y=fraction_delta,
             colour=other)) +
  geom_line(aes(group=location), alpha=0.3) +
  geom_smooth(se=FALSE, linetype=2) +
  scale_colour_viridis_d("Average prop.\nvaccinated", direction=-1) +
  scale_y_continuous(labels=scales::percent) +
  ylab("Delta penetration") +
  scale_x_continuous(breaks=seq(13, 23, 1)) +
  theme_classic() +
  xlab("Epi week")

ggsave("outputs/vaccination_prop_delta.pdf", width = 8,
       height = 4)
saveRDS(categorised_df, "data/processed/ranked_vaccine_2.rds")
