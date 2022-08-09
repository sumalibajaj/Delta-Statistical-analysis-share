# plots examining how the rate of Delta growth depends on relative self mobility

library(tidyverse)

dat_og <- readRDS("data/processed/England_weekly_processed.rds")

dat_og <- dat_og %>% 
  mutate(fraction_delta=delta_count / total_count)

av_self_mobility <- dat_og %>% 
  group_by(location) %>% 
  arrange(epiweek) %>% 
  summarise(relative_self_mobility_av=mean(relative_self_mobility))

totals <- dat_og %>% 
  group_by(location) %>% 
  summarise(total_count=sum(total_count))

locations_with_high <- totals %>% 
  filter(total_count >= 500) %>% 
  pull(location)

overall_av_self_mobility <- dat_og %>% 
  filter(location %in% locations_with_high) %>% 
  ungroup() %>% 
  summarise(relative_self_mobility_av=mean(relative_self_mobility)) %>% 
  pull(relative_self_mobility_av)

av_self_short <- av_self_mobility %>% 
  filter(location %in% locations_with_high) %>% 
  arrange(desc(relative_self_mobility_av)) %>% 
  mutate(top_index=seq_along(location)) %>% 
  arrange(relative_self_mobility_av) %>% 
  mutate(bottom_index=seq_along(location))

n_extreme <- 1:5
categorised_df <-  dat_og %>%
  filter(location %in% locations_with_high) %>% 
  left_join(av_self_short) %>% 
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
  scale_colour_viridis_d("Within UTLA\nmobility", direction=-1) +
  scale_y_continuous(labels=scales::percent) +
  ylab("Delta penetration") +
  scale_x_continuous(breaks=seq(13, 23, 1)) +
  theme_classic() +
  xlab("Epi week")

ggsave("outputs/relative_mobility_delta.pdf", width = 8,
       height = 4)
saveRDS(categorised_df, "data/processed/ranked_mobility.rds")

