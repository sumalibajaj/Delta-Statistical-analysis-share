library(tidyverse)

df <- readRDS("data/processed/England_weekly_processed.rds") %>% 
  mutate(proportion_delta=delta_count / total_count)

regions <- unique(df$region)

g <- df %>% 
  na.omit() %>% 
  filter(region %in% regions[1:16]) %>% 
  ggplot(aes(x=epiweek, y=proportion_delta)) +
  geom_point() +
  geom_line() +
  geom_line(aes(y=proportion_cum_2nd_dose_vaccinated),
            colour="orange") +
  facet_wrap(~location)

ggsave("outputs/cum_2_dose_vs_delta.pdf", g,
       width = 10, height = 5)  

g <- df %>% 
  na.omit() %>% 
  group_by(region) %>% 
  mutate(relative_introduction=number_introduction / max(number_introduction)) %>% 
  filter(region %in% regions[1:16]) %>% 
  ggplot(aes(x=epiweek, y=proportion_delta)) +
  geom_point() +
  geom_line() +
  geom_line(aes(y=relative_introduction),
            colour="orange") +
  facet_wrap(~location)

ggsave("outputs/number_introduction_vs_delta.pdf", g,
       width = 10, height = 5)
