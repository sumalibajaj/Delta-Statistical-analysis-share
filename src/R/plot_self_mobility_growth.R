library(tidyverse)

dat_og <- readRDS("data/processed/England_weekly_processed.rds")

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
  ggplot(aes(x=relative_self_mobility, y=growth)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~new_name, scales="free") +
  xlab("Within UTLA mobility") +
  ylab("Delta growth, weekly change in proportion")
ggsave("outputs/mobility_growth.pdf", g, width = 12, height = 8)

# look at models
coefs <- dat_growth %>% 
  group_by(new_name) %>% 
  group_map(~coef(lm(growth ~ relative_self_mobility, data = .))[2])
coefs <- unlist(coefs)
min(coefs)
max(coefs)

models <- dat_growth %>% 
  group_by(new_name) %>% 
  group_map(~lm(growth ~ relative_self_mobility, data = .))
