library(rstan)
library(tidyverse)

fit <- readRDS("data/processed/stan_fits/fit_best_shorter.rds")
beta <- rstan::extract(fit, "beta")[[1]]
beta_sum <- apply(beta, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
qs <- row.names(beta_sum)
colnames(beta_sum) <- c("epiweek", "within UTLA mobility")
beta_sum <- beta_sum %>% 
  as.data.frame() %>% 
  mutate(quantile=qs) %>% 
  mutate(type="epiweeks 13-21")
beta_short <- beta_sum
fit <- readRDS("data/processed/stan_fits/fit_best.rds")
beta <- rstan::extract(fit, "beta")[[1]]
beta_sum <- apply(beta, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
colnames(beta_sum) <- c("epiweek", "within UTLA mobility")
beta_sum <- beta_sum %>% 
  as.data.frame() %>% 
  mutate(quantile=qs) %>% 
  mutate(type="epiweeks 13-23")

beta_all <- beta_short %>% 
  bind_rows(beta_sum)
row.names(beta_all) <- seq_along(beta_all$epiweek)
beta_long <- beta_all %>% 
  pivot_longer(-c(quantile, type)) %>%
  mutate(value=round(value, 2)) %>% 
  group_by(type, name) %>% 
  summarise(estimate=paste0(value[2],
                            " (", value[1], ", ", value[3], ")")) %>% 
  rename(model=type)

write.csv(beta_long, "data/processed/beta_estimates.csv",
          row.names = FALSE)
