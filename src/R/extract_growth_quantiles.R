library(tidyverse)
source("src/R/helper.R")

dat_og <- readRDS("data/processed/England_weekly_processed.rds")

# pick the model with the best parameters to run
model_comparison <- readRDS("data/processed/model_comparison_results.rds")
covariates_df <- readRDS("data/processed/covariates_to_include.rds")

best_index <- model_comparison$id[which.max(model_comparison$log_likelihood)]
covariates_to_include <- get_covariates_to_include(
  best_index, covariates_df)

# get fit
stan_data <- prepare_stan_data_full(dat_og)
fit <- readRDS("data/processed/stan_fits/fit_best.rds")

# get rhos
convert_quantile <- function(q, rho, count_delta, dat_og) {
  rho_middle <- apply(rho, c(2, 3), function(x) quantile(x, q))
  count_delta <- stan_data$Y_delta
  epiweeks <- sort(unique(dat_og$epiweek))
  rho_middle <- rho_middle %>% 
    as.data.frame()
  colnames(rho_middle) <- colnames(count_delta)
  rho_middle <- rho_middle %>% 
    mutate(epiweek=epiweeks) %>% 
    pivot_longer(-epiweek)
  rho_middle
}
rho <- rstan::extract(fit, "rho")[[1]]
rho_middle <- convert_quantile(0.5, rho, count_delta, dat_og) %>% 
  rename(middle=value)
rho_lower <- convert_quantile(0.025, rho, count_delta, dat_og) %>% 
  rename(lower=value)
rho_upper <- convert_quantile(0.975, rho, count_delta, dat_og) %>% 
  rename(upper=value)
rho_all <- rho_middle %>% 
  left_join(rho_lower) %>% 
  left_join(rho_upper) %>% 
  rename(region=name) %>% 
  left_join(dat_og)

saveRDS(rho_all, "data/processed/growth_quantiles.rds")
