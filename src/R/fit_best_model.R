library(tidyverse)
library(rstan)
library(posterior)
source("src/R/helper.R")
options(mc.cores=4)

dat_og <- readRDS("data/processed/England_weekly_processed.rds")

# pick the model with the best parameters to run
model_comparison <- readRDS("data/processed/model_comparison_results.rds")
covariates_df <- readRDS("data/processed/covariates_to_include.rds")

best_index <- model_comparison$id[which.max(model_comparison$log_likelihood)]
covariates_to_include <- get_covariates_to_include(
  best_index, covariates_df)

# fit model
stan_data <- prepare_stan_data_full(dat_og)

model <- stan_model("src/stan/growth_full_simple.stan")

n_iters <- 60000
n_chains <- 4
n_thin <- 10
fit <- sampling(model, data=stan_data, iter=n_iters, chains=n_chains,
                thin=n_thin)

# only save if converged
sum_df <- summarise_draws(fit)

# it's ok to remove NAs since these are
# parameter values which are set to -99
rhats <- sum(sum_df$rhat>1.01, na.rm = TRUE)
ess_bulks <- sum(sum_df$ess_bulk<400, na.rm = TRUE)
ess_tails <- sum(sum_df$ess_tail<400, na.rm = TRUE)
is_converged <- rhats == 0 & ess_bulks == 0 & ess_tails == 0

if(is_converged) {
  saveRDS(fit, "data/processed/stan_fits/fit_best.rds")
} else {
  print("Not converged.")
  print(paste0("rhat = ", rhats))
  print(paste0("ess_bulk = ", ess_bulks))
  print(paste0("ess_tail = ", ess_tails))
}

