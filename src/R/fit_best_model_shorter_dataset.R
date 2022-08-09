library(tidyverse)
library(rstan)
library(posterior)
source("src/R/helper.R")
options(mc.cores=4)

dat_og <- readRDS("data/processed/England_weekly_processed.rds")

# determine training and testing weeks
epiweeks <- sort(unique(dat_og$epiweek))
number_testing_epiweeks <- 2
nweeks <- length(unique(dat_og$epiweek))
last_training_epiweek <- nweeks - number_testing_epiweeks
training_epiweeks <- epiweeks[1:last_training_epiweek]
testing_epiweeks <- epiweeks[(last_training_epiweek + 1):nweeks]

# select covariates to include
covariates_df <- readRDS("data/processed/covariates_to_include.rds")
covariates_to_include <- get_covariates_to_include(
  9,
  covariates_df)

stan_data <- prepare_stan_data(
  covariates_to_include, training_epiweeks,
  testing_epiweeks, last_training_epiweek,
  number_testing_epiweeks, dat_og
)

model <- stan_model("src/stan/growth_comparison.stan")
n_iters <- 60000
n_chains <- 4
n_thin <- 10
fit <- sampling(model, data=stan_data, iter=n_iters, chains=n_chains,
                thin=n_thin)

saveRDS(fit, "data/processed/stan_fits/fit_best_shorter.rds")


