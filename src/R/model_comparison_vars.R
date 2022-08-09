library(rstan)
library(tidyverse)
library(posterior)
library(loo)
source("src/R/helper.R")
options(mc.cores=4)
args <- commandArgs(trailingOnly = TRUE)

covariate_id <- as.numeric(args[1])
n_iters <- as.numeric(args[2])
n_chains <- as.numeric(args[3])
n_thins <- as.numeric(args[4])

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
  covariate_id,
  covariates_df)

stan_data <- prepare_stan_data(
  covariates_to_include, training_epiweeks,
  testing_epiweeks, last_training_epiweek,
  number_testing_epiweeks, dat_og
)

# fit Stan model
is_covariate_model <- length(covariates_to_include) > 0
if(is_covariate_model) {
  model <- stan_model("src/stan/growth_comparison.stan")
} else {
  model <- stan_model("src/stan/growth_comparison_no_covariates.stan")
}
  
fit <- sampling(model, data=stan_data,
                iter=n_iters, chains=n_chains,
                thin=n_thins)

# only save if converged
sum_df <- summarise_draws(fit)

# it's ok to remove NAs since these are
# parameter values which are set to -99
rhats <- sum(sum_df$rhat>1.01, na.rm = TRUE)
ess_bulks <- sum(sum_df$ess_bulk<400, na.rm = TRUE)
ess_tails <- sum(sum_df$ess_tail<400, na.rm = TRUE)
is_converged <- rhats == 0 & ess_bulks == 0 & ess_tails == 0
log_like <- get_loglikelihood(fit)
if(is_covariate_model) {
  beta <- rstan::extract(fit, "beta")[[1]]
} else {
  beta <- NA
}
result <- list(log_likelihood=log_like, beta=beta,
               covariates_included=covariates_to_include,
               is_converged=is_converged,
               rhats=rhats,
               ess_tails=ess_tails,
               ess_bulks=ess_bulks)

dir.create(file.path("data/processed", "stan_fits"), showWarnings = FALSE)
basefilename <- "data/processed/stan_fits/fit_model_comparison_"
filename <- paste0(basefilename, covariate_id, ".rds")

saveRDS(result, filename)
  
