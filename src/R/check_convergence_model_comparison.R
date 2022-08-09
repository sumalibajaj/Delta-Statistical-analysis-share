library(tidyverse)


get_converged <- function(
    covariate_id,
    basefilename="data/processed/stan_fits/fit_model_comparison_") {
  filename <- paste0(basefilename, covariate_id, ".rds")
  temp <- readRDS(filename)
  temp$ess_tails
}

ids <- seq(1, 15, 1)
sum(map_dbl(ids, ~get_converged(.))==0)
