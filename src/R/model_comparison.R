library(rstan)
library(loo)
library(tidyverse)
source("src/R/helper.R")

get_loglikelihood <- function(
    covariate_id, basefilename="data/processed/stan_fits/fit_model_comparison_",
    missing_data_code=-99) {
  filename <- paste0(basefilename, covariate_id, ".rds")
  temp <- readRDS(filename)
  log_like <- temp$log_likelihood
  overall_fit <- sum(colMeans(log_like))
  overall_fit
}

covariates_df <- readRDS("data/processed/covariates_to_include.rds")
log_likes <- map_dbl(covariates_df$id, get_loglikelihood)
covariates_df$log_likelihood <- log_likes
covars_list <- map(covariates_df$id, ~get_covariates_to_include(., covariates_df))
covars_list[[1]] <- "No covars"
covars_list <- map_chr(covars_list, ~paste(., collapse = ", "))
covariates_df$included <- covars_list
covariates_df <- covariates_df %>% 
  arrange(desc(log_likelihood))
saveRDS(covariates_df, "data/processed/model_comparison_results.rds")
short <- covariates_df %>% 
  select(included, log_likelihood) %>% 
  mutate(log_likelihood=round(log_likelihood)) %>% 
  rename(`log-likelihood`=log_likelihood) %>%
  rename(covariates=included) %>% 
  mutate(covariates=as.character(covariates)) %>% 
  mutate(covariates=gsub("relative_self_mobility",
                         "within UTLA mobility", covariates)) %>% 
  mutate(covariates=gsub("proportion_cum_2nd_dose_vaccinated",
                         "cumulative prop. 2nd dose vacc.", covariates)) %>%
  mutate(covariates=gsub("proportion_cum_1st_dose_vaccinated",
                         "cumulative prop. 1st dose vacc.", covariates)) %>% 
  mutate(covariates=gsub("number_introduction_lag2",
                         "introductions (2 week lag)", covariates)) %>% 
  mutate(covariates=gsub("number_introduction_lag1",
                         "introductions (1 week lag)", covariates)) %>% 
  mutate(covariates=gsub("number_introduction",
                         "introductions", covariates))
write.csv(short, "data/processed/model_comparison_results.csv",
          row.names = FALSE)

