################################################################################
################################################################################
# MODEL SETUP FOR ENGLAND

# MODELS AND THE VARIABLES USED:
# stan_dat_sm_time_plot USES MOBILITY, TIME
# fit_temp_blim_sm_plot_other USES MOBILITY, TIME AND NO.OF DELTA IMPORTS

# OUTPUT FROM stan_dat_sm_time_plot IS USED FOR FINAL PLOTS IN THE MANUSCRIPT

# FILE(S) USED:
# England_weekly_mob.rds
# England_weekly_intro.rds
# final.stan 

# OUTPUT:
# RESULTS FROM stan_dat_sm_time_plot IS USED FOR FINAL PLOTS IN THE MANUSCRIPT, 
# RUN IN plottingresults_England.R
# TO CHECK CONVERGENCE (check_convergence LIST OF THE FOLLOWING):
#   main_rhat_England, main_ess_bulk_England, main_ess_tail_England 
#   other_rhat_England, other_ess_bulk_England, other_ess_tail_England 

################################################################################
################################################################################
library(dplyr)
library(rstan)
library(reshape2)
library(loo)
library(posterior)
options(mc.cores=4)
rstan_options("auto_write")

rm(list = ls())

dat_og <- readRDS("data/processed/England_weekly_mob.rds")

n_iters = 40000
n_chains = 4
n_thins = 10

################################################################################
# HIERARCHICAL WITH VARIABLES
# NEED UPDATED MAIN REGION LEVEL DATA ON SELF MOBILITY AND IMPORTATION INTENSITY
################################################################################
  # VARIABLES (Sum is used as a function to replace NA's by zero)
  relative_self_mobility = reshape2::dcast(dat_og, epiweek ~ region, value.var="relative_self_mobility", fun.aggregate = sum) %>% select(-1)
  relative_self_mobility_long = reshape2::melt(relative_self_mobility)
  

  # LIKELIHOOD MATRIX
  include_lik = reshape2::dcast(dat_og, epiweek ~ region, value.var="delta_count") %>% select(-1)
  include_lik[!is.na(include_lik)] <- 1
  include_lik[is.na(include_lik)] <- 0
  
  # OTHER INPUTS FOR STAN
  region_vec <- dat_og %>% summarize(unique(region)) %>% select(1)
  nuts_region <- dat_og %>% group_by(nuts) %>% summarize(unique(region)) %>% rename(region = "unique(region)")
  nuts_vec <- merge(region_vec, nuts_region, by = "region")
  
  Y_delta = reshape2::dcast(dat_og, epiweek ~ region, value.var="delta_count", fun.aggregate = sum) %>% select(-1)
  
  Y_total = reshape2::dcast(dat_og, epiweek ~ region, value.var="total_count", fun.aggregate = sum) %>% select(-1)
  
  
  # relative_self_mobility and time
  model <- stan_model("src/stan/final.stan")
  time = reshape2::dcast(dat_og, epiweek ~ region, value.var="epiweek", fun.aggregate = sum) %>% select(-1)
  time_long = reshape2::melt(time)
  time_long$value <- time_long$value - mean(dat_og$epiweek)
  X_90 <- matrix(1, nrow = nrow(relative_self_mobility_long))
  stan_dat_sm_time_plot <-  list(T=length(unique(dat_og$epiweek)), I=length(unique(dat_og$region)), J=length(unique(dat_og$nuts)),
                            Y_delta=Y_delta, Y_total=Y_total, include_lik=include_lik, nuts = nuts_vec$nuts,
                            ncovs = 1, ncovs2 = 1,
                            X = as.matrix(relative_self_mobility_long$value),
                            X_90 = X_90,
                            X_time = as.matrix(time_long$value))
  fit_temp_sm_time_plot <- sampling(model, data = stan_dat_sm_time_plot, iter = n_iters, chains = n_chains, thin = n_thins)
  print(fit_temp_sm_time_plot, pars = c("beta", "beta_time", "sigma1", "sigma2"))
  sum_df <- summarise_draws(fit_temp_sm_time_plot)
  main_rhat_England <- table(sum_df$rhat<1.01)
  main_ess_bulk_England <- table(sum_df$ess_bulk>400)
  main_ess_tail_England <- table(sum_df$ess_tail>400)
  
  min((dat_og %>% group_by(region) %>%
         summarize(n=n()))$n)
  max((dat_og %>% group_by(region) %>%
         summarize(n=n()))$n)
  mean((dat_og %>% group_by(region) %>%
          summarize(n=n()))$n)

  min(dat_og$delta_count/dat_og$total_count*100)
  max(dat_og$delta_count/dat_og$total_count*100)
  mean(dat_og$delta_count/dat_og$total_count*100)  

  min(dat_og$delta_count)
  max(dat_og$delta_count)
  mean(dat_og$delta_count)
  
  min(dat_og$total_count)
  max(dat_og$total_count)
  mean(dat_og$total_count)
  
  min(dat_og$relative_self_mobility)
  max(dat_og$relative_self_mobility)
  mean(dat_og$relative_self_mobility)
  
  
  ################################################################################
  ################################################################################
  # DATA WITH IMPORTATION (REDUCES TO ALMOST HALF)
  ################################################################################
  ################################################################################
  
  dat_og <- readRDS("data/processed/England_weekly_intro.rds") 
  
  relative_self_mobility = reshape2::dcast(dat_og, epiweek ~ region, value.var="relative_self_mobility", fun.aggregate = sum) %>% select(-1)
  relative_self_mobility_long = reshape2::melt(relative_self_mobility)
  
  number_introduction = reshape2::dcast(dat_og, epiweek ~ region, value.var="number_introduction", fun.aggregate = sum) %>% select(-1)
  number_introduction_long = reshape2::melt(number_introduction)
  # number_introduction_long$value <- (number_introduction_long$value - mean(dat_og$number_introduction))/sd(dat_og$number_introduction)
  number_introduction_long$value <- sqrt(number_introduction_long$value)
  
  
  # LIKELIHOOD MATRIX
  include_lik = reshape2::dcast(dat_og, epiweek ~ region, value.var="delta_count") %>% select(-1)
  include_lik[!is.na(include_lik)] <- 1
  include_lik[is.na(include_lik)] <- 0
  
  # OTHER INPUTS FOR STAN
  region_vec <- dat_og %>% summarize(unique(region)) %>% select(1)
  nuts_region <- dat_og %>% group_by(nuts) %>% summarize(unique(region)) %>% rename(region = "unique(region)")
  nuts_vec <- merge(region_vec, nuts_region, by = "region")
  
  Y_delta = reshape2::dcast(dat_og, epiweek ~ region, value.var="delta_count", fun.aggregate = sum) %>% select(-1)
  
  Y_total = reshape2::dcast(dat_og, epiweek ~ region, value.var="total_count", fun.aggregate = sum) %>% select(-1)
  
  
  # relative_self_mobility, intro and time
  model <- stan_model("src/stan/final.stan")
  time = reshape2::dcast(dat_og, epiweek ~ region, value.var="epiweek", fun.aggregate = sum) %>% select(-1)
  time_long = reshape2::melt(time)
  time_long$value <- time_long$value - mean(dat_og$epiweek)
  X_90 <- matrix(1, nrow = nrow(relative_self_mobility_long))
  stan_dat_sm_time_plot_other <-  list(T=length(unique(dat_og$epiweek)), I=length(unique(dat_og$region)), J=length(unique(dat_og$nuts)),
                                 Y_delta=Y_delta, Y_total=Y_total, include_lik=include_lik, nuts = nuts_vec$nuts,
                                 ncovs = 2, ncovs2 = 1,
                                 X = as.matrix(cbind(relative_self_mobility_long$value, number_introduction_long$value)),
                                 X_90 = cbind(X_90, X_90),
                                 X_time = as.matrix(time_long$value))
  
  fit_temp_sm_time_plot_other <- sampling(model, data = stan_dat_sm_time_plot_other, iter = n_iters, chains = n_chains, thin = n_thins)
  print(fit_temp_sm_time_plot_other, pars = c("beta", "beta_time", "sigma1", "sigma2"))
  sum_df_other <- summarise_draws(fit_temp_sm_time_plot_other)
  other_rhat_England <- table(sum_df_other$rhat<1.01)
  other_ess_bulk_England <- table(sum_df_other$ess_bulk>400)
  other_ess_tail_England <- table(sum_df_other$ess_tail>400)
  
  t <- dat_og %>% group_by(region) %>% summarize(t = n())
  mean(t$t)
  length(unique(dat_og$region))
  
  # This is the final dataset used in the model
  dat_og <- readRDS("data/processed/England_weekly_mob.rds")
  
  check_convergence <- list(main_rhat_England,main_ess_tail_England,main_ess_bulk_England,
                            other_rhat_England,other_ess_tail_England,other_ess_bulk_England)
  
  ################################################################################
  save(check_convergence, file = "data/processed/check_convergence_England.RData")
  save.image("data/processed/England_final.RData")
  ################################################################################
  