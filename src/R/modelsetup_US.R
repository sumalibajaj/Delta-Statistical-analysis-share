################################################################################
################################################################################
# MODEL SETUP FOR US

# MODELS AND THE VARIABLES USED:
# fit_temp_blim_sm_plot USES BASELINE IMMUNITY (VAC+INFECTION), MOBILITY, TIME
# fit_temp_blim_sm_plot_other USES BASELINE IMMUNITY (INFECTION), MOBILITY, TIME
# fit_temp_blim_sm_plot_other1 USES BASELINE IMMUNITY (VAC+INFECTION), MOBILITY, 
#                                   TIME AND IMPORTATION INTENSITY OF DELTA

# FILE(S) USED:
# US_weekly.rds
# US_weekly_intro.rds
# final.stan 

# OUTPUT:
# RESULTS FROM fit_temp_blim_sm_plot IS USED FOR FINAL PLOTS IN THE MANUSCRIPT, 
# RUN IN plottingresults_US.R
# TO CHECK CONVERGENCE (check_convergence LIST OF THE FOLLOWING):
#   main_rhat_US, main_ess_bulk_US, main_ess_tail_US 
#   other_rhat_US, other_ess_bulk_US, other_ess_tail_US 
#   other1_rhat_US, other1_ess_bulk_US, other1_ess_tail_US 

################################################################################
################################################################################
library(dplyr)
library(rstan)
library(reshape2)
library(loo)
library(posterior)
options(mc.cores=4)

rm(list = ls())

dat_og <- readRDS("data/processed/US_weekly.rds")

n_iters = 20000
n_chains = 4
n_thins = 10

  # LIKELIHOOD MATRIX
  include_lik = reshape2::dcast(dat_og, epiweek ~ region, value.var="delta_count") %>% select(-1)
  include_lik[!is.na(include_lik)] <- 1
  include_lik[is.na(include_lik)] <- 0
  
  region_vec <- dat_og %>% summarize(unique(region)) %>% select(1)
  nuts_region <- dat_og %>% group_by(nuts) %>% summarize(unique(region)) %>% rename(region = "unique(region)")
  nuts_vec <- merge(region_vec, nuts_region, by = "region")
  
  Y_delta = reshape2::dcast(dat_og, epiweek ~ region, value.var="delta_count", fun.aggregate = sum) %>% select(-1)
  
  Y_total = reshape2::dcast(dat_og, epiweek ~ region, value.var="total_count", fun.aggregate = sum) %>% select(-1)
  
  model <- stan_model("src/stan/final.stan")
  blim = reshape2::dcast(dat_og, epiweek ~ region, value.var="bl_immunity", fun.aggregate = sum) %>% select(-1)
  blim_long = reshape2::melt(blim)
  sm = reshape2::dcast(dat_og, epiweek ~ region, value.var="relative_self_mobility", fun.aggregate = sum) %>% select(-1)
  sm_long = reshape2::melt(sm)
  time = reshape2::dcast(dat_og, epiweek ~ region, value.var="epiweek", fun.aggregate = sum) %>% select(-1)
  time_long = reshape2::melt(time)
  time_long$value <- time_long$value - mean(unique(dat_og$epiweek))
  blim_90 <- matrix(0.9, nrow = nrow(blim_long))
  stan_dat_blim_sm_plot <-  list(T=length(unique(dat_og$epiweek)), I=length(unique(dat_og$region)), J=length(unique(dat_og$nuts)),
                                     Y_delta=Y_delta, Y_total=Y_total, include_lik=include_lik, nuts = nuts_vec$nuts,
                                     ncovs = 1, ncovs2 = 2,
                                     X = as.matrix(blim_long$value), X_90 = blim_90,
                                     X_time = as.matrix(cbind(sm_long$value, time_long$value)))

  fit_temp_blim_sm_plot <- sampling(model, data = stan_dat_blim_sm_plot, iter = n_iters, chains = n_chains, thin = n_thins)
  print(fit_temp_blim_sm_plot, pars = c("beta", "beta_time", "sigma1", "sigma2"))
  sum_df <- summarise_draws(fit_temp_blim_sm_plot)
  main_rhat_US <- table(sum_df$rhat<1.01)
  main_ess_bulk_US <- table(sum_df$ess_bulk>400)
  main_ess_tail_US <- table(sum_df$ess_tail>400)
  
  
  blim = reshape2::dcast(dat_og, epiweek ~ region, value.var="bl_immunity1", fun.aggregate = sum) %>% select(-1)
  blim_long = reshape2::melt(blim)
  stan_dat_blim_sm_plot_other <-  list(T=length(unique(dat_og$epiweek)), I=length(unique(dat_og$region)), J=length(unique(dat_og$nuts)),
                                     Y_delta=Y_delta, Y_total=Y_total, include_lik=include_lik, nuts = nuts_vec$nuts,
                                     ncovs = 1, ncovs2 = 2,
                                     X = as.matrix(blim_long$value), X_90 = blim_90,
                                     X_time = as.matrix(cbind(sm_long$value, time_long$value)))
  fit_temp_blim_sm_plot_other <- sampling(model, data = stan_dat_blim_sm_plot_other, iter = n_iters, chains = n_chains, thin = n_thins)
  print(fit_temp_blim_sm_plot_other, pars = c("beta", "beta_time", "sigma1", "sigma2"))
  sum_df_other <- summarise_draws(fit_temp_blim_sm_plot_other)
  other_rhat_US <- table(sum_df_other$rhat<1.01)
  other_ess_bulk_US <- table(sum_df_other$ess_bulk>400)
  other_ess_tail_US <- table(sum_df_other$ess_tail>400)
  
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
  
  min(dat_og$bl_immunity)
  max(dat_og$bl_immunity)
  mean(dat_og$bl_immunity)
  
  ################################################################################
  ################################################################################
  # DATA WITH IMPORTATION (REDUCES TO ALMOST HALF)
  ################################################################################
  ################################################################################
  
  dat_og <- readRDS("data/processed/US_weekly_intro.rds")
  
  include_lik = reshape2::dcast(dat_og, epiweek ~ region, value.var="delta_count") %>% select(-1)
  include_lik[!is.na(include_lik)] <- 1
  include_lik[is.na(include_lik)] <- 0
  
  region_vec <- dat_og %>% summarize(unique(region)) %>% select(1)
  nuts_region <- dat_og %>% group_by(nuts) %>% summarize(unique(region)) %>% rename(region = "unique(region)")
  nuts_vec <- merge(region_vec, nuts_region, by = "region")
  
  Y_delta = reshape2::dcast(dat_og, epiweek ~ region, value.var="delta_count", fun.aggregate = sum) %>% select(-1)
  
  Y_total = reshape2::dcast(dat_og, epiweek ~ region, value.var="total_count", fun.aggregate = sum) %>% select(-1)
  
  model <- stan_model("src/stan/final.stan")
  blim = reshape2::dcast(dat_og, epiweek ~ region, value.var="bl_immunity", fun.aggregate = sum) %>% select(-1)
  blim_long = reshape2::melt(blim)
  sm = reshape2::dcast(dat_og, epiweek ~ region, value.var="relative_self_mobility", fun.aggregate = sum) %>% select(-1)
  sm_long = reshape2::melt(sm)
  time = reshape2::dcast(dat_og, epiweek ~ region, value.var="epiweek", fun.aggregate = sum) %>% select(-1)
  time_long = reshape2::melt(time)
  total_importation_intensity = reshape2::dcast(dat_og, epiweek ~ region, value.var="total_importation_intensity", fun.aggregate = sum) %>% select(-1)
  total_importation_intensity_long = reshape2::melt(total_importation_intensity)
  total_importation_intensity_long$value <- (total_importation_intensity_long$value - mean(dat_og$total_importation_intensity))/sd(dat_og$total_importation_intensity)
  
  time_long$value <- time_long$value - mean(unique(dat_og$epiweek))
  blim_90 <- matrix(0.9, nrow = nrow(blim_long))

  
  stan_dat_blim_sm_plot_other1 <-  list(T=length(unique(dat_og$epiweek)), I=length(unique(dat_og$region)), J=length(unique(dat_og$nuts)),
                                 Y_delta=Y_delta, Y_total=Y_total, include_lik=include_lik, nuts = nuts_vec$nuts,
                                 ncovs = 2, ncovs2 = 2,
                                 X = as.matrix(cbind(blim_long$value, total_importation_intensity_long$value)), X_90 = cbind(blim_90, blim_90),
                                 X_time = as.matrix(cbind(sm_long$value, time_long$value)))
  
  fit_temp_blim_sm_plot_other1 <- sampling(model, data = stan_dat_blim_sm_plot_other1, iter = n_iters, chains = n_chains, thin = n_thins)
  print(fit_temp_blim_sm_plot_other1, pars = c("beta", "beta_time", "sigma1", "sigma2"))
  sum_df_other1 <- summarise_draws(fit_temp_blim_sm_plot_other1)
  other1_rhat_US <- table(sum_df_other$rhat<1.01)
  other1_ess_bulk_US <- table(sum_df_other$ess_bulk>400)
  other1_ess_tail_US <- table(sum_df_other$ess_tail>400)
  
  t <- dat_og %>% group_by(region) %>% summarize(t = n())
  mean(t$t)
  length(unique(dat_og$region))
  
  # This is the final dataset used in the model
  dat_og <- readRDS("data/processed/US_weekly.rds")
  
  check_convergence <- list(main_rhat_US,main_ess_tail_US,main_ess_bulk_US,
                            other_rhat_US,other_ess_tail_US,other_ess_bulk_US,
                            other1_rhat_US,other1_ess_tail_US,other1_ess_bulk_US)
  
  ################################################################################ 
  save(check_convergence, file = "data/processed/check_convergence_US.RData")
  save.image("data/processed/US_final.RData")
  ################################################################################
  