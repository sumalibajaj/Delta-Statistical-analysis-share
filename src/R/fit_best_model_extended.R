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

dat_og <- readRDS("data/processed/England_weekly_processed.rds")

n_iters = 60000
n_chains = 4
n_thins = 100

################################################################################
# HIERARCHICAL WITH VARIABLES
# NEED UPDATED MAIN REGION LEVEL DATA ON SELF MOBILITY AND IMPORTATION INTENSITY
################################################################################

normalise <- function(X){ (X - mean(X))/(sd(X))}

# VARIABLES (Sum is used as a function to replace NA's by zero)
relative_self_mobility = reshape2::dcast(dat_og, epiweek ~ region, value.var="relative_self_mobility", fun.aggregate = sum) %>% select(-1)
relative_self_mobility_long = reshape2::melt(relative_self_mobility)
relative_self_mobility_long$value = normalise(relative_self_mobility_long$value)

time = reshape2::dcast(dat_og, epiweek ~ region, value.var="epiweek", fun.aggregate = sum) %>% select(-1)
time_long = reshape2::melt(time)
time_long$value <- time_long$value - mean(dat_og$epiweek)
time_long$value <- normalise(time_long$value)

# LIKELIHOOD MATRIX
include_lik = reshape2::dcast(dat_og, epiweek ~ region, value.var="delta_count") %>% select(-1)
include_lik[!is.na(include_lik)] <- 1
include_lik[is.na(include_lik)] <- 0

# OTHER INPUTS FOR STAN
region_vec <- dat_og %>% summarize(unique(region)) %>% select(1) %>% rename(region = "unique(region)")
nuts_region <- dat_og %>% group_by(nuts) %>% summarize(unique(region)) %>% rename(region = "unique(region)")
nuts_vec <- merge(region_vec, nuts_region, by = "region")

dat_og$delta_count = ifelse(is.na(dat_og$delta_count), -99, dat_og$delta_count)
Y_delta = reshape2::dcast(dat_og, epiweek ~ region, value.var="delta_count", fun.aggregate = sum) %>% select(-1)

dat_og$total_count = ifelse(is.na(dat_og$total_count), -99, dat_og$total_count)
Y_total = reshape2::dcast(dat_og, epiweek ~ region, value.var="total_count", fun.aggregate = sum) %>% select(-1)


# relative_self_mobility and time
model <- stan_model("src/stan/growth_full_extended.stan")

X_90 <- matrix(1, nrow = nrow(relative_self_mobility_long))

stan_dat_sm_time_plot <-  list(T=length(unique(dat_og$epiweek)), I=length(unique(dat_og$region)), J=length(unique(dat_og$nuts)),
                          Y_delta=Y_delta, Y_total=Y_total, include_lik=include_lik, nuts = nuts_vec$nuts,
                          ncovs = 1, ncovs2 = 1,
                          X = as.matrix(relative_self_mobility_long$value),
                          X_90 = X_90,
                          X_time = as.matrix(time_long$value))
fit <- sampling(model, data = stan_dat_sm_time_plot, iter = n_iters, chains = n_chains, thin = n_thins)
print(fit, pars = c("beta", "beta_time", "sigma1", "sigma2"))

# only save if converged
sum_df <- summarise_draws(fit)

# it's ok to remove NAs since these are
# parameter values which are set to -99
rhats <- sum(sum_df$rhat>1.01, na.rm = TRUE)
ess_bulks <- sum(sum_df$ess_bulk<400, na.rm = TRUE)
ess_tails <- sum(sum_df$ess_tail<400, na.rm = TRUE)
is_converged <- rhats == 0 & ess_bulks == 0 & ess_tails == 0

if(is_converged) {
  saveRDS(fit, "data/processed/stan_fits/fit_best_extended.rds")
} else {
  print("Not converged.")
  print(paste0("rhat = ", rhats))
  print(paste0("ess_bulk = ", ess_bulks))
  print(paste0("ess_tail = ", ess_tails))
}
