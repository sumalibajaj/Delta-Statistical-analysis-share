################################################################################
################################################################################
# MODEL COMPARISON FOR ENGLAND

# COMPARING OUT OF SAMPLE PREDICTION (2 WEEKS AT THE END) BETWEEN MODEL WITHOUT 
# ANY VARIABLES VS FINAL MODEL (WITH VARIABLES SPECIFIED BELOW)

# MODELS AND THE VARIABLES USED:
# fit_temp_novar USES NO VARIABLES
# fit_temp_v2 USES MOBILITY, TIME 
#                  (FINAL MODEL IN THE MANUSCRIPT)

# STAN FILE(S) USED:
# final_MC_novar.stan 
# final_MC_var.stan 

################################################################################
################################################################################
library(dplyr)
library(rstan)
library(reshape2)
library(loo)
library(posterior)
options(mc.cores=4)

rm(list = ls())

dat_og <- readRDS("data/processed/England_weekly_mob.rds")
n_iters = 40000
n_chains = 4
n_thins = 10
length(unique(dat_og$epiweek))
T1 = 2
T = length(unique(dat_og$epiweek)) -T1

################################################################################
# WITHOUT VARIABLES
################################################################################
# model_novar <- stan_model("src/stan/growthrate_hierarchical_utla_MC_extra.stan") # diff rho
model_novar <- stan_model("src/stan/final_MC_novar.stan")

# LIKELIHOOD MATRIX
include_lik = reshape2::dcast(dat_og, epiweek ~ region, value.var="delta_count") %>% select(-1)
include_lik[!is.na(include_lik)] <- 1
include_lik[is.na(include_lik)] <- 0
include_lik_T<- include_lik[1:T, ]
include_lik_extra <- include_lik[(T+1):nrow(include_lik), ]

# OTHER INPUTS FOR STAN
# region_vec <- dat_og %>% summarize(unique(region)) %>% select(1)%>% rename(region = "unique(region)")
# run the line below here instead if an error is thrown:
region_vec <- dat_og %>% summarize(unique(region)) %>% select(1)
nuts_region <- dat_og %>% group_by(nuts) %>% summarize(unique(region)) %>% rename(region = "unique(region)")
nuts_vec <- merge(region_vec, nuts_region, by = "region")

Y_delta = reshape2::dcast(dat_og, epiweek ~ region, value.var="delta_count", fun.aggregate = sum) %>% select(-1)
Y_delta_T = Y_delta[1:T,]
Y_delta_extra = Y_delta[(T+1):nrow(Y_delta),]
Y_delta_extra_long = reshape2::melt(Y_delta_extra) %>% rename(state = variable, observed = value)


Y_total = reshape2::dcast(dat_og, epiweek ~ region, value.var="total_count", fun.aggregate = sum) %>% select(-1)
Y_total_T = Y_total[1:T,]
Y_total_extra = Y_total[(T+1):nrow(Y_total),]

# no variable
stan_dat_novar <-  list(T=T, T1 = T1, I=length(unique(dat_og$region)), J=length(unique(dat_og$nuts)),
                        Y_delta=Y_delta_T, Y_total=Y_total_T, include_lik=include_lik_T, 
                        nuts = nuts_vec$nuts,
                        Y_delta_extra=Y_delta_extra, Y_total_extra=Y_total_extra, include_lik_extra=include_lik_extra)
fit_temp_novar <- sampling(model_novar, data = stan_dat_novar, iter = n_iters, chains = n_chains, thin = n_thins)
print(fit_temp_novar, pars = c("sigma1", "sigma2"))
print(fit_temp_novar, pars = c("rho_ov_reg"))
print(fit_temp_novar, pars = c("a", "b"))

sum_df_novar <- summarise_draws(fit_temp_novar)
main_rhat_England_novar <- table(sum_df_novar$rhat<1.01)
main_ess_bulk_England_novar <- table(sum_df_novar$ess_bulk>400)
main_ess_tail_England_novar <- table(sum_df_novar$ess_tail>400)

loglikelihood_novar <- extract_log_lik(fit_temp_novar, "loglikelihood_extra")
loglikelihood_novar_test <- loglikelihood_novar[ , apply(loglikelihood_novar, 2, function(x) !any(is.na(x)))]

means_novar <- colMeans(loglikelihood_novar_test)
sum_means_novar <- sum(colMeans(loglikelihood_novar_test))
mean_means_novar <- mean(colMeans(loglikelihood_novar_test))

Y_delta_extra_est <- as.data.frame(extract(fit_temp_novar, pars = "Y_delta_hat_post_extra"))
Y_delta_extra_est <- as.data.frame(cbind(c((T+1):length(unique(dat_og$epiweek))), t(Y_delta_extra_est))) 
for(j in 1:nrow(Y_delta_extra_est)){
  Y_delta_extra_est$mean[j] <- mean(as.numeric(Y_delta_extra_est[j,2:ncol(Y_delta_extra_est)]), na.rm = TRUE)
  temp <- quantile(Y_delta_extra_est[j,2:ncol(Y_delta_extra_est)],  probs = c(0.025, 0.975), na.rm = TRUE)
  Y_delta_extra_est$CI95_lower[j] <- temp[1]
  Y_delta_extra_est$CI95_upper[j] <- temp[2]
}

Y_delta_extra_est <- Y_delta_extra_est %>% select(V1, mean, CI95_lower, CI95_upper) %>% rename(epiweek = V1, estimated_mean = mean)

Y_delta_plot <- cbind(Y_delta_extra_long, Y_delta_extra_est)

################################################################################
# HIERARCHICAL WITH VARIABLES
################################################################################

model <- stan_model("src/stan/final_MC_var.stan")

relative_self_mobility = reshape2::dcast(dat_og, epiweek ~ region, value.var="relative_self_mobility", fun.aggregate = sum) %>% select(-1)
relative_self_mobility_T = relative_self_mobility[1:T,]
relative_self_mobility_extra = relative_self_mobility[(T+1):nrow(relative_self_mobility),]
relative_self_mobility_long_T = reshape2::melt(relative_self_mobility_T)
relative_self_mobility_long_extra = reshape2::melt(relative_self_mobility_extra)

epiweek = reshape2::dcast(dat_og, epiweek ~ region, value.var="epiweek", fun.aggregate = sum) %>% select(-1)
epiweek_T = epiweek[1:T,]
epiweek_extra = epiweek[(T+1):nrow(epiweek),]
epiweek_long_T = reshape2::melt(epiweek_T)
epiweek_long_extra = reshape2::melt(epiweek_extra)


# LIKELIHOOD MATRIX
include_lik = reshape2::dcast(dat_og, epiweek ~ region, value.var="delta_count") %>% select(-1)
include_lik[!is.na(include_lik)] <- 1
include_lik[is.na(include_lik)] <- 0
include_lik_T<- include_lik[1:T, ]
include_lik_extra <- include_lik[(T+1):nrow(include_lik), ]

# OTHER INPUTS FOR STAN
# region_vec <- dat_og %>% summarize(unique(region)) %>% select(1) %>% rename(region = "unique(region)")
# run the line below here instead if an error is thrown:
region_vec <- dat_og %>% summarize(unique(region)) %>% select(1)
nuts_region <- dat_og %>% group_by(nuts) %>% summarize(unique(region)) %>% rename(region = "unique(region)")
nuts_vec <- merge(region_vec, nuts_region, by = "region")

Y_delta = reshape2::dcast(dat_og, epiweek ~ region, value.var="delta_count", fun.aggregate = sum) %>% select(-1)
Y_delta_T = Y_delta[1:T,]
Y_delta_extra = Y_delta[(T+1):nrow(Y_delta),]
Y_delta_extra_long = reshape2::melt(Y_delta_extra) %>% rename(state = variable, observed = value)


Y_total = reshape2::dcast(dat_og, epiweek ~ region, value.var="total_count", fun.aggregate = sum) %>% select(-1)
Y_total_T = Y_total[1:T,]
Y_total_extra = Y_total[(T+1):nrow(Y_total),]


####### UK ####### 
# relative self mobility and time
stan_dat_v2 <-  list(T=T, T1 = T1, I=length(unique(dat_og$region)), J=length(unique(dat_og$nuts)),
                     Y_delta=Y_delta_T, Y_total=Y_total_T, include_lik=include_lik_T, 
                     nuts = nuts_vec$nuts, ncovs = 2, 
                     X = as.matrix(cbind(relative_self_mobility_long_T$value, epiweek_long_T$value)),
                     Y_delta_extra=Y_delta_extra, Y_total_extra=Y_total_extra, 
                     include_lik_extra=include_lik_extra,
                     X_extra = as.matrix(cbind(relative_self_mobility_long_extra$value, epiweek_long_extra$value)))

fit_temp_v2 <- sampling(model, data = stan_dat_v2, iter = n_iters, chains = n_chains, thin = n_thins)
print(fit_temp_v2, pars = c("beta", "sigma1", "sigma2"))

sum_df_v2 <- summarise_draws(fit_temp_v2)
main_rhat_England_v2 <- table(sum_df_v2$rhat<1.01)
main_ess_bulk_England_v2 <- table(sum_df_v2$ess_bulk>400)
main_ess_tail_England_v2 <- table(sum_df_v2$ess_tail>400)

loglikelihood_v2 <- extract_log_lik(fit_temp_v2, "loglikelihood_extra")
loglikelihood_v2_test <- loglikelihood_v2[ , apply(loglikelihood_v2, 2, function(x) !any(is.na(x)))]

means_v2 <- colMeans(loglikelihood_v2_test)
sum_means_v2 <- sum(colMeans(loglikelihood_v2_test))
mean_means_v2 <- mean(colMeans(loglikelihood_v2_test))

# NO VAR VS V2
# difference in sums (numerator)
difference_sum <- sum_means_novar - sum_means_v2
# denominator
se_sum <- sqrt(length(means_novar)) * sd(means_novar - means_v2)
pvalue_sum <- ifelse(difference_sum/se_sum <0, pnorm(difference_sum/se_sum),
                     1 - pnorm(difference_sum/se_sum))

# # difference in means (numerator)
difference_mean <- mean_means_novar - mean_means_v2
# denominator
se_mean <- sd(means_novar - means_v2)/sqrt(length(means_novar))
pvalue_mean <- ifelse(difference_mean/se_mean <0, pnorm(difference_mean/se_mean),
                      1 - pnorm(difference_mean/se_mean))



check_convergence <- list(main_rhat_England_novar,main_ess_tail_England_novar,main_ess_bulk_England_novar,
                          main_rhat_England_v2,main_ess_tail_England_v2,main_ess_bulk_England_v2)

################################################################################
save(check_convergence, file = "data/processed/check_convergence_England_MC.RData")
# 18 Oct 2021
save.image("data/processed/England_MC.RData")
################################################################################


