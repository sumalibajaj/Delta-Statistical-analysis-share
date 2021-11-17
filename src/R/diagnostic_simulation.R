################################################################################
################################################################################
# DIAGNOSTIC SIMULATION DONE FOR A HYPOTHETICAL POPULATION.

# TO TEST MODEL PARAMETER IDENTIFIABILITY WE FIX THE PARAMETERS BEFORE HAND, 
# GENERATE DATA BASED ON THE PARAMETERS AND ASSUMED MODEL STRUCTURE. 
# INFERENCE PERFORMED ON THIS DATA IS ABLE TO REASONABLY RECOVER THE PARAMETERS.

# STAN FILE(S) USED:
# final.stan 

# OUTPUT:
# PDF WITH A PLOT COMPARING THE OBSERVED AND ESTIMATED RELATIVE GROWTH OF DELTA
################################################################################
################################################################################
rm(list = ls())

library(boot)
library(reshape2)
library(rstan)
options(mc.cores=4)
library(dplyr)
library(posterior)
library(loo)

set.seed(123)

n_iters = 20000
n_chains = 4
n_thins = 10
T = 15 # number of weeks
I = 30 # number of areas
J = 5 # number of regions (areas are nested within regions)
nuts = rep(1:J, each = I/J) # repeat region values corresponding to areas
sigma1 = 0.6
sigma2 = 0.2
beta1 = 2 # parameter value for variable's relationship with rho (relative growth)
beta2_time = 0.2 # parameter value for time's relationship with rho (relative growth)
n_area <- 200 # number of total cases in each area, (keeping it same here)
a_og <- c(-1, -3, -4, -3, -2) # first mean value of phi for each region (5 regions)
a <- rep(a_og, each = I/J)
b_og <- c(2, 3, 4, 2, 3) # first sd value of phi for each region (5 regions)
b <- rep(b_og, each = I/J)


# area level variable data for each week
var1 <- matrix(seq(0.2, 0.5, length = I), nrow = 1) # vac data for utla 1
var <- matrix(NA, nrow = T-1, ncol = I)
var <- rbind(var1, var)
for(t in 2:T){
  for(i in 1:I){
    var[t,i] = var[t-1,i]# repeated value for each area (like baseline immunity)
  }
}

var_input <- melt(var)

X_time_og = rep(seq(1:T), I)
X_time <- X_time_og - mean(unique(X_time_og))

# for region level "relative growth"
rho_ov_reg <- matrix(NA, nrow = T, ncol = J)
for(t in 1:T){
  for(j in 1:J){
    if(t==1){
      rho_ov_reg[t,j] = rnorm(1, 0, 1) # starts with zero advantage
    }
    else{
      rho_ov_reg[t,j] = rnorm(1, rho_ov_reg[t-1,j], sigma2) # each region's growth rate depends on its previous weeks growth rate
    }
  }
}

# for area level relative growth (rho)
dev_rho = matrix(rnorm(T*I, 0, sigma1), nrow = T, ncol = I) # all deviations centered around zero
# dev_rho = matrix(0, nrow = T, ncol = I) # all deviations centered around zero
# 5 regions with 6 areas within
rho_ov_reg_rep <- cbind(rho_ov_reg[,1], rho_ov_reg[,1], rho_ov_reg[,1], rho_ov_reg[,1], rho_ov_reg[,1], rho_ov_reg[,1],
                        rho_ov_reg[,2], rho_ov_reg[,2], rho_ov_reg[,2], rho_ov_reg[,2], rho_ov_reg[,2], rho_ov_reg[,2],
                        rho_ov_reg[,3], rho_ov_reg[,3], rho_ov_reg[,3], rho_ov_reg[,3], rho_ov_reg[,3], rho_ov_reg[,3],
                        rho_ov_reg[,4], rho_ov_reg[,4], rho_ov_reg[,4], rho_ov_reg[,4], rho_ov_reg[,4], rho_ov_reg[,4],
                        rho_ov_reg[,5], rho_ov_reg[,5], rho_ov_reg[,5], rho_ov_reg[,5], rho_ov_reg[,5], rho_ov_reg[,5])
rho = dev_rho+rho_ov_reg_rep # area level relative growth is its regional relative growth + deviation
rho_var = dev_rho+rho_ov_reg_rep + beta1*var + beta2_time*X_time # additionally, area level variable affects the area level relative growth

# plot observed proportion of Delta over time by areas
obs_rho <- reshape2::melt(rho_var) %>% rename(time = Var1, area = Var2, observed_p = value)
# ggplot(data = obs_rho, aes(x = time, y = observed_p)) + geom_line() + facet_wrap(~area)

# phi is the logit of proportion of delta cases
# phi_t1 = matrix(rnorm(I,-3,3), nrow = 1, ncol = I) # first value of utla level phi comes from this distribution
phi_t1 = matrix(a, nrow = 1, ncol = I) # first value of area level phi comes from this distribution
phi <- matrix(NA, nrow = T-1, ncol = I)
phi <- rbind(phi_t1, phi)
for(t in 2:T){
  for(i in 1:I){
    phi[t,i] = phi[t-1,i] + rho_var[t-1,i]
  }
}

p <- inv.logit(phi)
Y_total = matrix(n_area, nrow = T, ncol = I)
Y_delta = round(p*Y_total)

ignore = which(Y_delta/n_area>0.95) # which Y_delta_var values are greater than 95% sample size

include_lik = matrix(1, nrow = T, ncol = I)
include_lik[ignore] = 0 # add zero to cells which have to be ignored from the likelihood

X_90 <- matrix(0.9, nrow = nrow(var_input))

################################################################################
model <- stan_model("src/stan/final.stan")
stan_dat <-  list(T=nrow(Y_delta), I=ncol(Y_delta), J=length(unique(nuts)),
                  nuts = nuts, Y_delta=Y_delta, Y_total=Y_total, 
                  include_lik=include_lik,
                  ncovs = 1, ncovs2 = 1, 
                  X = as.matrix(var_input$value), X_time = as.matrix(X_time),
                  X_90 = X_90)

set.seed(Sys.time())
fit_temp <- sampling(model, data = stan_dat, iter = n_iters, chains = n_chains, thin = n_thins)
print(fit_temp, pars = c("beta", "beta_time","sigma1", "sigma2"))
print(fit_temp, pars = c("a", "b"))
loglikelihood <- extract_log_lik(fit_temp, "loglikelihood_old")
loo <- loo(loglikelihood)

sum_df <- summarise_draws(fit_temp)
main_rhat <- table(sum_df$rhat<1.01)
main_ess_bulk <- table(sum_df$ess_bulk>400)
main_ess_tail <- table(sum_df$ess_tail>400)

sum_df <- summarise_draws(fit_temp)
# save.image("data/processed/final_simulation.RData")
# load("data/processed/final_simulation.RData")
################################################################################
fit <- fit_temp
delta <- melt(Y_delta) %>% rename(time = Var1, utla_code = Var2, delta_voc = value) %>% mutate(total = 200)
NUTS1 <- rep(1:J, each = T*I/J)
area <- NUTS1 <- rep(1:I, each = T)
delta <- cbind(delta, NUTS1)
nuts_rep <-   rep(unique(delta$NUTS1), each = length(unique(delta$time)))

# ALL relative growths
pred_rho <- as.data.frame(extract(fit, pars = c("rho")))
pred_rho <- as.data.frame(cbind(unique(X_time_og), t(pred_rho))) 
for(j in 1:nrow(pred_rho)){
  pred_rho$mean[j] <- mean(as.numeric(pred_rho[j,2:ncol(pred_rho)]), na.rm = TRUE)
  temp <- quantile(pred_rho[j,2:ncol(pred_rho)],  probs = c(0.025, 0.975), na.rm = TRUE)
  pred_rho$CI95_lower[j] <- temp[1]
  pred_rho$CI95_upper[j] <- temp[2]
}

dat_pred <- pred_rho %>% select(V1, mean, CI95_lower, CI95_upper) %>% rename(time = V1)
dat_pred <- cbind(dat_pred, area)

colors <- c("Est. relative growth (observed)" = "#ffd544", 
            "Raw relative growth (observed)" = "#018040")
sim_plot <- ggplot(data = obs_rho) + 
            geom_line(aes(x = time, y = observed_p, color = "Raw relative growth (observed)"), linetype = "dashed") +
            geom_line(data = dat_pred, aes(x = time, y = mean, color = "Est. relative growth (observed)")) +
            geom_ribbon(data = dat_pred, aes(x = time, ymin = CI95_lower, ymax = CI95_upper, fill = "Est. relative growth (observed)"), alpha = 0.2) +
            facet_wrap(~ area) +
            theme_bw() +
            scale_color_manual(values = colors)+ 
            scale_fill_manual(values = colors) +
            guides(fill=FALSE) +
            labs(y = "relative growth", x = "time",
                 color = "Legend")+
            theme(text = element_text(size=15),
                  legend.title = element_blank(),
                  legend.text = element_text(size = 10),
                  legend.position="bottom",
                  strip.background = element_rect(
                    color="white", fill="white", size=1.5, linetype="solid"
                  ))

pdf(file = "outputs/simulation.pdf", # The directory you want to save the file in
    width = 18, # The width of the plot in inches
    height = 8) # The height of the plot in inches

sim_plot

dev.off()  

check_convergence <- list(main_rhat,main_ess_tail,main_ess_bulk)
################################################################################
save(check_convergence, file = "data/processed/check_convergence_diagnostic_simulation.RData")
################################################################################





