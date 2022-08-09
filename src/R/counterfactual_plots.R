library(tidyverse)
library(rstan)
library(posterior)
source("src/R/helper.R")

dat_og <- readRDS("data/processed/England_weekly_processed.rds")
covariates_to_include <- c("epiweek", "relative_self_mobility")
epiweeks <- sort(unique(dat_og$epiweek))
number_testing_epiweeks <- 2
nweeks <- length(unique(dat_og$epiweek))
last_training_epiweek <- nweeks - number_testing_epiweeks
training_epiweeks <- epiweeks[1:last_training_epiweek]
testing_epiweeks <- epiweeks[(last_training_epiweek + 1):nweeks]
stan_data <- prepare_stan_data(
  covariates_to_include, training_epiweeks,
  testing_epiweeks, last_training_epiweek,
  number_testing_epiweeks, dat_og
)

fit <- readRDS("data/processed/stan_fits/fit_best_shorter.rds")
sum_df <- summarise_draws(fit)
rhats <- sum(sum_df$rhat>1.01, na.rm = TRUE)
ess_bulks <- sum(sum_df$ess_bulk<400, na.rm = TRUE)
ess_tails <- sum(sum_df$ess_tail<400, na.rm = TRUE)
is_converged <- rhats == 0 & ess_bulks == 0 & ess_tails == 0
conv_df <- list(rhat=rhats,
                ess_bulk=ess_bulks,
                ess_tail=ess_tails,
                is_converged=is_converged)
saveRDS(conv_df, "data/processed/convergence_check_best_shorter.rds")

beta <- rstan::extract(fit, "beta")[[1]]
sigma1 <- rstan::extract(fit, "sigma1")[[1]]
dev_rho_raw <- rstan::extract(fit, "dev_rho_raw")[[1]]
rho_ov_reg <- rstan::extract(fit, "rho_ov_reg")[[1]]
phi_1 <- rstan::extract(fit, "phi_1")[[1]]
iter <- 1

project_fractions <- function(iter, X, beta, sigma1, dev_rho_raw,
                              rho_ov_reg, phi_1, stan_data) {
  T <- stan_data$T
  I <- stan_data$I
  nuts <- stan_data$nuts
  dev_rho <- matrix(nrow = T, ncol = I)
  phi <- matrix(nrow = T, ncol = I)
  rho <- matrix(nrow = T, ncol = I)
  p <- matrix(nrow = T, ncol = I)
  beta <- beta[iter, ]
  sigma1 <- sigma1[iter]
  dev_rho_raw <- dev_rho_raw[iter, ,]
  rho_ov_reg <- rho_ov_reg[iter, ,]
  phi_1 <- phi_1[iter, 1,]
  
  for(i in 1:I){
    for(t in 1:T){
      dev_rho[t,i] = dev_rho_raw[t,i]*sigma1;
    }
  }
  
  c = 1;
  for(i in 1:I){
    for(t in 1:T){
      rho[t,i] = rho_ov_reg[t,nuts[i]] + dev_rho[t,i] + sum(X[c, ]*beta);
      c = c+1;
    }
  }
  
  
  for(i in 1:I){
    for(t in 1:T){
      if(t==1){
        phi[t,i] = phi_1[i];
      }
      else{
        phi[t,i] = phi[t-1,i] + rho[t-1,i];
      }
    }
  }
  
  inv_logit <- function(x) {1 / (1 + exp(-x))}
  inv_logit(phi)
}

X_true <- stan_data$X
means_df <- stan_data$means_df
sd_df <- stan_data$sd_df
X_lower <- X_true
X_upper <- X_true
X_lower[, 2] <- (0 - means_df$relative_self_mobility) / sd_df$relative_self_mobility
X_upper[, 2] <- (1 - means_df$relative_self_mobility) / sd_df$relative_self_mobility
epiweeks <- seq(13, 21, 1)
get_stacked <- function(X_lower, iters) {
  p_lower <- project_fractions(
    1, X_lower, beta, sigma1, dev_rho_raw, rho_ov_reg,
    phi_1, stan_data) %>% 
    as.data.frame()
  colnames(p_lower) <- colnames(stan_data$Y_delta)
  p_lower <- p_lower %>% 
    mutate(epiweek=epiweeks) %>% 
    relocate(epiweek) %>% 
    pivot_longer(-epiweek) %>% 
    mutate(iteration=1)
  for(i in iters) {
    temp <- project_fractions(
      i, X_lower, beta, sigma1, dev_rho_raw, rho_ov_reg,
      phi_1, stan_data) %>% 
      as.data.frame()
    colnames(temp) <- colnames(stan_data$Y_delta)
    temp <- temp %>% 
      mutate(epiweek=epiweeks) %>% 
      relocate(epiweek) %>% 
      pivot_longer(-epiweek) %>% 
      mutate(iteration=i)
    p_lower <- p_lower %>% 
      bind_rows(temp)
  }
  p_lower
}

max_iters <- dim(extract(fit, "phi")[[1]])[1]
p_lower <- get_stacked(X_lower, seq(2, max_iters, 100)) %>% 
  mutate(type="low mobility")
p_upper <- get_stacked(X_upper, seq(2, max_iters, 100)) %>% 
  mutate(type="high mobility")
p_factual <- get_stacked(X_true, seq(2, max_iters, 100)) %>% 
  mutate(type="observed mobility")
p_all <- p_lower %>% 
  bind_rows(p_upper) %>% 
  bind_rows(p_factual)
p_all_sum <- p_all %>% 
  group_by(type, name, epiweek) %>% 
  summarise(lower=quantile(value, 0.025),
            middle=quantile(value, 0.5),
            upper=quantile(value, 0.975))

df_actual <- dat_og %>% 
  mutate(middle=delta_count / total_count) %>% 
  select(region, location, epiweek, middle) %>% 
  mutate(type="actual")
lookup <- df_actual %>% 
  select(region, location) %>% 
  unique()
df_combined <- p_all_sum %>% 
  rename(region=name) %>%
  left_join(lookup) %>% 
  bind_rows(df_actual)
locations <- unique(df_combined$location)[1:4]

df_short <- df_combined %>% 
  mutate(type=as.factor(type)) %>% 
  mutate(type=fct_relevel(type, "low mobility", "observed mobility", "high mobility")) %>% 
  mutate(location=stringr::str_to_title(location))


g <- ggplot(df_short %>% filter(type != "actual"),
       aes(x=epiweek, y=middle, colour=type)) +
  geom_line() +
  scale_colour_viridis_d("") +
  scale_x_continuous(breaks=epiweeks) +
  scale_y_continuous(labels=scales::percent) +
  ylab("Delta penetration") +
  xlab("2021 epiweek") +
  facet_wrap(~location)

ggsave("outputs/plot_counterfactual.pdf",
       g, width = 12, height = 8)  
ggsave("outputs/plot_counterfactual.png",
       g, width = 12, height = 8)  