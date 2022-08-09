library(rstan)
library(tidyverse)
source("src/R/helper.R")
args <- commandArgs(trailingOnly = TRUE)
is_posterior_predictive <- if_else(as.numeric(args[1]) == 1, 1, 0)

# actual data
dat_og <- readRDS("data/processed/England_weekly_processed.rds")
covariates_to_include <- c("epiweek", "relative_self_mobility")
stan_data <- prepare_stan_data_full(dat_og)
Y_delta_actual <- stan_data$Y_delta
Y_total <- stan_data$Y_total

Y_delta_actual <- Y_delta_actual %>% 
  mutate(epiweek=seq_along(E10000003) + 12) %>% 
  pivot_longer(-epiweek) %>% 
  rename(x=value)

# fitted data
fit <- readRDS("data/processed/stan_fits/fit_best_extended.rds")
if(is_posterior_predictive) {
  Y_delta <- rstan::extract(fit, "Y_delta_hat_post")[[1]]
}else {
  Y_delta <- rstan::extract(fit, "Y_delta_hat_prior")[[1]]
}
lower <- apply(Y_delta, c(2, 3), function(x) quantile(x, 0.025))
lower[lower < 0] <- NA # handles data points that are missing
upper <- apply(Y_delta, c(2, 3), function(x) quantile(x, 0.975))
upper[upper < 0] <- NA
middle <- apply(Y_delta, c(2, 3), function(x) quantile(x, 0.5))
middle[middle < 0] <- NA

lower <- lower / Y_total
middle <- middle / Y_total
upper <- upper / Y_total
colnames(lower) <- colnames(Y_total)
colnames(middle) <- colnames(Y_total)
colnames(upper) <- colnames(Y_total)
lower <- lower %>% 
  mutate(epiweek=seq_along(E10000003) + 12) %>% 
  pivot_longer(-epiweek) %>% 
  rename(lower=value)
middle <- middle %>% 
  mutate(epiweek=seq_along(E10000003) + 12) %>% 
  pivot_longer(-epiweek) %>% 
  rename(middle=value)
upper <- upper %>% 
  mutate(epiweek=seq_along(E10000003) + 12) %>% 
  pivot_longer(-epiweek) %>% 
  rename(upper=value)
df_simulated <- lower %>% 
  left_join(middle) %>% 
  left_join(upper) %>% 
  mutate(type="simulated")

Y_total <- Y_total %>% 
  mutate(epiweek=seq_along(E10000003) + 12) %>% 
  pivot_longer(-epiweek) %>% 
  rename(n=value)

df_actual <- Y_delta_actual %>% 
  left_join(Y_total) %>% 
  filter(x != -99) %>% 
  mutate(middle=qbeta(0.5, 1 + x, 1 + n - x),
         lower=qbeta(0.025, 1 + x, 1 + n - x),
         upper=qbeta(0.975, 1 + x, 1 + n - x)) %>% 
  select(-one_of(c("x", "n"))) %>% 
  mutate(type="actual")

# combine all
df_both <- df_actual %>% 
  bind_rows(df_simulated) %>% 
  na.omit() %>% 
  rename(region=name) %>% 
  left_join(dat_og) %>% 
  mutate(location=stringr::str_to_title(location))

locations <- sort(unique(df_both$location))

plot_locations_subset <- function(ids, df_both) {
  df_short <- df_both %>% 
    filter(location %in% locations[ids])

  ggplot(data = filter(df_short, type=="actual"),
         aes(x=epiweek, y=middle)) +
  geom_ribbon(data=filter(df_short, type=="simulated"),
              aes(ymin=lower, ymax=upper),
              fill="blue", alpha=0.5) +
  geom_line(data=filter(df_short, type=="simulated"),
            linetype=1) +
    geom_pointrange(aes(ymin=lower, ymax=upper),
                    colour="orange", size=0.5) +
  scale_y_continuous(labels=scales::percent) +
  ylab("Delta penetration") +
  scale_x_continuous(breaks=seq(13, 23, 1)) +
  xlab("2021 epiweek") +
  facet_wrap(~location)
}

if(is_posterior_predictive) {
  basefilename <- "outputs/posterior_predictive"
  basedatafilename <- "data/processed/posterior_predictive.rds"
}else{
  basefilename <- "outputs/prior_predictive"
  basedatafilename <- "data/processed/prior_predictive.rds"
}
g <- plot_locations_subset(seq(1, 64), df_both)
ggsave(paste0(basefilename, ".pdf"), g,
       width = 12, height = 6)
ggsave(paste0(basefilename, ".png"), g,
       width = 12, height = 6)
saveRDS(df_both, basedatafilename)
