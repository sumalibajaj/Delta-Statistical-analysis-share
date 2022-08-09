library(rstan)
library(tidyverse)
source("src/R/helper.R")

dat_og <- readRDS("data/processed/England_weekly_processed.rds")
covariates_to_include <- c("epiweek", "relative_self_mobility")
stan_data <- prepare_stan_data_full(dat_og)
Y_delta_actual <- stan_data$Y_delta
Y_total <- stan_data$Y_total

Y_delta_actual <- Y_delta_actual %>% 
  mutate(epiweek=seq_along(E10000003) + 12) %>% 
  pivot_longer(-epiweek) %>% 
  rename(x=value)

fit <- readRDS("data/processed/stan_fits/fit_best_extended.rds")

rho <- rstan::extract(fit, "rho")[[1]]
lower <- apply(rho, c(2, 3), function(x) quantile(x, 0.025))
upper <- apply(rho, c(2, 3), function(x) quantile(x, 0.975))
middle <- apply(rho, c(2, 3), function(x) quantile(x, 0.5))
colnames(lower) <- colnames(Y_total)
colnames(middle) <- colnames(Y_total)
colnames(upper) <- colnames(Y_total)

lower <- lower %>% 
  as.data.frame() %>% 
  mutate(epiweek=seq_along(E10000003) + 12) %>% 
  pivot_longer(-epiweek) %>% 
  rename(lower=value)
middle <- middle %>% 
  as.data.frame() %>% 
  mutate(epiweek=seq_along(E10000003) + 12) %>% 
  pivot_longer(-epiweek) %>% 
  rename(middle=value)
upper <- upper %>% 
  as.data.frame() %>% 
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
  mutate(fraction=x/n) %>% 
  mutate(logit_f=log(fraction / (1 - fraction))) %>% 
  arrange(name, epiweek) %>% 
  group_by(name) %>% 
  mutate(growth=c(NA, diff(logit_f))) %>% 
  mutate(include=if_else(lag(x)==-99, 0, 1)) %>% 
  mutate(include=if_else(x==-99, 0, include)) %>% 
  filter(include==1) %>% 
  filter(growth!=-Inf) %>% 
  filter(growth!=Inf) %>% 
  select(-one_of(c("x", "n"))) %>% 
  mutate(type="actual") %>% 
  rename(middle=growth)

# combine all
df_both <- df_actual %>% 
  bind_rows(df_simulated) %>% 
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
    ylab("weekly Delta growth (log odds scale)") +
    scale_x_continuous(breaks=seq(13, 23, 1)) +
    xlab("2021 epiweek") +
    theme(axis.text.x = element_text(size=6)) +
    facet_wrap(~location)
}

basefilename <- "outputs/posterior_predictive_growth"
g <- plot_locations_subset(seq(1, 64), df_both)
ggsave(paste0(basefilename, ".pdf"), g,
       width = 12, height = 6)
ggsave(paste0(basefilename, ".png"), g,
       width = 12, height = 6)