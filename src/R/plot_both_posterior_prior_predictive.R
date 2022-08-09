library(tidyverse)
library(cowplot)

df_prior <- readRDS("data/processed/prior_predictive.rds")
df_posterior <- readRDS("data/processed/posterior_predictive.rds")

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
    theme(axis.text.x = element_text(size=6)) +
    theme(axis.text.y = element_text(size=6)) +
    facet_wrap(~location)
}

locations <- sort(unique(df_prior$location))
g1 <- plot_locations_subset(seq(1, 64), df_prior)
g2 <- plot_locations_subset(seq(1, 64), df_posterior)

g <- plot_grid(g2, g1, nrow=2, labels = c("A.", "B."))
save_plot("outputs/prior_and_posterior_predictive.png",
          g, base_width = 12, base_height = 14)
