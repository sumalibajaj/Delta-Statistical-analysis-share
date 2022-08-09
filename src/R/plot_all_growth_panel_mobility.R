library(tidyverse)
library(cowplot)

rho_all <- readRDS("data/processed/growth_quantiles.rds")
ranked_vaccines <- readRDS("data/processed/ranked_vaccine_2.rds") %>% 
  mutate(type="Vaccination")
ranked_mobility <- readRDS("data/processed/ranked_mobility.rds") %>% 
  mutate(type="Within-UTLA mobility")
ranked_both <- ranked_vaccines %>% 
  bind_rows(ranked_mobility) %>% 
  mutate(type=as.factor(type)) %>% 
  mutate(type=fct_relevel(type, "Within-UTLA mobility", "Vaccination"))

g1 <- ranked_both %>% 
  filter(type=="Within-UTLA mobility") %>% 
  ggplot(aes(x=epiweek, y=fraction_delta,
             colour=other)) +
  geom_line(aes(group=location), alpha=0.5) +
  geom_smooth(se=FALSE, linetype=2, size=1.3) +
  theme_bw() +
  scale_colour_viridis_d("", direction=-1) +
  scale_y_continuous(labels=scales::percent) +
  ylab("Delta penetration") +
  scale_x_continuous(breaks=seq(13, 23, 1)) +
  theme(legend.position = c(0.2, 0.6),
        legend.text = element_text(size=14),
        text=element_text(size=16),
        legend.background = element_blank(),
        strip.background = element_rect(
          color="white", fill="white", size=1.5, linetype="solid")) +
  xlab("2021 Epiweek")

examples <- ranked_both %>% 
  filter(location != "Greater London") %>% 
  group_by(type, other) %>% 
  summarise(location=last(location)) %>% 
  ungroup() %>% 
  filter(other != "medium")

rho_short <- rho_all %>% 
  left_join(examples) %>% 
  filter(!is.na(type))

g2 <- rho_short %>% 
  filter(type=="Within-UTLA mobility") %>%
  ggplot(aes(x=epiweek, y=middle)) +
  geom_ribbon(aes(ymin=lower, ymax=upper),
              fill="royalblue1") +
  geom_line() +
  theme_bw() + 
  geom_text(aes(x=18, y=-2, label=location), size=5) +
  scale_x_continuous(breaks=seq(13, 23, 1)) +
  ylab("Relative growth\n(log odds scale)") +
  xlab("2021 Epiweek") +
  theme(text=element_text(size=16), strip.background = element_rect(
    color="white", fill="white", size=1.5, linetype="solid"), 
    strip.text.x = element_blank(), strip.text.y = element_blank()) +
  facet_wrap(~other, nrow=1)

g <- plot_grid(g1, g2, nrow = 2, labels = c("A.", "B."))
save_plot("outputs/growth_panel_mobility_only.pdf", g,
          base_height = 8, base_width = 12)
save_plot("outputs/growth_panel_mobility_only.png", g,
          base_height = 8, base_width = 12)
