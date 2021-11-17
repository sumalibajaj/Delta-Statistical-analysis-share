library(dplyr)
library(rstan)
library(bayesplot) # for posterior predictive checks
library(ggplot2)
library(gridExtra)
library(data.table)
library(loo)
library(cowplot)

rm(list = ls())

load("data/processed/England_final.RData")


w1 <- 66-53
w2 <- max(dat_og$epiweek )
wks <- w1:w2
byweek <- expand.grid(regcd=unique(dat_og$region),time=wks,stringsAsFactors = F)
byweek$delta_count  <- apply(byweek,1,function(x)sum(subset(dat_og,epiweek==x[2]&region==x[1])$delta_count))
byweek$total <- apply(byweek,1,function(x)sum(subset(dat_og,epiweek==x[2]&region==x[1])$total_count))
byweek$epiweek <- byweek$time
dat_full <- byweek %>% select(-time) %>% rename(region_rep = regcd) %>% arrange(region_rep)
dat_full$p_obs <- dat_full$delta_count/dat_full$total
dat_full <- merge(dat_full, dat_og %>% group_by(region) %>% summarize(location = unique(location)),
                  by.x = c("region_rep"), by.y = c("region")) %>% rename(location_rep = location)

dat_new <- dat_og %>% rename(region_rep = region)
dat_new <- dat_og %>% rename(location_rep = location)

plot_fun <- function(fit_input){
  fit <- fit_input
  region_rep <- rep(unique(dat_og$region), each = length(unique(dat_og$epiweek)))
  location_rep <- rep(unique(dat_og$location), each = length(unique(dat_og$epiweek)))
  epiweek_rep <- rep(min(dat_og$epiweek):max(dat_og$epiweek), length(unique(dat_og$region)))
  nuts_rep <-   rep(unique(dat_og$NUTS1), each = length(unique(dat_og$epiweek)))
  nuts_region <- dat_og %>% group_by(NUTS1) %>% summarize(region = unique(region))
  
  pred_rho_ov <- as.data.frame(extract(fit, pars = c("rho_ov_reg")))
  pred_rho_ov <- as.data.frame(cbind(unique(dat_og$epiweek), t(pred_rho_ov))) 
  for(j in 1:nrow(pred_rho_ov)){
    pred_rho_ov$mean[j] <- mean(as.numeric(pred_rho_ov[j,2:ncol(pred_rho_ov)]), na.rm = TRUE)
    temp <- quantile(pred_rho_ov[j,2:ncol(pred_rho_ov)],  probs = c(0.025, 0.975), na.rm = TRUE)
    pred_rho_ov$CI95_lower[j] <- temp[1]
    pred_rho_ov$CI95_upper[j] <- temp[2]
  }
  
  dat_pred_ov <- pred_rho_ov %>% select(V1, mean, CI95_lower, CI95_upper) %>% rename(epiweek = V1)
  dat_pred_ov <- cbind(nuts_rep, dat_pred_ov)
  
  p_rho_ov <- ggplot(data = dat_pred_ov, aes(x = epiweek)) +
    geom_point(aes(y = mean), color = "#fe6362") + 
    geom_line(aes(y = mean), color = "#fe6362") +
    facet_wrap(~ nuts_rep) +
    geom_ribbon(aes(x = epiweek, ymin = CI95_lower, ymax = CI95_upper), fill = "#fe6362", alpha = 0.3) +
    # geom_vline(xintercept = imp_epiweek, color = "grey") + 
    theme_bw() +
    labs(title = "Estimated weekly growth in log odds of being delta variant",
         x = "2021 Epiweek", y = "Estimated growth") +
    theme(text = element_text(size=15))

  # State level growth
  states_imp <- c("BLACKPOOL", "GREATER LONDON", "GREATER MANCHESTER", "OXFORDSHIRE", "PORTSMOUTH")
  
  pred_rho <- as.data.frame(extract(fit, pars = c("rho")))
  pred_rho <- as.data.frame(cbind(unique(dat_og$epiweek), t(pred_rho))) 
  for(j in 1:nrow(pred_rho)){
    pred_rho$mean[j] <- mean(as.numeric(pred_rho[j,2:ncol(pred_rho)]), na.rm = TRUE)
    temp <- quantile(pred_rho[j,2:ncol(pred_rho)],  probs = c(0.025, 0.975), na.rm = TRUE)
    pred_rho$CI95_lower[j] <- temp[1]
    pred_rho$CI95_upper[j] <- temp[2]
  }
  
  dat_pred <- pred_rho %>% select(V1, mean, CI95_lower, CI95_upper) %>% rename(epiweek = V1)
  dat_pred <- cbind(region_rep, dat_pred)
  dat_pred <- cbind(location_rep, dat_pred)
  
  p_rho_all <- ggplot(data = dat_pred, aes(x = epiweek)) +
    geom_point(aes(y = mean), color = "grey") +
    geom_line(aes(y = mean), color = "grey") +
    # facet_wrap(~ region_rep, nrow = 1) +
    facet_wrap(~location_rep) +
    geom_ribbon(aes(x = epiweek, ymin = CI95_lower, ymax = CI95_upper), fill = "grey", alpha = 0.3) +
    # geom_vline(xintercept = imp_epiweek, color = "grey") +
    theme_bw() +
    labs(x = "2021 Epiweek", y = "weekly growth") +
    theme(text = element_text(size=15),
          strip.background = element_rect(
            color="white", fill="white", size=1.5, linetype="solid"),
          strip.text.x = element_text(size = 8)) +
    xlim(10,25)
  p_rho_all
  
  p_rho <- ggplot(data = dat_pred %>% filter(location_rep %in% states_imp), aes(x = epiweek)) +
    geom_point(aes(y = mean), color = "grey") +
    geom_line(aes(y = mean), color = "grey") +
    # facet_wrap(~ region_rep, nrow = 1) +
    facet_grid(cols = vars(location_rep)) +
    geom_ribbon(aes(x = epiweek, ymin = CI95_lower, ymax = CI95_upper), fill = "grey", alpha = 0.3) +
    # geom_vline(xintercept = imp_epiweek, color = "grey") +
    theme_bw() +
    labs(x = "", y = "relative growth") +
    theme(text = element_text(size=15),
          strip.background = element_rect(
            color="white", fill="white", size=1.5, linetype="solid"
          ),
          strip.text.x = element_blank(),
          axis.text.x = element_blank()) +
    xlim(10,25)
  p_rho

    p <- as.data.frame(extract(fit, pars = c("p")))
    p <- as.data.frame(cbind(unique(dat_og$epiweek), t(p))) 
    for(j in 1:nrow(p)){
      p$mean[j] <- mean(as.numeric(p[j,2:ncol(p)]), na.rm = TRUE)
      temp <- quantile(p[j,2:ncol(p)],  probs = c(0.025, 0.975), na.rm = TRUE)
      p$CI95_lower[j] <- temp[1]
      p$CI95_upper[j] <- temp[2]
    }
    p <- p %>% select(V1, mean, CI95_lower, CI95_upper) %>% rename(time = V1) %>% mutate(model = "variable")
    p <- cbind(region_rep, p)
    p <- cbind(location_rep, p)
    p <- merge(p, nuts_region, by.x = "region_rep", by.y = "region")
        
    p_wovar <- as.data.frame(extract(fit, pars = c("p_wovac")))
    p_wovar <- as.data.frame(cbind(unique(dat_og$epiweek), t(p_wovar))) 
    for(j in 1:nrow(p_wovar)){
      p_wovar$mean[j] <- mean(as.numeric(p_wovar[j,2:ncol(p_wovar)]), na.rm = TRUE)
      temp <- quantile(p_wovar[j,2:ncol(p_wovar)],  probs = c(0.025, 0.975), na.rm = TRUE)
      p_wovar$CI95_lower[j] <- temp[1]
      p_wovar$CI95_upper[j] <- temp[2]
    }
    p_wovar <- p_wovar %>% select(V1, mean, CI95_lower, CI95_upper) %>% rename(time = V1) %>% mutate(model = "w/o variable")
    p_wovar <- cbind(region_rep, p_wovar)
    p_wovar <- cbind(location_rep, p_wovar)
    p_wovar <- merge(p_wovar, nuts_region, by.x = "region_rep", by.y = "region")
    
    p_90var <- as.data.frame(extract(fit, pars = c("p_90vac")))
    p_90var <- as.data.frame(cbind(unique(dat_og$epiweek), t(p_90var))) 
    for(j in 1:nrow(p_90var)){
      p_90var$mean[j] <- mean(as.numeric(p_90var[j,2:ncol(p_90var)]), na.rm = TRUE)
      temp <- quantile(p_90var[j,2:ncol(p_90var)],  probs = c(0.025, 0.975), na.rm = TRUE)
      p_90var$CI95_lower[j] <- temp[1]
      p_90var$CI95_upper[j] <- temp[2]
    }
    p_90var <- p_90var %>% select(V1, mean, CI95_lower, CI95_upper) %>% rename(time = V1) %>% mutate(model = "90% variable")
    p_90var <- cbind(region_rep, p_90var)
    p_90var <- cbind(location_rep, p_90var)
    p_90var <- merge(p_90var, nuts_region, by.x = "region_rep", by.y = "region")

    colors <- c("Estimated Delta proportion" = "#ffd544", 
                "Est. Delta proportion with min self mobility" = "#484c7e", 
                "Est. Delta proportion with max self mobility" = "#fe6362", 
                "Observed Delta proportion" = "royalblue1", 
                "Relative self mobility" = "black",
                "Reported cases" = "black")
    
    p_freq_state_all <- ggplot(data = p) +
      geom_line(aes(x = time, y = mean, color = "Estimated Delta proportion")) +
      geom_ribbon(aes(x = time, ymin = CI95_lower, ymax = CI95_upper, fill = "Estimated Delta proportion"), alpha = 0.5) +
      geom_line(data = p_90var, aes(x = time, y = mean, color = "Est. Delta proportion with max self mobility")) +
      geom_ribbon(data = p_90var, aes(x = time, ymin = CI95_lower, ymax = CI95_upper, fill = "Est. Delta proportion with max self mobility"), alpha = 0.3) +
      geom_line(data = p_wovar, aes(x = time, y = mean, color = "Est. Delta proportion with min self mobility")) +
      geom_ribbon(data = p_wovar , aes(x = time, ymin = CI95_lower, ymax = CI95_upper, fill = "Est. Delta proportion with min self mobility"), alpha = 0.3) +
      geom_line(data = dat_full, aes(x = epiweek, y = p_obs, color = "Observed Delta proportion"), linetype = "dashed") +
      geom_line(data = dat_new , aes(x = epiweek, y = relative_self_mobility, color = "Relative self mobility"), linetype = "dashed") +
      # facet_wrap(~ region_rep, nrow = 1) +
      # facet_grid(cols = vars(region_rep)) +
      facet_wrap(~location_rep) +
      labs(y = "proportion",
           color = "Legend") +
      theme_bw() +
      scale_color_manual(values = colors)+
      scale_fill_manual(values = colors) +
      guides(fill=FALSE) +
      theme(text = element_text(size=15),
            legend.title = element_blank(),
            legend.text = element_text(size = 13),
            legend.position="bottom",
            strip.background = element_rect(
              color="white", fill="white", size=1.5, linetype="solid"),
            strip.text.x = element_text(size = 8))+
      xlim(10,25)+
      ylim(0,1)

    p_freq_state_all
    
    max_new_cases_reported_all <- max(dat_new$new_cases_reported)
    p_freq_state_all <- p_freq_state_all +
      geom_line(data = dat_new, aes(x = epiweek, y = new_cases_reported/max_new_cases_reported_all, color = "Reported cases"), linetype = "dotted") +
      scale_y_continuous(sec.axis = sec_axis(~.*max_new_cases_reported_all, name="cases"))
    p_freq_state_all
    
    colors_fig <- c("Estimated Delta proportion" = "#ffd544", 
                    "Est. Delta proportion without immunity" = "#484c7e", 
                    "Est. Delta proportion with 90% immunity" = "#fe6362", 
                    "Reported cases" = "royalblue1",
                    "Relative self mobility" = "black")

    p_freq_state <- ggplot(data = p %>% filter(location_rep %in% states_imp)) +
      geom_line(aes(x = time, y = mean, color = "Estimated Delta proportion")) +
      geom_ribbon(aes(x = time, ymin = CI95_lower, ymax = CI95_upper, fill = "Estimated Delta proportion"), alpha = 0.5) +
      geom_line(data = dat_new %>% filter(location_rep %in% states_imp), aes(x = epiweek, y = relative_self_mobility, color = "Relative self mobility"), linetype = "dashed") +
      # facet_wrap(~ region_rep, nrow = 1) +
      facet_grid(cols = vars(location_rep)) +

      labs(y = "proportion", x = "",
           color = "Legend") +
      theme_bw() +
      scale_color_manual(values = colors_fig)+
      scale_fill_manual(values = colors_fig) +
      guides(fill=FALSE) +
      theme(text = element_text(size=15),
            legend.title = element_blank(),
            legend.text = element_text(size = 13),
            strip.background = element_rect(
              color="white", fill="white", size=1.5, linetype="solid"
            ),
            axis.text.x = element_blank()) +
  xlim(10,25)+
      ylim(0,1)

    p_freq_state
    
    max_new_cases_reported <- max((dat_new %>% filter(location_rep %in% states_imp))$new_cases_reported)
    p_freq_state <- p_freq_state +
      geom_line(data = dat_new %>% filter(location_rep %in% states_imp), aes(x = epiweek, y = new_cases_reported/max_new_cases_reported, color = "Reported cases")) +
      scale_y_continuous(sec.axis = sec_axis(~.*max_new_cases_reported, name="cases"))
    
    # legend <- get_legend(p_freq_state)
    
    p_freq_state <- p_freq_state + theme(legend.position='bottom')
    p_freq_state
    


    # final <- ggdraw(plot_grid(plot_grid(p_freq_state, p_rho, labels = c('A', 'B'), ncol = 1, align = "v"),
    #                 plot_grid(legend, NULL, ncol = 1),
    #                 rel_widths=c(1, 0.2)))
    
    final <- plot_grid(p_freq_state, p_rho, labels = c('A', 'B'), ncol = 1, align = "v")
    final
  
  ## posterior predictive check
  po <- rstan::extract(fit)
  pp <- as.data.table(reshape2::melt( po$Y_delta_hat_post) )
  setnames(pp, 2:3, c('time','region'))
  pp <- pp[, list(value=quantile(value, p=c(.5,.025,.975)), stat=c('PP_M','PP_CL','PP_CU')), by=c('time','region')]
  pp <- dcast.data.table(pp, region+time~stat, value.var='value')
  tmp <- as.data.table(dat_full[, c("delta_count")])
  setnames(tmp, 1, c('observed'))
  pp <- cbind(pp,tmp)
  
  print( 'posterior predictive check' )
  print(
    round(nrow(subset(pp, !observed<PP_CL & !observed>PP_CU))/nrow(pp)*100)
  )
  
  ##  posterior predictive scatter plot
  p_post <- ggplot(pp, aes(x=observed, y=PP_M)) +
    geom_errorbar(aes(ymin=PP_CL, ymax=PP_CU), colour='grey80') +
    geom_point(aes(colour=factor(region))) +
    geom_abline(slope=1, intercept=0) +
    scale_x_log10(expand=c(0,0)) +
    scale_y_log10(expand=c(0,0)) +
    theme_bw() +
    theme(legend.position='') +
    labs(x='Observed',y='Posterior Predicted',colour='region',
         title = "Posterior predictive check")
  
  ## posterior predictive check (another)
  # calculating the uncertainty in observed proportion of Delta (with Beta(1,1) uninformative prior)
  dat_full$p_obs_lower <- NA
  dat_full$p_obs_upper <- NA
  
  for(i in 1:nrow(dat_full)){
    if(!is.na(dat_full$p_obs[i])){
      dat_full$p_obs_lower[i] = qbeta(0.025, 1+dat_full$delta_count[i], 1+dat_full$total[i]-dat_full$delta_count[i])
      dat_full$p_obs_upper[i] = qbeta(0.975, 1+dat_full$delta_count[i], 1+dat_full$total[i]-dat_full$delta_count[i])
    }
  }
  
  colors1 <- c("Estimated Delta proportion" = "grey50", 
              "Observed Delta proportion" = "royalblue1")
  
  p_post1 <- ggplot(data = p) +
    geom_line(aes(x = time, y = mean, color = "Estimated Delta proportion")) +
    geom_ribbon(aes(x = time, ymin = CI95_lower, ymax = CI95_upper, fill = "Estimated Delta proportion"), alpha = 0.5) +
    geom_point(data = dat_full, aes(x = epiweek, y = p_obs, color = "Observed Delta proportion"), size = 0.7) +
    geom_linerange(data = dat_full, aes(x = epiweek, y = p_obs, ymin = p_obs_lower, ymax = p_obs_upper, color = "Observed Delta proportion")) +
    facet_wrap(~location_rep) +
    labs(y = "proportion",
         color = "Legend") +
    theme_bw() +
    scale_color_manual(values = colors1)+
    scale_fill_manual(values = colors1) +
    guides(fill=FALSE) +
    theme(text = element_text(size=15),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.position="bottom",
          strip.background = element_rect(
            color="white", fill="white", size=1.5, linetype="solid"),
          strip.text.x = element_text(size = 8))+
    xlim(10,25)+
    ylim(0,1)
  
  p_post1
  

  ## prior predictive check
  po <- rstan::extract(fit)
  pp <- as.data.table(reshape2::melt( po$Y_delta_hat_prior) )
  setnames(pp, 2, c('time'))
  pp <- pp[, list(value=quantile(value, p=c(.5,.025,.975)), stat=c('PP_M','PP_CL','PP_CU')), by=c('time')]
  pp <- dcast.data.table(pp, time~stat, value.var='value')
  tmp <- as.data.table(dat_full[, c("delta_count")])
  setnames(tmp, 1, c('observed'))
  pp <- cbind(pp,tmp)

  print( 'prior predictive check' )
  print(
    round(nrow(subset(pp, !observed<PP_CL & !observed>PP_CU))/nrow(pp)*100)
  )
  
  
  
  ##  prior predictive scatter plot
  p_prior <- ggplot(pp, aes(x=observed, y=PP_M)) +
    geom_errorbar(aes(ymin=PP_CL, ymax=PP_CU), colour='grey80') +
    geom_jitter(width=0.05) +
    geom_abline(slope=1, intercept=0) +
    scale_x_log10(expand=c(0,0)) +
    scale_y_log10(expand=c(0,0)) +
    theme_bw() +
    theme(legend.position='') +
    labs(x='Observed',y='Prior Predicted',
         title = "Prior predictive check")
  
  
  ## prior predictive check (proportion over time)
  po <- rstan::extract(fit)
  pp <- as.data.table(reshape2::melt( po$p_prior) )
  setnames(pp, 2, c('time'))
  pp <- pp[, list(value=quantile(value, p=c(.5,.025,.975)), stat=c('PP_M','PP_CL','PP_CU')), by=c('time')]
  pp <- dcast.data.table(pp, time~stat, value.var='value')
  
  
  ##  prior predictive scatter plot (proportion over time)
  p_prior1 <- ggplot(pp, aes(x=time, y=PP_M)) +
    geom_errorbar(aes(ymin=PP_CL, ymax=PP_CU), colour='grey80') +
    geom_point() +
    theme_bw() +
    theme(legend.position='') +
    labs(x='Time',y='Proportion of Delta samples over time',
         title = "Prior predictive check")
  
  final1 <- plot_grid(p_post1, p_prior1, labels = c('A', 'B'), ncol = 1)

  return(list(final, final1, p_freq_state_all, p_rho_all))
}    
  
output <- plot_fun(fit_temp_sm_time_plot)


pdf(file = "outputs/results_England.pdf", # The directory you want to save the file in
    width = 17, # The width of the plot in inches
    height = 8) # The height of the plot in inches

output[[1]]

dev.off()

pdf(file = "outputs/results_England_supp.pdf", # The directory you want to save the file in
    width = 18, # The width of the plot in inches
    height = 8) # The height of the plot in inches

output[[3]]
output[[4]]

dev.off()

pdf(file = "outputs/results_England_supp1.pdf", # The directory you want to save the file in
    width = 18, # The width of the plot in inches
    height = 12) # The height of the plot in inches

output[[2]]

dev.off()

