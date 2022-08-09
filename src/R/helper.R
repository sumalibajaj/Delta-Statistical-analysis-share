
prepare_stan_data <- function(
    covariates_to_include,
    training_epiweeks, testing_epiweeks,
    last_training_epiweek, number_testing_epiweeks,
    dat_og, missing_value_code=-99) {
  
  # get all useful data: separate into training and testing sets
  # delta count
  count_delta <- dat_og %>% 
    select(delta_count, region, epiweek) %>% 
    pivot_wider(names_from=region, values_from=delta_count,
                id_cols=epiweek)
  count_delta[is.na(count_delta)] <- missing_value_code
  count_delta_training <- count_delta %>% 
    filter(epiweek %in% training_epiweeks) %>% 
    select(-epiweek)
  count_delta_testing <- count_delta %>% 
    filter(epiweek %in% testing_epiweeks) %>% 
    select(-epiweek)
  
  # total count
  count_total <- dat_og %>% 
    select(total_count, region, epiweek) %>% 
    pivot_wider(names_from=region, values_from=total_count,
                id_cols=epiweek)
  count_total[is.na(count_total)] <- missing_value_code
  count_total_training <- count_total %>% 
    filter(epiweek %in% training_epiweeks) %>% 
    select(-epiweek)
  count_total_testing <- count_total %>% 
    filter(epiweek %in% testing_epiweeks) %>% 
    select(-epiweek)
  
  # include in likelihood? highlight NAs
  include_lik_training <- count_total_training
  include_lik_training[include_lik_training!=missing_value_code] <- 1
  include_lik_training[include_lik_training==missing_value_code] <- 0
  include_lik_testing <- count_total_testing
  include_lik_testing[include_lik_testing!=missing_value_code] <- 1
  include_lik_testing[include_lik_testing==missing_value_code] <- 0
  
  # get nuts region
  nuts_vec <- dat_og %>% 
    select(region, nuts) %>% 
    unique()
  
  cnames_equal <- colnames(count_delta_testing)==nuts_vec$region
  if(mean(cnames_equal) != 1)
    stop("column names in nuts vector don't match those in counts.")
  
  # average vaccinations in each region since
  # there is so little temporal heterogeneity in these
  covars_df <- dat_og %>% 
    group_by(region) %>% 
    mutate(proportion_cum_1st_dose_vaccinated=mean(proportion_cum_1st_dose_vaccinated),
           proportion_cum_2nd_dose_vaccinated=mean(proportion_cum_2nd_dose_vaccinated)) %>%
    ungroup()
  
  # get covariates
  covars_df <- covars_df %>% 
    ungroup() %>% 
    select(epiweek, all_of(covariates_to_include))
  # and normalise
  normalise <- function(X){ (X - mean(X))/(sd(X))}
  X_training <- covars_df %>%
    filter(epiweek %in% training_epiweeks)
  means_df <- X_training %>% 
    summarise_all(mean)
  sd_df <- X_training %>% 
    summarise_all(sd)
  X_training <- X_training %>% 
    mutate_all(normalise)
  X_testing <- covars_df %>%
    filter(epiweek %in% testing_epiweeks) %>% 
    mutate_all(normalise)
  if(!("epiweek" %in% covariates_to_include)) {
    X_training <- X_training %>% 
      select(-epiweek)
    X_testing <- X_testing %>% 
      select(-epiweek)
  }

  stan_dat <-  list(T=last_training_epiweek,
                    T1 = number_testing_epiweeks,
                    I=length(unique(dat_og$region)),
                    J=length(unique(dat_og$nuts)),
                    Y_delta=count_delta_training,
                    Y_total=count_total_training,
                    include_lik=include_lik_training, 
                    nuts = nuts_vec$nuts,
                    ncovs = length(covariates_to_include), 
                    X = X_training,
                    Y_delta_extra=count_delta_testing, Y_total_extra=count_total_testing, 
                    include_lik_extra=include_lik_testing,
                    X_extra = X_testing,
                    means_df=means_df,
                    sd_df=sd_df)
  stan_dat
}

get_covariates_to_include <- function(covariate_id, covariates_df) {
  covars_indicator <- covariates_df %>% 
    filter(id==covariate_id) %>% 
    select(-id)
  covar_names <- colnames(covars_indicator)
  covars_indicator <- as.numeric(covars_indicator[1, ])
  covariates_to_include <- covar_names[covars_indicator==1]
  covariates_to_include
}

get_loglikelihood <- function(fit, missing_data_code=-99) {
  log_like <- extract_log_lik(fit, "loglikelihood_extra")
  keep <- which(log_like[1, ] != missing_data_code)
  log_like <- log_like[, keep]
  log_like
}

prepare_stan_data_full <- function(
    dat_og, missing_value_code=-99) {
  
  # get all useful data: separate into training and testing sets
  # delta count
  count_delta <- dat_og %>% 
    select(delta_count, region, epiweek) %>% 
    pivot_wider(names_from=region, values_from=delta_count,
                id_cols=epiweek) %>% 
    select(-epiweek)
  count_delta[is.na(count_delta)] <- missing_value_code
  
  # total count
  count_total <- dat_og %>% 
    select(total_count, region, epiweek) %>% 
    pivot_wider(names_from=region, values_from=total_count,
                id_cols=epiweek) %>% 
    select(-epiweek)
  count_total[is.na(count_total)] <- missing_value_code
  
  # include in likelihood? highlight NAs
  include_lik <- count_total
  include_lik[include_lik!=missing_value_code] <- 1
  include_lik[include_lik==missing_value_code] <- 0
 
  # get nuts region
  nuts_vec <- dat_og %>% 
    select(region, nuts) %>% 
    unique()
  
  cnames_equal <- colnames(count_delta)==nuts_vec$region
  if(mean(cnames_equal) != 1)
    stop("column names in nuts vector don't match those in counts.")
  
  # average vaccinations in each region since
  # there is so little temporal heterogeneity in these
  covars_df <- dat_og %>% 
    group_by(region) %>% 
    mutate(proportion_cum_1st_dose_vaccinated=mean(proportion_cum_1st_dose_vaccinated),
           proportion_cum_2nd_dose_vaccinated=mean(proportion_cum_2nd_dose_vaccinated)) %>%
    ungroup()
  
  # get covariates
  normalise <- function(X){ (X - mean(X))/(sd(X))}
  covars_df <- covars_df %>% 
    ungroup() %>% 
    select(epiweek, all_of(covariates_to_include)) %>% 
    mutate_all(normalise)
  
  if(!("epiweek" %in% covariates_to_include)) {
    covars_df <- covars_df %>% 
      select(-epiweek)
  }
  
  stan_dat <-  list(T=nrow(count_delta),
                    I=length(unique(dat_og$region)),
                    J=length(unique(dat_og$nuts)),
                    Y_delta=count_delta,
                    Y_total=count_total,
                    include_lik=include_lik, 
                    nuts = nuts_vec$nuts,
                    ncovs = length(covariates_to_include), 
                    X = covars_df)
  stan_dat
}
