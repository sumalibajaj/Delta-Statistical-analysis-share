.PHONY: all

MODELS := $(shell seq 1 15)
FITS := $(addsuffix .rds, $(addprefix data/processed/stan_fits/fit_model_comparison_, $(MODELS)))
# $(info VAR="$(FITS)")

all: data/processed/stan_fits/fit_best.rds\
	data/processed/data_summary.csv\
	data/processed/growth_quantiles.rds

plots: outputs/relative_mobility_delta.pdf\
	outputs/relative_mobility_delta_nothreshold.pdf\
	outputs/number_introduction_vs_delta.pdf\
	outputs/cum_2_dose_vs_delta.pdf\
	outputs/mobility_growth.pdf\
	outputs/vaccines_2_growth.pdf\
	outputs/vaccination_prop_delta.pdf\
	outputs/vaccination_prop_delta_no_threshold.pdf\
	outputs/growth_panel.pdf\
	outputs/growth_panel_mobility_only.pdf\
	outputs/prior_and_posterior_predictive.png\
	outputs/plot_counterfactual.pdf
	

data/processed/England_weekly_cleaned.rds: src/R/cleaning_England.R\
	data/raw/updated_delta.csv\
	data/raw/utla_2022-04-13.csv
	Rscript $<
data/processed/England_weekly_processed.rds: src/R/combine_with_vaccinations.R\
	data/processed/England_weekly_cleaned.rds\
	data/raw/utla_2022-04-13.csv
	Rscript $<
	
data/processed/data_summary.csv: src/R/data_summary.R\
	data/processed/England_weekly_processed.rds
	Rscript $<

data/processed/covariates_to_include.rds: src/R/prepare_covariates_for_comparison.R
	Rscript $<
	
$(FITS): data/processed/stan_fits/fit_model_comparison_%.rds: src/R/model_comparison_vars.R\
	data/processed/England_weekly_processed.rds\
	src/stan/growth_comparison.stan
	Rscript $< $* 50000 4 10

data/processed/model_comparison_results.rds: src/R/model_comparison.R\
	data/processed/covariates_to_include.rds\
	$(FITS)
	Rscript $<
data/processed/model_comparison_results.csv: data/processed/model_comparison_results.rds
	
data/processed/stan_fits/fit_best.rds: src/R/fit_best_model.R\
	data/processed/England_weekly_processed.rds\
	data/processed/model_comparison_results.rds
	Rscript $<

data/processed/growth_quantiles.rds: src/R/extract_growth_quantiles.R\
	data/processed/England_weekly_processed.rds\
	data/processed/stan_fits/fit_best.rds
	Rscript $<
	
outputs/growth_panel.pdf: src/R/plot_all_growth_panel.R\
	data/processed/ranked_vaccine_2.rds\
	data/processed/ranked_mobility.rds
	Rscript $<

outputs/growth_panel_mobility_only.pdf: src/R/plot_all_growth_panel_mobility.R\
	data/processed/ranked_mobility.rds
	Rscript $<
outputs/growth_panel_mobility_only.png: outputs/growth_panel_mobility_only.pdf

outputs/relative_mobility_delta.pdf: src/R/plot_self_mobility.R\
	data/processed/England_weekly_processed.rds
	Rscript $<
data/processed/ranked_mobility.rds: outputs/relative_mobility_delta.pdf

outputs/relative_mobility_delta_nothreshold.pdf: src/R/plot_self_mobility_nothreshold.R\
	data/processed/England_weekly_processed.rds
	Rscript $<
	
outputs/number_introduction_vs_delta.pdf: src/R/plot_covariates.R\
	data/processed/England_weekly_processed.rds
	Rscript $<
outputs/cum_2_dose_vs_delta.pdf: outputs/number_introduction_vs_delta.pdf

outputs/mobility_growth.pdf: src/R/plot_self_mobility_growth.R\
	data/processed/England_weekly_processed.rds
	Rscript $<

outputs/vaccines_2_growth.pdf: src/R/plot_vaccine_2_growth.R\
	data/processed/England_weekly_processed.rds
	Rscript $<

outputs/vaccination_prop_delta_no_threshold.pdf: src/R/plot_vaccination_2_nothreshold.R\
	data/processed/England_weekly_processed.rds
	Rscript $<

outputs/vaccination_prop_delta.pdf: src/R/plot_vaccination_2.R\
	data/processed/England_weekly_processed.rds
	Rscript $<
data/processed/ranked_vaccine_2.rds: outputs/vaccination_prop_delta.pdf
	
data/processed/stan_fits/fit_best_extended.rds: src/R/fit_best_model_extended.R\
	data/processed/England_weekly_processed.rds\
	src/stan/growth_full_extended.stan
	Rscript $<

outputs/posterior_predictive.pdf: src/R/plot_posterior_and_prior_predictive.R\
	data/processed/England_weekly_processed.rds\
	data/processed/stan_fits/fit_best_extended.rds
	Rscript $< 1
outputs/posterior_predictive.png: outputs/posterior_predictive.pdf
data/processed/posterior_predictive.rds: outputs/posterior_predictive.pdf

outputs/prior_predictive.pdf: src/R/plot_posterior_and_prior_predictive.R\
	data/processed/England_weekly_processed.rds\
	data/processed/stan_fits/fit_best_extended.rds
	Rscript $< 0
outputs/prior_predictive.png: outputs/prior_predictive.pdf 
data/processed/prior_predictive.rds: outputs/prior_predictive.pdf 

outputs/prior_and_posterior_predictive.png: src/R/plot_both_posterior_prior_predictive.R\
	data/processed/prior_predictive.rds\
	data/processed/posterior_predictive.rds
	Rscript $<

outputs/posterior_predictive_growth.pdf: src/R/plot_posterior_growth_location.R\
	data/processed/England_weekly_processed.rds\
	data/processed/stan_fits/fit_best_extended.rds
	Rscript $<
outputs/posterior_predictive_growth.png: outputs/posterior_predictive_growth.pdf

data/processed/stan_fits/fit_best_shorter.rds: src/R/fit_best_model_shorter_dataset.R\
	data/processed/England_weekly_processed.rds\
	src/stan/growth_comparison.stan
	Rscript $<
data/processed/convergence_check_best_shorter.rds: data/processed/stan_fits/fit_best_shorter.rds
outputs/plot_counterfactual.pdf: src/R/counterfactual_plots.R\
	data/processed/England_weekly_processed.rds\
	data/processed/stan_fits/fit_best_shorter.rds
	Rscript $<
	
