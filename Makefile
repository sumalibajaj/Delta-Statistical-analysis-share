.PHONY: all
	
all: outputs/results_England.pdf\
	outputs/results_US.pdf\
	data/processed/England_MC.RData\
	data/processed/US_MC.RData\
	outputs/simulation.pdf

data/processed/England_weekly_mob.rds: src/R/cleaning_England.R\
	data/raw/master_total_updated.csv
	Rscript $<

data/processed/England_weekly_vac.rds: data/processed/England_weekly_mob.rds

data/processed/England_weekly_intro.rds: data/processed/England_weekly_vac.rds

data/processed/England_final.RData: src/R/modelsetup_England.R\
	data/processed/England_weekly_mob.rds\
	data/processed/England_weekly_intro.rds\
	src/stan/final.stan
	Rscript $<
	
outputs/results_England.pdf: src/R/plottingresults_England.R\
	data/processed/England_final.RData
	Rscript $<

outputs/results_England_supp.pdf: outputs/results_England.pdf

outputs/results_England_supp1.pdf: outputs/results_England.pdf

data/processed/US_weekly.rds: src/R/cleaning_US.R\
	data/raw/master_total_updated.csv
	Rscript $<
	
data/processed/US_weekly_intro.rds: data/processed/US_weekly.rds

data/processed/US_final.RData: src/R/modelsetup_US.R\
	data/processed/US_weekly.rds\
	data/processed/US_weekly_intro.rds\
	src/stan/final.stan
	Rscript $<
	
outputs/results_US.pdf: src/R/plottingresults_US.R\
	data/processed/US_final.RData
	Rscript $<
	
outputs/results_US_supp.pdf: outputs/results_US.pdf

outputs/results_US_supp1.pdf: outputs/results_US.pdf

data/processed/England_MC.RData: src/R/modelsetup_England_MC.R\
	data/processed/England_weekly_mob.rds\
	src/stan/final_MC_novar.stan\
	src/stan/final_MC_var.stan
	Rscript $<
	
data/processed/US_MC.RData: src/R/modelsetup_US_MC.R\
	data/processed/US_weekly.rds\
	src/stan/final_MC_novar.stan\
	src/stan/final_MC_var.stan
	Rscript $<
	
outputs/simulation.pdf: src/R/diagnostic_simulation.R\
	src/stan/final.stan
	Rscript $<
