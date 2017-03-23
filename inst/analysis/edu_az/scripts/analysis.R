source("~/repo/mr_frailty/R/functions.R")

dat <- simulate_data(
	age_summary = age_summary,
	snp_beta = edu_snps$beta.exposure,
	snp_se = edu_snps$se.exposure,
	snp_eaf = edu_snps$eaf.exposure,
	exposure_mean = 14,
	exposure_sd = 3.51,
	exposure_lb = 5,
	exposure_ub = 30,
	outcome_prevalence_function = az_prevalence,
	survival_function = edu_survival,
	min_age = 40,
	max_age = 100,
	sample_size_multiplier=11
)

res <- analyse_data(dat, age_summary)

plot_quantiles(exposure = dat$dat$exposure, outcome = dat$dat$alive, ylab = "Survival", xlab="Education", exposure_breaks=5)

plot_quantiles(exposure = dat$dat$age, outcome = dat$dat$alive, ylab = "Survival", xlab="Age", exposure_breaks=10, mediator= dat$dat$exposure, mediator_breaks=3, mediator_name="Education")

plot_quantiles(exposure = dat$dat$age, outcome = dat$dat$alive, ylab = "Survival", xlab="Age", exposure_breaks=10, mediator= dat$dat$cc, mediator_name="AZ")

