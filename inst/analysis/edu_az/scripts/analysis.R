source("~/repo/mr_frailty/R/functions.R")

az_prevalence <- function(age, ...)
{
	prev <- read.table("~/repo/mr_frailty/data-raw/az_prevalence.txt", he=T)
	# Linear interpolation of prevalence for each age supplied
	pred <- approx(x = prev$age, y = prev$prevalence, xout=age, rule=2)
	return(pred$y)
}

edu_survival <- function(age, edu)
{
	edu_sd <- sd(edu, na.rm=TRUE)
	edu <- scale(edu)
	eff <- log(1/0.96) / edu_sd
	survival <- (1 - gompertz_makeham_cdf(age)) * exp(edu * eff)
	survival[survival > 1] <- 1
	survival[survival < 0] <- 0
	return(survival)
}


# Get age distributions for cases and controls

demog <- read.csv("~/repo/mr_frailty/data-raw/alz_gwas_ages.csv")
age_summary <- get_age_summary(demog, "cases", "controls", "case_age_mean", "control_age_mean", "case_age_sd", "control_age_sd")


# Get instruments for exposure - effect, se, effect allele frequency

edu_snps <- read.csv("~/repo/mr_frailty/data-raw/EduYears_Alz_ResultsBySNP.csv")

dat <- simulate_data(
	age_summary = age_summary,
	snp_beta = edu_snps$beta.exposure,
	snp_se = edu_snps$se.exposure,
	snp_eaf = edu_snps$eaf.exposure,
	exposure_mean = 12,
	exposure_sd = 3.5,
	exposure_lb = 5.01,
	exposure_ub = 20,
	outcome_prevalence_function = az_prevalence,
	survival_function = edu_survival,
	min_age = 40,
	max_age = 100,
	sample_size_multiplier=30
)

res <- analyse_data(dat, age_summary)

plot_quantiles(exposure = dat$dat$exposure, outcome = dat$dat$alive, ylab = "Survival", xlab="Education", exposure_breaks=5)

plot_quantiles(exposure = dat$dat$age, outcome = dat$dat$alive, ylab = "Survival", xlab="Age", exposure_breaks=10, mediator= dat$dat$exposure, mediator_breaks=3, mediator_name="Years of education")
ggsave("../results/example_survival.pdf")

plot_quantiles(exposure = dat$dat$age, outcome = dat$dat$alive, ylab = "Survival", xlab="Age", exposure_breaks=10, mediator= dat$dat$cc, mediator_name="AZ")

d <- dat$dat[new_sample_cases_controls(dat$dat, age_summary, 40, 100),]
plot_age_density(d)
ggsave("../results/sampled_age_density.pdf")


## ---- setup ----

library(ggplot2)
library(dplyr)
library(systemfit)
library(Runuran)

load("../results/model1.RData")
res1 <- res
res1$model <- "model1"
res1$modname <- "Edu ~ SNP; best estimate"
load("../results/model2.RData")
res2 <- res
res2$model <- "model2"
res2$modname <- "Edu ~ SNP; exaggerated effect"

res <- rbind(res1, res2)

res_plot <- subset(res, test %in% c("2SLS", "Observational", "MR Egger", "Inverse variance weighted"))

res_plot$test <- as.factor(res_plot$test)
levels(res_plot$test) <- c("2SLS", "IVW", "MR Egger", "Obs assoc")


## ---- method_comparisons ----

mc_dat <- group_by(res_plot, test, model) %>% 
	dplyr::summarise(
		nsim = n(),
		b=mean(beta), 
		se=sd(beta)/sqrt(n()), 
		modname=first(modname), 
		tval=b / se, 
		pval=2 * pt(abs(tval), n()-1, low=FALSE)
)

mc_dat <- group_by(res_plot, test, model) %>% do({print(t.test(.$beta))})

ggplot(subset(res_plot, test != "2SLS"), aes(x=(beta))) +
geom_density(aes(fill=modname), alpha=0.5) +
geom_vline(data=subset(mc_dat, model %in% c("model4", "model5") & test != "2SLS"), aes(xintercept=(b), colour=modname)) +
geom_vline(xintercept=0, linetype="dashed") +
facet_wrap(~ test, scale="free") +
scale_fill_brewer(type="qual") +
scale_colour_brewer(type="qual") +
labs(x="Effect estimate", fill="Model", colour="Model")

write.csv(subset(mc_dat, modname == modname[1], select=-c(modname, model)), file="../results/table.csv")

