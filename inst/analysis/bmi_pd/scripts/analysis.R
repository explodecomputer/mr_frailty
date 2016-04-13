## ---- setup ----

library(ggplot2)
library(dplyr)

load("~/repo/mr_frailty/inst/analysis/bmi_pd/results/model1.RData")
res1 <- res
res1$model <- "model1"
load("~/repo/mr_frailty/inst/analysis/bmi_pd/results/model2.RData")
res2 <- res
res2$model <- "model2"
res <- rbind(res1, res2)


res_plot <- subset(res, test %in% c("2sls", "obs", "MR Egger", "Inverse variance weighted"))

res_plot$test <- as.factor(res_plot$test)
levels(res_plot$test) <- c("2SLS", "IVW", "MR Egger", "Obs assoc")


## ---- method_comparisons ----

dat <- group_by(res_plot, test, model) %>% dplyr::summarise(b=mean(beta), se=sd(beta)/n(), sd=sd(beta))


ggplot(subset(res_plot, test != "grs"), aes(x=beta)) +
geom_density(aes(fill=model), alpha=0.2) +
geom_vline(data=subset(dat, test != "grs"), aes(xintercept=b, colour=model)) +
geom_vline(xintercept=0, linetype="dashed") +
facet_grid(test ~ ., scale="free_y") +
scale_fill_brewer(type="qual")
ggsave("~/repo/mr_frailty/inst/analysis/bmi_pd/images/method_comparison_both_models.pdf")


## ---- method_comparisons_model1 ----

dat <- group_by(res_plot, test, model) %>% dplyr::summarise(b=mean(beta), se=sd(beta)/n(), sd=sd(beta))

ggplot(subset(res_plot, test != "grs" & model == "model1"), aes(x=beta)) +
geom_density(alpha=0.2) +
geom_vline(data=subset(dat, test != "grs" & model == "model1"), aes(xintercept=b)) +
geom_vline(xintercept=0, linetype="dashed") +
facet_grid(test ~ ., scale="free_y") +
scale_fill_brewer(type="qual")
ggsave("~/repo/mr_frailty/inst/analysis/bmi_pd/images/method_comparison_model1.pdf")




## ---- empirical_results ----

# convert res to OR change per 5 kg/m2
res_plot$beta2 <- exp(res_plot$beta * 5)

dat_ci <- dplyr::group_by(res_plot, test, model) %>% dplyr::summarise(b=mean(beta2), lower_ci=quantile(beta2, 0.05), upper_ci=quantile(beta2, 0.95))

emp_dat <- data.frame(
	b = c(0.79, 0.79, 1.00),
	lower_ci = c(0.66, 0.38, 0.89),
	upper_ci = c(0.96, 1.27, 1.12),
	test = c("IVW", "MR Egger", "Obs assoc"),
	model = "Empirical"
)

dat_ci <- rbind(dat_ci, emp_dat)
dat_ci <- subset(dat_ci, model != "model2")
dat_ci$model <- as.factor(dat_ci$model)
levels(dat_ci$model) <- c("Empirical results", "Frailty simulation")

ggplot(dat_ci, aes(y=b, x=model)) +
geom_point() +
geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0) +
facet_grid(test ~ .) +
coord_flip() +
scale_colour_brewer(type="qual") +
theme_bw() +
labs(x="", y="Odds ratio and 95% confidence intervals")
ggsave("~/repo/mr_frailty/inst/analysis/bmi_pd/images/empircal.pdf")


## ---- example_simulation ----

source("~/repo/mr_frailty/R/functions.R")

pd_incidence <- function(age, ...)
{
	pd_inc <- read.table("~/repo/mr_frailty/data-raw/pd_incidence.txt", he=T)
	pd_inc$age <- (pd_inc$age_low + pd_inc$age_upp) / 2
	pd_inc$n <- pd_inc$personyears / ((pd_inc$age_upp + pd_inc$age_low) / 2)
	pd_inc$p <- pd_inc$cases / pd_inc$n
	mod <- lm(p ~ poly(age,4), weights=n, pd_inc)
	pred <- predict(mod, data.frame(age=age))
	pred[pred < 0] <- 0
	return(pred)
}


bmi_survival <- function(age, bmi)
{
	hr_summary <- read.table("~/repo/mr_frailty/data-raw/bmi_hr.txt", he=T)
	hr_summary$lhr <- log(hr_summary$hr)
	hr_summary$lhr_se <- (log(hr_summary$ci_upp) - log(hr_summary$ci_low)) / 3.92
	
	l <- list()
	for(i in 1:nrow(hr_summary))
	{
		l[[i]] <- data.frame(
			bmi = seq(hr_summary$bmi_low[i], hr_summary$bmi_upp[i], length.out=hr_summary$n[i]),
			hr = exp(rnorm(hr_summary$n[i], hr_summary$lhr[i], hr_summary$lhr_se[i]))
		)
	}
	hr_summary <- bind_rows(l)

	mod <- lm(hr ~ poly(bmi, 4), hr_summary)
	newdat <- data.frame(bmi = bmi)
	hr <- predict(mod, newdat)
	survival <- (1 - gompertz_makeham_cdf(age)) ^ hr
	return(survival)
}


set.seed(12)
demog <- read.csv("~/repo/mr_frailty/data-raw/pd_demographics.csv")
age_summary <- get_age_summary(demog, "Cases", "Controls", "Case_age_mean", "Control_age_mean", "Case_age_sd", "Control_age_sd")
bmi_snps <- read.table("~/repo/mr_frailty/data-raw/bmi_2015_clumped.txt", he=T)
bmi_snps$b <- bmi_snps$b
bmi_snps_mean <- 25
bmi_snps_sd <- 4.18

dat <- simulate_ages(age_summary$gn[3], age_summary$gm[3], age_summary$gs[3], max_age=100, min_age=40, sample_size_multiplier=4)
dat$cc <- simulate_events(dat$age, NULL, pd_incidence)
snps <- simulate_snps(nrow(dat), bmi_snps$Freq1.Hapmap)
dat$bmi <- simulate_exposure(nrow(dat), snps, bmi_snps$b * bmi_snps_sd, bmi_snps_mean, bmi_snps_sd, 15, 60)
dat$alive <- simulate_events(dat$age, dat$bmi, bmi_survival)
dat$grs <- snps %*% bmi_snps$b
index <- sample_cases_controls(dat, age_summary)


## ---- plot_sim_mort_bmi_alive ----

plot_quantiles(exposure = dat$bmi, outcome = dat$alive, ylab = "Mortality", xlab="BMI", exposure_breaks=10)

## ---- plot_sim_mort_age_bmi ----

plot_quantiles(exposure = dat$age, outcome = dat$alive, ylab = "Mortality", xlab="Age", exposure_breaks=10, mediator= dat$bmi, mediator_breaks=4, mediator_name="BMI")

## ---- plot_sim_mort_age_pd ----

plot_quantiles(exposure = dat$age, outcome = dat$alive, ylab = "Mortality", xlab="Age", exposure_breaks=10, mediator= dat$cc, mediator_name="PD")


## ---- plot_sim_pd_bmi ----

plot_quantiles(exposure = dat$bmi, outcome = dat$cc, ylab = "PD risk", xlab="BMI", exposure_breaks=10)





