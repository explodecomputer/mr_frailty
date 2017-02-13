## ---- setup ----

library(ggplot2)
library(dplyr)
library(systemfit)
library(Runuran)

load("../results/model1.RData")
res1 <- res
res1$model <- "model1"
res1$modname <- "BMI ~ SNP; nonlinear"
load("../results/model2.RData")
res2 <- res
res2$model <- "model2"
res2$modname <- "BMI ~ SNP; large frailty effect"
load("../results/model3.RData")
res3 <- res
res3$model <- "model3"
res3$modname <- "BMI ~ SNP; not sure"
load("../results/model4.RData")
res4 <- res
res4$model <- "model4"
res4$modname <- "BMI ~ SNP; linear"
load("../results/model5.RData")
res5 <- res
res5$model <- "model5"
res5$modname <- "BMI ~ SNP + age; linear"
res <- rbind(res1, res2, res3, res4, res5)


res_plot <- subset(res, test %in% c("2sls", "obs", "MR Egger", "Inverse variance weighted"))

res_plot$test <- as.factor(res_plot$test)
levels(res_plot$test) <- c("2SLS", "IVW", "MR Egger", "Obs assoc")


## ---- method_comparisons ----

mc_dat <- group_by(res_plot, test, model) %>% dplyr::summarise(b=mean(beta), se=sd(beta)/n(), sd=sd(beta), modname=first(modname))


ggplot(subset(res_plot), aes(x=beta)) +
geom_density(aes(fill=modname), alpha=0.5) +
geom_vline(data=subset(mc_dat), aes(xintercept=b, colour=modname)) +
geom_vline(xintercept=0, linetype="dashed") +
facet_grid(test ~ ., scale="free_y") +
scale_fill_brewer(type="qual") +
scale_colour_brewer(type="qual") 
ggsave("../images/method_comparison_models_with_age.pdf")


dc <- group_by(res_plot, model, test) %>%
dplyr::summarise(b=exp(quantile(beta * 5, 0.5)), lci=exp(quantile(beta * 5, 0.05)), uci=exp(quantile(beta * 5, 0.95)))

dc
head(res_plot)

ggplot(dc, aes(x=b, y=model)) +
geom_point() +
geom_errorbarh(aes(xmin=lci, xmax=uci), height=0) +
facet_grid(test ~ .)



## ---- method_comparisons_model1 ----

mc1_dat <- group_by(res_plot, test, model) %>% dplyr::summarise(b=mean(beta), se=sd(beta)/n(), sd=sd(beta))
mc1_dat$pval <- pnorm(abs(mc1_dat$b / mc1_dat$se), low=FALSE)

ggplot(subset(res_plot, test != "grs" & model == "model1"), aes(x=beta)) +
geom_density(alpha=0.2) +
geom_vline(data=subset(mc1_dat, test != "grs" & model == "model1"), aes(xintercept=b)) +
geom_vline(xintercept=0, linetype="dashed") +
facet_grid(test ~ ., scale="free_y") +
scale_fill_brewer(type="qual") +
labs(x="Estimate (log[OR] per kg/m2)")
# ggsave("../images/method_comparison_model1.pdf")


## ---- empirical_results ----

# convert res to OR change per 5 kg/m2

res_plot <- subset(res_plot, model != "model5")
res_plot$beta2 <- exp(res_plot$beta * 5)

dat_ci <- dplyr::group_by(res_plot, test, model) %>% dplyr::summarise(b=mean(beta2), lower_ci=quantile(beta2, 0.05), upper_ci=quantile(beta2, 0.95))

emp_dat <- data.frame(
	b = c(0.82, 0.76, 1.00),
	lower_ci = c(0.69, 0.51, 0.89),
	upper_ci = c(0.98, 1.14, 1.12),
	test = c("IVW", "MR Egger", "Obs assoc"),
	model = "Empirical"
)

dat_ci <- bind_rows(subset(dat_ci, model == "model4" & test != "2SLS"), emp_dat)
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
labs(x="", y="Odds ratio per 5 kg/m2 and 95% confidence intervals")
# ggsave("../images/empircal_model4.pdf")
ggsave("../images/empirical_model4.pdf")


## ---- bmi_hr ----

hr_summary <- read.table("../../../../data-raw/bmi_hr.txt", he=T)
hr_summary$x <- paste0(hr_summary$bmi_low, " - ", hr_summary$bmi_upp)

ggplot(hr_summary, aes(x=x, y=hr, ymin=ci_low, ymax=ci_upp)) +
geom_point() +
geom_errorbar(width=0) +
theme(axis.text.x=element_text(angle=90)) +
labs(y="Hazard ratio", x="BMI")


## ---- pd_hr ----

pd_inc <- read.table("../../../../data-raw/pd_incidence.txt", he=T)
pd_inc$age <- (pd_inc$age_low + pd_inc$age_upp) / 2
pd_inc$n <- pd_inc$personyears / ((pd_inc$age_upp + pd_inc$age_low) / 2)
pd_inc$p <- pd_inc$cases / pd_inc$n
pd_inc$x <- with(pd_inc, paste0(age_low, " - ", age_upp))

ggplot(pd_inc, aes(x=x, y=p)) +
geom_point() +
theme(axis.text.x=element_text(angle=90)) +
labs(y="p(PD diagnosis)", x="Age range")


## ---- example_simulation ----

source("../../../../R/functions.R")

pd_incidence <- function(age, ...)
{
	pd_inc <- read.table("../../../../data-raw/pd_incidence.txt", he=T)
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
	hr_summary <- read.table("../../../../data-raw/bmi_hr.txt", he=T)
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
demog <- read.csv("../../../../data-raw/pd_demographics.csv")
age_summary <- get_age_summary(demog, "Cases", "Controls", "Case_age_mean", "Control_age_mean", "Case_age_sd", "Control_age_sd")
bmi_snps <- read.table("../../../../data-raw/locke_clumped_0.001.txt", he=T)
bmi_snps$beta.exposure <- bmi_snps$beta.exposure
bmi_snps_mean <- 25
bmi_snps_sd <- 4.18

dat <- simulate_ages(age_summary$gn[3], age_summary$gm[3], age_summary$gs[3], max_age=100, min_age=40, sample_size_multiplier=4)
dat$cc <- simulate_events(dat$age, NULL, pd_incidence)
snps <- simulate_snps(nrow(dat), bmi_snps$eaf.exposure)
dat$bmi <- simulate_exposure(nrow(dat), snps, bmi_snps$beta.exposure * bmi_snps_sd, bmi_snps_mean, bmi_snps_sd, 15, 60)
dat$alive <- simulate_events(dat$age, dat$bmi, bmi_survival)
dat$grs <- snps %*% bmi_snps$beta.exposure
index <- sample_cases_controls(dat, age_summary)


## ---- plot_sim_mort_bmi_alive ----

plot_quantiles(exposure = dat$bmi, outcome = dat$alive, ylab = "Mortality", xlab="BMI", exposure_breaks=10)

## ---- plot_sim_mort_age_bmi ----

plot_quantiles(exposure = dat$age, outcome = dat$alive, ylab = "Mortality", xlab="Age", exposure_breaks=10, mediator= dat$bmi, mediator_breaks=4, mediator_name="BMI")

## ---- plot_sim_mort_age_pd ----

plot_quantiles(exposure = dat$age, outcome = dat$alive, ylab = "Mortality", xlab="Age", exposure_breaks=10, mediator= dat$cc, mediator_name="PD")


## ---- plot_sim_pd_bmi ----

plot_quantiles(exposure = dat$bmi, outcome = dat$cc, ylab = "PD risk", xlab="BMI", exposure_breaks=10)


## ---- baseline_gm ----

d <- data.frame(age=40:100)
d$survival <- survival <- (1 - gompertz_makeham_cdf(d$age))
qplot(data=d, x=age, y=survival, geom="line")



## ---- example_mr ----

a <- summary(glm(cc ~ bmi, dat[index,], family="binomial"))
b <- summary(systemfit(cc ~ bmi, "2SLS", inst = ~ grs, data = dat[index,]))

ex <- as.data.frame(rbind(
	coefficients(a)[2,c(1,2,4)],
	coefficients(b)[2,c(1,2,4)]
))

names(ex) <- c("beta", "se", "pval")
ex$method <- c("Obs assoc", "2SLS")


ss <- get_summary_stats(dat, snps, index)
mres <- do_mr(bmi_snps$beta.exposure, ss$b, bmi_snps$se.exposure, ss$se)
mres <- data.frame(beta = mres$b, se = mres$se, pval = mres$pval, method = mres$method)
ex <- rbind(ex, subset(mres, method %in% c("Inverse variance weighted", "MR Egger")))

ggplot(ex, aes(x = method, y = beta, ymin = beta+se*1.96, ymax = beta-se*1.96)) +
geom_point() +
geom_errorbar(width=0) +
coord_flip() +
labs(x="Method", y = "Estimate (log[OR] change per kg/m2")


## ---- nalls_age_dist ----

demog <- read.csv("../../../../data-raw/pd_demographics.csv")
age_summary <- get_age_summary(demog, "Cases", "Controls", "Case_age_mean", "Control_age_mean", "Case_age_sd", "Control_age_sd")
names(age_summary) <- c("Sample size", "Mean age", "SD", "Case/control")


## ---- simulated_age_dist ----

ggplot(dat[index,], aes(x=age)) +
geom_density(alpha=0.2, aes(y=..count.., fill=factor(cc))) +
labs(fill="PD status")
