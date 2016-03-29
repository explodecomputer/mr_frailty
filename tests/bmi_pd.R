# exposure, outcome and mortality

# exposure has no effect on outcome

# exposure has an effect on mortality

# outcome occurs based on age

# mortality is a function of age and exposure



# exposure = snp(s)

# mortality = age + exposure

# outcome = age


# simulate population with snp(s)

# give them ages from a specified age distribution


# total_n
# alive_n

# snp_freq
# snp_beta
# snp_effect_allele



pd_incidence <- function(age, ...)
{
	pd_inc <- read.table("data-raw/pd_incidence.txt", he=T)
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
	hr_summary <- read.table("data-raw/bmi_hr.txt", he=T)
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


# bmi_survival <- function(age, bmi)
# {
# 	hr <- rep(1, length(bmi))
# 	hr[bmi > 27] <- 5
# 	survival <- (1 - gompertz_makeham_cdf(age)) ^ hr
# 	return(survival)
# }



demog <- read.csv("data-raw/pd_demographics.csv")
ages <- get_age_summary(demog, "Cases", "Controls", "Case_age_mean", "Control_age_mean", "Case_age_sd", "Control_age_sd")
bmi_snps <- read.table("data-raw/bmi_2015_clumped.txt", he=T)
bmi_snps$b <- bmi_snps$b * 4.18
bmi_snps_mean <- 25
bmi_snps_sd <- 4.18


nsim <- 1000
l1 <- list()
l2 <- list()
for (i in 1:nsim)
{
	message(i)

	dat <- simulate_ages(ages$gn[3], ages$gm[3], ages$gs[3], max_age=100, min_age=40, sample_size_multiplier=4)
	dat$cc <- simulate_events(dat$age, NULL, pd_incidence)
	snps <- simulate_snps(nrow(dat), bmi_snps$Freq1.Hapmap)
	dat$bmi <- simulate_exposure(nrow(dat), snps, bmi_snps$b, bmi_snps_mean, bmi_snps_sd)
	dat$alive <- simulate_events(dat$age, dat$bmi, bmi_survival)
	dat$grs <- snps %*% bmi_snps$b
	s <- sample_cases_controls(dat, age_summary)

	a <- summary(glm(cc ~ bmi, subset(s), family="binomial"))
	b <- summary(glm(cc ~ grs, subset(s), family="binomial"))

	l1[[i]] <- coefficients(a)[2,]
	l2[[i]] <- coefficients(b)[2,]

}

save(l1, l2, file="tests/res.RData")
