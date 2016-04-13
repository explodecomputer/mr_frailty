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



source("~/repo/mr_frailty/R/functions.R")
library(dplyr)
library(systemfit)


main <- function()
{
	demog <- read.csv("~/repo/mr_frailty/data-raw/pd_demographics.csv")
	age_summary <- get_age_summary(demog, "Cases", "Controls", "Case_age_mean", "Control_age_mean", "Case_age_sd", "Control_age_sd")
	bmi_snps <- read.table("~/repo/mr_frailty/data-raw/bmi_2015_clumped.txt", he=T)
	bmi_snps$b <- bmi_snps$b
	bmi_snps_mean <- 25
	bmi_snps_sd <- 4.18


	parameters <- expand.grid(sim = 1:1000)

	# Parallel
	arguments <- commandArgs(T)
	jid <- "all"
	outdir <- "./"
	if(length(arguments) > 0)
	{
		jid <- as.numeric(arguments[1])
		splits <- as.numeric(arguments[2])
		outdir <- arguments[3]
		stopifnot(all(!is.na(jid), !is.na(splits), !is.na(outdir)))

		first <- (jid - 1) * splits + 1
		last <- min(nrow(parameters), jid * splits)
		parameters <- parameters[first:last, , drop=FALSE]
	}

	# Set output file
	outfile <- paste(outdir, "/results", jid, ".RData", sep="")
	message("Running ", jid, ": ", nrow(parameters), " rows")
	message("Saving in ", outfile)


	l1 <- list()
	l2 <- list()
	l3 <- list()
	l4 <- list()
	for (i in 1:nrow(parameters))
	{
		message(i)

		dat <- simulate_ages(age_summary$gn[3], age_summary$gm[3], age_summary$gs[3], max_age=100, min_age=40, sample_size_multiplier=4)
		dat$cc <- simulate_events(dat$age, NULL, pd_incidence)
		snps <- simulate_snps(nrow(dat), bmi_snps$Freq1.Hapmap)
		dat$bmi <- simulate_exposure(nrow(dat), snps, bmi_snps$b * bmi_snps_sd, bmi_snps_mean, bmi_snps_sd, lb=15,ub=60)
		dat$alive <- simulate_events(dat$age, dat$bmi, bmi_survival)
		dat$grs <- snps %*% bmi_snps$b
		index <- sample_cases_controls(dat, age_summary)

		a <- summary(glm(cc ~ bmi, dat[index,], family="binomial"))
		b <- summary(glm(cc ~ grs, dat[index,], family="binomial"))
		c <- summary(systemfit(cc ~ bmi, "2SLS", inst = ~ grs, data = dat[index,]))

		l1[[i]] <- coefficients(a)[2,]
		l2[[i]] <- coefficients(b)[2,]
		l3[[i]] <- coefficients(c)[2,]

		ss <- get_summary_stats(dat, snps, index)
		mres <- do_mr(bmi_snps$b, ss$b, bmi_snps$se, ss$se)
		l4[[i]] <- data.frame(beta = mres$b, se = mres$se, tval = NA, pval = mres$pval, sim = parameters$sim[i], test = mres$method)
	}

	l1 <- as.data.frame(do.call(rbind, l1))
	names(l1) <- c("beta", "se", "tval", "pval")
	l1$sim <- parameters$sim
	l1$test <- "obs"

	l2 <- as.data.frame(do.call(rbind, l2))
	names(l2) <- c("beta", "se", "tval", "pval")
	l2$sim <- parameters$sim
	l2$test <- "grs"

	l3 <- as.data.frame(do.call(rbind, l3))
	names(l3) <- c("beta", "se", "tval", "pval")
	l3$sim <- parameters$sim
	l3$test <- "2sls"

	l4 <- bind_rows(l4)

	res <- rbind(l1, l2, l3, l4)

	save(res, file=outfile)

}




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


# bmi_survival <- function(age, bmi)
# {
# 	hr <- rep(1, length(bmi))
# 	hr[bmi > 27] <- 5
# 	survival <- (1 - gompertz_makeham_cdf(age)) ^ hr
# 	return(survival)
# }

main()