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



source("../../../../R/functions.R")
library(dplyr)
library(systemfit)
library(TwoSampleMR)
library(survival)

main <- function()
{
	demog <- read.table("../../../../data-raw/az_demographics.txt", sep="\t", header=TRUE)
	age_summary <- get_age_summary(demog, "n_case", "n_control", "mean_age_case", "mean_age_control", "sd_age_case", "sd_age_control")
	sbp_snps <- read_exposure_data("../../../../data-raw/sbp_inst.txt", sep="\t")
	sbp_mean <- 130
	sbp_sd <- 5.6


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

		dat <- simulate_ages(age_summary$gn[3], age_summary$gm[3], age_summary$gs[3], max_age=110, min_age=60, sample_size_multiplier=5)
		dat$cc <- simulate_events_old(dat$age, NULL, ad_incidence)
		snps <- simulate_snps(nrow(dat), sbp_snps$eaf.exposure)
		dat$sbp <- simulate_exposure(nrow(dat), snps, sbp_snps$beta.exposure, sbp_mean, sbp_sd)
		dat$alive <- simulate_events_old(dat$age, dat$sbp, sbp_survival)
		dat$dead <- as.numeric(!dat$alive)
		dat$grs <- snps %*% sbp_snps$beta.exposure
		coxph(Surv(age, dead) ~ sbp, dat)

		table(dat$cc, dat$alive)
		index <- sample_cases_controls(dat, age_summary)

		a <- summary(glm(cc ~ bmi, dat[index,], family="binomial"))
		b <- summary(glm(cc ~ grs, dat[index,], family="binomial"))
		c <- summary(systemfit(cc ~ bmi, "2SLS", inst = ~ grs, data = dat[index,]))

		l1[[i]] <- coefficients(a)[2,]
		l2[[i]] <- coefficients(b)[2,]
		l3[[i]] <- coefficients(c)[2,]

		ss <- get_summary_stats(dat, snps, index)
		mres <- do_mr(bmi_snps$beta.exposure, ss$b, bmi_snps$se.exposure, ss$se)
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


ad_incidence <- function(age, ...)
{
	ad_inc <- data.frame(
		age = c(57, 62, 67, 72, 77, 82, 87, 92),
		p = c(0.001, 0.005, 0.01, 0.025, 0.035, 0.11, 0.235, 0.30) * 2
	)
	d <- approx(x = ad_inc$age, y = ad_inc$p, xout=age, yleft=min(ad_inc$p), yright=max(ad_inc$p))
	return(d$y)
}


sbp_survival <- function(age, sbp)
{
	# hr of 1.06 per sd increase in sbp
	survival <- (1 - gompertz_makeham_cdf(age)) ^ (sbp * log(1.28) * sd(sbp, na.rm=TRUE))
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



library(survival)
library(KMsurv)
