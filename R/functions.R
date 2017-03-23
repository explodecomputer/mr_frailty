#' gompertz makeham cdf
#'
#' <full description>
#'
#' @param t <what param does>
#' @param alpha=7.478359e-05 <what param does>
#' @param beta=8.604875e-02 <what param does>
#' @param lambda=-1.846973e-03 <what param does>
#'
#' @export
#' @return vector
gompertz_makeham_cdf <- function(t, alpha=7.478359e-05, beta=8.604875e-02, lambda=-1.846973e-03)
{
	pmin(pmax(1 - exp( -lambda * t - alpha / beta * (exp(beta*t) - 1) ), 0), 1)
}


#' simulate snps
#'
#' <full description>
#'
#' @param n <what param does>
#' @param snp_freq <what param does>
#'
#' @export
#' @return matrix
simulate_snps <- function(n, snp_freq)
{
	snps <- sapply(snp_freq, function(x) rbinom(n, 2, x))
	colnames(snps) <- paste0("snp", 1:ncol(snps))
	return(snps)
}


#' simulate exposure
#'
#' <full description>
#'
#' @param n <what param does>
#' @param snps <what param does>
#' @param snp_beta <what param does>
#' @param exposure_mean <what param does>
#' @param exposure_sd <what param does>
#'
#' @export
#' @return vector
simulate_exposure <- function(n, snps, snp_beta, exposure_mean, exposure_sd, lb=-Inf, ub=Inf)
{
	grs <- as.numeric(snps %*% snp_beta)
	v_grs <- var(grs)
	v_resid <- exposure_sd^2 - v_grs
	grs <- scale(grs) * sqrt(v_grs / (v_grs + v_resid))
	resid <- rnorm(n, 0, sqrt(v_resid / (v_grs + v_resid)))
	exposure <- as.numeric(scale(grs + resid) * exposure_sd + exposure_mean)

	dist <- sort(urnorm(n, exposure_mean, exposure_sd, lb, ub))
	ord <- order(exposure)
	exposure[ord] <- dist

	return(exposure)
}


#' simulate events
#'
#' <full description>
#'
#' @param age <what param does>
#' @param exposure <what param does>
#' @param survival_function <what param does>
#'
#' @export
#' @return vector
simulate_events <- function(age, exposure, survival_function)
{
	n <- length(age)
	survival <- survival_function(age, exposure)
	event <- rbinom(n, 1, survival)
	return(event)
}


#' plot quantiles
#'
#' <full description>
#'
#' @param exposure <what param does>
#' @param outcome <what param does>
#' @param exposure_breaks=10 <what param does>
#' @param mediator=NULL <what param does>
#' @param mediator_breaks=NULL <what param does>
#' @param mediator_name=NULL <what param does>
#' @param ylab="outcome" <what param does>
#' @param xlab="exposure" <what param does>
#'
#' @export
#' @return plot
plot_quantiles <- function(exposure, outcome, exposure_breaks=10, mediator=NULL, mediator_breaks=NULL, mediator_name=NULL, ylab="outcome", xlab="exposure")
{
	require(dplyr)
	require(ggplot2)
	d <- data.frame(exposure=exposure, outcome=outcome, cuts=cut(exposure, breaks=exposure_breaks))

	if(!is.null(mediator))
	{
		if(is.numeric(mediator_breaks))
		{
			d$mediator <- cut(mediator, mediator_breaks)
		} else {
			d$mediator <- mediator
		}

		res <- d %>% dplyr::group_by(cuts, mediator) %>%
			dplyr::summarise(prob = mean(outcome), prob_sd = sd(outcome)/n())


		p <- ggplot(res, aes(x=cuts, y=prob, colour=as.factor(mediator), group=as.factor(mediator))) +
			geom_errorbar(aes(ymin=pmax(prob - prob_sd * 1.96, 0), ymax=pmin(prob + prob_sd * 1.96, 1)), width=0) +	
			stat_summary(fun.y=identity, geom="line") +
			geom_point() +
			theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
			labs(x=xlab, y=ylab, colour=mediator_name) +
			scale_colour_brewer(type="qual")

	} else {

		res <- d %>% dplyr::group_by(cuts) %>%
			dplyr::summarise(prob = mean(outcome), prob_sd = sd(outcome)/n())

		p <- ggplot(res, aes(x=cuts, y=prob)) +
			geom_point() +
			geom_errorbar(aes(ymin=pmax(prob - prob_sd * 1.96, 0), ymax=pmin(prob + prob_sd * 1.96, 1)), width=0) +
			theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
			labs(x=xlab, y=ylab)

	}
	return(p)
}




#' pool mean sd
#'
#' <full description>
#'
#' @param n <what param does>
#' @param m <what param does>
#' @param s <what param does>
#'
#' @export
#' @return data frame
pool_mean_sd <- function(n, m, s)
{
	gn <- sum(n, na.rm=TRUE)
	gm <- sum(m * n, na.rm=TRUE) / gn
	gs <- sqrt(sum(n * (s^2 + (m - gm)^2), na.rm=TRUE) / gn)
	return(data.frame(gn=gn, gm=gm, gs=gs))
}


#' simulate outcome binary
#'
#' DEPRECATED
#'
#' @param age_summary <what param does>
#' @param max_age=100 <what param does>
#' @param min_age=18 <what param does>
#' @param sample_size_multiplier=3 <what param does>
#'
#' @export
#' @return data frame
simulate_outcome_binary <- function(age_summary, max_age=100, min_age=18, sample_size_multiplier=3)
{
	require(dplyr)
	dat <- group_by(age_summary, cc) %>%
		do({
			data.frame(age = rnorm(.$gn * sample_size_multiplier, .$gm, .$gs), cc = .$cc, stringsAsFactors=FALSE)
		}) %>%
		as.data.frame()
	dat$age <- pmin(max_age, dat$age) %>% pmax(min_age)
	return(dat)
}


#' simulate ages
#'
#' <full description>
#'
#' @param n <what param does>
#' @param age_mean <what param does>
#' @param age_sd <what param does>
#' @param max_age=100 <what param does>
#' @param min_age=18 <what param does>
#' @param sample_size_multiplier=3 <what param does>
#'
#' @export
#' @return data frame
simulate_ages <- function(n, age_mean, age_sd, max_age=100, min_age=18, sample_size_multiplier=3)
{
	require(Runuran)
	dat <- data.frame(age = urnorm(n * sample_size_multiplier, age_mean, age_sd, ub=max_age, lb=min_age), stringsAsFactors=FALSE)
	dat$age <- pmin(max_age, dat$age) %>% pmax(min_age)
	return(dat)
}


#' get age summary
#'
#' <full description>
#'
#' @param demog <what param does>
#' @param ncases <what param does>
#' @param ncontrol <what param does>
#' @param case_mean <what param does>
#' @param control_mean <what param does>
#' @param case_sd <what param does>
#' @param control_sd <what param does>
#'
#' @export
#' @return data frame
get_age_summary <- function(demog, ncases, ncontrol, case_mean, control_mean, case_sd, control_sd)
{
	ages <- rbind(
		pool_mean_sd(demog[[ncases]], demog[[case_mean]], demog[[case_sd]]),
		pool_mean_sd(demog[[ncontrol]], demog[[control_mean]], demog[[control_sd]])
	)
	ages <- rbind(ages, pool_mean_sd(ages$gn, ages$gm, ages$gs))
	ages$cc <- c(1,0,NA)
	return(ages)
}


#' sample cases controls
#'
#' @param dat <what param does>
#' @param age_summary <what param does>
#'
#' @export
#' @return array
sample_cases_controls <- function(dat, age_summary)
{
	dat$id <- 1:nrow(dat)
	temp1 <- subset(dat, alive==1 & cc==1)
	temp1 <- temp1[sample(1:nrow(temp1), age_summary$gn[age_summary$cc==1]), ]
	temp0 <- subset(dat, alive==1 & cc==0)
	temp0 <- temp0[sample(1:nrow(temp0), age_summary$gn[age_summary$cc==0]), ]
	return(sort(rbind(temp1, temp0)$id))
}



#' Get summary stats from dat and snps
#'
#' @param dat <what param does>
#' @param snps <what param does>
#' @param index <what param does>
#'
#' @export
#' @return data frame
get_summary_stats <- function(dat, snps, index)
{
	d <- dat[index, ]
	s <- snps[index, ]

	b <- array(0, ncol(s))
	se <- array(0, ncol(s))
	pval <- array(0, ncol(s))
	for(i in 1:ncol(s))
	{
		message(i)
		mod <- summary(glm(d$cc ~ s[,i], family="binomial"))
		b[i] <- coefficients(mod)[2,1]
		se[i] <- coefficients(mod)[2,2]
		pval[i] <- coefficients(mod)[2,4]
	}
	return(data.frame(b=b, se=se, pval=pval))
}


do_mr <- function(exposure_effects, outcome_effects, exposure_se, outcome_se)
{
	require(TwoSampleMR)

	dat <- data.frame(
		id.exposure = "exposure",
		id.outcome = "outcome",
		exposure = "exposure",
		outcome = "outcome",
		beta.exposure = exposure_effects,
		beta.outcome = outcome_effects,
		se.exposure = exposure_se,
		se.outcome = outcome_se,
		mr_keep = TRUE
	)

	return(mr(dat))
}



#' Simulation of the survival probability
#'
#' The simulation requires each individual has some probability of death.
#' This is a function of age and some predictor, related by the Gompertz-Makeham
#' cumulative density function.
#'
#' @param age Array of ages
#' @param predictor Predicted values of independent influence on death
#'
#' @export
#' @return Numeric array of survival probabilities
simulate_survival_probability <- function(age, predictor)
{
	survival <- (1 - gompertz_makeham_cdf(age)) ^ predictor
	return(survival)
}



#' Simulate dataset
#'
#' Simulates a dataset with SNPs, case/control status, exposure trait, age, alive/dead status, polygenic risk score.
#'
#' @param age_summary Data frame output from \code{get_age_summary}
#' @param snp_beta Array of SNP effect sizes
#' @param snp_se Array of SNP standard errors
#' @param snp_eaf Array of SNP effect allele frequencies
#' @param exposure_mean Mean of the exposure
#' @param exposure_sd Standard deviation of the exposure
#' @param exposure_lb Lower bound of the exposure
#' @param exposure_ub Upper bound of the exposure
#' @param outcome_prevalence_function Function that relates age to case/control prevalence. Expects age as an input and returns a prevalence for every age value
#' @param survival_function Function that relates the age and exposure to the survival probability. Takes a vector of ages and exposures (same lengths) and returns a vector of survival probabilities.
#' @param min_age=40 Minimum age to simulate
#' @param max_age=100 Maximum age to simulate
#' @param sample_size_multiplier=5 How many times larger than age summary should the dataset be
#'
#' @export
#' @return List of data frames - SNPs, Phenotype data, SNP effects
simulate_data <- function(age_summary, snp_beta, snp_se, snp_eaf, exposure_mean, exposure_sd, exposure_lb=-Inf, exposure_ub=Inf, outcome_prevalence_function, survival_function, min_age=40, max_age=100, sample_size_multiplier=5)
{
	dat <- simulate_ages(age_summary$gn[3], age_summary$gm[3], age_summary$gs[3], max_age=max_age, min_age=min_age, sample_size_multiplier=sample_size_multiplier)
	dat$cc <- simulate_events(dat$age, NULL, outcome_prevalence_function)
	snps <- simulate_snps(nrow(dat), snp_eaf)
	dat$exposure <- simulate_exposure(nrow(dat), snps, snp_beta * exposure_sd, exposure_mean, exposure_sd, lb=exposure_lb, ub=exposure_ub)

	dat$alive <- simulate_events(dat$age, dat$exposure, survival_function)
	dat$dead <- as.numeric(!dat$alive)
	dat$grs <- snps %*% snp_beta
	eff <- data.frame(b=snp_beta, se=snp_se, eaf=snp_eaf)
	return(list(dat=dat, snps=snps, eff=eff))
}


#' Perform MR and observational analysis on simulated data
#'
#' @param dat Output from \code{simulate_data}
#' @param age_summary Output from \code{get_age_summary}
#'
#' @export
#' @return Data frame of analysis results
analyse_data <- function(dat, age_summary)
{
	require(systemfit)
	snps <- dat$snps
	eff <- dat$eff
	dat <- dat$dat
	index <- sample_cases_controls(dat, age_summary)

	a <- summary(glm(cc ~ exposure, dat[index,], family="binomial"))
	b <- summary(glm(cc ~ grs, dat[index,], family="binomial"))
	c <- summary(systemfit(cc ~ exposure, "2SLS", inst = ~ grs, data = dat[index,]))

	l1 <- coefficients(a)[2,]
	l2 <- coefficients(b)[2,]
	l3 <- coefficients(c)[2,]
	l <- as.data.frame(rbind(l1, l2, l3), stringsAsFactors=FALSE)
	names(l) <- c("beta", "se", "tval", "pval")
	l$test <- c("Observational", "GRS", "2SLS")

	ss <- get_summary_stats(dat, snps, index)
	mres <- do_mr(eff$b, ss$b, eff$se, ss$se)
	l4 <- data.frame(beta = mres$b, se = mres$se, tval = NA, pval = mres$pval, test = mres$method)

	l <- rbind(l, l4)
	l <- subset(l, select=c(test, beta, se, pval))
	return(l)
}


#' Run simulations
#'
#' Simulates data and performs analysis a number of times to return results from multiple runs
#'
#' @param sim_start First simulation number
#' @param sim_end Last simulation number
#' @param age_summary Data frame output from \code{get_age_summary}
#' @param snp_beta Array of SNP effect sizes
#' @param snp_se Array of SNP standard errors
#' @param snp_eaf Array of SNP effect allele frequencies
#' @param exposure_mean Mean of the exposure
#' @param exposure_sd Standard deviation of the exposure
#' @param exposure_lb Lower bound of the exposure
#' @param exposure_ub Upper bound of the exposure
#' @param outcome_prevalence_function Function that relates age to case/control prevalence. Expects age as an input and returns a prevalence for every age value
#' @param survival_function Function that relates the age and exposure to the survival probability. Takes a vector of ages and exposures (same lengths) and returns a vector of survival probabilities.
#' @param min_age=40 Minimum age to simulate
#' @param max_age=100 Maximum age to simulate
#' @param sample_size_multiplier=5 How many times larger than age summary should the dataset be
#'
#' @export
#' @return Data frame of analysis results from multiple simulations
run_simulations <- function(sim_start, sim_end, age_summary, snp_beta, snp_se, snp_eaf, exposure_mean, exposure_sd, exposure_lb, exposure_ub, outcome_prevalence_function, survival_function, min_age=40, max_age=100, sample_size_multiplier=5)
{
	require(dplyr)
	res <- list()
	sim <- sim_start:sim_end
	nsim <- length(sim)
	for (i in 1:nsim)
	{
		message(sim[i], " of ", sim_end)

		dat <- simulate_data(age_summary, snp_beta, snp_se, snp_eaf, exposure_mean, exposure_sd, exposure_lb, exposure_ub, outcome_prevalence_function, survival_function, min_age, max_age, sample_size_multiplier)

		res[[i]] <- analyse_data(dat, age_summary)
		res[[i]]$sim <- sim[i]

	}
	res <- bind_rows(res)
	return(res)
}

