source("~/repo/mr_frailty/R/functions.R")

# Required user defined functions

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
	eff <- log(1/0.5) / edu_sd
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


# Run simulations

parameters <- expand.grid(sim = 1:1000)


# Parallel
arguments <- commandArgs(T)
jid <- "all"
outdir <- "./"
if(length(arguments) > 0)
{
	jid <- as.numeric(arguments[1])
	chunksize <- as.numeric(arguments[2])
	outdir <- arguments[3]
	stopifnot(all(!is.na(jid), !is.na(chunksize), !is.na(outdir)))

	first <- (jid - 1) * chunksize + 1
	last <- min(nrow(parameters), jid * chunksize)
	parameters <- parameters[first:last, , drop=FALSE]
}

# Set output file
outfile <- paste(outdir, "/results_", jid, ".RData", sep="")
message("Running ", jid, ": ", nrow(parameters), " rows")
message("Saving in ", outfile)


res <- run_simulations(
	sim_start = parameters$sim[1], 
	sim_end = parameters$sim[nrow(parameters)],
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

save(res, file=outfile)
