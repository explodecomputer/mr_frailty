library(dplyr)

arguments <- commandArgs(T)
splits <- 10
outdir <- "~/repo/mr_frailty/inst/analysis/edu_az/scratch2"
savefile <- "~/repo/mr_frailty/inst/analysis/edu_az/results/model2.RData"

l <- list()
for(i in 1:splits)
{
	cat(i, "\n")
	filename <- paste0(outdir, "/results_", i, ".RData")
	if(file.exists(filename))
	{
		load(filename)
		l[[i]] <- res
	} else {
		message("Missing ", filename)
	}
}

res <- bind_rows(l)
save(res, file = savefile)
