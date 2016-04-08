library(dplyr)

arguments <- commandArgs(T)
splits <- 100
outdir <- "~/repo/mr_frailty/tests/bmi_pd/scratch"
savefile <- "~/repo/mr_frailty/tests/bmi_pd/results/model1"

l <- list()
for(i in 1:splits)
{
	cat(i, "\n")
	filename <- paste0(outdir, "/results", i, ".RData")
	if(file.exists(filename))
	{
		load(filename)
		l[[i]] <- res
	} else {
		message("Missing ", filename)
	}
}

res <- rbind_all(l)
save(res, file = savefile)
