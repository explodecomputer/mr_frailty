library(dplyr)


splits <- 100
outdir <- "../scratch"
savefile <- "../results/model2.RData"

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
