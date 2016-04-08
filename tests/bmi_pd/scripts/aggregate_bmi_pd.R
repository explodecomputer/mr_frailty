library(dplyr)

arguments <- commandArgs(T)
splits <- as.numeric(arguments[1])
outdir <- arguments[2]
savefile <- arguments[3]

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

parameters <- rbind_all(l)
save(parameters, file = savefile)
