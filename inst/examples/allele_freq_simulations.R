library(ggplot2)

d1 <- data.frame(age = 20:100, af = 0.5, cat="never smokers") 
d2 <- data.frame(age = 20:100, af = 0.5, cat="ever smokers") 
d2$af <- d2$af - ((d2$age - 20) / 150)^2

d <- rbind(d1, d2) 
d$n <- (1/(d$age))^2 * 1000000

ggplot(d, aes(x=age, y=af)) + 
geom_point(aes(colour=cat)) + 
geom_line(aes(colour=cat)) + 
labs(x="Age group / years", y="Allele frequency of smoking heaviness variant",colour="") + 
ylim(c(0, 0.7))
ggsave("images/smoking_af.pdf", width=10)

library(dplyr)

dat <- group_by(d, age, cat) %>%
	do({
		x <- as.data.frame(.)
		aa <- x$af[1]^2 * x$n[1]
		Aa <- 2 * x$af[1] * (1 - x$af[1]) * x$n[1]
		AA <- (1 - x$af[1])^2 * x$n[1]
		return(data.frame(g = c("aa", "Aa", "AA"), ng = c(aa, Aa, AA), af=x$af[1]))
	})


ggplot(subset(dat, age > 60), aes(x=age, y=ng)) +
geom_bar(position="stack", aes(fill=g), stat="identity") +
facet_grid(cat ~ .) +
geom_text(data = subset(d, age > 60), aes(label=round(af,2), y=250), angle=90, hjust=0) +
scale_fill_brewer(type="qual") +
labs(fill="Genotype", x="Age group / years", y="Individuals still alive")
ggsave("images/smoking_af_n.pdf", width=10)


ggplot(subset(dat, age > 60 & cat == "ever smokers"), aes(x=age, y=ng)) +
geom_bar(position="stack", aes(fill=g), stat="identity") +
geom_text(data = subset(d, age > 60 & cat == "ever smokers"), aes(label=round(af,2), y=250), angle=90, hjust=0) +
scale_fill_brewer(type="qual") +
labs(fill="Genotype", x="Age group / years", y="Individuals still alive")
ggsave("images/smoking_af_n_ever.pdf", width=10)
