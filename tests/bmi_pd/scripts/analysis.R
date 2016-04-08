library(ggplot2)
library(dplyr)


load("~/repo/mr_frailty/tests/bmi_pd/results/model1.RData")
res1 <- res
res1$model <- "model1"
load("~/repo/mr_frailty/tests/bmi_pd/results/model2.RData")
res2 <- res
res2$model <- "model2"
res <- rbind(res1, res2)



res_plot <- subset(res, test %in% c("2sls", "obs", "MR Egger", "Inverse variance weighted"))

res_plot$test <- as.factor(res_plot$test)
levels(res_plot$test) <- c("2SLS", "IVW", "MR Egger", "Obs assoc")

dat <- group_by(res_plot, test, model) %>% summarise(b=mean(beta), se=sd(beta)/n(), sd=sd(beta))


ggplot(subset(res_plot, test != "grs"), aes(x=beta)) +
geom_density(aes(fill=model), alpha=0.2) +
geom_vline(data=subset(dat, test != "grs"), aes(xintercept=b, colour=model)) +
geom_vline(xintercept=0, linetype="dashed") +
facet_grid(test ~ ., scale="free_y") +
scale_fill_brewer(type="qual")
ggsave("~/repo/mr_frailty/tests/bmi_pd/images/method_comparison.pdf")

