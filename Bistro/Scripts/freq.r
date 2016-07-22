## R code to compare the frequency of trees

library(ape)
dat <- read.table("run15.txt", header=TRUE)
str(dat)
logweight = dat$logWt
logw = logweight - max(logweight)
w = exp(logw) / sum(exp(logw))
dat$w <- w

means=with(dat, sapply( split( w,tree ),sum ) )
dat$mean = means[8]
head(dat)
dat$indicator = dat$tree == levels(dat$tree)[8]
v = sum((dat$w*dat$indicator - dat$mean)^2)/100

1/sum(w^2)
