## R script to study the weights in logw.txt
## after ./bl ....
## Claudia May 2016

seed = 646
dat = read.table(paste0("logw",seed,".txt"))
head(dat)
summary(dat$V1)
mean(dat$V1)

logw = dat$V1 - mean(dat$V1)
summary(logw)
hist(logw)
