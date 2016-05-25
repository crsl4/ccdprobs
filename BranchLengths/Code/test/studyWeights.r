## R script to study the weights in logw.txt
## after ./bl ....
## Claudia May 2016

dat = read.table("logw.txt")
head(dat)
summary(dat$V1)
mean(dat$V1)

logw = dat$V1 - mean(dat$V1)
summary(logw)
