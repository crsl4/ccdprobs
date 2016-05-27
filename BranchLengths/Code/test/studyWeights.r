## R script to study the weights in logw.txt
## after ./bl ....
## Claudia May 2016

seed = 646
dat = read.table(paste0("logw",seed,"_norm3.txt"))
head(dat)
summary(dat$V1)
mean(dat$V1)

logw = dat$V1 - mean(dat$V1)
summary(logw)
hist(logw)
w = exp(logw) / sum(exp(logw))
summary(w)
length(w[w>0.01])
w[w>0.01]
hist(w)
plot(1:length(w),cumsum(rev(sort(w))))
(1/sum(w^2))/1000
