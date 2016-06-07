## R script to study the weights in logw.txt
## after ./bl ....
## Claudia May 2016
## Modified by Bret on June 1, 2016

dat = read.csv("logw.txt",header=TRUE)
#dat = read.csv("logw.txt",header=FALSE)
logweight = dat[,ncol(dat)]
#logweight = logweight[!is.nan(logweight)]
logw = logweight - max(logweight)
w = exp(logw) / sum(exp(logw))
print(summary(w))
print(length(w[w>0.01]))
print(w[w>0.01])
plot(1:length(w),cumsum(rev(sort(w))),type="l")
print( (1/(length(w)*sum(w^2)) ) )
