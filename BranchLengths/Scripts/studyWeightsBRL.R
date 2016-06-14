## R script to study the weights in logw.txt
## after ./bl ....
## Claudia May 2016
## Modified by Bret on June 1, 2016

dat = read.csv("logw.txt",header=TRUE)
#dat = read.csv("logw.txt",header=FALSE)
logweight = dat[,ncol(dat)]
logweight <- logweight[!is.nan(logweight)]
logweight <- logweight[is.finite(logweight)]
#logweight = logweight[!is.nan(logweight)]
logw = logweight - max(logweight)
w = exp(logw) / sum(exp(logw))
print(summary(w))
cat("# weights larger than 1%\n")
print(length(w[w>0.01]))
print(w[w>0.01])
plot(1:length(w),cumsum(rev(sort(w))),type="l")
cat("ESS proportion:\n")
print( (1/(length(w)*sum(w^2)) ) )
cat("ESS:\n")
print( (1/sum(w^2) ) )
