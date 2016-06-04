## R script to study the weights in logw.txt
## after ./bl ....
## Claudia May 2016

data = read.table("logw.txt", header=TRUE, sep=",")
head(data)
nc = ncol(data)

mylogw = data[,nc] - mean(data[,nc])
data$w = exp(mylogw)/sum(exp(mylogw))
data[data$w>0.01,]$w
length(data[data$w>0.01,]$w)
hist(data$w)
plot(1:length(data$w),cumsum(rev(sort(data$w))))
(1/sum(data$w^2))/1000
(1/sum(data$w^2))/907
##summary(data$bl2)
