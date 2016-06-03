## R script to study the weights in logw.txt
## after ./bl ....
## Claudia May 2016

data = read.table("logw.txt", header=TRUE, sep=",")
head(data)

mylogw = data$logweight - mean(data$logweight)
data$w = exp(mylogw)/sum(exp(mylogw))
data[data$w>0.01,]
length(data[data$w>0.01,]$w)
hist(data$w)
plot(1:length(data$w),cumsum(rev(sort(data$w))))
(1/sum(data$w^2))/1000
summary(data$bl2)
