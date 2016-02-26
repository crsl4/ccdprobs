# r script to sample tree from ccdprobs, and sample branch lengths with TN,
# and calculate importance weight
# Claudia February 2016

library(ape)
# library(bigvis)
source('branch-length.r')
source('4taxa_cats_dogs_functions.r')

who = "cats"
d=read.dna("../datasets/4taxa-cats.phy") #needs to be 4 taxa
dat.tre=read.table("../datasets/4taxa-cats_ccdprobs.out", header=FALSE)

who = "birds"
d=read.dna("../datasets/birds4-clean.phy") #needs to be 4 taxa
dat.tre=read.table("../datasets/birds4-clean_ccdprobs.out", header=FALSE)


print(dat.tre)
t=sampleTopQuartet(dat.tre)
print(t)
b=sampleBLQuartet(d,t$tre)
print(b)
t$tre$edge.length <- b$bl

logw = b$logprior+b$loglik-log(t$prob)-b$logdensity
print(logw)

nreps = 1000
trees = rep(NA,nreps)
logwv = rep(0,nreps)
branch.lengths = matrix(0,nreps,5)
for(i in 1:nreps){
    t=sampleTopQuartet(dat.tre)
    b=sampleBLQuartet(d,t$tre)
    branch.lengths[i,]=b$bl
    logw = b$logprior+b$loglik-log(t$prob)-b$logdensity
    #print(logw)
    trees[i] = write.tree(t$tre)
    logwv[i] = logw
}
data = data.frame(trees,logwv,branch.lengths)
head(data)
summary(data)
my.logw = logwv - mean(logwv)
data$w = exp(my.logw)/sum(exp(my.logw))
save(data,file=paste0("data_",who,".Rda"))
#load("data_birds.Rda")
data[data$w>0.01,]

tree1 = "((1,2),3,4);"
data1=which(data$trees== tree1)
tree2 = "((1,3),2,4);"
data2=which(data$trees== tree2)
tree3 = "(1,(2,3),4);"
data3=which(data$trees== tree3)

df=data.frame(tree=c(tree1,tree2,tree3,"combined"))
df$probTop = round(c(sum(data$w[data1]), sum(data$w[data2]), sum(data$w[data3]), 1.0),6)
# mean bl
df$meanX1 = round(c(weighted.mean(data$X1[data1],data$w[data1]),
    weighted.mean(data$X1[data2],data$w[data2]),
    weighted.mean(data$X1[data3],data$w[data3]),
    weighted.mean(data$X1,data$w)),4)
df$meanX2 = round(c(weighted.mean(data$X2[data1],data$w[data1]),
    weighted.mean(data$X2[data2],data$w[data2]),
    weighted.mean(data$X2[data3],data$w[data3]),
    weighted.mean(data$X2,data$w)),4)
df$meanX3 = round(c(weighted.mean(data$X3[data1],data$w[data1]),
    weighted.mean(data$X3[data2],data$w[data2]),
    weighted.mean(data$X3[data3],data$w[data3]),
    weighted.mean(data$X3,data$w)),4)
df$meanX4 = round(c(weighted.mean(data$X4[data1],data$w[data1]),
    weighted.mean(data$X4[data2],data$w[data2]),
    weighted.mean(data$X4[data3],data$w[data3]),
    weighted.mean(data$X4,data$w)),4)
df$meanX5 = round(c(weighted.mean(data$X5[data1],data$w[data1]),
    weighted.mean(data$X5[data2],data$w[data2]),
    weighted.mean(data$X5[data3],data$w[data3]),
    weighted.mean(data$X5,data$w)),4)

# 2.5% quantile bl
df$q025.X1 = round(c(weighted.quantile(data$X1[data1],data$w[data1], probs=0.025),
    weighted.quantile(data$X1[data2],data$w[data2], probs=0.025),
    weighted.quantile(data$X1[data3],data$w[data3], probs=0.025),
    weighted.quantile(data$X1,data$w, probs=0.025)),4)
df$q025.X2 = round(c(weighted.quantile(data$X2[data1],data$w[data1], probs=0.025),
    weighted.quantile(data$X2[data2],data$w[data2], probs=0.025),
    weighted.quantile(data$X2[data3],data$w[data3], probs=0.025),
    weighted.quantile(data$X2,data$w, probs=0.025)),4)
df$q025.X3 = round(c(weighted.quantile(data$X3[data1],data$w[data1], probs=0.025),
    weighted.quantile(data$X3[data2],data$w[data2], probs=0.025),
    weighted.quantile(data$X3[data3],data$w[data3], probs=0.025),
    weighted.quantile(data$X3,data$w, probs=0.025)),4)
df$q025.X4 = round(c(weighted.quantile(data$X4[data1],data$w[data1], probs=0.025),
    weighted.quantile(data$X4[data2],data$w[data2], probs=0.025),
    weighted.quantile(data$X4[data3],data$w[data3], probs=0.025),
    weighted.quantile(data$X4,data$w, probs=0.025)),4)
df$q025.X5 = round(c(weighted.quantile(data$X5[data1],data$w[data1], probs=0.025),
    weighted.quantile(data$X5[data2],data$w[data2], probs=0.025),
    weighted.quantile(data$X5[data3],data$w[data3], probs=0.025),
    weighted.quantile(data$X5,data$w, probs=0.025)),4)

# 97.5% quantile bl
df$q975.X1 = round(c(weighted.quantile(data$X1[data1],data$w[data1], probs=0.975),
    weighted.quantile(data$X1[data2],data$w[data2], probs=0.975),
    weighted.quantile(data$X1[data3],data$w[data3], probs=0.975),
    weighted.quantile(data$X1,data$w, probs=0.975)),4)
df$q975.X2 = round(c(weighted.quantile(data$X2[data1],data$w[data1], probs=0.975),
    weighted.quantile(data$X2[data2],data$w[data2], probs=0.975),
    weighted.quantile(data$X2[data3],data$w[data3], probs=0.975),
    weighted.quantile(data$X2,data$w, probs=0.975)),4)
df$q975.X3 = round(c(weighted.quantile(data$X3[data1],data$w[data1], probs=0.975),
    weighted.quantile(data$X3[data2],data$w[data2], probs=0.975),
    weighted.quantile(data$X3[data3],data$w[data3], probs=0.975),
    weighted.quantile(data$X3,data$w, probs=0.975)),4)
df$q975.X4 = round(c(weighted.quantile(data$X4[data1],data$w[data1], probs=0.975),
    weighted.quantile(data$X4[data2],data$w[data2], probs=0.975),
    weighted.quantile(data$X4[data3],data$w[data3], probs=0.975),
    weighted.quantile(data$X4,data$w, probs=0.975)),4)
df$q975.X5 = round(c(weighted.quantile(data$X5[data1],data$w[data1], probs=0.975),
    weighted.quantile(data$X5[data2],data$w[data2], probs=0.975),
    weighted.quantile(data$X5[data3],data$w[data3], probs=0.975),
    weighted.quantile(data$X5,data$w, probs=0.975)),4)

save(df,file=paste0("df_",who,".Rda"))
write.table(df,file=paste0("df_",who,".txt"), sep=",", row.names=FALSE)

## #things to check:
## # weight well-distributed:
#round(data$w,4)
#data[data$w>0.01,] # we don't want all weight in few trees

## #want to compare sample branch lengths, to true branch lengths (posterior dist)
## library(ggplot2)
## tre=read.tree(text=tree1)
## tre$edge #which branch length corresponds each branch
## ggplot(data[data1,],aes(x=X3)) + geom_density()
## trueBL = c(0.087513,0.079521,0.048142,0.062346,0.020288);
## ## (Felis_catus___domestic_cat:0.087513,Neofelis_nebulosa___clouded_leopard:0.079521,(Panthera_pardus___leopard:0.048142,Panthera_tigris___tiger:0.062346):0.020288);

## ## want to compare not only the true mean, but 95% CI for each branch length
