# r script to sample tree from ccdprobs, and sample branch lengths with TN,
# and calculate importance weight
# Claudia February 2016

library(ape)
source('branch-length.r')
source('4taxa_cats_dogs_functions.r')


d=read.dna("../datasets/4taxa-cats.phy") #needs to be 4 taxa
dat.tre=read.table("../datasets/4taxa-cats_ccdprobs.out", header=FALSE)
print(dat.tre)
t=sampleTopQuartet(dat.tre)
print(t)
b=sampleBLQuartet(d,t$tre)
print(b)
t$tre$edge.length <- b$bl

prior = 1/3 #fixit: indep exponential mean=0.1
logw = log(prior)+b$loglik-log(t$prob)-b$logdensity #posterior = prior*lik, not normalized, logdensity not normalized
print(logw)


reps = 100
trees = rep(NA,nreps)
logwv = rep(0,nreps)
branch.lengths = matrix(0,nreps,5)
for(i in 1:nreps){
    t=sampleTopQuartet(dat.tre)
    b=sampleBLQuartet(d,t$tre)
    branch.lengths[i,]=b$bl
    logw = log(prior)+b$loglik-log(t$prob)-b$logdensity
    print(logw)
    trees[i] = write.tree(t$tre)
    logwv[i] = logw
}
data = data.frame(trees,logwv,branch.lengths)
head(data)
summary(data)
my.logw = logwv - mean(logwv)
data$w = exp(my.logw)/sum(exp(my.logw))
save(data,file="data.Rda")
#load("data.Rda")

tree1 = "((1,2),3,4);"
data1=which(data$trees== tree1)
tree2 = "((1,3),2,4);"
data2=which(data$trees== tree2)
tree3 = "(1,(2,3),4);"
data3=which(data$trees== tree3)

#things to check:
# weight well-distributed:
round(data$w,4)
data[data$w>0.01,] # we don't want all weight in few trees


#want to compare sample branch lengths, to true branch lengths (posterior dist)
apply(branch.lengths[data1,],2,mean) # we actually want weighted mean by w
tre=read.tree(text=tree1)
tre$edge #which branch length corresponds each branch
ggplot(data[data1,],aes(x=X1)) + geom_density()
trueBL = c(0.087513,0.079521,0.048142,0.062346,0.020288);
## (Felis_catus___domestic_cat:0.087513,Neofelis_nebulosa___clouded_leopard:0.079521,(Panthera_pardus___leopard:0.048142,Panthera_tigris___tiger:0.062346):0.020288);

## want to compare not only the true mean, but 95% CI for each branch length
