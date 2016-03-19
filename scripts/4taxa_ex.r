## r script to sample tree from ccdprobs, and sample branch lengths with TN,
## and calculate importance weight
## Claudia February 2016
## for cats, birds and sim data

library(ape)
source('branch-length_lik.r')
source('4taxa_functions.r')

who = "cats"
d=read.dna("../datasets/4taxa-cats.phy") #needs to be 4 taxa
dat.tre=read.table("../datasets/4taxa-cats_ccdprobs.out", header=FALSE)
r=c(2.815,51.982,1.903,1.275,65.402,1.000)
p=c(0.2590,0.2379,0.1900,0.3131);
den = r[6]
r = r/den
Q = makeQ(r,p,4, rescale=TRUE)

who = "birds"
d=read.dna("../datasets/birds4-clean.phy") #needs to be 4 taxa
dat.tre=read.table("../datasets/birds4-clean_ccdprobs.out", header=FALSE)
r = c(0.2463,0.1764,0.1231,0.0187,0.4185,0.0170)
p = c(0.2776,0.2937,0.1612,0.2675)
den = r[6]
r = r/den
Q = makeQ(r,p,4, rescale=TRUE)


who = "simulated" #tree (1,2),3,4
nsites=1500
branch.length = c(0.03,0.11,0.078,0.091,0.098) #bl as birds mb: dxy,d1x,d2x,d3y,d4y
Q = randomQ(4,rescale=TRUE)
r=Q$r
p=Q$p
d = simulateData(Q,branch.length, nsites,filename="simSeq1_Qrand.txt")
dat.tre=read.table("../datasets/birds4-clean_ccdprobs.out", header=FALSE)
dat.tre$V2 = c(0.0,1.0,0.0)

print(dat.tre)
t=sampleTopQuartet(dat.tre)
print(t)
b=sampleBLQuartet(d,t$tre, trueT0=TRUE, verbose=TRUE)
print(b)
t$tre$edge.length <- b$bl

logw = b$logprior+b$loglik-log(t$prob)-b$logdensity
print(logw)

nreps = 1000
trees = rep(NA,nreps)
logwv = rep(0,nreps)
branch.lengths = matrix(0,nreps,5)
err = 0
for(i in 1:nreps){
    t=sampleTopQuartet(dat.tre)
    b=try(sampleBLQuartet(d,t$tre, trueT0=FALSE, Q=Q))
    if(class(b) == "try-error"){
        err = err + 1
    } else {
        branch.lengths[i,]=b$bl
        logw = b$logprior+b$loglik-b$logdensity
        ##print(logw)
        trees[i] = write.tree(t$tre)
        logwv[i] = logw
    }
}

data = data.frame(trees,logwv,branch.lengths)
lines=which(!is.na(data$trees))
data=data[lines,]
head(data)
summary(data)
my.logw = logwv[lines] - mean(logwv[lines])
data$w = exp(my.logw)/sum(exp(my.logw))
save(data,file=paste0("data_",who,".Rda"))

load("data_birds1.Rda")
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

