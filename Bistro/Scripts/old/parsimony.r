## R script to compute the parsominy score of each tree in a sample
## to try to see if we can use the parsimony score to weight sample of NJ trees
## Claudia July 2016

## in Examples/Whales
library(ape)
library(phangorn)
NJtrees = read.tree("whales.in.trees")
dna = read.FASTA("../../Data/whales.fasta")
pdat = phyDat(dna)
mbtree = read.tree(text="((((((1,2),(3,(((4,5),14),((6,7),((((8,9),11),13),(10,12)))))),15),((16,17),((18,19),(20,21)))),(22,23)),(((24,(27,28)),(25,26)),(29,30)),31);")
mbtree$tip.label <- names(pdat)[as.numeric(mbtree$tip.label)]
mbscore = parsimony(mbtree,pdat)

data = read.table("whales.in.trees", sep=" ")
names(data) <- c("tree", "count")
data$score = rep(0,length(mbtrees))


data = data.frame( tree = rep(NA,length(NJtrees)), score = rep(0,length(NJtrees)))

for(i in 1:length(NJtrees)){
    tree = NJtrees[[i]]
    data$tree[i] = write.tree(tree)
    tree$tip.label <- names(pdat)[as.numeric(tree$tip.label)]
    data$score[i] = parsimony(tree,pdat)
}

hist(data$score, main="Parsimony score NJ trees (red=MrBayes consensus tree)")
abline(v=mbscore,col="red")


## mrbayes trees ----------------------------------------------------
library(ape)
library(phangorn)
mbtrees = read.tree("whales-mb.in.trees")
dna = read.FASTA("../../Data/whales.fasta")
pdat = phyDat(dna)

datamb = read.table("whales-mb.in.trees", sep=" ")
names(datamb) <- c("tree", "count")
datamb$score = rep(0,length(mbtrees))

for(i in 1:length(mbtrees)){
    tree = mbtrees[[i]]
    tree$tip.label <- names(pdat)[as.numeric(tree$tip.label)]
    datamb$score[i] = parsimony(tree,pdat)
}

head(datamb)
plot(datamb$score,datamb$count)
abline(v=mbscore,col="red")

mi = min(datamb$score)
k = 0.5
datamb$weight = exp(k*(mi - datamb$score))
datamb$weight2 = datamb$weight * datamb$count
s = sum(datamb$weight2)
datamb$weight2 = datamb$weight2 / s

plot(datamb$count, datamb$weight2)



plot(datamb$count)
plot(datamb$score)

bestpar = sum(datamb[datamb$score<=mbscore,]$count)
worsepar = sum(datamb[datamb$score>mbscore,]$count)
bestpar / sum(datamb$count)
worsepar / sum(datamb$count)

cor(datamb$count,datamb$score) ## -0.212
cor(datamb$count,exp(datamb$score-max(datamb$score))) ## -0.0355
cor(datamb$count,exp(min(datamb$score)-datamb$score)) ## -0.0104


m = max(datamb$score)
m2 = min(datamb$score)
datamb$weight = exp(datamb$score-m) ## does not look good
datamb$weight2 = exp(m2-datamb$score) ## does not look good
datamb$weight = datamb$weight / min(datamb$weight)
min(datamb$weight2)
max(datamb$weight2)


## need to: decide weights function, write whales.in and run ccdprobs
## ccdprobs needs integers?
## weight
data$weight = (1/data$score)
plot(data$weight)
data$weight = data$weight / min(data$weight)
max(data$weight)
sum(data$weight)

m = max(data$score)
m2 = min(data$score)
data$weight = exp(data$score-m) ## does not look good
data$weight2 = exp(m2-data$score) ## does not look good
data$weight = data$weight / min(data$weight)
min(data$weight)
max(data$weight)


