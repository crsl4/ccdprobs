# r script to macth branches in mb trees and branches in IS
# Claudia February 2016

library(ape)
library(phangorn)
who="birds"
load("data_birds_mb.Rda")
data_mb=data
load("data_birds.Rda")

who="cats"
load("data_cats_mb.Rda")
data_mb=data
load("data_cats.Rda")

head(data)
head(data_mb)

trees = levels(data$trees)[1:3]
trees
trees_mb= levels(data_mb$V1)[1:3]
trees_mb

# trees in IS
t1=read.tree(text=trees[1])
t2=read.tree(text=trees[2])
t3=read.tree(text=trees[3])

# trees in mb
t1_mb=read.tree(text=trees_mb[1])
t2_mb=read.tree(text=trees_mb[2])
t3_mb=read.tree(text=trees_mb[3])
t1_mb<- unroot(t1_mb)
t2_mb<- unroot(t2_mb)
t3_mb<- unroot(t3_mb)
trees_mb = c(t1_mb,t2_mb, t3_mb)
# need to match order of trees in mb to IS
ind1 = which(c(RF.dist(t1,t1_mb), RF.dist(t1,t2_mb), RF.dist(t1,t3_mb)) == 0)
ind2 = which(c(RF.dist(t2,t1_mb), RF.dist(t2,t2_mb), RF.dist(t2,t3_mb)) == 0)
ind3 = which(c(RF.dist(t3,t1_mb), RF.dist(t3,t2_mb), RF.dist(t3,t3_mb)) == 0)
print(ind1)
print(ind2)
print(ind3)

data=data.frame(tree=c(rep(trees[1],5), rep(trees[2],5), rep(trees[3],5)), is=rep(0,15), mb=rep(0,15))

# first tree
e <- t1$edge
o <- order(t1$tip.label)
emb <- trees_mb[[ind1]]$edge
omb <- order(trees_mb[[ind1]]$tip.label)
# internal branch
data[1,2] = which(rowSums(e) == 11)
data[1,3] = which(rowSums(emb) == 11)
# branch leading to 1
data[2,2] = which(e[,2] == o[1])
data[2,3] = which(emb[,2] == omb[1])
# branch leading to 2
data[3,2] = which(e[,2] == o[2])
data[3,3] = which(emb[,2] == omb[2])
# branch leading to 3
data[4,2] = which(e[,2] == o[3])
data[4,3] = which(emb[,2] == omb[3])
# branch leading to 4
data[5,2] = which(e[,2] == o[4])
data[5,3] = which(emb[,2] == omb[4])

# second tree
e <- t2$edge
o <- order(t2$tip.label)
emb <- trees_mb[[ind2]]$edge
omb <- order(trees_mb[[ind2]]$tip.label)
# internal branch
data[6,2] = which(rowSums(e) == 11)
data[6,3] = which(rowSums(emb) == 11)
# branch leading to 1
data[7,2] = which(e[,2] == o[1])
data[7,3] = which(emb[,2] == omb[1])
# branch leading to 2
data[8,2] = which(e[,2] == o[2])
data[8,3] = which(emb[,2] == omb[2])
# branch leading to 3
data[9,2] = which(e[,2] == o[3])
data[9,3] = which(emb[,2] == omb[3])
# branch leading to 4
data[10,2] = which(e[,2] == o[4])
data[10,3] = which(emb[,2] == omb[4])

# third tree
e <- t3$edge
o <- order(t3$tip.label)
emb <- trees_mb[[ind3]]$edge
omb <- order(trees_mb[[ind3]]$tip.label)
# internal branch
data[11,2] = which(rowSums(e) == 11)
data[11,3] = which(rowSums(emb) == 11)
# branch leading to 1
data[12,2] = which(e[,2] == o[1])
data[12,3] = which(emb[,2] == omb[1])
# branch leading to 2
data[13,2] = which(e[,2] == o[2])
data[13,3] = which(emb[,2] == omb[2])
# branch leading to 3
data[14,2] = which(e[,2] == o[3])
data[14,3] = which(emb[,2] == omb[3])
# branch leading to 4
data[15,2] = which(e[,2] == o[4])
data[15,3] = which(emb[,2] == omb[4])

data
save(data,file=paste0("bl_",who,".Rda"))
write.table(data,file=paste0("bl_",who,".csv"), sep=",", row.names=FALSE)
