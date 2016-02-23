# r script to read trees in *.tre, *.top files (after mb2badger and summarize, see README4cats)
# returns a df*.csv file with a summary of branch lengths
# Claudia February 2016

library(ape)
who="cats"
trees= read.tree("../datasets/4taxa-cats-clean.tre")
top=read.table("../datasets/4taxa-cats-clean.top")

who="birds"
trees= read.tree("../datasets/birds4.tre")
top=read.table("../datasets/birds4.top")

n=length(top$V1)
print(n)
length(trees) #11001
trees = trees[1002:length(trees)]

branch.lengths = matrix(0,length(trees),5)
i = 1
for(t in trees){
    branch.lengths[i,] = t$edge.length[2:6] #first bl is 0
    i = i+1
}
head(branch.lengths)
nrow(branch.lengths)
data <- data.frame(V1=top[1002:n,], branch.lengths)
head(data)

t = unique(data$V1)
data1=which(data$V1== t[1])
data2=which(data$V1== t[2])
data3=which(data$V1== t[3])

summary(data)
save(data,file=paste0("data_",who,"_mb.Rda"))

w=table(data$V1)
d=as.data.frame(w)
d

df=data.frame(tree=c(levels(data$V1)[data$V1[data1[1]]],levels(data$V1)[data$V1[data2[1]]],levels(data$V1)[data$V1[data3[1]]],"combined"))
ind1 = which(d$Var1 == levels(data$V1)[data$V1[data1[1]]])
ind2 = which(d$Var1 == levels(data$V1)[data$V1[data2[1]]])
ind3 = which(d$Var1 == levels(data$V1)[data$V1[data3[1]]])
n=length(data$V1)
df$probTop = round(c(d[ind1,2]/n, d[ind2,2]/n, d[ind3,2]/n, 1.0),4)


df$meanX1 = round(c(mean(data$X1[data1]),
    mean(data$X1[data2]),
    mean(data$X1[data3]),
    mean(data$X1)),4)
df$meanX2 = round(c(mean(data$X2[data1]),
    mean(data$X2[data2]),
    mean(data$X2[data3]),
    mean(data$X2)),4)
df$meanX3 = round(c(mean(data$X3[data1]),
    mean(data$X3[data2]),
    mean(data$X3[data3]),
    mean(data$X3)),4)
df$meanX4 = round(c(mean(data$X4[data1]),
    mean(data$X4[data2]),
    mean(data$X4[data3]),
    mean(data$X4)),4)
df$meanX5 = round(c(mean(data$X5[data1]),
    mean(data$X5[data2]),
    mean(data$X5[data3]),
    mean(data$X5)),4)

df$q025.X1 = round(c(quantile(data$X1[data1], probs=0.025),
    quantile(data$X1[data2],probs=0.025),
    quantile(data$X1[data3], probs=0.025),
    quantile(data$X1, probs=0.025)),4)
df$q025.X2 = round(c(quantile(data$X2[data1], probs=0.025),
    quantile(data$X2[data2],probs=0.025),
    quantile(data$X2[data3], probs=0.025),
    quantile(data$X2, probs=0.025)),4)
df$q025.X3 = round(c(quantile(data$X3[data1], probs=0.025),
    quantile(data$X3[data2],probs=0.025),
    quantile(data$X3[data3], probs=0.025),
    quantile(data$X3, probs=0.025)),4)
df$q025.X4 = round(c(quantile(data$X4[data1], probs=0.025),
    quantile(data$X4[data2],probs=0.025),
    quantile(data$X4[data3], probs=0.025),
    quantile(data$X4, probs=0.025)),4)
df$q025.X5 = round(c(quantile(data$X5[data1], probs=0.025),
    quantile(data$X5[data2],probs=0.025),
    quantile(data$X5[data3], probs=0.025),
    quantile(data$X5, probs=0.025)),4)

df$q975.X1 = round(c(quantile(data$X1[data1], probs=0.975),
    quantile(data$X1[data2],probs=0.975),
    quantile(data$X1[data3], probs=0.975),
    quantile(data$X1, probs=0.975)),4)
df$q975.X2 = round(c(quantile(data$X2[data1], probs=0.975),
    quantile(data$X2[data2],probs=0.975),
    quantile(data$X2[data3], probs=0.975),
    quantile(data$X2, probs=0.975)),4)
df$q975.X3 = round(c(quantile(data$X3[data1], probs=0.975),
    quantile(data$X3[data2],probs=0.975),
    quantile(data$X3[data3], probs=0.975),
    quantile(data$X3, probs=0.975)),4)
df$q975.X4 = round(c(quantile(data$X4[data1], probs=0.975),
    quantile(data$X4[data2],probs=0.975),
    quantile(data$X4[data3], probs=0.975),
    quantile(data$X4, probs=0.975)),4)
df$q975.X5 = round(c(quantile(data$X5[data1], probs=0.975),
    quantile(data$X5[data2],probs=0.975),
    quantile(data$X5[data3], probs=0.975),
    quantile(data$X5, probs=0.975)),4)

df

save(df,file=paste0("df_",who,"_mb.Rda"))
write.table(df,file=paste0("df_",who,"_mb.csv"), sep=",", row.names=FALSE)
