## r script to check branch-lengths_lik.r
## Claudia March 2016

library(ape)
source("branch-length_lik.r")

who = "birds"
d=read.dna("../datasets/birds4-clean.phy") #needs to be 4 taxa
seq1 = as.vector(unname(as.character(d[1,])))
seq2 = as.vector(unname(as.character(d[2,])))
seq3 = as.vector(unname(as.character(d[3,])))
seq4 = as.vector(unname(as.character(d[4,])))

## remove missing:
s1 <-seq1!="-"
s2 <- seq2!="-"
s3 <- seq3!="-"
s4 <- seq4!="-"
seq1 <- seq1[s1&s2&s3&s4]
seq2 <- seq2[s1&s2&s3&s4]
seq3 <- seq3[s1&s2&s3&s4]
seq4 <- seq4[s1&s2&s3&s4]

d1x = 0.11
d2x = 0.078
d3y = 0.091
d4y = 0.098
dxy = 0.026

seq1.dist = seqMatrix(seq1)
seq2.dist = seqMatrix(seq2)
seq3.dist = seqMatrix(seq3)
seq4.dist = seqMatrix(seq4)

out12 = countsMatrix(seq1,seq2)
r = rep(1,6)
Q = optim.gtr(out12,r)


## 1------2 (d12)
d12=d1x+d2x
P = matrixExp(Q$Q,d12)
suma=0
for(i in 1:ncol(seq1.dist)){
    suma = suma + log(Q$Q$p[which(seq1.dist[,i]==1)]*P[which(seq1.dist[,i]==1),which(seq2.dist[,i]==1)])
}
suma   # -2972.42

suma=0
for(i in 1:ncol(seq1.dist)){
    A=seq1.dist[,i] %*% t(seq2.dist[,i])
    suma = suma + log(sum(Q$Q$p * P * A))
}
suma # -2972.42

suma=0
for(i in 1:ncol(seq1.dist)){
    f = fk(seq1.dist[,i],seq2.dist[,i],Q,d12)
    suma = suma + log(f$fk)
}
suma # -2972.42


siteLik(d1x,d2x,seq1.dist[,1],seq2.dist[,1],Q$Q)
P1 = matrixExp(Q$Q,d1x)
P2 = matrixExp(Q$Q,d2x)
P1[,3]*P2[,3]


## (1,2),3,4
l = gtr.log.lik.all(d1x,d2x,dxy,d3y,d4y,seq1.dist,seq2.dist, seq3.dist,seq4.dist,Q)
l ## -4405.753



seqx.dist = sequenceDist(d1x,d2x,seq1.dist, seq2.dist, Q)
seqy.dist = sequenceDist(d3y,d4y,seq3.dist, seq4.dist, Q)
l2 = loglik(seqx.dist, seqy.dist, Q, dxy)
l2$ll ## -4405.753

