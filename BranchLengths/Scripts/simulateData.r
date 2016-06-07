## R script to simulate data with 4 taxa to analyze with bl
## Claudia May 2016

library(ape)
source('../../scripts/branch-length_lik.r')
source('../../scripts/4taxa_functions.r')
library(ggplot2)
library(weights)
library(mvtnorm)

seed = 0646
set.seed(seed)
d1x0=0.11
d2x0=0.078
dxy0 = 0.03
d3y0 = 0.091
d4y0 = 0.098
## d1x0=0.5
## d2x0=0.5
## dxy0 = 0.03
## d3y0 = 0.5
## d4y0 = 0.5
eta = 0.5
nsites=1500
nuc <- c('a','c','g','t')
Q = randomQ(4,rescale=TRUE)
r=Q$r
r
p=Q$p
p
## simulate seqx
seqx = sample(nuc,size=nsites,prob=Q$p,replace=TRUE)

## simulate seq1
P = matrixExp(Q,d1x0)
seq1 = numeric(nsites)
for ( i in 1:nsites )
    seq1[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])
## simulate seq2
P = matrixExp(Q,d2x0)
seq2 = numeric(nsites)
for ( i in 1:nsites )
    seq2[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])

## simulate seqy
P = matrixExp(Q,dxy0)
seqy = numeric(nsites)
for ( i in 1:nsites )
    seqy[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])

## simulate seq3
P = matrixExp(Q,d3y0)
seq3 = numeric(nsites)
for ( i in 1:nsites )
    seq3[i] = sample(nuc,size=1,prob=P[which(nuc==seqy[i]),])
## simulate seq4
P = matrixExp(Q,d4y0)
seq4 = numeric(nsites)
for ( i in 1:nsites )
    seq4[i] = sample(nuc,size=1,prob=P[which(nuc==seqy[i]),])

## write to file:
rootname = paste0("simSeq",seed)
filename = paste0(rootname,".fasta")
l2 = paste0(seq1,collapse="")
l3 = paste0(seq2,collapse="")
l4 = paste0(seq3,collapse="")
l5 = paste0(seq4,collapse="")

write(">1",file=filename)
write(l2,file=filename,append=TRUE)
write("",file=filename,append=TRUE)
write(">2",file=filename, append=TRUE)
write(l3,file=filename, append=TRUE)
write("",file=filename,append=TRUE)
write(">3",file=filename, append=TRUE)
write(l4,file=filename, append=TRUE)
write("",file=filename,append=TRUE)
write(">4",file=filename, append=TRUE)
write(l5,file=filename, append=TRUE)

