## r script to compare the TN approach to the likelihood
## on the internal branch length
## Claudia February 2016

## x-----y, bl=t
source('branch-length.r')
nsites=500
branch.length=0.15
eta.jc=0.5
eta.tn=0.5

p1 = doit(nsites=nsites,branch.length=branch.length,eta.jc=eta.jc, eta.tn=eta.tn)
plot(p1)

Q = randomQ(4)
branch.length=0.02
nsites=1500
x = simulateSequenceSummary(nsites,Q,branch.length)
print(x)
Qlist.gtr = optim.gtr(x,Q$r)
log.like2 = gtr.log.like(x,branch.length,Qlist.gtr$Q)

P = matrixExp(Q,branch.length)
P
diag(Q$p)
diag(Q$p) %*% P
log(diag(Q$p) %*% P)
x * log(diag(Q$p) %*% P)


## A>x----y, bl=t
who = "birds"
d=read.dna("../datasets/birds4-clean.phy") #needs to be 4 taxa
seq1 = as.vector(unname(as.character(d[1,])))
seq2 = as.vector(unname(as.character(d[2,])))
seq3 = as.vector(unname(as.character(d[3,])))

## remove missing:
s1 <-seq1!="-"
s2 <- seq2!="-"
s3 <- seq3!="-"
seq1 <- seq1[s1&s2&s3]
seq2 <- seq2[s1&s2&s3]
seq3 <- seq3[s1&s2&s3]

nsites= length(seq3)

##need seq1.dist and seq2.dist which are matrices
n=4
seq1.dist = matrix(0,n,nsites)
seq2.dist = matrix(0,n,nsites)
seq3.dist = matrix(0,n,nsites)
for(i in 1:nsites){
    if(seq1[i] == 'a'){
        seq1.dist[1,i] = 1
    } else if(seq1[i] == 'c'){
        seq1.dist[2,i] = 1
    } else if(seq1[i] == 'g'){
        seq1.dist[3,i] = 1
    } else if(seq1[i] == 't'){
        seq1.dist[4,i] = 1
    }
    if(seq2[i] == 'a'){
        seq2.dist[1,i] = 1
    } else if(seq2[i] == 'c'){
        seq2.dist[2,i] = 1
    } else if(seq2[i] == 'g'){
        seq2.dist[3,i] = 1
    } else if(seq2[i] == 't'){
        seq2.dist[4,i] = 1
    }
    if(seq3[i] == 'a'){
        seq3.dist[1,i] = 1
    } else if(seq3[i] == 'c'){
        seq3.dist[2,i] = 1
    } else if(seq3[i] == 'g'){
        seq3.dist[3,i] = 1
    } else if(seq3[i] == 't'){
        seq3.dist[4,i] = 1
    }
}



nuc = c("a","c","g","t")
out12 = matrix(0,n,n) # distance between 1 and 2
for(i in 1:n)
    for(j in 1:n)
        out12[i,j] = sum(seq1==nuc[i] & seq2==nuc[j])
print(out12)
checkMatCounts(out12)

r = rep(1,6)
Q = optim.gtr(out12,r)

d1x = 0.11
d2x = 0.078
dxy=0.1 # need to vary this

suma = 0
for(site in 1:nsites){
    print(site)
    pA.given.x = siteLik(d1x,d2x,seq1.dist[,site], seq2.dist[,site], Q$Q)
    print(pA.given.x)
    Q$Q$p * pA.given.x
    col=which(seq3.dist[,site] == 1)
    Pxy = matrixExp(Q$Q,dxy)
    print(Pxy)
    print(Pxy[,col])
    print(log(sum(Q$Q$p * pA.given.x * Pxy[,col])))
    suma = suma + log(sum(Q$Q$p * pA.given.x * Pxy[,col]))
}

## to do: need to vary dxy, and create a function to comput the lik for a vector of dxy
## similar to gtr.log.like, and compare it to density TN
## follow steps in branch-length.r
## plot densities
## see if TN is centered in 0.1, and the other one is centered in 0.02, as it should
