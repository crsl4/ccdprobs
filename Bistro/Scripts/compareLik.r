## R script to compare lik---?-?.out files
## by checking the two columns of logLik and logLik0
## Claudia September 2016


file = "../Examples/Artiodactyl/lik---0-249.out"
file = "../Examples/Artiodactyl/lik---250-499.out"
file = "../Examples/Artiodactyl/lik---500-749.out"
file = "../Examples/Artiodactyl/lik---750-999.out"

dat1 = read.table(file, sep=" ", header=TRUE)
str(dat1)

dat1 <- within(dat1, diff <- logl-logliknew)
diffInd <- which(dat1$diff != 0)
length(diffInd)
dat1[5,]
dat1[diffInd[4],]
hist(dat1$diff)

which(dat1$diff < -7)
dat1[223,]

dat1[diffInd,1:2]
dat1[-diffInd,1:2]
