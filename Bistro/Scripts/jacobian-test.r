## starting with multivariate normals

d1 = 8
seed = 118945
set.seed(seed)
rho = runif(d1*(d1-1)/2) ## correlations
C = matrix(0,d1,d1)
C[row(C) > col(C)] = rho
C = C + t(C)
diag(C) = rep(1,d1)
C == t(C)
library(matrixcalc)
is.positive.definite(C) ## FALSE

## genPositiveDefMat {clusterGeneration}
## genPositiveDefMat(dim, covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"),
##                   alphad=1, eta=1, rangeVar=c(1,10), lambdaLow=1, ratioLambda=10)
library(clusterGeneration)
seed = 118945
set.seed(seed)
C = genPositiveDefMat(dim=d1)

library(mvtnorm)
##dmvnorm(x, mean, sigma, log=FALSE)
x = rmvnorm(n=1, mean=rep(0,d1), sigma=C)
