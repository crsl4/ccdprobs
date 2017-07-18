tab0 = read.table("test-0.par")
tab1 = read.table("test-1.par")
tab2 = read.table("test-2.par")
tab3 = read.table("test-3.par")
head(tab1)

n = length(tab0$V1)
start = n/2
logl0 = tab0$V1[start:n]
logl1 = tab1$V1[start:n]
logl2 = tab2$V1[start:n]
logl3 = tab3$V1[start:n]

m0 = mean(logl0)
m1 = mean(logl1)
m2 = mean(logl2)
m3 = mean(logl3)

v0 = var(logl0)
v1 = var(logl1)
v2 = var(logl2)
v3 = var(logl3)


meanVar = mean(c(v0,v1,v2,v3)) ## w
varMean = (n-start)*var(c(m0,m1,m2,m3)) ## b

v = (1-1/(n-start))*meanVar + (1/(n-start))*varMean
rstat = sqrt(v/meanVar)
rstat
