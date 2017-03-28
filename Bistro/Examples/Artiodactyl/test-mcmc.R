# test MCMC p and s

p = read.table("p.txt")
names(p) = paste0("pi",c("A","C","G","T"))
nrow(p)
sapply(p,mean)
sapply(p,var)
s = read.table("s.txt")
names(s) = paste0("s",c("AC","AG","AT","CG","CT","GT"))
sapply(s,mean)
sapply(s,sd)
