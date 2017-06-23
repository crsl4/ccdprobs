# Checking output of MCMC within bistro with the sampled pi and s using generalized Dirichlet

require(dplyr)
require(ggplot2)

#file_head = "birge-bistroFixT-3-024"
file_head = "bistroFixT-3-024"
mcmc_out = read.table(paste0(file_head,".par"),header=FALSE)
n = 1000
burn = nrow(mcmc_out) - n
# remove burnin
#mcmc_out = mcmc_out[101:1100,]
names(mcmc_out) = c("logl",paste0("pi",1:4),paste0("s",1:6))

# find mean and variance for each column
apply(mcmc_out,2,mean)
apply(mcmc_out,2,var)
good.mean = apply(mcmc_out[-(1:350),-1],2,mean)
good.sd = apply(mcmc_out[-(1:350),-1],2,sd)
bad.mean = apply(mcmc_out[-(1:burn),-1],2,mean)
bad.sd = apply(mcmc_out[-(1:burn),-1],2,sd)

## read in IS sample
source("../../Scripts/readBistro.R")
bistro = readBistro(file_head)
b = bistro %>%
  select(num_range("pi",1:4),num_range("s",1:6))
bistro.mean = apply(b,2,mean)
bistro.sd = apply(b,2,sd)

print( rbind(good.mean,bad.mean,bistro.mean) )
print( rbind(good.sd,bad.sd,bistro.sd) )

plot( ggplot(mcmc_out[-(1:burn),],aes(x=1:n,y=logl)) + geom_point() )

plot(
  ggplot(mcmc_out[-(1:burn),],aes(x=1:n)) +
#  ggplot(mcmc_out,aes(x=-(burn-1):n)) +
  geom_point(aes(y=pi1),col="blue") +
  geom_point(aes(y=pi2),col="red") +
  geom_point(aes(y=pi3),col="green") +
  geom_point(aes(y=pi4),col="yellow")   
)
