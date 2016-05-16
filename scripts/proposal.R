## Bret Larget
## May 11, 2016
## Code to test the derivation

## Example from fitted dog/wolf/coyote example
n = 10000

mu = c(0.00194982,0.0429325,0.00262431)
vc = matrix(c(1.39396e-06,-1.20382e-08,-1.10568e-07,-1.20382e-08,3.0445e-05,5.98742e-09,-1.10568e-07,5.98742e-09,1.82874e-06),3,3)
L = t(chol(vc))

a1 = mu[1]^2/L[1,1]^2
b1 = mu[1]/L[1,1]^2

t1 = rgamma(n,a1,b1)

z1 = (t1-mu[1])/L[1,1]

a2 = (mu[2] + L[2,1]*z1)^2/L[2,2]^2
b2 = (mu[2] + L[2,1]*z1)/L[2,2]^2

t2 = rgamma(n,a2,b2)

z2 = (t2 - mu[2] - L[2,1]*z1)/L[2,2]

a3 = (mu[3] + L[3,1]*z1 + L[3,2]*z2)^2/L[3,3]^2
b3 = (mu[3] + L[3,1]*z1 + L[3,2]*z2)/L[3,3]^2

t3 = rgamma(n,a3,b3)

t.all = cbind(t1,t2,t3)
m = apply(t.all,2,mean)
v = cov(t.all)




