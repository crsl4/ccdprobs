## Code to get the bistro input for Q matrix from MrBayes posterior sample
## Bret Larget
## June 8, 2016

r = c(2.118317e-01,2.565228e-01,6.828224e-02,3.311636e-02,4.182256e-01,1.202127e-02)
p = c(2.865626e-01,2.718173e-01,1.456956e-01,2.959245e-01)

p = p/sum(p)

temp = p %o% p

pp = temp[row(temp) > col(temp)]

s = r*pp

s = s/sum(s)


