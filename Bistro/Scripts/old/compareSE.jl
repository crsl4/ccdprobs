## julia script to compare the SE from the crude version: p(1-p)/ESS to the
## one from importance sampling formula
## August 2016

using PhyloNetworks
using DataFrames

trees = readMultiTopology("../Examples/Artiodactyl/pars10k.treeBL");
trees = readMultiTopology("../Examples/Artiodactyl/new.txt.treeBL");
length(trees)

dat = readtable("../Examples/Artiodactyl/pars10k.out", header=true, separator=' ')
dat = readtable("../Examples/Artiodactyl/new.txt", header=true, separator=' ')
logwt = convert(Array,dat[:logWt])
logwt = logwt - maximum(logwt)
wt = exp(logwt)
wt = wt / sum(wt)
sum(wt)

tree1 = readTopology("(1,2,(3,(4,(5,6))));") ##0.5388
tree2 = readTopology("(1,2,(3,((4,6),5)));") ##0.4291

wt1 = Float64[]
wt2 = Float64[]
id1 = Int[]
id2 = Int[]

i = 1
for(t in trees)
    if(hardwiredClusterDistance(tree1,t,false) == 0)
        push!(wt1,wt[i])
        push!(id1,1)
        push!(id2,0)
    elseif(hardwiredClusterDistance(tree2,t,false) == 0)
        push!(wt2,wt[i])
        push!(id2,1)
        push!(id1,0)
    else
        push!(id1,0)
        push!(id2,0)
    end
    i += 1
end

mu1 = sum(wt1) ## 0.5388
mu2 = sum(wt2) ## 0.4291

ess = 1/sum(wt .* wt) ## 5296.71

crudeSE1 = sqrt(mu1*(1-mu1)/ess) ## 0.006849
crudeSE2 = sqrt(mu2*(1-mu2)/ess) ## 0.006801

SE1 = sqrt(sum((id1 - mu1) .* (id1-mu1) .* wt .* wt)) ## 0.006945
SE2 = sqrt(sum((id2 - mu2) .* (id2-mu2) .* wt .* wt)) ## 0.007004
