## julia script to summarize branch lengths in *.treeBL
## Claudia August 2016

using PhyloNetworks

trees = readMultiTopology("../Examples/Artiodactyl/pars1k.treeBL");

## Now we want to find the mean BL (without weights) for edge above {4,6} in trees2 and {5-6} in trees1
## 7	...***
## 8	..****
## 9	....**
## 10	...*.*

taxa = tipLabels(tree1) ## 2,5,6,4,3,1
split7 = [false,true,true,true,false,false]
split8 = [false,true,true,true,true,false]
split9 = [false,true,true,false,false,false]
split10 = [false,false,true,true,false,false]

## WEIGHTS
using DataFrames
dat = readtable("../Examples/Artiodactyl/pars1k", header=true, separator=' ')
logwt = convert(Array,dat[:logWt])
logwt = logwt - maximum(logwt)
wt = exp(logwt)
wt = wt / sum(wt)
sum(wt)

wt7 = Float64[]
wt8 = Float64[]
wt9 = Float64[]
wt10 = Float64[]

bl7 = Float64[]
bl8 = Float64[]
bl9 = Float64[]
bl10 = Float64[]
i = 1
for(t in trees)
    for(e in t.edge)
        if(hardwiredCluster(e,taxa) == split9)
            push!(bl9,e.length)
            push!(wt9, wt[i])
        elseif(hardwiredCluster(e,taxa) == split10)
            push!(bl10,e.length)
            push!(wt10, wt[i])
        elseif(hardwiredCluster(e,taxa) == split7)
            push!(bl7,e.length)
            push!(wt7, wt[i])
        elseif(hardwiredCluster(e,taxa) == split8)
            push!(bl8,e.length)
            push!(wt8, wt[i])
        end
    end
    i += 1
end

## length[7]	5.078913e-02
## length[8]	2.752901e-02
## length[9]	1.871516e-02
## length[10]	1.797566e-02

mean(bl7) ## 0.04911991193415638
mean(bl8) ## 0.02660738650793651
mean(bl9) ## 0.0194147
mean(bl10) ## 0.0166257

std(bl7) ## 0.008629417071962504
std(bl8) ## 0.00801165488110768
std(bl9) ## 0.006019359625570532
std(bl10) ## 0.006074630766935639

dot(bl7,wt7) ## 0.02284941008057233
dot(bl8,wt8) ## 0.005876769453781266
dot(bl9,wt9) ## 0.011324144495792485
dot(bl10,wt10) ## 0.0025176030897723295




#############################################################################
tree1 = readTopology("(2,(((6,5),4),3),1);")
tree2 = readTopology("(2,((5,(6,4)),3),1);")

equal1 = repeat([-1],inner = [length(trees)])
equal2 = repeat([-1],inner = [length(trees)])
for(i in 1: length(trees))
    equal1[i] = hardwiredClusterDistance(trees[i], tree1, false)
    equal2[i] = hardwiredClusterDistance(trees[i], tree2, false)
end

f(x) = x == 0
trees1 = trees[find(f,equal1)]
length(trees1) ## 365
trees2 = trees[find(f,equal2)]
length(trees2) ## 197

## Now we want to find the mean BL (without weights) for edge above {4,6} in trees2 and {5-6} in trees1
taxa = tipLabels(tree1) ## 2,5,6,4,3,1

bl56 = 0
for(t in trees1)
    for(e in t.edge)
        if(hardwiredCluster(e,taxa) == [false,true,true,false,false,false])
            bl56 += e.length
        end
    end
end
bl56/length(trees1) ## 0.01315549 (mb 0.0187152)

bl46 = 0
for(t in trees2)
    for(e in t.edge)
        if(hardwiredCluster(e,taxa) == [false,false,true,true,false,false])
            bl46 += e.length
        end
    end
end
bl46/length(trees2) ## 0.0132982 (mb 0.017977)

