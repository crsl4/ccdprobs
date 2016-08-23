## julia script to summarize branch lengths in *.treeBL
## Claudia August 2016

using PhyloNetworks

##################### primates ###############################################################
trees = readMultiTopology("../Examples/Primates/pars10k.treeBL");
length(trees)

taxa = tipLabels(trees[1])
taxa = ["11", "4", "5", "6", "3", "1", "2", "12", "10", "9", "8", "7"]
split20 = [false,true,true,true,true,false,false,false,false,false,false,false]
split21 = [false,true,false,false,true,false,false,false,false,false,false,false]


## WEIGHTS
using DataFrames
dat = readtable("../Examples/Primates/pars10k.out", header=true, separator=' ')
logwt = convert(Array,dat[:logWt])
logwt = logwt - maximum(logwt)
wt = exp(logwt)
wt = wt / sum(wt)
sum(wt)


wt20 = Float64[]
wt21 = Float64[]

bl20 = Float64[]
bl21 = Float64[]



i = 1
for(t in trees)
    for(e in t.edge)
        if(hardwiredCluster(e,taxa) == split20)
            push!(bl20,e.length)
            push!(wt20, wt[i])
        elseif(hardwiredCluster(e,taxa) == split21)
            push!(bl21,e.length)
            push!(wt21, wt[i])
        end
    end
    i += 1
end

mean(bl20)
mean(bl21)


std(bl20)
std(bl21)


wt20 = wt20 / sum(wt20)
wt21 = wt21 / sum(wt21)

dot(bl20,wt20)
dot(bl21,wt21)


##################### cats ###############################################################
trees = readMultiTopology("../Examples/cats-dogs/pars10k.treeBL");
trees = readMultiTopology("../Examples/cats-dogs/pars1k-scale.treeBL");
length(trees)

taxa = tipLabels(trees[1])
taxa = ["11", "4", "5", "6", "3", "1", "2", "12", "10", "9", "8", "7"]
split13 = [false,false,false,false,false,false,false,false,false,false,true,true]
split14 = [true,false,false,false,false,false,false,true,true,true,true,true]
split15 = [false,false,false,false,false,false,false,false,false,true,true,true]
split16 = [false,true,true,false,false,false,false,false,false,false,false,false]
split17 = [false,true,true,true,true,false,false,false,false,false,false,false]
split18 = [false,false,false,false,false,false,false,false,true,true,true,true]
split19 = [false,false,false,false,false,false,false,true,true,true,true,true]
split20 = [false,true,true,true,false,false,false,false,false,false,false,false]
split21 = [true,false,false,false,false,false,true,true,true,true,true,true]
split22 = [true,true,true,true,true,false,false,true,true,true,true,true]

## WEIGHTS
using DataFrames
dat = readtable("../Examples/cats-dogs/pars10k.out", header=true, separator=' ')
dat = readtable("../Examples/cats-dogs/pars1k-scale.out", header=true, separator=' ')
logwt = convert(Array,dat[:logWt])
logwt = logwt - maximum(logwt)
wt = exp(logwt)
wt = wt / sum(wt)
sum(wt)


wt13 = Float64[]
wt14 = Float64[]
wt15 = Float64[]
wt16 = Float64[]
wt17 = Float64[]
wt18 = Float64[]
wt19 = Float64[]
wt20 = Float64[]
wt21 = Float64[]
wt22 = Float64[]

bl13 = Float64[]
bl14 = Float64[]
bl15 = Float64[]
bl16 = Float64[]
bl17 = Float64[]
bl18 = Float64[]
bl19 = Float64[]
bl20 = Float64[]
bl21 = Float64[]
bl22 = Float64[]


i = 1
for(t in trees)
    for(e in t.edge)
        if(hardwiredCluster(e,taxa) == split13)
            push!(bl13,e.length)
            push!(wt13, wt[i])
        elseif(hardwiredCluster(e,taxa) == split14)
            push!(bl14,e.length)
            push!(wt14, wt[i])
        elseif(hardwiredCluster(e,taxa) == split15)
            push!(bl15,e.length)
            push!(wt15, wt[i])
        elseif(hardwiredCluster(e,taxa) == split16)
            push!(bl16,e.length)
            push!(wt16, wt[i])
        elseif(hardwiredCluster(e,taxa) == split17)
            push!(bl17,e.length)
            push!(wt17, wt[i])
        elseif(hardwiredCluster(e,taxa) == split18)
            push!(bl18,e.length)
            push!(wt18, wt[i])
        elseif(hardwiredCluster(e,taxa) == split19)
            push!(bl19,e.length)
            push!(wt19, wt[i])
        elseif(hardwiredCluster(e,taxa) == split20)
            push!(bl20,e.length)
            push!(wt20, wt[i])
        elseif(hardwiredCluster(e,taxa) == split21)
            push!(bl21,e.length)
            push!(wt21, wt[i])
        elseif(hardwiredCluster(e,taxa) == split22)
            push!(bl22,e.length)
            push!(wt22, wt[i])
        end
    end
    i += 1
end

mean(bl13)
mean(bl14)
mean(bl15)
mean(bl16)
mean(bl17)
mean(bl18)
mean(bl19)
mean(bl20)
mean(bl21)
mean(bl22)

std(bl13)
std(bl14)
std(bl15)
std(bl16)
std(bl17)
std(bl18)
std(bl19)
std(bl20)
std(bl21)
std(bl22)



wt13 = wt13 / sum(wt13)
wt14 = wt14 / sum(wt14)
wt15 = wt15 / sum(wt15)
wt16 = wt16 / sum(wt16)
wt17 = wt17 / sum(wt17)
wt18 = wt18 / sum(wt18)
wt19 = wt19 / sum(wt19)
wt20 = wt20 / sum(wt20)
wt21 = wt21 / sum(wt21)
wt22 = wt22 / sum(wt22)

dot(bl13,wt13)
dot(bl14,wt14)
dot(bl15,wt15)
dot(bl16,wt16)
dot(bl17,wt17)
dot(bl18,wt18)
dot(bl19,wt19)
dot(bl20,wt20)
dot(bl21,wt21)
dot(bl22,wt22)


##################### artiodactyl ###############################################################
trees = readMultiTopology("../Examples/Artiodactyl/pars1k.treeBL"); ## before fixing error
trees = readMultiTopology("../Examples/Artiodactyl/new.txt.treeBL"); ## after fixing error
trees = readMultiTopology("../Examples/Artiodactyl/pars10k.treeBL"); ## after fixing error

## Now we want to find the mean BL (without weights) for edge above {4,6} in trees2 and {5-6} in trees1
## 7	...***
## 8	..****
## 9	....**
## 10	...*.*

taxa = tipLabels(trees[1]) ## 2,5,6,4,3,1
taxa = ["2","5","6","4","3","1"]
split7 = [false,true,true,true,false,false]
split8 = [false,true,true,true,true,false]
split9 = [false,true,true,false,false,false]
split10 = [false,false,true,true,false,false]

## WEIGHTS
using DataFrames
dat = readtable("../Examples/Artiodactyl/pars1k", header=true, separator=' ') ## with error
dat = readtable("../Examples/Artiodactyl/new.txt", header=true, separator=' ') ## without error
dat = readtable("../Examples/Artiodactyl/pars10k.out", header=true, separator=' ') ## without error
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

## length[7]	5.083729e-02
## length[8]	2.749753e-02
## length[9]	1.871160e-02
## length[10]	1.802578e-02

mean(bl7)
mean(bl8)
mean(bl9)
mean(bl10)

std(bl7)
std(bl8)
std(bl9)
std(bl10)

wt7 = wt7 / sum(wt7)
wt8 = wt8 / sum(wt8)
wt9 = wt9 / sum(wt9)
wt10 = wt10 / sum(wt10)

dot(bl7,wt7)
dot(bl8,wt8)
dot(bl9,wt9)
dot(bl10,wt10)




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

