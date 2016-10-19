## julia script to read root---*.treeBL files after bistro
## to compare with mrbayes.nex.vstat
## needs to read mrbayes.nex.parts to identify the important splits
## based on old/sumamrizeBL.jl
## creates bistroroot.vstat to compare to mrbayes.vstat
## Claudia September 2016

## folder = "../Examples/Artiodactyl/"
## mbroot = "artiodactyl-6"
## folder = "../Examples/cats-dogs/"
## mbroot = "cats-dogs-no-fixed-q"
folder = "../Examples/Simulations/"
mbroot = "sim-cats-dogs"
bistroroot = "eta11"
verbose = false

using PhyloNetworks
using DataFrames

searchdirTreeBL(path,key) = filter(x->(contains(x,key) && contains(x,".treeBL")), readdir(path))
searchdirOut(path,key) = filter(x->(contains(x,key) && contains(x,".out")), readdir(path))

## read all trees with branch lengths
files = searchdirTreeBL(folder,"$bistroroot---")
trees = readMultiTopology("$folder$(files[1])");
for j=2:length(files)
    t = readMultiTopology("$folder$(files[j])");
    trees = vcat(trees,t)
end
taxa = map(x->parse(Int,x),tipLabels(trees[1]))
taxa = sort(taxa)
taxa = map(x->string(x),taxa)


## read weights
files = searchdirOut(folder,"$bistroroot---")
dat = readtable("$folder$(files[1])", header=true, separator=' ');
for j=2:length(files)
    t = readtable("$folder$(files[j])", header=true, separator=' ')
    dat = vcat(dat,t)
end

logwt = convert(Array,dat[:logWt]);
logwt = logwt - maximum(logwt);
wt = exp(logwt);
wt = wt / sum(wt);

PhyloNetworks.approxEq(sum(wt),1.0) || error("weights do not sum up to 1")
length(trees) == length(wt) || error("different length of trees vector and weights vector")

mbfile = "$folder$mbroot.nex.parts"
f = open(mbfile);
lines = readlines(f);
splits = AbstractString[]
column1 = Float64[]
column2 = Float64[]
column3 = Float64[]
column4 = Float64[]
column5 = Float64[]
column6 = Float64[]
column7 = Float64[]
column8 = Float64[]

i = 1
for l in lines
    if(i > 2)
        sb = split(l)
        split1 = falses(length(sb[2])) ## cannot call split bc of function split
        for j=1:length(sb[2])
            if(sb[2][j] == '*')
                split1[j] = true
            end
        end
        verbose && println("$(sb[2])")
        verbose && println("$(split1)")

        wtsplit = Float64[]
        blsplit = Float64[]

        k = 1
        for t in trees
            verbose && println("$t")
            for e in t.edge
                clus = hardwiredCluster(e,taxa)
                if(clus == split1 || clus == !split1)
                    verbose && println("found $(clus)")
                    push!(blsplit,e.length)
                    push!(wtsplit, wt[k])
                end
            end
            k += 1
        end

        wtsplit = wtsplit / sum(wtsplit)
        push!(splits,sb[2])
        if(!isempty(blsplit))
            push!(column7,mean(blsplit))
            push!(column8,std(blsplit))
            push!(column1,dot(blsplit,wtsplit)) ## same as mean(blsplit, weights(wtsplit))
            push!(column2,var(blsplit,weights(wtsplit)))
            push!(column3,std(blsplit,weights(wtsplit)))
            push!(column4,quantile(blsplit, weights(wtsplit), 0.025))
            push!(column5,quantile(blsplit, weights(wtsplit), 0.975))
            push!(column6,median(blsplit, weights(wtsplit)))
        else
            push!(column7,0.0)
            push!(column8,0.0)
            push!(column1,0.0)
            push!(column2,0.0)
            push!(column3,0.0)
            push!(column4,0.0)
            push!(column5,0.0)
            push!(column6,0.0)
        end

    end
    i += 1
end

dig = 6
dataBL = DataFrame(bl=splits, mean=round(column1,dig), var=round(column2,dig), std=round(column3,dig), CILower=round(column4,dig), CIUpper=round(column5,dig), median=round(column6,dig), sampleMean=round(column7,dig), sampleStd=round(column8,dig))
println("$dataBL")
writetable("$folder$bistroroot.vstat",dataBL, separator='\t')

