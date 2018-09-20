## julia script to run bistro,
## based on runningBistro.jl
## running inside Bistro/Simulations
## run: julia runningBistro-paper.jl 024 r b 0
## the last argument is 0=not fixed tree, 1=fixed tree
## Claudia September 2018

cd("/Users/Clauberry/Documents/github/CladeCondProb/ccdprobs/Bistro/Simulations")
include("../Scripts/runningBistro-fns.jl")

r = 1000
b = 1000
dat = "024"
fix = 0
if(length(ARGS) > 0)
    dat = ARGS[1]
    r = parse(Int,ARGS[2])
    b = parse(Int,ARGS[3])
    fix = parse(Int,ARGS[4])
end


stem = string("sim-",dat,"-nsites-1500")

if(fix == 0)
    runBistro(stem,r,b)
elseif(fix == 1)
    treefile = string("../Data/datasets/ExaBayes_ConsensusExtendedMajorityRuleNewick.",dat,".cons")
    f = open(treefile)
    tree = readlines(f)[1]
    close(f)
    ## remove "1.00:"
    tree2 = replace(tree,r"[0-9].[0-9]+:" => ":")
    runBistroFixT(stem,tree2,r,b)
end
