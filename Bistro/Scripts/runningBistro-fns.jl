## functions for runningBistro-paper.jl
## running inside Bistro/Simulations
## Claudia September 2018

function runBistro(stem,r,b)
    println("+++++++++++ running bistro for $stem +++++++++++++")
    try
        out = readstring(`../Code/bistro/bistro -f ../Data/datasets/$stem.fasta -o bistro-2-$stem -r $r -b $b`);
        f = open("bistro-2-$stem.log","w")
        write(f,out)
        close(f)
    catch
        warn("Error in Bistro for $stem")
    end
end

function runBistroFixT(stem,tree,r,b)
    println("+++++++++++ running bistro for $stem for fixed tree +++++++++++++")
    try
        out = readstring(`../Code/bistro/bistro -f ../Data/datasets/$stem.fasta -o bistroFixT-2-$stem -r 1000 -b 1000 -t "$tree"`);
        f = open("bistroFixT-2-$stem.log","w")
        write(f,out)
        close(f)
    catch
        warn("Error in Bistro fixed tree for $stem")
    end
end

