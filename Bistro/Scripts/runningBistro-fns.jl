## functions for runningBistro-paper.jl
## running inside Bistro/Simulations
## Claudia September 2018

function runBistro(stem,r,b)
    println("+++++++++++ running bistro for $stem +++++++++++++")
    try
        out = readstring(`../Code/bistro/bistro -f ../Data/datasets/$stem.fasta -o bistro-$stem-$r-$b -r $r -b $b`);
        f = open("bistro-$stem-$r-$b.log","w")
        write(f,out)
        close(f)
    catch
        warn("Error in Bistro for $stem")
    end
end

function runBistroFixT(stem,tree,r,b)
    println("+++++++++++ running bistro for $stem for fixed tree +++++++++++++")
    try
        out = readstring(`../Code/bistro/bistro -f ../Data/datasets/$stem.fasta -o bistroFixT-$stem-$r-$b -r $r -b $b -t "$tree"`);
        f = open("bistroFixT-$stem-$r-$b.log","w")
        write(f,out)
        close(f)
    catch
        warn("Error in Bistro fixed tree for $stem")
    end
end

