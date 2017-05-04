## julia script to run several bistro runs sequentially
## Claudia May 2017
## cd("/Users/Clauberry/Documents/phylo/projects/present/CladeCondProb/ccdprobs/Bistro/Examples/Artiodactyl")
## stem="artiodactyl-6"
## out = readstring(`../../Code/bistro/bistro -f ../../Data/$stem.fasta -o test -r 100 -b 100`);
## f = open("test.log","w")
## write(f,out)
## close(f)

function runBistro(stem)
    try
        out = readstring(`../../Code/bistro/bistro -f ../../Data/datasets/$stem.fasta -o bistro$stem -r 1000 -b 1000`);
        f = open("bistro$stem.log","w")
        write(f,out)
        close(f)
    catch
        warn("Error in Bistro for $stem")
    end
end

function runBistroFixT(stem,tree)
    try
        out = readstring(`../../Code/bistro/bistro -f ../../Data/datasets/$stem.fasta -o bistro$stem -r 1000 -b 1000 -t $tree`);
        f = open("bistroFixT$stem.log","w")
        write(f,out)
        close(f)
    catch
        warn("Error in Bistro fixed tree for $stem")
    end
end


cd("/Users/Clauberry/Documents/phylo/projects/present/CladeCondProb/ccdprobs/Bistro/Examples/datasets")
runBistro("024")
runBistro("027")
runBistro("036")
runBistro("041")
runBistro("043")
