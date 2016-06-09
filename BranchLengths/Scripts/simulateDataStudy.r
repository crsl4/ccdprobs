## R script to simulate data with simData() in genetrees.R
## for ntaxa=3,4,...,12 to run bl
## Claudia June 2016
## run as R CMD BATCH simulateDataStudy.r rep seed
## and it will create all simulated datasets by selecting a subset of ntax species
## from the cats-dogs example

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
    rep = 1
    seed = 1234
}else{
    rep = args[1]
    seed = args[2]
}

seed = 1234 # to debug first
source("genetrees.R")
source("../../scripts/branch-length_lik.R")
ntaxa = c(3,4,5,6,7,8,9,10,11,12)
nsites = 1500
alpha = 1000
tree.text = "(cheetah:0.11850973637208513101,((((snow_leopard:0.04020488777776567990,leopard:0.03846672365840048818):0.01445254156731264929,tiger:0.07079306712878565000):0.01190623639760595223,clouded_leopard:0.10461902411036745619):0.04344639957238824457,(red_fox:0.11974055940327851810,(((coyote:0.00840558068745050208,(gray_wolf:0.00206050882985083861,dog:0.00185256446369396789):0.03205946058703370433):0.02609285257533808938,dhole:0.07049077201732806275):0.13276609809571754406,raccoon_dog:0.15542990325076813662):0.07955504846187926027):0.79869116234474835103):0.03995629108638096977,cat:0.03751335233479641956):0.0;"
tree = read.tree(text=tree.text)
set.seed(seed)

## p = c(0.267,0.235,0.181,0.317)
## r = c(0.0108929,0.3428420,0.0398413,0.0106884,0.5949126,0.0008228)
## den = r[6]
## r = r/den
## Q = makeQ(r,p,4, rescale=TRUE)

for(j in 1:10){
    ntax = ntaxa[j]
    species <- sample(tree$tip.label,size=ntax)
    pruned.tree<-drop.tip(tree, setdiff(tree$tip.label, species))
    nl = length(pruned.tree$edge.length)
    pruned.tree$edge.length <- rep(0.15,nl) ## cheating to make bl comparable
    tr2 <- pruned.tree # copy tree
    tr2$edge.length <- NULL
    tr2 <- unroot(tr2)

    Q = randomQ(4,rescale=TRUE)
    r = Q$r
    p = Q$p

    ## reparametrize
    m = p %*% t(p)
    m2 = m[row(m)>col(m)]
    r = m2*r*2
    r = r/sum(r)


    setfile = paste0("../Examples/Simulations/settings",rep,".txt")
    if(j == 1){
        write(paste("ntaxa =",ntax),file=setfile)
    }else {
        write(paste("ntaxa =",ntax),file=setfile,append=TRUE)
    }
    write(paste("treeBL =",write.tree(pruned.tree)),file=setfile,append=TRUE)
    write(paste("tree =",write.tree(tr2)),file=setfile,append=TRUE)
    write(paste("rates =",paste(r[1],r[2],r[3],r[4],r[5],r[6],sep=",")),file=setfile,append=TRUE)
    write(paste("probs =",paste(p[1],p[2],p[3],p[4],sep=",")),file=setfile,append=TRUE)

    ## save list of all trees
    write(write.tree(pruned.tree),file="trees.txt", append=TRUE)

    ## x has nrows=number of nodes, the first n rows correspond to the tips, n+1 corresponds to the root
    x = simData(pruned.tree,nsites,Q,alpha)
    ## dat only has nrows=number of tips, tips ordered by node number, the actual taxa can be obtained with
    ## pruned.tree$tip.label, that is,
    ## dat[[1]] is the sequence for taxa tree$tip.label[1]
    dat = convertData(x,ntax)

    filename = paste0("../Data/sim",ntax,".fasta")
    for(i in 1:ntax){
        if(i == 1){
            write(paste0(">",pruned.tree$tip.label[i]), file=filename)
        } else{
            write(paste0(">",pruned.tree$tip.label[i]), file=filename, append=TRUE)
        }
        write(paste0(dat[[i]],collapse=""), file=filename, append=TRUE)
        write("", file=filename, append=TRUE)
    }

}
