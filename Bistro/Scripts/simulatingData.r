## r script to simulate data from the cats-dogs tree
## Claudia September 2016

source("old/genetrees.R")
source("old/branch-length_lik.R")
source("old/transformRates.r")

## ===============================================
## artiodactyl

seed = 0558
set.seed(seed)

fileRoot = "sim-16-500"
##nsites = 1500
nsites = 500
alpha = 1000
#tree.text = "(Camel:1.460221e-01,Dolphin:1.214757e-01,(Pig:1.072255e-01,(Cow:8.941711e-02,(Sheep:8.610459e-02,Giraffe:7.642583e-02):1.846791e-02):5.146787e-02):2.783294e-02);"
tree.text = "((((A:0.0025,B:0.0025):0.0025,(C:0.0025,D:0.0025):0.0025):0.0325,((E:0.0125,F:0.0125):0.0125,(G:0.0125,H:0.0125):0.0125):0.0125):0.0625,(((((((I:0.0125,J:0.0125):0.0125,K:0.025):0.0125,L:0.0375):0.0125,M:0.05):0.0125,N:0.0625):0.0125,O:0.075):0.0125,P:0.0875):0.0125);"

tree = read.tree(text=tree.text)
ntax = length(tree$tip.label)

## from artiodactyl-6.nex.pstat
#p = c(0.2870113,0.2715106,0.1457293,0.2957489)
#r = c(0.2116111, 0.2557697, 0.6806020, 0.3338476, 0.4193366, 0.1183771)

#s = getS(c(r,p))
#s = s[1:6]
p = c(0.3,0.3,0.1,0.3)
s = c(0.1,0.3,0.1,0.1,0.3,0.1)
Q = makeQfromS(s,p,4)

##r=r/r[6]
##Q = makeQ(r,p,4)

x = simData(tree,nsites,Q,alpha)
dat = convertData(x,ntax)

## writing fasta file
##filename = paste0("../Data/sim-artiodactyl-6.fasta")
##filename = paste0("../Data/sim-16.fasta")
filename = paste0("../Data/",fileRoot,".fasta")
for(i in 1:ntax){
    if(i == 1){
        write(paste0(">",tree$tip.label[i]), file=filename)
    } else{
        write(paste0(">",tree$tip.label[i]), file=filename, append=TRUE)
    }
    write(paste0(dat[[i]],collapse=""), file=filename, append=TRUE)
    write("", file=filename, append=TRUE)
}

## writing nex file
##filename = paste0("../Examples/Simulations/sim-artiodactyl-6.nex")
filename = paste0("../Data/",fileRoot,".nex")

## need to write a nex file as well
write("#NEXUS", file=filename)
write("", file=filename,append=TRUE)
write("begin taxa;", file=filename, append=TRUE)
write(paste0("dimensions ntax=",ntax,";"), file=filename, append=TRUE)
write("taxlabels", file=filename, append=TRUE)
for(i in 1:ntax){
    write(tree$tip.label[i],file=filename,append=TRUE)
}
write(";", file=filename, append=TRUE)
write("end;", file=filename, append=TRUE)
write("", file=filename,append=TRUE)

write("begin characters;", file=filename, append=TRUE)
write(paste0("dimensions nchar=",nsites,";"), file=filename, append=TRUE)
write("format datatype=dna gap=-;", file=filename, append=TRUE)
write("matrix", file=filename, append=TRUE)
for(i in 1:ntax){
    seq = paste0(dat[[i]], collapse="")
    write(paste0(sprintf("%-20s",paste0(tree$tip.label[i])),seq), file=filename, append=TRUE)
}
write(";", file=filename, append=TRUE)
write("end;", file=filename, append=TRUE)
write("begin mrbayes;", file=filename, append=TRUE)
write("      set autoclose=yes nowarn=yes;", file=filename, append=TRUE)
write("      lset nst=6;", file=filename, append=TRUE)
write("      outgroup 1;", file=filename, append=TRUE)
write("", file=filename, append=TRUE)
write("      mcmc ngen=1100000 printfreq=10000 samplefreq=50;", file=filename, append=TRUE)
write("      sumt burnin=100000;", file=filename, append=TRUE)
write("      sump burnin=100000;", file=filename, append=TRUE)
write("end;", file=filename, append=TRUE)



if ( FALSE )
    {
## =======================================================================
## cats and dogs
seed = 1234
set.seed(seed)

nsites = 1500
alpha = 1000
## mle tree with mle bl
tree.text = "(cheetah:0.11850973637208513101,((((snow_leopard:0.04020488777776567990,leopard:0.03846672365840048818):0.01445254156731264929,tiger:0.07079306712878565000):0.01190623639760595223,clouded_leopard:0.10461902411036745619):0.04344639957238824457,(red_fox:0.11974055940327851810,(((coyote:0.00840558068745050208,(gray_wolf:0.00206050882985083861,dog:0.00185256446369396789):0.03205946058703370433):0.02609285257533808938,dhole:0.07049077201732806275):0.13276609809571754406,raccoon_dog:0.15542990325076813662):0.07955504846187926027):0.79869116234474835103):0.03995629108638096977,cat:0.03751335233479641956):0.0;"

## mle tree with artificially bigger bl for dog and wolf
## no branch<0.03
tree.text = "(cheetah:0.11850973637208513101,((((snow_leopard:0.04020488777776567990,leopard:0.03846672365840048818):0.01445254156731264929,tiger:0.07079306712878565000):0.03190623639760595223,clouded_leopard:0.10461902411036745619):0.04344639957238824457,(red_fox:0.11974055940327851810,(((coyote:0.0840558068745050208,(gray_wolf:0.0306050882985083861,dog:0.0385256446369396789):0.03205946058703370433):0.03609285257533808938,dhole:0.07049077201732806275):0.13276609809571754406,raccoon_dog:0.15542990325076813662):0.07955504846187926027):0.79869116234474835103):0.03995629108638096977,cat:0.03751335233479641956):0.0;"

## mle tree with artificially bigger bl for dog and wolf
## no branch<0.03, but also small internal bl
tree.text = "(cheetah:0.11850973637208513101,((((snow_leopard:0.04020488777776567990,leopard:0.03846672365840048818):0.01445254156731264929,tiger:0.07079306712878565000):0.03190623639760595223,clouded_leopard:0.10461902411036745619):0.04344639957238824457,(red_fox:0.11974055940327851810,(((coyote:0.0840558068745050208,(gray_wolf:0.0306050882985083861,dog:0.0385256446369396789):0.03205946058703370433):0.03609285257533808938,dhole:0.07049077201732806275):0.13276609809571754406,raccoon_dog:0.15542990325076813662):0.07955504846187926027):0.179869116234474835103):0.03995629108638096977,cat:0.03751335233479641956):0.0;"

tree = read.tree(text=tree.text)
ntax = length(tree$tip.label)

## from cats-dogs-no-fixed-q.nex.pstat
p = c(0.2429295,0.2356097,0.2067433,0.3147175)
r = c(0.05274190, 0.3742855, 0.05746096, 0.01611567, 0.4916158, 0.007780191)

s = getS(c(r,p))
s = s[1:6]
Q = makeQfromS(s,p,4) ## Error in eigen(S, symmetric = TRUE) (from branch-length_lik.R#26) : infinite or missing values in 'x'

##r=r/r[6]
##Q = makeQ(r,p,4)

x = simData(tree,nsites,Q,alpha)
dat = convertData(x,ntax)

## writing fasta file
filename = paste0("../Data/sim-cats-dogs.fasta")
filename = paste0("../Data/sim-cats-dogs-long.fasta")
filename = paste0("../Data/sim-cats-dogs-long-short.fasta")
for(i in 1:ntax){
    if(i == 1){
        write(paste0(">",tree$tip.label[i]), file=filename)
    } else{
        write(paste0(">",tree$tip.label[i]), file=filename, append=TRUE)
    }
    write(paste0(dat[[i]],collapse=""), file=filename, append=TRUE)
    write("", file=filename, append=TRUE)
}

## writing nex file
filename = paste0("../Examples/Simulations/sim-cats-dogs.nex")

## need to write a nex file as well
write("#NEXUS", file=filename)
write("", file=filename,append=TRUE)
write("begin taxa;", file=filename, append=TRUE)
write(paste0("dimensions ntax=",ntax,";"), file=filename, append=TRUE)
write("taxlabels", file=filename, append=TRUE)
for(i in 1:ntax){
    write(tree$tip.label[i],file=filename,append=TRUE)
}
write(";", file=filename, append=TRUE)
write("end;", file=filename, append=TRUE)
write("", file=filename,append=TRUE)

write("begin characters;", file=filename, append=TRUE)
write(paste0("dimensions nchar=",nsites,";"), file=filename, append=TRUE)
write("format datatype=dna gap=-;", file=filename, append=TRUE)
write("matrix", file=filename, append=TRUE)
for(i in 1:ntax){
    seq = paste0(dat[[i]], collapse="")
    write(paste0(sprintf("%-20s",paste0(tree$tip.label[i])),seq), file=filename, append=TRUE)
}
write(";", file=filename, append=TRUE)
write("end;", file=filename, append=TRUE)
write("begin mrbayes;", file=filename, append=TRUE)
write("      set autoclose=yes nowarn=yes;", file=filename, append=TRUE)
write("      lset nst=6;", file=filename, append=TRUE)
write("      outgroup 1;", file=filename, append=TRUE)
write("", file=filename, append=TRUE)
write("      mcmc ngen=1100000 printfreq=10000 samplefreq=50;", file=filename, append=TRUE)
write("      sumt burnin=100000;", file=filename, append=TRUE)
write("      sump burnin=100000;", file=filename, append=TRUE)
write("end;", file=filename, append=TRUE)

## making new names for the real data
sim.cd = read.dna("../Data/sim-cats-dogs.fasta",format="fasta")
real.cd = read.dna("../Data/cats-dogs.fasta",format="fasta")
## rename with short names to match
cd.names = c("cat","cheetah","clouded_leopard","snow_leopard","leopard","tiger","dog","gray_wolf","coyote","dhole","red_fox","raccoon_dog")
dimnames(real.cd)[[1]] = cd.names

## part of code
require(phangorn)
parsimony(tree,as.phyDat(real.cd))
parsimony(tree,as.phyDat(sim.cd))
round(as.matrix(dist.dna(sim.cd,model="TN93")),3)
round(as.matrix(dist.dna(real.cd,model="TN93")),3)
round(cophenetic.phylo(tree),3)

}
