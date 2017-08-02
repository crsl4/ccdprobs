## r script to simulate data from the cats-dogs tree
## Claudia September 2016

## if we want to use consensus tree from MrBayes:
## sed 's/\[[^]]*\]//g' artiodactylMBTree.txt > cleanMBTree.txt

source("old/genetrees.R")
source("old/branch-length_lik.R")
source("old/transformRates.r")

seed = 0558
set.seed(seed)
## which = "artiodactyl"
## which = "cats-dogs"
##which = "sim-16"
##which = "whales"
which = "024"
nsites = 1500
alpha = 1000

if( which == "artiodactyl" ){
    fileRoot = "sim-artiodactyl-6"
    tree.text = "(Camel:1.460221e-01,Dolphin:1.214757e-01,(Pig:1.072255e-01,(Cow:8.941711e-02,(Sheep:8.610459e-02,Giraffe:7.642583e-02):1.846791e-02):5.146787e-02):2.783294e-02);"
    ## from artiodactyl-6.nex.pstat
    p = c(0.2870113,0.2715106,0.1457293,0.2957489)
    r = c(0.2116111, 0.2557697, 0.6806020, 0.3338476, 0.4193366, 0.1183771)
    s = getS(c(r,p))
    s = s[1:6]
}else if( which == "cats-dogs"){
    fileRoot = "sim-cats-dogs"
    tree.text = "(cheetah:0.11850973637208513101,((((snow_leopard:0.04020488777776567990,leopard:0.03846672365840048818):0.01445254156731264929,tiger:0.07079306712878565000):0.01190623639760595223,clouded_leopard:0.10461902411036745619):0.04344639957238824457,(red_fox:0.11974055940327851810,(((coyote:0.00840558068745050208,(gray_wolf:0.00206050882985083861,dog:0.00185256446369396789):0.03205946058703370433):0.02609285257533808938,dhole:0.07049077201732806275):0.13276609809571754406,raccoon_dog:0.15542990325076813662):0.07955504846187926027):0.79869116234474835103):0.03995629108638096977,cat:0.03751335233479641956):0.0;"
    ## from cats-dogs-no-fixed-q.nex.pstat
    p = c(0.2429295,0.2356097,0.2067433,0.3147175)
    r = c(0.05274190, 0.3742855, 0.05746096, 0.01611567, 0.4916158, 0.007780191)
    s = getS(c(r,p))
    s = s[1:6]
}else if( which == "sim-16"){
    fileRoot = "sim-16"
    tree.text = "((((A:0.0025,B:0.0025):0.0025,(C:0.0025,D:0.0025):0.0025):0.0325,((E:0.0125,F:0.0125):0.0125,(G:0.0125,H:0.0125):0.0125):0.0125):0.0625,(((((((I:0.0125,J:0.0125):0.0125,K:0.025):0.0125,L:0.0375):0.0125,M:0.05):0.0125,N:0.0625):0.0125,O:0.075):0.0125,P:0.0875):0.0125);"
    p = c(0.3,0.3,0.1,0.3)
    s = c(0.1,0.3,0.1,0.1,0.3,0.1)
}else if( which == "whales"){
    fileRoot = "sim-whales"
    ##From mr bayes (has polytomy)
    ## tree.text = "(1:2.217398e-02,2:3.718874e-02,((3:1.005738e-01,(((4:2.672636e-02,5:2.075132e-02):2.363831e-02,14:7.496660e-02):1.007081e-02,((6:2.627892e-02,7:1.792107e-02):2.262076e-02,(((8:1.532786e-02,9:2.047032e-02):2.290625e-02,11:3.914722e-02):9.123712e-03,10:3.152600e-02,12:3.750166e-02,13:3.419180e-02):9.589063e-03):1.632728e-02):1.814481e-02):1.585957e-02,(15:1.028530e-01,(((16:5.223536e-02,17:4.210148e-02):4.477917e-02,((18:1.089527e-02,19:5.694007e-03):2.426046e-02,(20:1.025593e-02,21:9.416916e-03):2.295827e-02):5.475090e-02):7.715475e-02,((22:9.064491e-02,23:8.766028e-02):2.203347e-02,((((24:8.141102e-02,(27:7.664552e-02,28:6.969899e-02):1.887018e-02):8.089943e-03,(25:5.734373e-02,26:4.017261e-02):3.613051e-02):1.383720e-02,(29:5.281744e-02,30:1.907180e-02):3.809643e-02):1.603342e-02,31:1.351928e-01):1.913932e-02):1.270959e-02):1.763299e-02):3.786344e-02):6.752775e-02);"
    ##From bistro
    tree.text = "(1:0.0231136,2:0.0372603,((3:0.094114,(((4:0.025742,5:0.022279):0.0249355,14:0.0685793):0.00162364,((6:0.027438,7:0.0173734):0.0210568,((((8:0.0170381,9:0.0189656):0.0173591,11:0.0375513):0.00117174,(10:0.0324844,12:0.0337327):0.00110449):0.000772143,13:0.0327726):0.00476205):0.00772679):0.014146):0.00534651,(15:0.101767,((((16:0.0503164,17:0.0470883):0.0368309,((18:0.0110323,19:0.00505243):0.02207,(20:0.0116344,21:0.0071754):0.0214576):0.0500723):0.0548426,((24:0.0799893,((25:0.0553201,26:0.0455578):0.0283923,((27:0.0766089,28:0.0689195):0.00240495,(29:0.0489559,30:0.0222969):0.0290598):0.00387214):0.00247192):0.0113158,31:0.12489):0.0064916):0.00406787,(22:0.0938027,23:0.0901243):0.00790296):0.00718762):0.0331146):0.064614);"

    ##pi and s to use for simulation
    p = c(0.231333,0.225334,0.173454,0.369879)
    s = c(0.151508,0.198957,0.0548603,0.0118846,0.578106,0.00468519)
}else if( which == "024"){
    fileRoot = "sim-024"
    ##From bistro
    tree.text = "(1:0.0402746,((((2:0.0396321,3:0.0391927):0.0436553,(4:0.0398195,5:0.0425747):0.0376978):0.0116913,(((6:0.0423926,7:0.0417774):0.0499896,(8:0.0411028,(9:0.0377416,10:0.037722):0.0831132):0.00265454):0.093231,(((11:0.0361578,12:0.0371425):0.0611912,(((13:0.0392235,14:0.039483):0.0647671,15:0.0438121):0.023461,(16:0.0390756,17:0.0409884):0.0732874):0.0990023):0.0546981,((18:0.0422706,19:0.0387413):0.0472085,(20:0.0410891,21:0.0416439):0.0615112):0.021677):0.0578341):0.0368816):0.0173533,(22:0.0415874,23:0.0397845):0.0472142):0.0651573,24:0.0397275);"

    ##pi and s to use for simulation
    p = c(0.242788,0.256868,0.186394,0.313951)
    s = c(0.364941,0.159708,0.11059,0.0356572,0.304613,0.0244902)
}

fileRoot = paste0(fileRoot,"-nsites",nsites)
tree = read.tree(text=tree.text)
ntax = length(tree$tip.label)

Q = makeQfromS(s,p,4)
##r=r/r[6]
##Q = makeQ(r,p,4)
x = simData(tree,nsites,Q,alpha)
dat = convertData(x,ntax)

## writing fasta file
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
## ===============================================
## artiodactyl

seed = 0558
set.seed(seed)

nsites = 1500
alpha = 1000
tree.text = "(Camel:1.460221e-01,Dolphin:1.214757e-01,(Pig:1.072255e-01,(Cow:8.941711e-02,(Sheep:8.610459e-02,Giraffe:7.642583e-02):1.846791e-02):5.146787e-02):2.783294e-02);"

tree = read.tree(text=tree.text)
ntax = length(tree$tip.label)

## from artiodactyl-6.nex.pstat
p = c(0.2870113,0.2715106,0.1457293,0.2957489)
r = c(0.2116111, 0.2557697, 0.6806020, 0.3338476, 0.4193366, 0.1183771)

s = getS(c(r,p))
s = s[1:6]
Q = makeQfromS(s,p,4)

##r=r/r[6]
##Q = makeQ(r,p,4)

x = simData(tree,nsites,Q,alpha)
dat = convertData(x,ntax)

## writing fasta file
filename = paste0("../Data/sim-artiodactyl-6.fasta")
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
filename = paste0("../Examples/Simulations/sim-artiodactyl-6.nex")

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
