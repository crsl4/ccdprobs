## r script to simulate data
## based on simulatingData.r
## run Rscript simulatingData-paper.r 024 500 1234
## where 024 is the dataset, 500 is the number of sites, 1234 is the seed
## Claudia September 2018

## if we want to use consensus tree from MrBayes:
## sed 's/\[[^]]*\]//g' artiodactylMBTree.txt > cleanMBTree.txt

source("old/genetrees.R")
source("old/branch-length_lik.R")
source("old/transformRates.r")

datafolder = "/Users/Clauberry/Documents/github/CladeCondProb/ccdprobs/Bistro/Data/datasets/"

# test if there is at least one argument: if not, return an error
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  data = "024"
  nsites = 500
  seed = 1234
} else if(length(args) == 3){
  data = args[1]
  nsites = strtoi(args[2])
  seed = strtoi(args[3])
} else {
    stop("need three arguments: data,nsites,seed")
}

set.seed(seed)
fileRoot = paste0("sim-",data,"-nsites-",nsites)

## Reading tree file
treefile = paste0(datafolder,"ExaBayes_ConsensusExtendedMajorityRuleNewick.",data,".cons")
tree = read.tree(file=treefile)
ntax = length(tree$tip.label)


## Reading parameters file
parsfile = paste0(datafolder,"ExaBayes_parameterStatistics.",data,".pars")
pars = read.table(parsfile, header=TRUE)
r = pars$mean[c(7,3,12,9,2,1)]
p = pars$mean[c(13,8,5,6)]
s = getS(c(r,p))
s = s[1:6]
##alpha = pars$mean[4]
alpha = 1000
Q = makeQfromS(s,p,4)


## Simulating data
x = simData(tree,nsites,Q,alpha)
dat = convertData(x,ntax)

## Writing fasta file
filename = paste0(datafolder,fileRoot,".fasta")
for(i in 1:ntax){
    if(i == 1){
        write(paste0(">",tree$tip.label[i]), file=filename)
    } else{
        write(paste0(">",tree$tip.label[i]), file=filename, append=TRUE)
    }
    write(paste0(dat[[i]],collapse=""), file=filename, append=TRUE)
}


## Writing nex file
filename2 = paste0(datafolder,fileRoot,".nex")
write("#NEXUS", file=filename2)
write("", file=filename2,append=TRUE)
write("begin taxa;", file=filename2, append=TRUE)
write(paste0("dimensions ntax=",ntax,";"), file=filename2, append=TRUE)
write("taxlabels", file=filename2, append=TRUE)
for(i in 1:ntax){
    write(tree$tip.label[i],file=filename2,append=TRUE)
}
write(";", file=filename2, append=TRUE)
write("end;", file=filename2, append=TRUE)
write("", file=filename2,append=TRUE)

write("begin characters;", file=filename2, append=TRUE)
write(paste0("dimensions nchar=",nsites,";"), file=filename2, append=TRUE)
write("format datatype=dna gap=-;", file=filename2, append=TRUE)
write("matrix", file=filename2, append=TRUE)
for(i in 1:ntax){
    seq = paste0(dat[[i]], collapse="")
    write(paste0(sprintf("%-20s",paste0(tree$tip.label[i])),seq), file=filename2, append=TRUE)
}
write(";", file=filename2, append=TRUE)
write("end;", file=filename2, append=TRUE)
write("begin mrbayes;", file=filename2, append=TRUE)
write("      set autoclose=yes nowarn=yes;", file=filename2, append=TRUE)
write("      lset nst=6;", file=filename2, append=TRUE)
write("      outgroup 1;", file=filename2, append=TRUE)
write("", file=filename2, append=TRUE)
write("      mcmc ngen=1100000 printfreq=10000 samplefreq=50;", file=filename2, append=TRUE)
write("      sumt burnin=100000;", file=filename2, append=TRUE)
write("      sump burnin=100000;", file=filename2, append=TRUE)
write("end;", file=filename2, append=TRUE)
