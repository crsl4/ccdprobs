# Code to do bootstrap NJ and write files out to be read by ccdprobs
require(ape)
f = function(x) return(nj(dist.dna(x)))

do.nj = function(seq.file="",B=10000,out.file="")
{
  x = read.dna(seq.file,format="fasta")
  phy = f(x)
  out = boot.phylo(phy,x,f,B=B)
  write.nexus(out,file=out.file)
}