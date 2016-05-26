## Bret Larget
## May 19, 2016
## Code to do bootstrap NJ and write files out to be read by ccdprobs
## Meant to be run as a script

## Arguments expected to be sequenceFile, outFile, #bootstrapReplicates
args = commandArgs(trailingOnly = TRUE)
print(args)

require(ape)
f = function(x) return(nj(dist.dna(x)))

do.nj = function(seq.file,out.file,B)
{
  x = read.dna(seq.file,format="fasta")
  phy = f(x)
  names = phy$tip.label
  n = length(names)
  rownames(x) = 1:n
  phy$tip.label = 1:n
  out = boot.phylo(phy,x,f,B=B,trees=TRUE)
  sink(out.file)
  cat("translate\n")
  for ( i in 1:n )
  {
    cat(formatC(i,width=8))
    cat(' ')
    cat(names[i])
    if ( i < n )
      cat(',')
    else
      cat(';')
    cat("\n")
  }
  sink()
  write.tree(out$trees,file=out.file,digits=1,append=TRUE)
}

do.nj(args[1],args[2],as.integer(args[3]))