library(ape)
library(phangorn)
tree = read.tree(text="(1,(2,((((((3,13),((((8,9),11),12),10)),(6,7)),(4,5)),14),(15,(((((16,17),((18,19),(20,21))),31),(23,((24,(25,26)),(27,(28,(29,30)))))),22)))));")
dat = read.FASTA("../../Data/whales.fasta")
pdat = phyDat(dat)
parsimony(tree,pdat)

