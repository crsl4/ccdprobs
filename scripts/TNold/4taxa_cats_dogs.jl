# julia script to sample branch lengths in a given 4-taxon topology
# we use the cats-dogs dataset to test
# Claudia February 2016

using PhyloNetworks


# input: 4-taxon topology
tree="(A,(B,(C,D)));"
top = readTopologyLevel1(tree)

# need read phylip
