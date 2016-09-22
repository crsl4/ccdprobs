# Code to generate random gene trees from a species tree (not general)
# Bret Larget
# March 21-23, 2012
# April 26-27, 2012     ---  add hybridization simulation
#                       ---  make data generation functions

library(ape)

############################################################
#
# Functions to simulate coalecent trees
#
############################################################

### coal.edge ###############################################################
# input:
#   nodes = array of node indices
#   heights = array of lengths (at edge beginning) for each node of edge so far (0 for tips)
#   len = length of population edge (NA if ancestral)
#   nextNode = number of next node at coalescent event (decreasing so last at root is one more than number of tips)
#   edge = matrix of edges for the entire tree, passed along to each edge in preorder traversal
#   edge.length = array of edge lengths corresponding to edge
#
# return value:
#   list of:
#   nodes
#   heights
#   nextNode
#   edge
#   edge.length
#

coal.edge = function(nodes,heights,len,nextNode,edge,edge.length) {
  n = length(nodes)
  rate = n*(n-1)/2
  x = len # remaining edge length
  
  while ( n > 1 && (is.na(x) || x > 0) ) {
    z = rexp(1,rate=rate)
    if ( !(is.na(x)) && (z > x) ) {
      heights = heights+x
      x = 0
    }
    else {
      if ( !is.na(x) )
        x = x-z
      ab = sample(n,size=2)
      heights = heights+z
      edge = rbind(c(nextNode,nodes[ab[1]]),edge)
      edge.length = c(heights[ab[1]],edge.length)
      edge = rbind(c(nextNode,nodes[ab[2]]),edge)
      edge.length = c(heights[ab[2]],edge.length)
      nodes = c(nodes[-ab],nextNode)
      heights = c(heights[-ab],0)
      n = n-1
      rate = n*(n-1)/2
      nextNode = nextNode - 1
    }
  }

  if ( !is.na(x) )
    heights = heights+x
  
  return( list(nodes=nodes,heights=heights, nextNode=nextNode, edge=edge, edge.length=edge.length) )
}

################################################################################

### coal.tree ############################################################
#
# input: a population tree (species tree)
#
# output: a gene tree simulated from the coalescent (one individual per species)
#

coal.tree = function(tree) {
  tree = reorder(tree,"pruningwise")
  n = length(tree$tip.label)
  nextNode = 2*n-1
  edge = matrix(NA,0,2)
  edge.length = numeric()
  nodes.list = list()
  for ( i in 1:nrow(tree$edge) ) {
    child = tree$edge[i,2]
    if ( child <= n ) { # this is a terminal edge
      out = coal.edge(nodes=child,heights=0,len=tree$edge.length[i],nextNode=nextNode,edge=edge,edge.length=edge.length)
    }
    else { # this is an internal edge, children already processed
      offspring = which( tree$edge[,1] == child )
      nodes = numeric()
      heights = numeric()
      for ( k in offspring ) {
        nodes = c(nodes,nodes.list[[k]]$nodes)
        heights = c(heights,nodes.list[[k]]$heights)
      }
      out = coal.edge(nodes=nodes,heights=heights,len=tree$edge.length[i],nextNode=nextNode,edge=edge,edge.length=edge.length)
    }      
    nextNode = out$nextNode
    edge = out$edge
    edge.length = out$edge.length
    nodes.list[[i]] = list(nodes=out$nodes,heights=out$heights)
  }
  # now do the ancestral edge
  offspring = which( tree$edge[,1] == n+1 )
  nodes = numeric()
  heights = numeric()
  for ( k in offspring ) {
    nodes = c(nodes,nodes.list[[k]]$nodes)
    heights = c(heights,nodes.list[[k]]$heights)
  }
  out = coal.edge(nodes=nodes,heights=heights,len=NA,nextNode=nextNode,edge=edge,edge.length=edge.length)
  genetree = list(edge=out$edge,tip.label=tree$tip.label,Nnode=n-1,edge.length=out$edge.length)
  class(genetree) = "phylo"
  return( genetree )
}

################################################################################

### function to generate ils gene trees
generate.trees = function(tree,nTrees,TOTAL=0.2) {
  genetrees = list()
  for ( i in 1:nTrees )
    genetrees[[i]] = coal.tree(tree)
    ### rescale gene trees from coalescent units to substitution units
    ### here, divide each tree by the average height and then multiply by TOTAL rate (try 0.2)
  simtrees = genetrees
  ht = rep(0,nTrees)

  for ( i in 1:nTrees )
    ht[i] = as.numeric(branching.times(genetrees[[i]])[1])
  ht.mean = mean(ht)
  for ( i in 1:nTrees )
    simtrees[[i]]$edge.length = simtrees[[i]]$edge.length / ht.mean * TOTAL

  return(list(genetrees=genetrees,simtrees=simtrees,ht=ht))
}

################################################################################
#
# Functions to simulate sequence data on trees
#
################################################################################

rdirichlet = function(n,a) {
  len = length(a)
  x = matrix(rgamma(len * n, a), ncol = len, byrow = TRUE)
  sm = apply(x,1,sum)
  return( x/sm )
}

############################################################

makeQ = function(p,r,full=T) {
  p = p/sum(p)           # normalize to sum to 1
  r = r/r[6]             # normalize so r[6] = 1
  R = matrix(0,4,4)      # initialize R=0 in all entries
  R[row(R) > col(R)] = r # set R for the lower triangle
  R = R + t(R)           # set R for upper triangle
  Q = R %*% diag(p)      # Q[i,j] = R[i,j] * p[j] for i != j
  diag(Q) = -1 * apply(Q,1,sum) # set diag(Q) to the negative row sum of the off-diagonal
  mu = -1 * sum(p*diag(Q))      # Compute expected number of substitutions per unit time
  Q = Q/mu                      # Rescale
  if(!full)
    return(Q)
  e = eigen(Q)
  U = e$vectors
  Uinv = solve(U)
  v = e$values
  return(list(Q=Q,p=p,U=U,v=v,Uinv=Uinv))
}

# TN93() makes a Tamura-Nei 1993 Q matrix

TN93 = function(p,kappaR,kappaY,full=T) {
  Q = makeQ(p,c(1,kappaR,1,1,kappaY,1),full=full)
  return(Q)
}

# matrixExp(x,a) computes a matrix exponential and returns the matrix exp(a*x)
#   where a is a scalar and x is a matrix.
# If full=T, then the input x is a list that includes a decomposition of the matrix.
# Otherwise, the decomposition is performed and then the matrix exponential is calculated.
# If x is a Q-matrix and a is a time, the result is the probability transition matrix for an edge of length a.
# More efficient algorithms are possible, but this is good enough for the assignment.

matrixExp = function(x,a,full=T) {
  if(full)
    return( x$U %*% diag(exp(a*x$v)) %*% x$Uinv )
  e = eigen(x)
  out = e$vectors %*% diag(exp(a*e$vectors)) %*% solve(e$vectors)
  if ( any(out < 0 ) ) {
    print(x)
    print(a)
    stop("No negative probabilities!!!")
  }
  return( out )
}

### function to simulate entire data set for a genetree by recursing from root to leaves
### internally, A=1, C=2, G=3, T=4

simData = function(tree,nSites,Q,alpha) {
  ### reorder tree cladewise to ensure that parents precede children
  ###  tree = reorder(tree,"clade")
  ### create local variables for less dereferencing in loops
  edge = tree$edge
  edge.length = tree$edge.length
  ### stationary distribution
  p = Q$p
  ### number of taxa
  n = length(tree$tip.label)
  ### total number of nodes
  nNodes = n + tree$Nnode
  ### number of edges
  nEdges = nrow(tree$edge)
  ### root is defined to have label one more than number of tips
  root = n+1
  ### matrix for all data, including internal nodes
  x = matrix(0,nNodes,nSites)
  ### simulate root data sequence
  x[root,] = sample(1:4,size=nSites,prob=p,replace=T)
  ### simulate rate for each site, protect against rates too close to zero
  rate = rgamma(nSites,alpha,alpha) + .Machine$double.eps^0.5
  ### loop through edges, loop through sites, generate random sequences
  for ( i in 1:nEdges ) {
    parent = edge[i,1]
    child = edge[i,2]
    for ( k in 1:nSites ) {
      P = matrixExp(Q,tree$edge.length[i]*rate[k])
      prob = P[x[parent,k],]
      if ( any(prob < 0) ) {
        print(prob)
        print(paste("i =",i))
        print(paste("rate =",rate[k]))
        print(Q)
        stop("Bad!!!")
      }
      x[child,k] = sample(1:4,size=1,prob=prob)
    }
  }
  return(x)
}

### function to change matrix of integer data to list of characters
convertData = function(x,n) {
  ### n is the number of tips and must be the first n rows of the matrix
  x = x[1:n,]
  ### out is an empty list which gets character arrays added
  out = list()
  ### nSites is number of sites
  nSites = ncol(x)
  ### loop through rows of x
  for ( i in 1:n ) {
    foo = rep("N",nSites)
    foo[ x[i,]==1 ] = "A"
    foo[ x[i,]==2 ] = "C"
    foo[ x[i,]==3 ] = "G"
    foo[ x[i,]==4 ] = "T"
    out[[i]] = foo
  }
  ### return list of sequences
  return( out )
}

################################################################################

# incomplete lineage sorting
sim.1 = function() {
  tree1.text = "(((((1:0.3,2:0.3):0.3,3:0.6):4.4,4:5.0):4.0,((5:1.0,6:1.0):4.0,7:5.0):4.0):4.0,8:13.0);"
  tree1 = read.tree(text=tree1.text)
  nTrees = 100
### generate gene trees for ils
  trees = generate.trees(tree1,nTrees)
### define file prefix for all output
  file.root = "ils-"
### Write genetrees and simtrees to file
  file.genetrees = paste(file.root,"genetrees",sep="")
  file.simtrees = paste(file.root,"simtrees",sep="")
  digits = 8 # so, 6 after decimal when 2 before; grr!
  write.tree(trees$genetrees[[1]],file.genetrees,append=F,digits=digits)
  write.tree(trees$simtrees[[1]],file.simtrees,append=F,digits=digits)
  for ( i in 2:nTrees ) {
    write.tree(trees$genetrees[[i]],file.genetrees,append=T,digits=digits)
    write.tree(trees$simtrees[[i]],file.simtrees,append=T,digits=digits)
  }

# Compute stationary distribution for each simtree in a single giant matrix with a row for each tree
# Use rdirichlet() with scale = rep(10,4) [ 95% of pi_i between 0.13 and 0.40 ]

### all stationary distributions ###
  p = rdirichlet(nTrees,rep(10,4))
# all r values
  r = rdirichlet(nTrees,c(1,3,1,1,1,4)*10) #[ random parameters for GTR model with transitions more frequent than transversions with high probability]
### all Q matrices
  Q.list = list()
  for ( i in 1:nTrees )
    Q.list[[i]] = makeQ(p[i,],r[i,])

### all gene lengths (3 * random negative binomial with mean 100, sd 33 )
  nSites = rnbinom(nTrees,size=10,mu=100)*3

### choose alpha parameters for rate variation among sites for each gene, uniform from 0.1 to 0.5
  alpha = runif(nTrees,0.1,0.5)

### Write p, r, alpha, and nSites to file
  file.params = paste(file.root,"params",sep="")
  params = data.frame(p,r,alpha,nSites)
  names(params) = c("pi.a","pi.c","pi.g","pi.t","r.ac","r.ag","r.at","r.cg","r.ct","r.gt","alpha","nSites")
  write.table(params,file.params,row.names=F)

### Loop over trees to produce data sets

  for ( i in 1:nTrees ) {
    n = length(trees$simtrees[[i]]$tip.label) # number of taxa
    dna = simData(trees$simtrees[[i]],nSites[i],Q.list[[i]],alpha[i])
    out = convertData(dna,n)
    if ( i < 10 )
      foo = "00"
    else if ( i < 100 )
      foo = "0"
    else
      foo = ""
    write.nexus.data(out,paste(file.root,foo,i,".nex",sep=""),interleaved=F)
  }

  ### Return list of output
  return( list(species.tree=tree1,trees=trees,nTrees=nTrees,p=p,r=r,alpha=alpha,nSites=nSites,Q.list=Q.list) )
}

################################################################################
# hybridization simulation

sim.2 = function() {
  tree2a.text = "(((((1:1.0,2:1.0):2.0,3:3.0):2.0,4:5.0):4.0,((5:1.0,6:1.0):4.0,7:5.0):4.0):4.0,8:13.0);"
  tree2a = read.tree(text=tree2a.text)
  tree2b.text = "(((((1:1.0,2:1.0):2.0,3:3.0):2.0,(4:1.0,5:1.0):4.0):4.0,(6:5.0,7:5.0):4.0):4.0,8:13.0);"
  tree2b = read.tree(text=tree2b.text)

### generate gene trees for hybridization
  nTrees=100
  nTreesA = rbinom(1,nTrees,0.5)
  nTreesB = nTrees - nTreesA
  hybridA.trees = generate.trees(tree2a,nTreesA)
  hybridB.trees = generate.trees(tree2b,nTreesB)
  trees = list(genetrees=c(hybridA.trees$genetrees,hybridB.trees$genetrees),simtrees=c(hybridA.trees$simtrees,hybridB.trees$simtrees),ht=c(hybridA.trees$ht,hybridB.trees$ht))

### define file prefix for all output
  file.root = "hyb-"
### Write genetrees and simtrees to file
  file.genetrees = paste(file.root,"genetrees",sep="")
  file.simtrees = paste(file.root,"simtrees",sep="")
  digits = 8 # so, 6 after decimal when 2 before; grr!
  write.tree(trees$genetrees[[1]],file.genetrees,append=F,digits=digits)
  write.tree(trees$simtrees[[1]],file.simtrees,append=F,digits=digits)
  for ( i in 2:nTrees ) {
    write.tree(trees$genetrees[[i]],file.genetrees,append=T,digits=digits)
    write.tree(trees$simtrees[[i]],file.simtrees,append=T,digits=digits)
  }

# Compute stationary distribution for each simtree in a single giant matrix with a row for each tree
# Use rdirichlet() with scale = rep(10,4) [ 95% of pi_i between 0.13 and 0.40 ]

### all stationary distributions ###
  p = rdirichlet(nTrees,rep(10,4))
# all r values
  r = rdirichlet(nTrees,c(1,3,1,1,1,4)*10) #[ random parameters for GTR model with transitions more frequent than transversions with high probability]
### all Q matrices
  Q.list = list()
  for ( i in 1:nTrees )
    Q.list[[i]] = makeQ(p[i,],r[i,])

### all gene lengths (3 * random negative binomial with mean 100, sd 33 )
  nSites = rnbinom(nTrees,size=10,mu=100)*3

### choose alpha parameters for rate variation among sites for each gene, uniform from 0.1 to 0.5
  alpha = runif(nTrees,0.1,0.5)

### Write p, r, alpha, and nSites to file
  file.params = paste(file.root,"params",sep="")
  params = data.frame(p,r,alpha,nSites)
  names(params) = c("pi.a","pi.c","pi.g","pi.t","r.ac","r.ag","r.at","r.cg","r.ct","r.gt","alpha","nSites")
  write.table(params,file.params,row.names=F)

### Loop over trees to produce data sets

  for ( i in 1:nTrees ) {
    n = length(trees$simtrees[[i]]$tip.label) # number of taxa
    dna = simData(trees$simtrees[[i]],nSites[i],Q.list[[i]],alpha[i])
    out = convertData(dna,n)
    if ( i < 10 )
      foo = "00"
    else if ( i < 100 )
      foo = "0"
    else
      foo = ""
    write.nexus.data(out,paste(file.root,foo,i,".nex",sep=""),interleaved=F)
  }

  ### Return list of output
  return( list(hybridA.tree=tree2a,hybridB.tree=tree2b,trees=trees,nTrees=nTrees,nA=nTreesA,nB=nTreesB,p=p,r=r,alpha=alpha,nSites=nSites,Q.list=Q.list) )
}

################################################################################
    
  
  






