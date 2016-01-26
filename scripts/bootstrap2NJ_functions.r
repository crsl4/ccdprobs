# functions for the step sequences -> bootstrap -> NJ trees
# needs ape package, and x = read.dna(phylip)

# returns bootstrap datasets
bootstrapDNA <- function(x,B=100){
    # function based on ape boot.phylo
    boot.samp <- vector("list", B)
    y <- nc <- ncol(x)
    for (i in 1:B){
        index <- unlist(sample(y, replace = TRUE))
        boot.samp[[i]] <- x[, index]
    }
    ans <- boot.samp
    ans
}

# returns list of bootstrap trees obtained with nj
bootstrapNJTree <- function(x,model="JC69",B=100){
    # function based on ape boot.phylo
    boot.tree <- vector("list", B)
    y <- nc <- ncol(x)
    for (i in 1:B){
        boot.samp <- unlist(sample(y, replace = TRUE))
        boot.tree[[i]] <- nj(dist.dna(x[, boot.samp],model=model))
    }
    ans <- boot.tree
    ans
}
