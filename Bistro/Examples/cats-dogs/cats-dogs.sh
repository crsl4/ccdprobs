#!/bin/bash

# Run R to do neighbor-joining on the cats-dogs data

echo 'Running R script'
Rscript ../../Code/R/nj.R ../../Data/cats-dogs.fasta cats-dogs.tre 100
