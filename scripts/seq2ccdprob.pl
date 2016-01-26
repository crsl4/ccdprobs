#!/usr/bin/perl

# INPUT: phylip or nexus sequence data
# RUN: perl seq2ccdprob.pl -phylip file.phy (or -nexus file.nex)
# 1) run in R bootstrap2NJ.r to get list of trees from
# sequence data: trees.out
# 2) run mbsum: ../src/mbsum trees.out
# (inside datasets). creates trees.in
# 3) run ccdprobs: ../src/ccdprobs --in trees.in --out trees
# (inside datasets)
# Claudia January 2016

use Getopt::Long;
use File::Path qw( make_path );
use strict;
use warnings;
use Carp;
use File::Basename;



