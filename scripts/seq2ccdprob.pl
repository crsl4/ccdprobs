#!/usr/bin/perl

# INPUT: phylip or nexus sequence data
# RUN: perl seq2ccdprob.pl -phylip file.phy (or -nexus file.nex)
# 0) create translate table and get phylip file with taxon numbers
# 1) run in R bootstrap2NJ.r to get list of trees
#  careful here to give input and output files to R, and to check that output file is empty
# 2) create input file for mbsum by pasting the translate table to the trees.out
# modify mbsum so that it will accept format from bucky (with tree name =) and without that (just list of trees after translate table)
# 3) run mbsum: ../src/mbsum trees.out
# (inside datasets). creates trees.in
# check that trees.in has the desired format for ccdprobs (it should)
# 4) run ccdprobs: ../src/ccdprobs --in trees.in --out trees
# (inside datasets)
# Claudia January 2016

use Getopt::Long;
use File::Path qw( make_path );
use strict;
use warnings;
use Carp;
use File::Basename;



