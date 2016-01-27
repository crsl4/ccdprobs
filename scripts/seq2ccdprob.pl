#!/usr/bin/perl

# INPUT: phylip or nexus sequence data
# RUN: perl seq2ccdprob.pl -phylip file.phy (or -nexus file.nex)
# WARNING: assumes it is in the same path with all the scripts

# 0) create translate table and get phylip file with taxon numbers
# 1) run in R bootstrap2NJ.r to get list of trees
#  careful here to give input and output files to R, and to check that output file is empty
# 2) create input file for mbsum by pasting the translate table to the trees.out
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

my $nexus = "none";
my $phylip = "none";
my $cmd;

# -------------- read arguments from command-line -----------------------
GetOptions('nexus=s' => \$nexus, #you can leave this comma
	   'phylip=s' => \$phylip, #you can leave this comma
    );

# -------------- converting to phylip ----------------------------------
if($phylip eq "none" && $nexus eq "none"){
    die "you need to input either a nexus or a phylip file";
} elsif($phylip eq "none"){ #entered nexus file
    print "converting nexus file to phylip...\n";
    $cmd = "perl nexus2phylip.pl -nexus $nexus";
    print "$cmd\n";
    system($cmd);
    my($filename, $dirs, $suffix) = fileparse($nexus, qr/\.[^.]*/);
    $phylip = $dirs.$filename.".phy";
    print "created phylip file: $phylip\n";
}

# -------------- translate table and taxon numbers ----------------------------------
print "creating translate table...\n";
$cmd = "perl createTranslateTable.pl -phylip $phylip";
print "$cmd\n";
system($cmd);

my($filename, $dirs, $suffix) = fileparse($phylip, qr/\.[^.]*/);
my $newphylip = $dirs.$filename."_numbers.phy"; #names need to match those in createTranslateTable.pl
my $table = $dirs.$filename."_table.out";
print "created new phylip file $newphylip and translate table $table\n";


my $Rout = $dirs.$filename."_listtrees.out";
my $mbsumIN = $dirs.$filename."_trees.out";


if(-f $Rout){
    print "warning: file $Rout exists and will be overwritten\n";
    unlink("$Rout");
}

# -------------- bootstrap + NJ in R ----------------------------------
print "running bootstrap + NJ...\n";
$cmd = "Rscript --vanilla bootstrap2NJ.r $newphylip $Rout";
print "$cmd\n";
system($cmd);
print "created list of bootstrap trees: $Rout\n";

# -------------- create mbsum input file ----------------------------------

open my $FHi, $Rout or die "can't open list of trees file";
open my $FHtab, $table or die "can't open translate table";
open my $FHo, ">", $mbsumIN or die "can't open new list of trees file";
print $FHo "begin trees;\n"; # need this line at the beginning
while (<$FHtab>){
    print $FHo "$_";
}
close $FHtab;
print $FHo "\n";
while (<$FHi>){
    print $FHo "$_";
}
close $FHi;

print "created mbsum input file: $mbsumIN\n";

# -------------- run mbsum ----------------------------------
$cmd = "../src/mbsum $mbsumIN";
print "$cmd\n";
#system($cmd);
my $mbsumOUT = $dirs.$filename."_trees.in"; 
print "created mbsum output file: $mbsumOUT\n";

# -------------- run mbsum ----------------------------------
my $ccdprobsOUT = $dirs.$filename."_ccdprobs";
$cmd = "../src/ccdprobs --in $mbsumOUT --out $ccdprobsOUT";
print "$cmd\n";
#system($cmd);

# ccdprobs

# - need to modify mbsum to accept two types of file: from bucky as it
# - is (with tree name =), from list of trees with artificial
# - translate table on top run mbsum in M336_trees.out and check what
# - file name it produces (update seq2ccdprobs) run ccdprobs in the
# - output of mbsum

# - finish seq2ccdprobs
