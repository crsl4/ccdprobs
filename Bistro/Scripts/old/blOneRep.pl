#!/usr/bin/perl

# perl script to run one replicate of bl analysis:
# - simulate fasta data for ntax=3,..,12
# for each ntax:
# - run bl on each data (time it)
# - run studyWeights.r on created logw.txt
# - rename logw.txt to current rep and ntax
# this will create a file weights.txt with one row per replicate
# warning: need to run inside BranchLength/Scripts
# modified to use bistro instead of bl,
# need to be inside Bistro/Scripts


use Getopt::Long;
use File::Path qw( make_path );
use strict;
use warnings;
use Carp;


# ================= parameters ======================
# Number of replicate of the whole process
my $irep = 1;
my $longseed = time;
my $seed = substr($longseed, -5);

# -------------- read arguments from command-line -----------------------
GetOptions( 'irep=i' => \$irep,
    );

my $logfile = "blOneRep.txt";
open FHlog, ">> $logfile";
print FHlog "============================================\n";
print FHlog "blOneRep for replicate $irep\n";
close FHlog;
system("date >> $logfile");
system("hostname >> $logfile");


# 1) Simulate the data for each ntax
open FHlog, ">> $logfile";
print FHlog "running simulateDataStudy.r\n";
my $Rcmd = "Rscript simulateDataStudy.r $irep $seed";
print FHlog "$Rcmd\n";
close FHlog;
system($Rcmd);

# 2) read settings from settings$irep.txt
my $setFile = "../Examples/Simulations/settings$irep.txt";
open FHlog, ">> $logfile";
print FHlog "reading file $setFile\n";
my @trees = `grep "tree =" $setFile`;
my @rates = `grep "rates =" $setFile`;
my @probs = `grep "probs =" $setFile`;

my @splittree;
my @splitrates;
my @splitprobs;
my $ntaxtree;
my $ntaxrates;
my $ntaxprobs;
my $blcmd;
my $start;
my $duration;
my $newfile;

# one replicate will run analysis for each ntaxa
foreach my $ntax (3,4,5,6,7,8,9,10,11,12){
    print FHlog "-----------------ntaxa = $ntax------------------\n";
    print "ntaxa = $ntax, ";
    # the actual tree, rates and probs is the third element after splitting by space
    # maybe there is a better way to do this step
    @splittree = split / /, $trees[$ntax-3];
    @splitrates = split / /, $rates[$ntax-3];
    @splitprobs = split / /, $probs[$ntax-3];
    $ntaxtree = $splittree[2];
    $ntaxrates = $splitrates[2];
    $ntaxprobs = $splitprobs[2];
    chomp $ntaxtree;
    chomp $ntaxrates;
    chomp $ntaxprobs;

    $seed = $seed +($ntax-3)*int(rand(1000));

    $blcmd = "../Code/bistro/bistro -f ../Data/sim$ntax.fasta -q $ntaxrates -p $ntaxprobs -s $seed -b 10000 -r 1000 >& out$ntax";
    print "running bistro\n";
    print FHlog "$blcmd\n";
    $start = time;
    system($blcmd);
    $duration = time - $start;

    $Rcmd = "Rscript studyWeights.r $ntax $irep $duration $seed";
    print FHlog "$Rcmd\n";
    system($Rcmd);

    $newfile = "logw_${ntax}_${irep}.txt";
    `mv out.txt $newfile`;
}

close FHlog;
