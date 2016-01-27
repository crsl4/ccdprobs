#!/usr/bin/perl

# perl script to read a phylip file and create
# a translate table, and transform the phylip file to
# have taxon numbers instead of names.
# run: perl createTranslateTable.pl -phylip name.phy
# Claudia January 2016

use Getopt::Long;
use File::Path qw( make_path );
use strict;
use warnings;
use Carp;
use File::Basename;

my $phylip = "test.phy";

# -------------- read arguments from command-line -----------------------
GetOptions('phylip=s' => \$phylip, #you can leave this comma
    );

my($filename, $dirs, $suffix) = fileparse($phylip, qr/\.[^.]*/);

my $oufn = $dirs.$filename."_numbers.phy";
my $table = $dirs.$filename."_table.out";

open my $FHi, $phylip or die "can't open PHYLIP gene sequence file";
open my $FHo, ">", $oufn or die "can't open new PHYLIP gene sequence file";
open my $FHtab, ">", $table or die "can't open new table file";
print $FHtab " translate\n"; # need space at the beginning
my $numLine = 0;
my $nReadNames = 1;
my $ntax = 0;
while (<$FHi>){
    $numLine++;
    if($numLine < 3){
	print $FHo $_;
	if($numLine == 1){
	    my @line = split /\s+/, $_;
	    $ntax = $line[1];
	}
	next
    }
    my @line = split /\s+/, $_;
    print $FHo "$nReadNames";
    # cannot use \t and needs to be aligned
    if($nReadNames < 10){
	print $FHo "       "; #7 spaces => sequence starts in 9th
    }elsif($nReadNames < 100){
	print $FHo "      "; #6 spaces => sequence starts in 9th
    }elsif($nReadNames < 1000){
	print $FHo "     "; #5 spaces => sequence starts in 9th
    }elsif($nReadNames < 10000){
	print $FHo "    "; #4 spaces => sequence starts in 9th
    }elsif($nReadNames < 100000){
	print $FHo "   "; #3 spaces => sequence starts in 9th
    } else{
	die "cannot handle more than 100000 taxa in createTranslateTable";
    }
    print $FHo "$line[-1]\n";
    print $FHtab "       $nReadNames $line[0]";
    #print "readnames $nReadNames, ntax $ntax \n";
    if($nReadNames == $ntax){
	print $FHtab ";";
    }else{
	print $FHtab ",\n";
    }
    $nReadNames++;
}

close $FHo;
close $FHi;
close $FHtab;

