#!/usr/bin/perl

# perl script to convert nexus file to phylip file
# run: perl nexus2phylip.pl -nexus name.nex
# Claudia January 2016

use Getopt::Long;
use File::Path qw( make_path );
use strict;
use warnings;
use Carp;
use File::Basename;

#-----------------------------------------------#
#  convert nexus to phylip files                #
#  interleaved format: do *not* repeat taxon names
#-----------------------------------------------#

my $nexus = "test.nex";

# -------------- read arguments from command-line -----------------------
GetOptions('nexus=s' => \$nexus, #you can leave this comma
    );

my($filename, $dirs, $suffix) = fileparse($nexus, qr/\.[^.]*/);

my $oufn = $dirs.$filename.".phy";

my $read = 0;
my $removeNames = 0; my $nReadNames = 0;
my $ntax = 0;
my $nchar = 0;
my $ncopy = 0;
open my $FHi, $nexus or die "can't open NEXUS gene sequence file";
open my $FHo, ">", $oufn or die "can't open PHYLIP gene sequence file";
while (<$FHi>){
    if ($read){
	# if (/^\s;$/){ last; } # for old data
	if (/^\s*;/){
	    print "found last in $nReadNames";
	    last;
	} # old data: /^;/
	#if (/copy/){ $ncopy++; next; }
	else {
	    #if ($removeNames){
		#print "inside remove names";
		#if (/^[^\s]+\s+(.*)/) { print $FHo "$1\n"; }
	    #} else {
	    print "$_ \n";
	    print $FHo $_;
	    $nReadNames++;
	    #if ($nReadNames == $ntax){ $removeNames = 1; }
	}
    #}
    }
    if ($read==0){
	if (/ntax=(\d+)/i){
	    $ntax = $1;
	    print $FHo " $ntax ";
	}
	if (/nchar=(\d+)/i){
	    $nchar = $1;
	    if ($ntax==0){ print "problem in gene: found nchar before ntax\n"}
	    print $FHo "$nchar\n";
	}
    }
    if (/^\s*matrix/i){
	$read=1;
	if ($ntax==0 or $nchar==0){
	    print "problem in gene: was unable to find ntax ($ntax) or nchar ($nchar)\n";
	}
    }
}
close $FHi;
close $FHo;

