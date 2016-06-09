#!/usr/bin/perl

# perl script to repeat blOneRep.pl
# for many replicates
# run as perl blAllReps.pl -from 1 -to 10

use Getopt::Long;
use File::Path qw( make_path );
use strict;
use warnings;
use Carp;


# ================= parameters ======================
# Number of replicate of the whole process
my $from = 1;
my $to = 2;

# -------------- read arguments from command-line -----------------------
GetOptions( 'from=i' => \$from,
	    'to=i' => \$to,
    );

for my $i ($from..$to){
    print "REPLICATE $i =================================\n";
    system("perl blOneRep.pl -irep $i");
}
