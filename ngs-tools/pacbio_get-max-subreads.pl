#!/usr/bin/env perl

# 24.07.2012 18:21:57 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
GetOptions("help!"   => \$sHelp);

# PRINT HELP
$sHelp = 1 unless(@ARGV>0);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP
   
   Filter longest subreads from a PacBio filtered_subreads file.
   
   Usage: $sScriptName <filtered_subreads.fasta>
    
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

my @asIDs;
open IN, $ARGV[0] or die "Error: can't open '$ARGV[0]': $!\n";
while (<IN>){
   if(/^>/){
      s/[\n\r]+$//;
      s/^>//;
      my ($sSMRTid, $nZMW, $sSegment) = split /\//, $_, -1;
      die "Error: no subread segment found in fasta header on line $.\n" unless ($sSegment);
      my ($nStart, $nEnd) = split /_/, $sSegment;
      die "Error: segment end missing in fasta header on line $.\n" unless ($nEnd);
      push @asIDs, [($sSMRTid, $nZMW, $sSegment, $nEnd-$nStart)];
   }
}
close IN;

@asIDs = sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $b->[3] <=> $a->[3]} @asIDs;
my $sLastID = "";
foreach my $rSegment (@asIDs){
   if ("$rSegment->[0]/$rSegment->[1]" ne $sLastID){
      print "$rSegment->[0]/$rSegment->[1]/$rSegment->[2]\n";
   }
   $sLastID = "$rSegment->[0]/$rSegment->[1]";
}

