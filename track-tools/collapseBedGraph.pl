#!/usr/bin/env perl

# 03.06.2013 11:14:15 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp         = 0;
my $nGap          = 0;
GetOptions("help!"   => \$sHelp,
           "gap:i"   => \$nGap);
my $sBedGraphFile = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sBedGraphFile);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [-g] <bedGraphFile>
   
   Merges consecutive regions in a bedgraph file together if they
   are within a specified distance. Values in merged regions are averaged.
   The default setting will only merge directly adjacent regions.
   
   Arguments:
    -gap <integer>
      Maximum gap between adjacent regions for merging. Default: 0
    -help
      This help message
   
HELP
}

##########
## MAIN ##
##########

# Open the bed file, make sure that start>end position and pass through sort prior to perl processing
my ($sLastChr, $nLastStart, $nLastEnd, $nRegionSum) = ('',0,0,0);
open INPUT, $sBedGraphFile or die "Can't open $sBedGraphFile: $!\n";
while (<INPUT>){
   next if /^\s*$/;
   next if /^\s*#/;
   s/[\n\r]//g;
   my ($sChr,$nStart,$nEnd,$nVal) = split /\t/, $_, -1;
   die "Error: data is not sorted by position\n" if ($sLastChr eq $sChr and $nStart<$nLastStart);
   if ( ($sChr ne $sLastChr) or ($nStart > $nLastEnd+$nGap) or eof){
      if (eof){
         $nRegionSum += ($nEnd - $nStart) * $nVal;
         $nLastEnd    = $nEnd;
      }
      if ($sLastChr){
         my $nRegionVal  = $nRegionSum / ($nLastEnd - $nLastStart);
         print join("\t", $sLastChr, $nLastStart, $nLastEnd, $nRegionVal, $nRegionSum, $nLastEnd - $nLastStart), "\n";
         
      }
      $sLastChr   = $sChr;
      $nLastStart = $nStart;
      $nLastEnd   = $nEnd;
      $nRegionSum = ($nEnd - $nStart) * $nVal;
   }
   else{
      $nRegionSum += ($nEnd - $nStart) * $nVal;
      $nLastEnd = $nEnd;
   }
}
close INPUT;
