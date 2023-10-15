#!/usr/bin/env perl

# 15.12.2015 16:49:24 EST
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

   Usage: $sScriptName <bed-file>

   Extract regions with uniform coverage that are flanked
   by empty space on either side.

   Arugments:    
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Get singleton regions
my ($sLastChr, $nLastEnd, $sLastRegion) = ("", 0, "");
open IN, $ARGV[0] or die "Error: can't open $ARGV[0]: $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sChr, $nStart, $nEnd, $nCov) = split /\t/;
   if ( ($sChr eq $sLastChr) and ($nLastEnd == $nStart) ){
      $sLastRegion = "";
      $nLastEnd    = $nEnd;
   }
   else{
      print $sLastRegion, "\n" if ($sLastRegion);
      $sLastRegion = join("\t", $sChr, $nStart, $nEnd, $nCov);
      $nLastEnd    = $nEnd;
      $sLastChr    = $sChr;
   }
}
close IN;
