#!/usr/bin/env perl

# 14.09.2015 20:59:41 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp    = 0;
my $sStrand  = '+';
GetOptions("help!"    => \$sHelp,
           "strand:s" => \$sStrand);

# PRINT HELP
$sHelp = 1 unless(@ARGV>0);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName
   
   Arguments:
    -strand <string>
      Strand to extract from the wig file. Default: $sStrand 
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Check arguments
die "Error: unknown strand '$sStrand': $!\n" unless(($sStrand eq "+") or ($sStrand eq "-"));

my ($sChr, $nSegStart, $nSegEnd, $nSegVal, $nSpan) = ("", 0, 0, 0, 0);
open IN, $ARGV[0] or die "Error: can't open '$ARGV[0]': $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   next if (/^track/);
   s/[\n\r]+$//;
   if (/^variableStep/){
      print join("\t", $sChr, $nSegStart, $nSegEnd, $nSegVal), "\n" if ($nSegVal >0);
      ($sChr, $nSpan) = $_ =~ /variableStep chrom=(\S+) span=(\d+)/;
      ($nSegStart, $nSegEnd, $nSegVal) = (0,0,0);
   }
   else{
      my (@anLine) = split /\t/;
      die "Error: not enough wig fields on line $.\n" if ( ($sStrand eq "+") and (@anLine<2) );
      die "Error: not enough wig fields on line $.\n" if ( ($sStrand eq "-") and (@anLine<3) );
      my $nCurPos = $anLine[0];
      my $nCurVal = $sStrand eq "+" ? $anLine[1] : $anLine[2];
      if ( ($nCurPos == $nSegEnd + $nSpan) and $nCurVal == $nSegVal){
         $nSegEnd = $nCurPos;
      }
      else{
         print join("\t", $sChr, $nSegStart, $nSegEnd, $nSegVal), "\n" if ($nSegVal >0);
         $nSegStart = $nCurPos - 1;
         $nSegEnd   = $nCurPos;
         $nSegVal   = $nCurVal;
      }
   }
}
close IN;

# Print value of last segment
print join("\t", $sChr, $nSegStart, $nSegEnd, $nSegVal), "\n" if ($nSegStart);
