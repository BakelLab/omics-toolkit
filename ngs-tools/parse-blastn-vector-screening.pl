#!/usr/bin/env perl

# 11.10.2011 20:24:28 EDT
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

   Usage: $sScriptName <filterin>
    
   Parse blastn output of a set of sequences against UniVec.
   The blastn output must have an additional column added at
   the end with the length of the query sequence
    
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

my $nBoundarySize = 100;
my $nBoundaryEval = 300;
my $nMiddleEval   = 1;

my %hCoordinates;
open IN, $ARGV[0] or die "Error: can't open $ARGV[0]: $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sQid, $sSid, $nIdent, $nAlnLen, $nMM, $nGO, $nQstart, $nQend, $nSstart, $nSend, $nE, $nB, $nQlen) = split /\t/;
   die "Error: query length missing on line $., did you add the extra query length column?\n" unless ($nQlen =~ /^\d+$/);
   
   if ($nQlen >= 2*$nBoundarySize){
      # Process the hits
      if ( ($nQstart < $nBoundarySize) and ($nQend < $nBoundarySize)){  # First boundary
         if ($nE <= $nBoundaryEval){
            my $nNewStart = $nQend > $nQstart ? $nQend : $nQstart;
            if (exists $hCoordinates{$sQid}){
               $hCoordinates{$sQid}{start} = $nNewStart if($nNewStart > $hCoordinates{$sQid}{start});
            }
            else{
               $hCoordinates{$sQid}{start} = $nNewStart;
               $hCoordinates{$sQid}{end}   = $nQlen;
               $hCoordinates{$sQid}{size}  = $nQlen;
            }
         }
      }
      elsif ( (($nQlen - $nQstart) < $nBoundarySize) and (($nQlen - $nQend) < $nBoundarySize)){  # Last boundary
         if ($nE <= $nBoundaryEval){
            my $nNewEnd   = $nQstart < $nQend ? $nQstart : $nQend;
            if (exists $hCoordinates{$sQid}){
               $hCoordinates{$sQid}{end} = $nNewEnd if($nNewEnd < $hCoordinates{$sQid}{end});
            }
            else{
               $hCoordinates{$sQid}{end}   = $nNewEnd;
               $hCoordinates{$sQid}{start} = 0;
               $hCoordinates{$sQid}{size}  = $nQlen;
            }
         }
      }
      else{  # Somewhere in the middle
         # Only consider when E-value < 1
         if ($nE < $nMiddleEval){
            warn("$sQid has a significant match outside boundary region\n");
         }
      }
   }
   else{
      warn("Skipping $sQid, too small\n");
   }
   
}
close IN;

# Print the new coordinates
my $nTotalCut = 0;
foreach my $sQid (keys %hCoordinates){
   print join("\t", $sQid, $hCoordinates{$sQid}{start}, $hCoordinates{$sQid}{end}, $sQid), "\n";
   $nTotalCut += $hCoordinates{$sQid}{start};
   $nTotalCut += $hCoordinates{$sQid}{size} - $hCoordinates{$sQid}{end};
}
warn ("Removed $nTotalCut bases of contamination\n");

