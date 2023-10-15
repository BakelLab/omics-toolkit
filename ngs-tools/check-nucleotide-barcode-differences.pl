#!/usr/bin/env perl

# 10.11.2011 09:19:43 EST
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

   Usage: $sScriptName <barcodefile>
   
   Count the minimum amount of base changes between each barcode 
   and any other barcode in a list. Used to assess whether a set of
   barcodes is distinct enough to be used for multiplexing. The input
   file should contain one barcode per line.
    
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Read the barcodes
my @asBarcodes;
open IN, $ARGV[0] or die "Error: can't open $ARGV[0]: $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   push @asBarcodes, $_;
}
close IN;

# Check differences
for (my $i=0 ; $i<@asBarcodes ; $i++){
   my $nMinDiff = -1;
   my $sMinComp = "";
   for (my $j=0 ; $j<@asBarcodes ; $j++){
      next if ($i==$j);
      my $nMinLen  = length($asBarcodes[$i]) < length($asBarcodes[$j]) ? length($asBarcodes[$i]) : length($asBarcodes[$j]);
      
      # Count mismatches between barcode pair
      my $nMismatches = 0;
      for (my $nPos=0 ; $nPos<$nMinLen ; $nPos++){
         my $nBaseI = substr($asBarcodes[$i], $nPos, 1);
         my $nBaseJ = substr($asBarcodes[$j], $nPos, 1);
         $nMismatches++ unless ($nBaseI eq $nBaseJ);
      }
      
      # Update minimum difference count
      if ($nMinDiff == -1){
         $nMinDiff = $nMismatches;
         $sMinComp = $asBarcodes[$j];
      }
      else{
         if ($nMismatches < $nMinDiff){
            $nMinDiff = $nMismatches;
            $sMinComp = $asBarcodes[$j];
         }
      }
   }
   print join("\t", $asBarcodes[$i], $nMinDiff, $sMinComp), "\n";
}
