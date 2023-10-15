#!/usr/bin/env perl

# 26.08.2015 15:50:47 EDT
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

   Usage: $sScriptName <interproscan.tsv>
   
   Convert interproscan output into a format that can be read
   by the bioconductor TopGO package to provide a custom GO
   reference file.
   
   Arguments: 
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Parse output
my %hGeneIDs;
my %hLUT;
open IN, $ARGV[0] or die "Error: can't open $ARGV[0]: $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sAcc, $sMd5, $nLength, $sAnalysis, $sSigAcc, $sSigDesc, $nDomStart, $nDomEnd, $nScore, $sStatus, $nDate, $sIprAcc, $sIprDesc, $sGO) = split /\t/, $_, -1;
   if ($sGO){
      my @asGO = split /\|/, $sGO;
      foreach my $sGOid (@asGO){
         $hLUT{$sAcc}{$sGOid}++;
      }
   }
   $hGeneIDs{$sAcc}++;
}
close IN;

# Write topGO input
foreach my $sID ( sort(keys %hGeneIDs) ){
   if (exists $hLUT{$sID}){
      my $sGOIDs = join(", ", keys( %{$hLUT{$sID}} ) );
      print "$sID\t$sGOIDs\n";
   }
   else{
      print "$sID\t\n";
   }
}
