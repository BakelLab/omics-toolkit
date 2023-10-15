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
   
   Compress interproscan results into a single line of interpro and GO 
   annotations to facilitate unique file joins.
   
   Arguments: 
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Parse output
my %hLUT;
open IN, $ARGV[0] or die "Error: can't open $ARGV[0]: $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sAcc, $sMd5, $nLength, $sAnalysis, $sSigAcc, $sSigDesc, $nDomStart, $nDomEnd, $nScore, $sStatus, $nDate, $sIprAcc, $sIprDesc, $sGO) = split /\t/, $_, -1;
   if ($sIprAcc and $sIprDesc){
      $hLUT{$sAcc}{IPR}{$sIprAcc} = $sIprDesc;
   }
   if ($sGO){
      my @asGO = split /\|/, $sGO;
      foreach my $sGOid (@asGO){
         $hLUT{$sAcc}{GO}{$sGOid}++;
      }
   }
}
close IN;

# Write condensed annotations
print "#Accession\tInterPro\tGO\n";
foreach my $sAcc (sort(keys %hLUT)){
   print "$sAcc";
   
   # Print interpro matches
   my @asInterpro;
   if (exists $hLUT{$sAcc}{IPR}){
      foreach my $sIprAcc (keys %{$hLUT{$sAcc}{IPR}}){
         push @asInterpro, join(":", $sIprAcc, $hLUT{$sAcc}{IPR}{$sIprAcc});
      }
      print "\t", join("|", @asInterpro);
   }
   else{
      print "\t";
   }
   
   # Print GO matches
   if (exists $hLUT{$sAcc}{GO}){
      print "\t", join("|", sort(keys %{$hLUT{$sAcc}{GO}}));
   }
   else{
      print "\t";
   }
   print "\n";
}
