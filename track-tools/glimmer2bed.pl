#!/usr/bin/env perl

# 09.10.2012 11:32:43 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
my $sGlimmerFile = "";
my $sInfoFile  = "";
GetOptions("help!"     => \$sHelp,
           "glimmer:s" => \$sGlimmerFile,
           "info:s"    => \$sInfoFile);

# PRINT HELP
$sHelp = 1 unless($sGlimmerFile and $sInfoFile);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Convert a glimmer gene prediction file to bed format

   Usage: $sScriptName -g <glimmer-predictions.txt> -i <reference-lengths>
    
   Options
    -g <string>
      Name of file with glimmer predictions
    -i <string>
      Name of file with the lengths of each sequence with
      glimmer predictions
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Read reference sequence lengths
my %hLengths;
open INFO, $sInfoFile or die "Error: can't open '$sInfoFile': \n";
while (<INFO>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sID, $nLength, @asRest) = split /\t/;
   $hLengths{$sID} = $nLength;
}
close INFO;


# Process glimmer output
my $sContigID = "";
open IN, $sGlimmerFile or die "Error: can't open '$sGlimmerFile': $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   if (/^>(.*)$/){
      $sContigID = $1;
   }
   else{
      if (exists $hLengths{$sContigID}){
         my ($sOrfID, $nGlStart, $nGlEnd, $nFrame, $nScore, @asRest) = split /\s+/;
         $nFrame =~ s/^\+//;
         ($nGlStart, $nGlEnd) = sort {$a <=> $b} ($nGlStart, $nGlEnd);
         $nGlStart--; # Glimmer coordinates are 1-based
         my $nBedStart  = $nGlStart <0 ? 0 : $nGlStart;
         my $nBedEnd    = $nGlEnd > $hLengths{$sContigID} ? $hLengths{$sContigID} : $nGlEnd;
         my $sBedStrand = $nFrame < 0 ? "-" : "+";
         my $sGlStrand  = $nFrame < 0 ? 1 : 0;
         my $sBedID     = join("_", $sContigID, "$nGlStart-$nGlEnd", $sGlStrand);
         print join("\t", $sContigID, $nBedStart, $nBedEnd, $sBedID, $nScore, $sBedStrand), "\n";
      }
      else{
         die "Error: reference length missing for contig '$sContigID'\n";
      }
   }
}
close IN;
