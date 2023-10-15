#!/usr/bin/env perl

use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use File::Basename;
use ValidateFiletypes qw(check_bed);

# GLOBALS
$ENV{SORT}          ||= 'sort';     # Unix sort binary

# ARGUMENTS
my $sHelp             = 0;
my $sFeaturesFile     = '';
my $sRegionFile       = '';
GetOptions("help!"            => \$sHelp,
           "features-file:s"  => \$sFeaturesFile,
           "region-file:s"    => \$sRegionFile);
           
# PRINT HELP
my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
$sHelp = 1 unless($sFeaturesFile and $sRegionFile);
if ($sHelp) {
    die <<HELP

    $sScriptName -features-file <file> -region-file <file>

   Shuffle feature coordinates within a specified set of regions (e.g. whole
   chromosomes or intergenic regions). Feature size is maintained and shuffling
   occurs randomly between and within chromosomes. NOTE: only the first three
   fields of a bed file are shuffled, strand and relative exon positions
   are maintained!

    Required:
    -features-file <string>
      Bed file (full or basic) with genomic features for shuffling
    -region-file <string>
      Bed file (full or basic) with genomic regions to sample from
    
    -help
      This help message
      
HELP
}

#########
# START #
#########

# Quick input checking
die "Error: features file '$sFeaturesFile' does not exist\n" unless (-e $sFeaturesFile);
die "Error: region file '$sRegionFile' does not exist\n" unless (-e $sRegionFile);

# Check the input bed file
my @asBedErrors = check_bed($sFeaturesFile, 4);
if (@asBedErrors){
   unshift @asBedErrors, "The following errors were found in the input bed file:";
   die join("\n", @asBedErrors), "\n";
}

# Read the region file
my $rRegions = read_regions($sRegionFile);

# Shuffle and print regions
shuffle_and_print_features($rRegions, $sFeaturesFile);


###############
# SUBROUTINES #
###############

# read_regions
#
# Read the region file into an array and keep track of the cumulative
# size
sub read_regions {
   my $sRegionFile = shift @_;
   my @aaRegions;
   my $nCumulativeSize = 0;
   
   open REGIONS, $sRegionFile or die "Error: can't open region file: $!\n";
   while (<REGIONS>){
      next if /^\s*$/;
      next if /^\s*#/;
      next if /^\s*[Tt]rack/;
      s/[\n\r]$//g;
      my ($sChr, $nStart, $nEnd) = (split /\t/)[0,1,2];
      die "Error: insufficient number of fields in '$sRegionFile', line $.\n"  unless ($nEnd);
      die "Error: start position must be numeric in '$sRegionFile', line $.\n" unless ($nStart =~ /^\d+$/);
      die "Error: end position must be numeric in '$sRegionFile', line $.\n"   unless ($nEnd =~ /^\d+$/);
      die "Error: end < start in '$sRegionFile', line $.\n"                    unless ($nEnd >= $nStart);
      my $nSegmentStart = $nCumulativeSize + 1;
      $nCumulativeSize += $nEnd - $nStart + 1;
      push @aaRegions, [$nSegmentStart, $nCumulativeSize, $sChr, $nStart];
   }
   close REGIONS;
   
   return \@aaRegions;
}


# shuffle_and_print_features
#
# Shuffle features within the supplied areas
sub shuffle_and_print_features{
   my ($rRegions, $sFeaturesFile) = @_;

   # Get the maximum region size for the random shuffle
   my @aaRegions  = @$rRegions;
   my $nMaxLength = $aaRegions[$#aaRegions][1];

   # Get an array of range-end positions for the region list
   my @aRegionEnds;
   foreach my $rRegion (@aaRegions){
      push @aRegionEnds, $rRegion->[1];
   }

   # Open the region file and start shuffling
   open IN, $sFeaturesFile or die "Error: can't open feature file for shuffling: $!\n";
   while (<IN>){
      next if /^\s*$/;
      next if /^\s*#/;
      next if /^\s*[Tt]rack/;
      s/[\n\r]$//g;
      my ($sFtChr, $nFtStart, $nFtEnd, @asRest) = split /\t/;
      die "Error: insufficient number of fields in '$sFeaturesFile', line $.\n"  unless ($nFtEnd);
      die "Error: start position must be numeric in '$sFeaturesFile', line $.\n" unless ($nFtStart =~ /^\d+$/);
      die "Error: end position must be numeric in '$sFeaturesFile', line $.\n"   unless ($nFtEnd =~ /^\d+$/);
      die "Error: end < start in '$sFeaturesFile', line $.\n"                    unless ($nFtEnd >= $nFtStart);
      my $nFtSize = $nFtEnd - $nFtStart + 1;
      
      # Pick a random position in a region that will accomodate the full feature
      # Make sure to give up after a certain number of tries, otherwise features that
      # are too large to fit in any regions will basically hang the script forever
      my ($sNewChr, $nNewStart, $nNewEnd) = ('',0,0);
      my $nCounter = 1000;
      while ($nCounter > 0){
         my $nShufStart = int(rand($nMaxLength));
         my $nShufEnd   = $nShufStart + $nFtSize - 1;
         
         # Get the array element ID in which this shuffle falls
         my $nRegionID  = bsearch_range_of_ends($nShufStart, \@aRegionEnds);
         my ($nCmlStart, $nCmlEnd, $sRegChr, $nRegStart) = @{$aaRegions[$nRegionID]};
         
         if ( ($nShufStart > $nCmlStart) and ($nShufEnd < $nCmlEnd) ){  # We found a region that fits
            my $nOffset = $nShufStart - $nCmlStart;
            $sNewChr    = $sRegChr;
            $nNewStart  = $nRegStart + $nOffset;
            $nNewEnd    = $nRegStart + $nOffset + $nFtSize - 1;
            $nCounter   = 0;
            last;
         }

         $nCounter--;
      }
      
      # Check if we could successfully reshuffle
      die "Error: could not reassign feature on line $. to a random region position\n" unless ($sNewChr);
      
      # Make sure to also update the thick start and thick end, maintaining the proper offset relative to Ft start and end
      if (@asRest >=5){
         $asRest[3] = $nNewStart + ($asRest[3] - $nFtStart);
         $asRest[4] = $nNewEnd - ($nFtEnd - $asRest[4]);
      }
      
      # Print the reshuffled position and continue with the next feature
      print join("\t", $sNewChr, $nNewStart, $nNewEnd, @asRest), "\n";
   }
   close IN;
}



# bsearch_range_of_ends
#
# Given a value and an array of sorted range ends, this function will 
# return the array element ID corresponding to the range in which the
# value falls
sub bsearch_range_of_ends {
   my ($x, $a) = @_;         # search for x in array a
   my ($l, $u) = (1, scalar(@$a)-1);   # lower, upper end of search interval
   my $i;                    # index of probe
   while ($l <= $u) {
      $i = int(($l + $u)/2);
      if ( ($a->[$i] < $x) and ($a->[$i-1] < $x) ) {
         $l = $i+1;
      }
      elsif ( ($a->[$i] > $x) and ($a->[$i-1] > $x) ) {
         $u = $i-1;
      } 
      else {
         return $i; # found
      }
    }
    return -1; # not found
}


