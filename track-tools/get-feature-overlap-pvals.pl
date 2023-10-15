#!/usr/bin/env perl

# 01.03.2010 13:30:33 EST
# Harm van Bakel <hvbakel@gmail.com>

# GLOBALS
$ENV{TMPDIR}         ||= "/tmp";     # location for tmp file storage
$ENV{FJOIN}          ||= "fjoin";    # Fjoin binary
$SIG{'INT'}=$SIG{'HUP'}=$SIG{'ABRT'}=$SIG{'QUIT'}=$SIG{'TRAP'}=$SIG{'STOP'}=\&INTERRUPT;

# MODULES
use strict;
use warnings;
use Getopt::Long;
use File::Temp qw(tempfile tempdir);
use ValidateFiletypes qw(check_gff);
use Carp qw(croak);

# GET PARAMETERS
my $sHelp            = 0;
my $sTestFile        = '';
my $sReferenceFile   = '';
my $sRegionFile      = '';
my $nIterations      = 1000;
my $flStrandSpecific = 0;
my $flShuffleStrands = 0;
GetOptions("help!"            => \$sHelp,
           "test:s"           => \$sTestFile,
           "reference:s"      => \$sReferenceFile,
           "regions:s"        => \$sRegionFile,
           "iterations:i"     => \$nIterations,
           "strand-specific!" => \$flStrandSpecific,
           "shuffle-strands!" => \$flShuffleStrands);

# PRINT HELP
$sHelp = 1 unless($sTestFile and $sRegionFile and $sReferenceFile);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

    $sScriptName -features-file <file> -region-file <file>

   Calculate pvalues for the overlap between two feature sets

    Required:
    -test <string>
      Gff file with features to test
    -reference <string>
      Gff file with features to compare to
    -regions <string>
      Bed file (full or basic) with genomic regions to sample from
    -iterations <integer>
      Number of iterations. Default: $nIterations
    -strand-specific
      Take strand into account when doing overlaps
    -shuffle-strands
      Do strand shuffling in the test file (default is position only)
    -help
      This help message
      
HELP
}


##########
## MAIN ##
##########

# Quick input checking
die "Error: test feature file '$sTestFile' does not exist\n" unless (-e $sTestFile);
die "Error: reference feature file '$sReferenceFile' does not exist\n" unless (-e $sReferenceFile);
die "Error: region file '$sRegionFile' does not exist\n" unless (-e $sRegionFile);

# Check the input files
for my $sFile ($sTestFile, $sReferenceFile){
   my @asGffErrors = check_gff($sFile);
   if (@asGffErrors){
      unshift @asGffErrors, "The following errors were found in '$sFile':";
      die join("\n", @asGffErrors), "\n";
   }
}

# Read the region file
my $rRegions  = read_regions($sRegionFile);

# Write the reference to a temporary file
my ($sTmpRef, $nRefLines)   = gff_to_tempfile($sReferenceFile);

# Get the observed overlap count
my ($sTmpTest, $nTestLines) = gff_to_tempfile($sTestFile);
my $nObsCount = get_feature_overlap_count(test=>$sTmpTest, reference=>$sTmpRef, strandspecific=>$flStrandSpecific);
print join("\t", "observed", $nTestLines, $nObsCount), "\n";

# Start shuffling the test file and count overlaps
for (my $i=1 ; $i<=$nIterations ; $i++){
   my ($sTmpShuf, $nFeatCount) = shuffle_and_print_features(test=>$sTmpTest, regions=>$rRegions, shufflestrands=>$flShuffleStrands);
   my $nExpCount = get_feature_overlap_count(test=>$sTmpShuf, reference=>$sTmpRef, strandspecific=>$flStrandSpecific);
   unlink($sTmpShuf) or warn("Warning: could not remove temporary file '$sTmpShuf'\n");
   print join("\t", $i, $nFeatCount, $nExpCount), "\n";
}


#################
## SUBROUTINES ##
#################


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



# gff_to_tempfile
#
# Convert an input gff file into a temporary file for fjoin comparisons
# Tab-delimited temp file contains: chr, start, end, strand, id
sub gff_to_tempfile {
   my ($sFile) = @_;
   my @aaOutput;
   my $nLineCount = 0;
   open IN, $sFile or croak("Error: can't open '$sFile'");
   while (<IN>){
      next if (/^\s*$/);
      next if (/^ *#/);
      next if /^\s*[Tt]rack/;
      s/[\n\r]+$//;
      my @asLine = (split /\t/)[0,3,4,6,8];
      push @aaOutput, [@asLine];
      $nLineCount++;
   }
   close IN;
   
   # Print the output, sorted by start position (for fjoin)
   my ($fhTmp, $sTmp) = tempfile('getmatrix-XXXXX', DIR=>$ENV{TMPDIR}, UNLINK=>1);
   foreach my $rLine (sort {$a->[1] <=> $b->[1]} @aaOutput){
      $fhTmp->print(join("\t", @$rLine), "\n");
   }
   $fhTmp->close();
   return ($sTmp, $nLineCount);
}



# get_feature_overlap_count
#
# Count the number of test features overlapping reference features
sub get_feature_overlap_count {
   my %args = (test           => "",
               reference      => "",
               strandspecific => 1,
               @_);
   
   my %hIDs;
   my $sCols = $args{strandspecific} ? '1,2,3,4' : '1,2,3';
   open FJOIN, "$ENV{FJOIN} --columns1=$sCols --columns2=$sCols -1 $args{test} -2 $args{reference} 2>/dev/null |";
   while (<FJOIN>){
      next if (/^\s*$/);
      next if (/^ *#/);
      s/[\n\r]+$//;
      my $sID = (split /\t/)[5];
      $hIDs{$sID} ='' unless(exists $hIDs{$sID});
   }
   close FJOIN;
   return scalar(keys %hIDs);
}



# shuffle_and_print_features
#
# Shuffle features within the supplied areas
sub shuffle_and_print_features{
   my %args = (test           => '',
               regions        => '',
               shufflestrands => 1,
               @_);

   # Get the maximum region size for the random shuffle
   my @aaRegions  = @{$args{regions}};
   my $nMaxLength = $aaRegions[$#aaRegions][1];

   # Get an array of range-end positions for the region list
   my @aRegionEnds;
   foreach my $rRegion (@aaRegions){
      push @aRegionEnds, $rRegion->[1];
   }

   # Open the region file and start shuffling
   my $nShuffleCount;
   my @aaOutput;
   open IN, $args{test} or die "Error: can't open feature file for shuffling: $!\n";
   while (<IN>){
      next if /^\s*$/;
      next if /^\s*#/;
      next if /^\s*[Tt]rack/;
      s/[\n\r]$//g;
      my ($sFtChr, $nFtStart, $nFtEnd, $sFtStrand, $sFtID) = split /\t/;
      my $nFtSize = $nFtEnd - $nFtStart + 1;
      
      # Pick a random position in a region that will accomodate the full feature
      # Make sure to give up after a certain number of tries, otherwise features that
      # are too large to fit in any regions will basically hang the script forever
      my ($sNewChr, $nNewStart, $nNewEnd) = ('',0,0);
      my $nCounter = 500;
      while ($nCounter > 0){
         my $nShufStart = int(rand($nMaxLength));
         my $nShufEnd   = $nShufStart + $nFtSize - 1;
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
      if ($sNewChr){
         my $sNewStrand = $sFtStrand;
         if($args{shufflestrands}){
            my $nRand = int(rand(2));
            $sNewStrand = $nRand ? '+' : '-';
         }
         push @aaOutput, [$sNewChr, $nNewStart, $nNewEnd, $sNewStrand, $sFtID];
         $nShuffleCount++;
      }
   }
   close IN;
   
   
   # Sort the shuffle file
   my ($fhTmp, $sTmp) = tempfile('getmatrix-XXXXX', DIR=>$ENV{TMPDIR}, UNLINK=>1);
   foreach my $rLine (sort {$a->[1] <=> $b->[1]} @aaOutput){
      $fhTmp->print(join("\t", @$rLine), "\n");
   }
   $fhTmp->close();

   return ($sTmp, $nShuffleCount);
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



# INTERRUPT
#
# Interrupt routine, make sure we exit gracefully for tmp file cleanup
sub INTERRUPT{
   exit(1); # This will call END
}

