#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use ValidateFiletypes qw(check_gff);

# GLOBALS
$ENV{SORT}        = 'sort';
$ENV{SORT_BUFFER} = '2G';

# ARGUMENTS
my $sHelp       = 0;
GetOptions("help!"      => \$sHelp);
my $sGffFile = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sGffFile);
if ($sHelp) {
    my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
    die <<HELP

    $sScriptName <gff-file>

    Calculate the processed/spliced length of each feature in a gff file
    Also returns outermost genomic coordinates of each feature.

    -help
      This help message
      
HELP
}


# Check the gff file format
my @asGffErrors = check_gff($sGffFile);
if (@asGffErrors){
   unshift @asGffErrors, "The following errors were found in the input gff file:";
   die join("\n", @asGffErrors), "\n";
}


# Sort the gff file by name and process the entries
my %hFeatures;
open GFF, "$ENV{SORT} -S $ENV{SORT_BUFFER} -s -t '\t' -k9,9 $sGffFile |" or die "Error sorting gff file: $!\n";
while (<GFF>){
   next if /^\s*$/;
   next if /^\s*[Tt]rack/;
   next if /^\s*#/;
   s/[\n\r]$//g;
   my ($sChr, $sSource, $sType, $nStart, $nEnd, $nScore, $sStrand, $nFrame, $sID) = split /\t/;
   $sType = lc($sType);
   if (exists($hFeatures{$sID})){
      die "Error: feature '$sID' was found on multiple chromosomes\n" unless($hFeatures{$sID}{'chr'} eq $sChr);
      die "Error: feature '$sID' was found on multiple strands\n"     unless($hFeatures{$sID}{'strand'} eq $sStrand);
      $hFeatures{$sID}{$sType}{'start'} ||= $nStart;
      $hFeatures{$sID}{$sType}{'end'}   ||= $nEnd;
      $hFeatures{$sID}{$sType}{'start'} = $nStart if ($nStart < $hFeatures{$sID}{$sType}{'start'});
      $hFeatures{$sID}{$sType}{'end'}   = $nEnd   if ($nEnd   > $hFeatures{$sID}{$sType}{'end'});
      push @{$hFeatures{$sID}{$sType}{'features'}}, [$nStart, $nEnd];
   }
   else{
      print_length(\%hFeatures) if (keys(%hFeatures));
      %hFeatures = ();
      $hFeatures{$sID}{'strand'}        = $sStrand;
      $hFeatures{$sID}{'chr'}           = $sChr;
      $hFeatures{$sID}{$sType}{'start'} = $nStart;
      $hFeatures{$sID}{$sType}{'end'}   = $nEnd;
      push @{$hFeatures{$sID}{$sType}{'features'}}, [$nStart, $nEnd];
   }   
}
print_length(\%hFeatures) if (keys(%hFeatures));


###############
# SUBROUTINES #
###############

# print_length
#
# Prints a gff feature as a bedline
sub print_length{
   my $rhFeatures = shift @_;
   
   foreach my $sID (keys(%$rhFeatures)){
      if (exists($rhFeatures->{$sID}{'exon'})){
         my ($nStart, $nEnd, $nLength) = get_length($rhFeatures->{$sID}{'exon'}{'features'});
         print join ("\t", $rhFeatures->{$sID}{'chr'}, $nStart, $nEnd, $sID, $nLength), "\n";
      }
      else{
         warn "Warning: no exon information was found for feature '$sID', skipping\n";
      }
      
   }
}

# get_blocks
#
# Get block count, sizes and starts
sub get_length {
   my ($raBlocks) = @_;
   my @aaBlocks = @$raBlocks;
   
   # Now sort by start then end and process each block
   my $nLength;
   my $nLastStart = 0;
   my $nLastEnd   = 0;
   @aaBlocks = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @aaBlocks;
   foreach my $rBlock (@aaBlocks){
      my ($nStart, $nEnd) = @{$rBlock};
      if ($nLastEnd){
         if ($nStart > $nLastEnd){
            $nLength += $nLastEnd-$nLastStart+1;
            ($nLastStart,$nLastEnd) = ($nStart, $nEnd);
         }
         else{
            $nLastEnd = $nEnd if ($nEnd > $nLastEnd);
         }
      }
      else{
         ($nLastStart,$nLastEnd) = ($nStart, $nEnd);
      }
   }
   $nLength += $nLastEnd-$nLastStart+1; # still need to add the last exon block here
   return($aaBlocks[0][0], $nLastEnd, $nLength);
}
