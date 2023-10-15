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

    Takes a gff file and converts it into a bed file

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
      print_bedline(\%hFeatures) if (keys(%hFeatures));
      %hFeatures = ();
      $hFeatures{$sID}{'strand'}        = $sStrand;
      $hFeatures{$sID}{'chr'}           = $sChr;
      $hFeatures{$sID}{$sType}{'start'} = $nStart;
      $hFeatures{$sID}{$sType}{'end'}   = $nEnd;
      push @{$hFeatures{$sID}{$sType}{'features'}}, [$nStart, $nEnd];
   }   
}
print_bedline(\%hFeatures) if (keys(%hFeatures));


###############
# SUBROUTINES #
###############

# print_bedline
#
# Prints a gff feature as a bedline
sub print_bedline{
   my $rhFeatures = shift @_;
   
   foreach my $sID (keys(%$rhFeatures)){
      if (exists($rhFeatures->{$sID}{'exon'}) and exists($rhFeatures->{$sID}{'cds'})){
         my $nStart      = $rhFeatures->{$sID}{'exon'}{'start'} - 1;
         my $nThickStart = $rhFeatures->{$sID}{'cds'}{'start'} - 1;
         print join ("\t", $rhFeatures->{$sID}{'chr'}, $nStart, $rhFeatures->{$sID}{'exon'}{'end'},
                     $sID, '0', $rhFeatures->{$sID}{'strand'}, $nThickStart, $rhFeatures->{$sID}{'cds'}{'end'},
                     '0', get_blocks($rhFeatures->{$sID}{'exon'}{'start'}, $rhFeatures->{$sID}{'exon'}{'features'}));
      }
      elsif (exists($rhFeatures->{$sID}{'exon'})){
         my $nStart      = $rhFeatures->{$sID}{'exon'}{'start'} - 1;
         my $nThickStart = $nStart;
         print join ("\t", $rhFeatures->{$sID}{'chr'}, $nStart, $rhFeatures->{$sID}{'exon'}{'end'},
                     $sID, '0', $rhFeatures->{$sID}{'strand'}, $nThickStart, $rhFeatures->{$sID}{'exon'}{'end'},
                     '0', get_blocks($rhFeatures->{$sID}{'exon'}{'start'}, $rhFeatures->{$sID}{'exon'}{'features'}));
      }
      elsif (exists($rhFeatures->{$sID}{'cds'})){
         my $nStart      = $rhFeatures->{$sID}{'cds'}{'start'} - 1;
         my $nThickStart = $nStart;
         print join ("\t", $rhFeatures->{$sID}{'chr'}, $nStart, $rhFeatures->{$sID}{'cds'}{'end'},
                     $sID, '0', $rhFeatures->{$sID}{'strand'}, $nThickStart, $rhFeatures->{$sID}{'cds'}{'end'},
                     '0', get_blocks($rhFeatures->{$sID}{'cds'}{'start'}, $rhFeatures->{$sID}{'cds'}{'features'}));
      }
      else{
         die "Error: no exon or CDS information was found for feature '$sID'\n";
      }
      print "\n";
   }
}


# get_blocks
#
# Get block count, sizes and starts
sub get_blocks {
   my ($nGenomeStart, $raBlocks) = @_;
   my @aaBlocks = @$raBlocks;
   
   if (scalar(@aaBlocks)==1){
      my ($nStart,$nEnd) = @{$aaBlocks[0]};
      my $blocksize = $nEnd-$nStart+1;
      return join("\t", 1,"$blocksize,","0,");
   }
   else{
      
      # Now sort by start then end and process each block
      @aaBlocks = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @aaBlocks;
      my @anBlockStarts;
      my @anBlockSizes;
      foreach my $rBlock (@aaBlocks){
         my ($nStart, $nEnd) = @{$rBlock};
         push @anBlockStarts, $nStart-$nGenomeStart;
         push @anBlockSizes,  $nEnd-$nStart+1;
      }
      my $sBlockStarts = join(',', @anBlockStarts);
      my $sBlockSizes  = join(',', @anBlockSizes);
      return join("\t", scalar(@aaBlocks), "$sBlockSizes,","$sBlockStarts,");
   }
}
