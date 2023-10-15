#!/usr/bin/perl

# 14.01.2020 12:17:24 EST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp                = 0;
my $sBamFile             = "";
my $sLengthFile          = "";
my $nPoffset             = 0;
my $nMinUTRlength        = 0;
my $nMinTranscriptLength = 200;
my $nAverage             = 0;
GetOptions("help!"       => \$sHelp,
           "bedgraph:s"  => \$sBamFile,
           "lengths:s"   => \$sLengthFile,
           "poffset:i"   => \$nPoffset,
           "U:i"         => \$nMinUTRlength,
           "T:i"         => \$nMinTranscriptLength,
           "average:i"   => \$nAverage);

# PRINT HELP
$sHelp = 1 unless($sBamFile and $sLengthFile and $nPoffset);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName
   
   Arguments:
    -b --bam <string>
      Name of bedgraph file with 5'-read-end (init) counts per transcript.
    -l --lengths <string>
      Name of tab-delimited file with AUG position and transcript lengths
    -p --poffset <integer>
      P-site offset to add to track. Default: 0
    -U <integer>
      Minimum UTR length
    -T <integer>
      Minimum transcript length
    -a <integer>
      Average fragment length. Reads with length greater
      than the average fragment length will be corrected
      by an additional nucleotide.
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Read transcript and UTR lengths
my %hTranscriptLengths;
my %hUTRlengths;
open LENGTHS, $sLengthFile or die "Error: can't open '$sLengthFile': $!\n";
while (<LENGTHS>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sTranscriptID, $nUTRlength, $nTranscriptLength) = split /\t/;
   if ( ($nUTRlength >= $nMinUTRlength) and ($nTranscriptLength >= $nMinTranscriptLength) ){
      $hTranscriptLengths{$sTranscriptID} = $nTranscriptLength;
      $hUTRlengths{$sTranscriptID} = $nUTRlength;
   }
}
close LENGTHS;

# Gather transcript psite densities
my %hTranscriptPsites;
my %hTranscriptReadCounts;
open IN, "bamToBed -bed12 -i $sBamFile |" or die "Error: can't open bam file: $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sTranscriptID, $nStart, $nEnd, $sReadID, $nScore, $sStrand, $nTstart, $nTend, $sRGB, $nBlocks, $nBlockLengths, $nBlockStarts) = split /\t/;
   
   # Process unspliced reads mapping to the sense strand only
   if ( exists($hTranscriptLengths{$sTranscriptID}) and ($nBlocks== 1) and ($sStrand eq "+") ){
      my $nFragLength = $nBlockLengths;
      my $nFragOffset = ($nAverage and $nFragLength>$nAverage) ? $nPoffset+1 : $nPoffset;
      my $nPosition   = $nStart - $hUTRlengths{$sTranscriptID} + $nFragOffset;
      $hTranscriptPsites{$sTranscriptID}{$nPosition}++;
      $hTranscriptReadCounts{$sTranscriptID}++;
   }
}
close IN;

# Write raw and normalized p-site tracks
foreach my $sTranscriptID (sort keys %hTranscriptPsites){
   foreach my $nPosition ( sort {$a <=>$b} keys %{$hTranscriptPsites{$sTranscriptID}}){
      my $nAveragePsiteDensity = $hTranscriptReadCounts{$sTranscriptID} / $hTranscriptLengths{$sTranscriptID};
      print join("\t", $sTranscriptID, $nPosition, $hTranscriptPsites{$sTranscriptID}{$nPosition}, $hTranscriptPsites{$sTranscriptID}{$nPosition}/$nAveragePsiteDensity), "\n";
   }
}
