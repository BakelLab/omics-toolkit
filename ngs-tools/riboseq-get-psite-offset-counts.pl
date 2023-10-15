#!/usr/bin/perl

# 12.12.2019 12:20:00 EST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp                = 0;
my $sBedFile             = "";
my $sLengthFile          = "";
my $nPoffset             = 0;
my $nMinUTRlength        = 20;
my $nMinTranscriptLength = 200;
my $nPsiteBandwidth      = 3;
my $nMaxCDSlength        = 100;
GetOptions("help!"       => \$sHelp,
           "bedgraph:s"  => \$sBedFile,
           "lengths:s"   => \$sLengthFile,
           "poffset:i"   => \$nPoffset,
           "U:i"         => \$nMinUTRlength,
           "T:i"         => \$nMinTranscriptLength);

# PRINT HELP
$sHelp = 1 unless($sBedFile and $sLengthFile);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName
   
   Arguments:
    -b --bedgraph <string>
      Name of bedgraph file with 5'-read-end (init) counts per transcript.
    -l --lengths <string>
      Name of tab-delimited file with AUG position and transcript lengths
    -p --poffset <integer>
      P-site offset to add to track. Default: 0
    -U <integer>
      Minimum UTR length
    -T <integer>
      Minimum transcript length
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

# Read bed file and convert to offset track
my %hOut;
open BEDGRAPH, $sBedFile or die "Error: can't open '$sBedFile': $!\n";
while (<BEDGRAPH>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sTranscriptID, $nStart, $nEnd, $nCount) = split /\t/;
   if (exists $hTranscriptLengths{$sTranscriptID}){
      if ( ($hUTRlengths{$sTranscriptID} + $nPsiteBandwidth + $nMaxCDSlength) < $hTranscriptLengths{$sTranscriptID} ){
         my $nCorrPos = $nStart - $hUTRlengths{$sTranscriptID} - $nPoffset;
      
         if ($nCorrPos < -$nPsiteBandwidth ){
            $hOut{$sTranscriptID}{utr} += $nCount;
         }
         elsif ( ($nCorrPos >= -$nPsiteBandwidth) and ($nCorrPos <= $nPsiteBandwidth) ){
            $hOut{$sTranscriptID}{psite} += $nCount;
         }
         elsif ( ($nCorrPos > $nPsiteBandwidth) and ($nCorrPos <= ($nPsiteBandwidth + $nMaxCDSlength)) ){
            $hOut{$sTranscriptID}{cds} += $nCount;
         }
         else{
            # Skip the rest
         }
      }
   }
}
close BEDGRAPH;

# Print output
foreach my $sTranscriptID (keys %hOut){
   my $nUcount = exists $hOut{$sTranscriptID}{utr} ? $hOut{$sTranscriptID}{utr} : 0;
   my $nPcount = exists $hOut{$sTranscriptID}{psite} ? $hOut{$sTranscriptID}{psite} : 0;
   my $nCcount = exists $hOut{$sTranscriptID}{cds} ? $hOut{$sTranscriptID}{cds} : 0;
   my $nUfpkm = ($nUcount / ($hUTRlengths{$sTranscriptID} - $nPsiteBandwidth + 1)) * 1000;
   my $nPfpkm = ($nPcount / ((2 * $nPsiteBandwidth) + 1)) * 1000;
   my $nCfpkm = ($nCcount / $nMaxCDSlength) * 1000;
   my $nUratio = $nPfpkm > 0 ? $nUfpkm/$nPfpkm : "NA";
   my $nCratio = $nPfpkm > 0 ? $nCfpkm/$nPfpkm : "NA";
   print join("\t", $sTranscriptID, $nUfpkm, $nPfpkm, $nCfpkm, $nUratio, $nCratio), "\n";
}
