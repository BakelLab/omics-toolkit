#!/usr/bin/env perl

# 08.09.2010 11:03:11 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
my $sInput       = '';
my $nEcut        = 0.0001;
my $nLengthCut   = 110;
GetOptions("help!"     => \$sHelp,
           "input:s"   => \$sInput,
           "ethresh:n" => \$nEcut,
           "lthresh:i" => \$nLengthCut);

# PRINT HELP
$sHelp = 1 unless($sInput);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName
    
    -i --input <string>
      Blast hit input file (tabular format)
    -e --ethresh <number>
      E-value threshold. Default: $nEcut
    -l --lthresh <integer>
      Length threshold. Default: $nLengthCut
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Check args
die "Error: length cutoff must be greater than zero\n"  unless ($nLengthCut > 0);
die "Error: e-value cutoff must be greater than zero\n" unless ($nEcut > 0);

# Process the input file
my %hHits;
open IN, "<$sInput" or die "Error: can't open '$sInput': $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sQID, $sSID, $nPctIdent, $nAlignLength, $nMismatches, $nGaps, $nQstart, $nQend, $nSstart, $nSend, $nEval, $nBit) = split /\t/;
   
   # Set alignment length
   $hHits{$sQID}{$sSID}{length} += $nAlignLength;
   
   # Set eval
   if (exists $hHits{$sQID}{$sSID}{eval}){
      $hHits{$sQID}{$sSID}{eval} = $nEval if ($nEval < $hHits{$sQID}{$sSID}{eval});
   }
   else{
      $hHits{$sQID}{$sSID}{eval} = $nEval;
   }
   
   # Set orientation
   my $nQor = $nQstart < $nQend ? 0 : 1;
   my $nSor = $nSstart < $nSend ? 0 : 1;
   my $nRevSeq = $nQor == $nSor ? 0 : 1;
   if (exists $hHits{$sQID}{$sSID}{revseq}){
      warn("Warning: hit '$sSID' matches '$sQID' in multiple orientations\n") unless ($hHits{$sQID}{$sSID}{revseq} == $nRevSeq);
   }
   else{
      $hHits{$sQID}{$sSID}{revseq} = $nRevSeq;
   }
}
close IN;

# Process the hits
foreach my $sQID (sort keys %hHits){
   foreach my $sSID (sort {$a <=> $b} keys %{$hHits{$sQID}}){
      if ( ($hHits{$sQID}{$sSID}{length} >= $nLengthCut) and ($hHits{$sQID}{$sSID}{eval} <= $nEcut) ){
         print join("\t", $sQID, $sSID, $hHits{$sQID}{$sSID}{revseq}, $hHits{$sQID}{$sSID}{length}, $hHits{$sQID}{$sSID}{eval}), "\n";
      }
   }
}
