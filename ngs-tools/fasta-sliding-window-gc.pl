#!/usr/bin/env perl

# 20.01.2012 13:03:00 EST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp     = 0;
my $nWinSize  = 200;
my $nStepSize = 50;
GetOptions("help!"     => \$sHelp,
           "window:i"  => \$nWinSize,
           "step:i"    => \$nStepSize);
my $sInput = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sInput);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [options] <fasta-file>
  
   Calculate the GC content of sequences within a multi-fasta 
   file across a sliding window applied to each sequence. 
    
   Options:
    -w --window <integer>
      Size of sliding window
      Default: $nWinSize
    -s --step <integer>
      Step size for sliding window
      Default: $nStepSize
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Process the INPUT
my %hData;
my %hBins;
my $sFastaHeader = '';
my $sFastaSeq    = '';
open INPUT, "<$sInput" or die "Error: can't read the fasta file\n";
while (<INPUT>){
   if (/^>/ or eof){
      if (eof){
         die "Error: file ends in fasta header without sequence\n" if (/^>/);
         $sFastaSeq .= $_;
      }
      if ($sFastaHeader){
         $sFastaSeq    =~ s/\s//g;
         $sFastaSeq    =~ s/[\n\r]+//g;
         $sFastaHeader =~ s/[\n\r]+//g;
         for (my $i=0 ; $i <= length($sFastaSeq)-$nWinSize ; $i+=$nStepSize) { 
            my $sSegment   = substr($sFastaSeq,$i,$nWinSize);
            my $nGCcount   = $sSegment =~ tr/GCgc/GCgc/;
            my $nGCcontent = 100 * $nGCcount/length($sSegment);
            $hData{$sFastaHeader}{$i+1} = $nGCcontent;
            $hBins{$i+1}++;
         }
      }
      $sFastaHeader = $_;
      $sFastaSeq    = "";
   }
   else{
      next if (/^\s*$/);
      next if (/^ *#/);
      $sFastaSeq .= $_ if ($sFastaHeader);
   }
}
close INPUT;

# Print output
my @anBins = sort {$a<=>$b} keys(%hBins);
print join("\t", "#ID", @anBins), "\n";
foreach my $sRead (keys %hData){
   print $sRead;
   foreach my $nBin (@anBins){
      if (exists $hData{$sRead}{$nBin}){
         print "\t$hData{$sRead}{$nBin}";
      }
      else{
         print "\tNA";
      }
   }
   print "\n";
}



