#!/usr/bin/env perl

# 20.01.2012 13:03:00 EST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp         = 0;
my $nBasesPerLine = 100;
GetOptions("help!"     => \$sHelp,
           "length:i"  => \$nBasesPerLine);
my $sInput = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sInput);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [options] <fasta-file>
   
   Reflow sequences in a fasta file so that each sequence
   line is set to the specified length.
    
   Options:
    -l --length <integer>
      Number of bases per line in the fasta file
      Default: $nBasesPerLine
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Process the INPUT
my $sFastaHeader = '';
my $sFastaSeq    = '';
my $flRevSeq     = 0;
open INPUT, "<$sInput" or die "Error: can't read the fasta file\n";
while (<INPUT>){
   if (/^>/ or eof){
      if (eof){
         die "Error: file ends in fasta header without sequence\n" if (/^>/);
         $sFastaSeq .= $_;
      }
      if ($sFastaHeader){
         $sFastaSeq  =~ s/\s//g;
         $sFastaSeq  =~ s/[\n\r]+//g;
         $sFastaSeq =~ s/.{$nBasesPerLine}/$&\n/sg;
         $sFastaSeq =~ s/\n+$//;
         $sFastaSeq .= "\n";
         print join("", $sFastaHeader, $sFastaSeq);
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
