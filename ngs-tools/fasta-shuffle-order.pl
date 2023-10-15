#!/usr/bin/env perl

# 07.09.2010 20:00:48 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use List::Util 'shuffle';

# GET PARAMETERS
my $sHelp         = 0;
GetOptions("help!"       => \$sHelp);
my $sInput = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sInput);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName <fasta-file>

   Randomly shuffle the order of sequences in a multi-fasta file.
    
   Options:
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Process the INPUT
my @asFasta;
my $sFastaHeader = '';
my $sFastaSeq    = '';
open INPUT, "<$sInput" or die "Error: can't read the fasta file\n";
while (<INPUT>){
   if (/^>/ or eof){
      if (eof){
         die "Error: file ends in fasta header without sequence\n" if (/^>/);
         $sFastaSeq .= $_;
      }
      push @asFasta, [($sFastaHeader, $sFastaSeq)] if ($sFastaHeader);
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

# Print shuffled array
foreach my $sSequence (shuffle(@asFasta)){
   print @$sSequence;
}

