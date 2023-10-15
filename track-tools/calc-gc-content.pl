#!/usr/bin/env perl

# 16.06.2010 16:29:01 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $nWinSize     = 25;
my $nStepSize    = 4;
my $sHelp        = 0;
GetOptions("help!"      => \$sHelp,
           "stepsize=i" => \$nStepSize,
           "winsize=i"  => \$nWinSize);
my $sSeqFile = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sSeqFile);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [-w <windowsize> -s <stepsize>] <fasta-file>
    
    -window <integer>
      Window size in bp. Default: $nWinSize
    -step <integer>
      Step size in bp. Default: $nStepSize
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Check arguments
die "Error: window size must be an integer\n" unless ($nWinSize  =~ /^[1-9]\d*$/);
die "Error: step size must be an integer\n"   unless ($nStepSize =~ /^[1-9]\d*$/);

# Work the magic
my $rhSeq      = getSeq($sSeqFile);
foreach my $sKey (sort keys(%$rhSeq)){
   my $sSeq       = $rhSeq->{$sKey};
   die "Error: length of sequence '$sKey' is smaller than the window size\n" unless(length($sSeq) > $nWinSize);
   my $nSeqLength = length($sSeq);
   for (my $i = 0; ($i <= ($nSeqLength - $nWinSize - 1)); $i += $nStepSize){
      my $sWindowSeq = substr($sSeq, $i, $nWinSize);
      my $nGCcontent = calcGC($sWindowSeq);
      my $nMidpoint   = $i + int($nWinSize/2);
      print join("\t", $sKey, $nMidpoint, sprintf("%1.5f", $nGCcontent)), "\n";
   }
}


#################
## SUBROUTINES ##
#################

# calcGC
#
# Calculate the GC content
sub calcGC{
   my $sSeq   = shift @_;
   my $nGCcnt = $sSeq =~ tr/CG/CG/;
   return $nGCcnt/length($sSeq);
}

# getSeq
#
# Read sequences from FASTA file
sub getSeq{
   my $sInput = shift @_;
   my $sCurrentSeqID = "";
   my %hSeq;
   
   open INPUT, $sInput or die "Could not open $sInput: $!\n";
   while (<INPUT>){
      next if /^\s*$/;
      next if (/^ *#/);
      s/[\n\r]+$//;
      if (/^\>/){
         $sCurrentSeqID = $_;
         $sCurrentSeqID =~ s/>//g;
         die "Error: duplicate fasta header found on line $.\n" if (exists $hSeq{$sCurrentSeqID});
      }
      else {
         my $sSeq = uc($_);
         die "Error: sequence file does not appear to have a fasta header\n" unless ($sCurrentSeqID);
         die "Error: illegal characters on line $.\n" unless ($sSeq =~ /^[ACGTUN]+$/);
         $hSeq{$sCurrentSeqID} .= $sSeq;
      }
   }
   close INPUT;
   return \%hSeq;
}
