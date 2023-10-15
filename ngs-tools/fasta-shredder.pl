#!/usr/bin/env perl

# 07.09.2010 20:00:48 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp         = 0;
my $nSegmentsize  = 900;
my $nOverlap      = 300;
my $flSplitN      = 0;
GetOptions("help!"         => \$sHelp,
           "Nsplit!"       => \$flSplitN,
           "segmentsize:s" => \$nSegmentsize,
           "overlap:s"     => \$nOverlap);
my $sInput = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sInput);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [options] <fasta-file>
   
   Shred a (multi) fasta file into smaller segments with a specified 
   amount of overlap between segments.
   
   Options:
    -s --segmentsize <integer>
      The size of the shredded fragments
    -o --overlap <integer>
      The overlap between segments
    -N --Nsplit
      Split input sequences on 'N' nucleotides
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Process the fasta file
my $sFastaHeader = '';
my $sFastaSeq    = '';
open FASTA, "<$sInput" or die "Error: can't read '$sInput': $!\n";
while (<FASTA>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   if (/^>/){
      # Print the previous match if there is any
      if ($sFastaHeader){
         shredder($sFastaHeader, $sFastaSeq, $nSegmentsize, $nOverlap, $flSplitN);
      }
   
      # Reset everything for the new sequence
      $sFastaHeader = $_;
      $sFastaSeq    = '';
   }
   else{
      $sFastaSeq .= $_ if ($sFastaHeader);
   }
}
close FASTA;

# Print the last sequence
if ($sFastaHeader){
   shredder($sFastaHeader, $sFastaSeq, $nSegmentsize, $nOverlap, $flSplitN);
}


#################
## SUBROUTINES ##
#################


# shredder
#
# Shreds fasta seqs
sub shredder {
   my ($sFastaHeader, $sFastaSeq, $nSegmentsize, $nOverlap, $flSplitN) = @_;
   
   # Split fastaheader
   $sFastaHeader =~ s/^>\s+/>/;
   my @asFastaHeader = split /\s/, $sFastaHeader;
   my $sPrimaryHeader = $asFastaHeader[0];
   
   # Split fastaseq on 'N'
   my @asFastaSeq;
   if ($flSplitN){
       @asFastaSeq    = split /[nN]+/, $sFastaSeq;
   }
   else{
       push @asFastaSeq, $sFastaSeq;
   }
   
   # Process the individual sequences
   my $nScaffFragment = 0;
   foreach my $sSeq (@asFastaSeq){
      my $nOffset  = 0;
      my $nLastEnd = 0;
      my $nSeqlen  = length $sSeq;
      while($nLastEnd < $nSeqlen){
         my $nChunk   = $nOffset + $nSegmentsize < $nSeqlen ? $nSegmentsize : $nSeqlen - $nOffset;
         my $sFragSeq = substr ($sSeq, $nOffset, $nSegmentsize);
         $asFastaHeader[0] = join('', $sPrimaryHeader, ":", $nScaffFragment, ":", $nOffset+1, "-", $nOffset+$nChunk);
         print join(' ', @asFastaHeader), "\n";
         print "$sFragSeq\n";
         $nLastEnd = $nOffset+$nChunk;
         $nOffset += $nSegmentsize - $nOverlap;
      }
      $nScaffFragment++;
   }
}
