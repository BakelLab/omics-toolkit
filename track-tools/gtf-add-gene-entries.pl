#!/usr/bin/perl

# 22.08.2019 18:31:20 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
GetOptions("help!"   => \$sHelp);

# PRINT HELP
$sHelp = 1 unless(@ARGV>0);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Script to add 'gene' lines to a GTF file that does not have them.

   Usage: $sScriptName <gtf-file>
    
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# We buffer data per chromosome and make the assumption that the data is at least sorted by chromosome
my @asPrintOrder;
my %hBuffer;
my $sLastChr = "";

open IN, "sort -t '\t' -k1,1 -k4,4n -k5,5n $ARGV[0] | gffread -F --keep-exon-attrs -T -o- |" or die "Error: can't open 'gffread $ARGV[0]': $!\n";
while (<IN>){
   
   # Split the line
   my ($sChr, $sSrc, $sType, $nStart, $nEnd, $nScore, $sStrand, $sOffset, $sName) = split /\t/;

   # Print the buffer
   if ( ($sChr ne $sLastChr) or eof ){
      foreach my $sGeneID (@asPrintOrder){
         print join("\t", $hBuffer{$sGeneID}{GeneChr}, $hBuffer{$sGeneID}{GeneSrc}, "gene", $hBuffer{$sGeneID}{GeneStart}, $hBuffer{$sGeneID}{GeneEnd}, ".", $hBuffer{$sGeneID}{GeneStrand}, ".", "gene_id \"$sGeneID\";"), "\n";
         print $hBuffer{$sGeneID}{buffer};
      }
      @asPrintOrder = ();
      %hBuffer = ();
   }
   
   # Fill the buffer
   if (/gene_id "(\S+)";/){
      my $sGeneID = $1;
      if (exists $hBuffer{$sGeneID}){
         $hBuffer{$sGeneID}{GeneStart}  = $nStart if ($nStart < $hBuffer{$1}{GeneStart} );
         $hBuffer{$sGeneID}{GeneEnd}    = $nEnd   if ($nEnd   > $hBuffer{$1}{GeneEnd} );
      }
      else{
         $hBuffer{$sGeneID}{GeneStart}  = $nStart;
         $hBuffer{$sGeneID}{GeneEnd}    = $nEnd;
         $hBuffer{$sGeneID}{GeneStrand} = $sStrand;
         $hBuffer{$sGeneID}{GeneChr}    = $sChr;
         $hBuffer{$sGeneID}{GeneSrc}    = $sSrc;
         push @asPrintOrder, $sGeneID;
      }
      $hBuffer{$sGeneID}{buffer} .= $_;
   }
   $sLastChr = $sChr;
}
close IN;
