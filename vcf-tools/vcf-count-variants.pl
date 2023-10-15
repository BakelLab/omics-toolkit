#!/usr/bin/env perl

# 09.09.2011 11:15:48 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
my $sHeader      = '';
GetOptions("help!"    => \$sHelp,
           "header:s" => \$sHeader);
my $sInput = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sInput and $sHeader);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName -h header <vcf-file>
    
    -header <string>
      Comma-separated header info for the multi-sample calls
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Parse the header
my @asHeader = split /,/, $sHeader;

# Now process the vcf file
my %hCalls;
open VCF, $sInput or die "Error: can't open '$sInput': $!\n";
while (<VCF>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sChr, $nPos, $sID, $sRef, $sAlt, $sQual, $sFilter, $sInfo, $sFormat, @asGenotypes) = split /\t/, $_, -1;
   die "Error: number of header fields don't match the number of genotypes\n" unless(@asHeader == @asGenotypes);
   
   # Clean up the genotype calls
   for (my $i=0 ; $i<@asGenotypes ; $i++){
      $asGenotypes[$i] =~ s/\:.*$//;
   }
   
   # Get the heterozygosity calls in each sample
   for (my $i=0 ; $i<@asGenotypes ; $i++){
      my ($a, $b) = split /\/|\|/, $asGenotypes[$i];
      $hCalls{$asHeader[$i]}{$asHeader[$i]}++ if ($a ne $b);
   }
   
   # Get the differences between each sample
   for (my $i=0 ; $i<@asGenotypes ; $i++){
      for (my $j=1; $j<@asGenotypes ; $j++){
         $hCalls{$asHeader[$i]}{$asHeader[$j]}++ if ($asGenotypes[$i] ne $asGenotypes[$j]);
      }
   }
   
}
close VCF;


foreach my $sRow (sort keys %hCalls){
   foreach my $sCol (sort keys %{$hCalls{$sRow}}){
      print join("\t", $sRow, $sCol, $hCalls{$sRow}{$sCol}), "\n";
   }
}
