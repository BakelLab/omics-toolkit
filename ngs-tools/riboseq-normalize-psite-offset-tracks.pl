#!/usr/bin/perl

# 13.01.2020 16:30:30 EST
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

   Usage: $sScriptName
    
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Get sum of read counts per transcript
my %hSums;
open IN, $ARGV[0] or die "Error: can't open '$ARGV[0]:': $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sID, $nPos, $nCount) = split /\t/;
   if ($nPos >=-15 and $nPos <= 30){
      $hSums{$sID} += $nCount;
   }
}
close IN;

# Normalize each transcript by total read count per transcript
open IN, $ARGV[0] or die "Error: can't open '$ARGV[0]:': $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sID, $nPos, $nCount) = split /\t/;
   if ($nPos >=-15 and $nPos <= 30 and $hSums{$sID} > 50){
      print join("\t", $sID, $nPos, ($nCount/$hSums{$sID})*100), "\n"; 
   }
}
close IN;





