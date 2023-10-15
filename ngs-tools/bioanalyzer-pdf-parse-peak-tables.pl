#!/usr/bin/env perl

# 15.08.2018 14:10:06 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

# GET PARAMETERS
my $sHelp        = 0;
GetOptions("help!"   => \$sHelp);

# PRINT HELP
$sHelp = 1 unless(@ARGV>0);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName <bioanalyzer-output.pdf>
    
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

my $sChipName = basename($ARGV[0]);
$sChipName =~ s/^.*DE13805199_//;
$sChipName =~ s/\.pdf$//;
my $sSampleID = "";
open IN, "less \"$ARGV[0]\" | " or die "Error: can't open $ARGV[0]: $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   if (/^Peak table for sample .*:\s+(\S+.*)$/){
      $sSampleID = $1;
   }
   if (/^(\d+)\s+(\d+)\s+(\d+\.?\d*)\s+(\d+\.?\d*)$/){
      if ($sSampleID){
         print "$sChipName\t$sSampleID\t$1\t$2\t$3\t$4\n"
      }
   }
}
close IN;
