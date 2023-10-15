#!/usr/bin/env perl

# 29.05.2018 19:48:10 EDT
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

my ($nMin, $nMax, $sBlockScaff) = (0,0,"");
open IN, $ARGV[0] or die "Error: can't open file\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   next if (/\*+/);
   s/[\n\r]+$//;
   if (/^BLOCK/){
      if ($nMin > 0 or $nMax > 0){
         print "$sBlockScaff\t$nMin\t$nMax\n";
      }
      ($nMin, $nMax, $sBlockScaff) = (0,0,"");
   }
   if (/^\d+\t/){
      my ($nID, $nA, $nB, $sCID, $nPos, @rest) = split /\t/;
      if ($nMin == 0 and $nMax == 0){
         $nMin = $nPos;
         $nMax = $nPos;
         $sBlockScaff = $sCID;
      }
      else{
         $nMax = $nPos if ($nPos > $nMax);
      }
   }
   if (eof){
      print "$sBlockScaff\t$nMin\t$nMax\n";
   }
}
close IN;
