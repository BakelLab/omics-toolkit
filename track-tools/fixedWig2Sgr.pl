#!/usr/bin/env perl

# 14.09.2015 20:59:41 EDT
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

my ($sChr, $nPos, $nStep) = ("", 0, 0);
open IN, $ARGV[0] or die "Error: can't open '$ARGV[0]': $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   if (/^fixedStep/){
      ($sChr, $nPos, $nStep) = $_ =~ /fixedStep chrom=(\S+) start=(\d+) step=(\d+)/;
   }
   else{
      print join("\t", $sChr, $nPos, $_), "\n";
      $nPos += $nStep;
   }
}
close IN;
