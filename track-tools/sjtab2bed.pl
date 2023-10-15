#!/usr/bin/env perl

# 06.03.2016 10:10:36 EST
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

   Usage: $sScriptName <SJ.out.tab>
   
   Convert a STAR SJ.out.tab file to a bed-formatted junction file
    
   Arguments:
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

my $nSJcount = 1;
open SJTAB, $ARGV[0] or die "Error: can't open '$ARGV[0]': $!\n";
while (<SJTAB>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sChr, $nIntronStart, $nIntronEnd, $nStr, $nIntronMotif, $flAnnotated, $nUniqueReads, $nMultiReads, $nMaxOverhang) = split /\t/;
   $nIntronStart--;
   
   # Format strand
   my $sStr = '.';
   $sStr = '+' if ($nStr == 1);
   $sStr = '-' if ($nStr == 2);
   
   # Feature color
   my $sCol = '0,255,0';
   $sCol = '0,0,255' if ($nStr == 1);
   $sCol = '255,0,0' if ($nStr == 2);
   
   # Feature start and end
   my $nFtStart = $nIntronStart - $nMaxOverhang;
   my $nFtEnd   = $nIntronEnd   + $nMaxOverhang;
   
   # Block end
   my $nBlockEnd = $nIntronEnd - $nIntronStart + $nMaxOverhang;
   
   print join("\t",
              $sChr,
              $nFtStart,
              $nFtEnd,
              "SJ$nSJcount|$nIntronMotif",
              $nUniqueReads + $nMultiReads,
              $sStr,
              $nFtStart,
              $nFtEnd,
              $sCol,
              2,
              "$nMaxOverhang,$nMaxOverhang",
              "0,$nBlockEnd"), "\n";
   $nSJcount++;
}
close SJTAB;

