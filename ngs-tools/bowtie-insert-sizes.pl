#!/usr/bin/env perl

# 17.10.2010 13:33:47 EDT
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

   Usage: $sScriptName <bowtie-paired-output>
   
   Print insert sizes for paired reads in a bowtie output file. 
   The file must be sorted according to read pairs.
    
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

open IN, $ARGV[0] or die "Error: can't open $ARGV[0]: $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sName1,$sStrand1,$sTarget1,$nPos1,$sSeq1,$sQual1,$sMM1,$sMMlist1) = split /\t/;
   $_ = <IN>;
   if ($_){
      my ($sName2,$sStrand2,$sTarget2,$nPos2,$sSeq2,$sQual2,$sMM2,$sMMlist2) = split /\t/;
      $sName1 = $sName1 =~ / / ? (split(/\s/,$sName1))[0] : substr($sName1, 0, -2);
      $sName2 = $sName2 =~ / / ? (split(/\s/,$sName2))[0] : substr($sName2, 0, -2);
      if ($sName1 eq $sName2){
	 my $nFragmentSize = $nPos2>$nPos1 ? $nPos2 - $nPos1 + length($sSeq2) : $nPos1 - $nPos2 + length($sSeq1);
	 print join("\t",$sName1, $nFragmentSize), "\n";
      }
      else{
	 die "Error: Paired read mismatch ($sName1, $sName2) on line $.\n";
      }
   }
   else{
      die "Error: Missing paired read at end of file\n";
   }
}
close IN;
