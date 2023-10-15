#!/usr/bin/env perl

# 08.11.2011 13:22:29 EST
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

   Usage: $sScriptName <bam file>
    
   Print insert sizes for paired reads in a bam file. The
   file must be sorted according to read pairs.
    
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

open SAM, "samtools view $ARGV[0] |" or die "Error: can't run samtools to view '$ARGV[0]': $!\n";
while (<SAM>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sQnameA, $nFlagA, $sRnameA, $nPosA, $nMapQA, $sCigarA, $sRnextA, $nPnextA, $nTlenA, $sSeqA, $sQualA, @asRestA) = split /\t/;
   my ($sQnameB, $nFlagB, $sRnameB, $nPosB, $nMapQB, $sCigarB, $sRnextB, $nPnextB, $nTlenB, $sSeqB, $sQualB, @asRestB) = split /\t/, <SAM>;
   ($nTlenA, $nTlenB) = (abs($nTlenA), abs($nTlenB));
   if ( ($sQnameA eq $sQnameB) and ($nTlenA == $nTlenB) ){
      print join("\t", $sQnameA, $nTlenA), "\n";
   }
   else{
      die "Error: did not find read pairs in consecutive lines, please resort bam file with '-n' argument\n" ;
   }
}
close SAM;
