#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use ValidateFiletypes qw(check_gff);

# GLOBALS
$ENV{SORT}        = 'sort';
$ENV{SORT_BUFFER} = '2G';

# ARGUMENTS
my $sHelp       = 0;
GetOptions("help!"      => \$sHelp);
my $sGffFile = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sGffFile);
if ($sHelp) {
    my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
    die <<HELP

    $sScriptName <gff-file>

    Takes a gff file and converts it into a file that can be loaded as a 
    custom annotation track in gbrowse

    -help
      This help message
      
HELP
}


# Check the gff file format
my @asGffErrors = check_gff($sGffFile);
if (@asGffErrors){
   unshift @asGffErrors, "The following errors were found in the input gff file:";
   die join("\n", @asGffErrors), "\n";
}


# Sort the gff file by name and process the entries
my $sLastChr = '';
open GFF, "$ENV{SORT} -S $ENV{SORT_BUFFER} -s -t '\t' -k1,1 -k4,4n $sGffFile |" or die "Error sorting gff file: $!\n";
while (<GFF>){
   next if /^\s*$/;
   next if /^\s*[Tt]rack/;
   next if /^\s*#/;
   s/[\n\r]$//g;
   my ($sChr, $sSource, $sType, $nStart, $nEnd, $nScore, $sStrand, $nFrame, $sID) = split /\t/;
   if ($sChr ne $sLastChr){
      print "\n" if ($sLastChr);
      print "reference = $sChr\n";
      $sLastChr = $sChr;
   }
   print join("\t", "\"$sSource\"", "\"$sID\"", "$nStart-$nEnd"), "\n";
}

