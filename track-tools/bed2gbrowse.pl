#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use ValidateFiletypes qw(check_bed);

# GLOBALS
$ENV{SORT}        = 'sort';
$ENV{SORT_BUFFER} = '2G';

# ARGUMENTS
my $sHelp       = 0;
my $sTrackName  = "CustomBedTrack";
GetOptions("help!"       => \$sHelp,
           "trackname=s" => \$sTrackName);
my $sBedFile = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sBedFile);
if ($sHelp) {
    my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
    die <<HELP

    $sScriptName <bed-file>

    Takes a bed file and converts it into a file that can be 
    loaded as a custom annotation track in gbrowse.
    The current version only supports bed basic (4-column) formats.

    Arguments:
    -t <string>
      The track name. Default: $sTrackName
    -help
      This help message
      
HELP
}

# Check the bed file format
my @asBedErrors = check_bed($sBedFile);
if (@asBedErrors){
   unshift @asBedErrors, "The following errors were found in the input bed file:";
   die join("\n", @asBedErrors), "\n";
}


# Sort the Bed file by name and process the entries
my $sLastChr = '';
open BED, "$ENV{SORT} -S $ENV{SORT_BUFFER} -s -t '\t' -k1,1 -k2,2n $sBedFile |" or die "Error sorting bed file: $!\n";
while (<BED>){
   next if /^\s*$/;
   next if /^\s*[Tt]rack/;
   next if /^\s*#/;
   s/[\n\r]$//g;
   my ($sChr, $nStart, $nEnd, $sID, @asRest) = split /\t/;
   if ($sChr ne $sLastChr){
      print "\n" if ($sLastChr);
      print "reference = $sChr\n";
      $sLastChr = $sChr;
   }
   print join("\t", "\"$sTrackName\"", "\"$sID\"", "$nStart-$nEnd"), "\n";
}
close BED;
