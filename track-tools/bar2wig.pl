#!/usr/bin/env perl

# 03.08.2011 17:35:54 EDT
# Harm van Bakel <hvbakel@gmail.com>

# GLOBALS
$ENV{BAR2TXT} ||= 'bar2txt';

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
my $sDrop        = '';
GetOptions("help!"   => \$sHelp,
           "drop:s"  => \$sDrop);
my $sInput = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sInput);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName <bar-file>
   
   Arguments:
    -drop <string>
      Optional comma-separated list of chromosomes that you want to
      exclude from the wig file
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Build drop list
my %hDrop;
if ($sDrop){
   my @asDrop = split /,/, $sDrop;
   foreach my $sKey (@asDrop) {$hDrop{$sKey}++}
}

# Process the bar file
my $sLastChr = '';
open IN, "bar2txt -bar $sInput |" or die "Error: can't run bar2txt: $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sChr, $nPos, $nVal) = split /\t/, $_, -1;
   next if (exists $hDrop{$sChr});
   if ($sChr ne $sLastChr){
      print "variableStep chrom=$sChr\n";
      $sLastChr = $sChr;
   }
   print join("\t", $nPos, $nVal), "\n";
}
close IN;
