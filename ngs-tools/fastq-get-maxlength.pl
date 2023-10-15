#!/usr/bin/env perl

# 07.09.2010 20:00:48 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use IO::File;
use IO::Zlib qw(:gzip_external 1);

# GET PARAMETERS
my $sHelp    = 0;
my $nLines   = 10000;
GetOptions("help!"  => \$sHelp,
           "n:i"    => \$nLines);

# PRINT HELP
$sHelp = 1 unless(@ARGV);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [options] <fastq-file>
   
   Returns the maximum read length of an input fastq file.
    
   Options:
    -n <integer>
      Number of fastq lines to read from input.
      Default: $nLines
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

foreach my $sInput (@ARGV){
   my ($nMaxLength, $nCounter) = (0,0);

   # Open the input files
   my ($fhInput);
   if($sInput =~ /\.gz$/) {
      $fhInput = new IO::Zlib($sInput, "rb") or die "Error: can't open '$sInput': $!\n";
   }
   else{
      $fhInput = new IO::File($sInput, "r") or die "Error: can't open '$sInput': $!\n";
   }
      
   # Process the reads
   while(<$fhInput>) {
      my $sName = $_;
      my $sSeq  = <$fhInput>;
      my $sTmp  = <$fhInput>;
      my $sQual = <$fhInput>;
      $sSeq =~ s/[\n\r]+$//;
      my $nSeqLength = length($sSeq);
      $nMaxLength = $nSeqLength if ($nSeqLength > $nMaxLength);
      $nCounter++;
      last if $nCounter >= $nLines;
   }
   close($fhInput);
   
   print $nMaxLength, "\n";
}
