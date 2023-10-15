#!/usr/bin/env perl

# Modified from a bowtie script written by Ben Langmead

# MODULES
use strict;
use warnings;
use Getopt::Long;
use IO::File;
use IO::Zlib qw(:gzip_external 1);

# GET PARAMETERS
my $sHelp      = 0;
my $nMaxLength = 32;
GetOptions("help!"   => \$sHelp);

# PRINT HELP
$sHelp = 1 unless(@ARGV);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName <fastq-file>
    
   Gets read and base count stats from a fastq file
   
   Arguments:
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

foreach my $sInput (@ARGV){
   
   # Open the input files
   my ($fhInput);
   if($sInput =~ /\.gz$/) {
      $fhInput = new IO::Zlib($sInput, "rb") or die "Error: can't open '$sInput': $!\n";
   }
   else{
      $fhInput = new IO::File($sInput, "r") or die "Error: can't open '$sInput': $!\n";
   }
      
   # Process the reads
   my $nReadCount = 0;
   my $nBaseCount = 0;
   while(<$fhInput>) {
      my $sName = $_;
      my $sSeq  = <$fhInput>;
      my $sTmp  = <$fhInput>;
      my $sQual = <$fhInput>;
      $nReadCount++;
      $sSeq =~ s/[\n\r]+$//;
      $nBaseCount += length($sSeq);
   }
   close($fhInput);
   
   print "$sInput\treads:\t$nReadCount\n";
   print "$sInput\tbases:\t$nBaseCount\n";
}
