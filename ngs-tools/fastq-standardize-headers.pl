#!/usr/bin/env perl

# 26.07.2012 14:08:43 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use IO::File;
use IO::Zlib qw(:gzip_external 1);

# GET PARAMETERS
my $sHelp        = 0;
my $nSuffix      = 0;
GetOptions("help!"    => \$sHelp,
           "suffix:n" => \$nSuffix);

# PRINT HELP
$sHelp = 1 unless($nSuffix and @ARGV);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName -s <suffix> <fastq-file1> .. [fastq-fileN]
    
   Converts any fastq header to a standard format where each
   read ends in /1 or /2 to denote the orientation
   
   Arguments:
    -suffix <integer>
      The read suffix integer (required)
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
   my $nSkipped = 0;
   while(<$fhInput>) {
      my $sName = $_;
      my $sSeq  = <$fhInput>;
      my $sTmp  = <$fhInput>;
      my $sQual = <$fhInput>;
      $sName =~ s/[\n\r]+$//;
      $sName =~ s/ .*$//;
      $sName =~ s/[\/:]\d$//;
      print "$sName/$nSuffix\n" . $sSeq . "+\n" . $sQual;
   }
   close $fhInput;
}
