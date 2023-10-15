#!/usr/bin/env perl

# Modified from a bowtie script written by Ben Langmead

# MODULES
use strict;
use warnings;
use Getopt::Long;
use IO::File;
use IO::Zlib qw(:gzip_external 1);

# GET PARAMETERS
my $sHelp    = 0;
my $sMate1File    = "";
my $sMate2File    = "";
GetOptions ("1=s"      => \$sMate1File,
            "2=s"      => \$sMate2File) || die "Bad option";

# PRINT HELP
$sHelp = 1 unless($sMate1File and $sMate2File);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName -1 <fastq1> -2 <fastq2>
    
   Check whether all reads are properly paired
   
   Arguments:
    -1 <string>
      Name of file containing the first set of paired reads.
    -2 <string>
      Name of file containing the matched paired reads.
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Open the input files
my ($fhM1, $fhM2);
if($sMate1File =~ /\.gz$/) {
   $fhM1 = new IO::Zlib($sMate1File, "rb") or die "Error: can't open '$sMate1File': $!\n";
}
else{
   $fhM1 = new IO::File($sMate1File, "r") or die "Error: can't open '$sMate1File': $!\n";
}
if($sMate2File =~ /\.gz$/) {
   $fhM2 = new IO::Zlib($sMate2File, "rb") or die "Error: can't open '$sMate2File': $!\n";
}
else{
   $fhM2 = new IO::File($sMate2File, "r") or die "Error: can't open '$sMate2File': $!\n";
}

# Process the reads
my ($nMatched, $nUnmatched) = (0,0);
while(<$fhM1>) {
   my $sName1 = $_;
   my $sName2 = <$fhM2>;
   my $sSeq1  = <$fhM1>;
   my $sSeq2  = <$fhM2>;
   my $sTmp1  = <$fhM1>;
   my $sTmp2  = <$fhM2>;
   my $sQual1 = <$fhM1>;
   my $sQual2 = <$fhM2>;
   if ($sQual1 and $sQual2){
      die "Error: read name '$sName1' in file '$sMate1File' does not conform to the fastq standard\n" unless ($sName1 =~ /^\@/);
      die "Error: read name '$sName2' in file '$sMate2File' does not conform to the fastq standard\n" unless ($sName2 =~ /^\@/);
      my ($sMatch1, $sMatch2) = ($sName1, $sName2);
      $sMatch1 =~ s/[\n\r]+$//;
      $sMatch2 =~ s/[\n\r]+$//;
      $sMatch1 = substr $sMatch1, 1;
      $sMatch2 = substr $sMatch2, 1;
      if( (substr($sMatch1,0,-2) eq substr($sMatch2,0,-2)) and (length($sSeq1) == length($sQual1)) and (length($sSeq2) == length($sQual2)) ){
         $nMatched++;
	 print "$sMatch1\n";
      }
      else{
         $nUnmatched++
      }
   }
   else{
      die "Error: Missing paired reads at end of file\n";
   }
}
close($fhM1);
close($fhM2);
warn "$nMatched matched reads\n$nUnmatched unmatched reads\n";
