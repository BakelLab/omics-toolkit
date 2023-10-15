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
my $sOut1File     = "";
my $sMate2File    = "";
my $sOut2File     = "";
my $sIDfile       = "";
my $flStrip       = 0;
my $flInvertMatch = 0;
GetOptions ("1=s"      => \$sMate1File,
            "2=s"      => \$sMate2File,
            "o1=s"     => \$sOut1File,
            "o2=s"     => \$sOut2File,
            "ids=s"    => \$sIDfile,
            "strip!"   => \$flStrip,
            "help!"    => \$sHelp,
            "v|invert-match!" => \$flInvertMatch) || die "Bad option";

# PRINT HELP
$sHelp = 1 unless($sMate1File and $sMate2File and $sOut1File and $sOut2File and ($sIDfile or @ARGV));
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName -1 <fastq1> -2 <fastq2> -o1 <out1> -o2 <out2> [-i <id file>] [<id1> ... <idN>]
    
   Remove reads from a set of paired-end fastq files based on the common part of the
   read identifier. IDs can be read from a file or specified on the command line.
   The script reads/writes uncompressed or gzip compressed files, depending
   on whether a .gz extension is present in the input or output filename
   
   Arguments:
    -1 <string>
      Name of file containing the first set of paired reads.
    -2 <string>
      Name of file containing the matched paired reads.
    -o1 <string>
      Name of first read output file
    -o2 <string>
      Name of second read output file
    -ids <string>
      Name of the file with read identifiers
    -strip
      Strip off pair identifier (e.g. /1 and /2) before matching IDs
    -v
      Invert the sense of matching, to select non-matching reads.
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Read the IDs to filter
my %hIDs;
if ($sIDfile){
   open IN, $sIDfile or die "Error: can't open '$sIDfile': $!\n";
   while (<IN>){
      next if (/^\s*$/);
      next if (/^ *#/);
      s/[\n\r]+$//;
      s/^@//;
      s/[\/:]\d$// if($flStrip);
      $hIDs{$_} = "";
   }
   close IN;
}

# Read IDs from the command line
foreach my $sID (@ARGV){
   $sID =~ s/^@//;
   $sID =~ s/[\/:]\d$// if($flStrip);
   $hIDs{$sID} = "";
}

# Make sure we have IDs at this point
die "Error: please specify IDs to extract\n" unless(%hIDs);

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

# Open the output files
my ($fhO1, $fhO2);
if($sOut1File =~ /\.gz$/) {
   $fhO1 = new IO::Zlib($sOut1File, "wb6") or die "Error: can't open '$sOut1File': $!\n";
} 
else {
   $fhO1 = new IO::File($sOut1File, "w") or die "Error: can't open '$sOut1File': $!\n";
}
if($sOut2File =~ /\.gz$/) {
   $fhO2 = new IO::Zlib($sOut2File, "wb6") or die "Error: can't open '$sOut2File': $!\n";
} 
else {
   $fhO2 = new IO::File($sOut2File, "w") or die "Error: can't open '$sOut2File': $!\n";
}

# Process the reads
my $nExtracted = 0;
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
      $sMatch1 =~ s/ .*$//;
      $sMatch2 =~ s/ .*$//;
      die "Error: mismatch between paired reads ($sMatch1 vs $sMatch2)\n" unless(substr($sMatch1,0,-2) eq substr($sMatch2,0,-2));
      $sMatch1 =~ s/[\/:]\d$// if($flStrip);
      $sMatch2 =~ s/[\/:]\d$// if($flStrip);
      my $flExtract = exists($hIDs{$sMatch1}) and exists($hIDs{$sMatch2}) ? 1 : 0;
      $flExtract = 1 - $flExtract if ($flInvertMatch);
      if ($flExtract){
         print $fhO1 $sName1,$sSeq1,$sTmp1,$sQual1;
         print $fhO2 $sName2,$sSeq2,$sTmp2,$sQual2;
         $nExtracted++;
      }
   }
   else{
      die "Error: Missing paired reads at end of file\n";
   }
}
close($fhM1);
close($fhM2);
close($fhO1);
close($fhO2);

# Report the number of sequences dropped
my $sMessage = $nExtracted == 1 ? "Extracted $nExtracted read pair\n" : "Extracted $nExtracted read pairs\n";
warn ($sMessage);
