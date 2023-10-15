#!/usr/bin/env perl

# 11.06.2013 12:07:46 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use IO::File;
use IO::Zlib qw(:gzip_external 1);

# GET PARAMETERS
my $sHelp         = 0;
my $sFastq        = "";
my $sMatchFile    = "";
my $sPrefix       = "";
my $sPostfix      = "";
my $flInvertMatch = 0;
GetOptions ("fastq=s"         => \$sFastq,
            "match=s"         => \$sMatchFile,
            "help!"           => \$sHelp,
            "prefix:s"        => \$sPrefix,
            "o|postfix:s"     => \$sPostfix,
            "v|invert-match!" => \$flInvertMatch) || die "Bad option";

# PRINT HELP
$sHelp = 1 unless($sFastq and ($sMatchFile or @ARGV));
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName -f <fasta file> [-m <id file>] [<id1> ... <idN>]
    
   Filter reads in a fastq file based on ID matches
   
   Arguments:
    -f --fastq <string>
      Name of the file with fastq data
    -m --match <string>
      Name of the file with read identifiers
    -p --prefix <string>
      Optional prefix to be added to the reference IDs before matching
    -o --postfix <string>
      Optional postfix to be added to the reference IDs before matching
    -v --invert-match
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
if ($sMatchFile){
   open IN, $sMatchFile or die "Error\n";
   while (<IN>){
      next if (/^\s*$/);
      next if (/^ *#/);
      my $sID = $_;
      $sID =~ s/[\n\r]+$//;
      $sID =~ s/^@//;
      $sID =~ s/\s.*$//;
      $sID = join('', $sPrefix, $sID, $sPostfix);
      $hIDs{$sID} = '';
   }
   close IN;
}

# Read IDs from the command line
foreach my $sID (@ARGV){
   $sID =~ s/^@//;
   $sID = join('', $sPrefix, $sID, $sPostfix);
   $hIDs{$sID} = '';
}

# Open the input files
my ($fhFastq);
if($sFastq =~ /\.gz$/) {
   $fhFastq = new IO::Zlib($sFastq, "rb") or die "Error: can't open '$sFastq': $!\n";
}
else{
   $fhFastq = new IO::File($sFastq, "r") or die "Error: can't open '$sFastq': $!\n";
}

# Process the reads
my $nExtracted = 0;
while(<$fhFastq>) {
   my $name = $_;
   my $seq  = <$fhFastq>;
   my $tmp  = <$fhFastq>;
   my $qual = <$fhFastq>;
   if ($qual){
      die "Error: read name '$name' in file '$sFastq' does not conform to the fastq standard\n" unless ($name =~ /^\@/);
      my $match = $name;
      $match =~ s/^@//;
      $match =~ s/[\n\r]+$//;
      $match =~ s/\s.*$//;
      my $flExtract = exists $hIDs{$match} ? 1 : 0;
      $flExtract = 1 - $flExtract if ($flInvertMatch);
      if ($flExtract){
          print $name,$seq,$tmp,$qual;
          $nExtracted++;
      }
   }
   else{
      die "Error: Missing reads at end of file\n";
   }
}
close($fhFastq);

# Report the number of sequences extracted
my $sMessage = $nExtracted == 1 ? "Extracted $nExtracted sequence\n" : "Extracted $nExtracted sequences\n";
warn ($sMessage);
