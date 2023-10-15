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
my $sHelp         = 0;
my $nMinLength    = 0;
my $nMaxLength    = 0;
GetOptions("help!"       => \$sHelp,
           "minlength:i" => \$nMinLength,
           "maxlength:i" => \$nMaxLength);

# PRINT HELP
$sHelp = 1 unless(@ARGV);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [options] <fastq-file1> .. <fastq-fileN>
   
   Extract sequences from a fastq file that meet
   certain length thresholds.
    
   Options:
    -minlength <integer>
      Minumum length of fasta sequence
      Default: 0 (no filtering)
    -maxlength <integer>
      Maximum length of fasta sequence
      Default: 0 (no filtering)
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
   while(<$fhInput>) {
      my $sName = $_;
      my $sSeq  = <$fhInput>;
      my $sTmp  = <$fhInput>;
      my $sQual = <$fhInput>;
      $sSeq =~ s/[\n\r]+$//;
      my $flPrint = 1;
      $flPrint = 0 if ($nMinLength and (length($sSeq) < $nMinLength));
      $flPrint = 0 if ($nMaxLength and (length($sSeq) > $nMaxLength));
      print "$sName$sSeq\n$sTmp$sQual" if $flPrint;
   }
   close($fhInput);
}
