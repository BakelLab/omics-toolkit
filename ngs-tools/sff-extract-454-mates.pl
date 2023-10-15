#!/usr/bin/env perl

# 14.02.2011 19:28:53 EST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
my $sPrefix      = '';
GetOptions("help!"    => \$sHelp,
           "prefix:s" => \$sPrefix);

# PRINT HELP
$sHelp = 1 unless(@ARGV>0);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName -p <prefix> <fastq-file>
    
    Extract 454 mate pairs from a fastq file generated by
    the ssf-extract utility.
    
    -p --prefix <string>
      Output file prefix
    -help
      This help message
   
HELP
}

##########
## MAIN ##
##########

foreach my $sFile (@ARGV){
   my ($sLastID, $nLastCount) = ('', 2);
   open INPUT, $sFile or die "Error: can't open '$sFile': $!\n";
   open READ1, ">$sPrefix-read1.fastq" or die "Error: can't write to '$sPrefix-read1.fastq': $!\n";
   open READ2, ">$sPrefix-read2.fastq" or die "Error: can't write to '$sPrefix-read2.fastq': $!\n";
   open SINGLE, ">$sPrefix-single.fastq" or die "Error: can't write to '$sPrefix-single.fastq': $!\n";
   while (<INPUT>){
      next if (/^\s*$/);
      next if (/^ *#/);
      my $sID   = $_;
      my $sSeq  = <INPUT>;
      my $sTmp  = <INPUT>;
      my $sQual = <INPUT>;
      chomp $sID;
      chomp $sSeq;
      chomp $sTmp;
      chomp $sQual;
      die "Error: '$sFile' does not appear to be a fastq file\n" unless ($sID =~ /^\@/ and $sTmp =~ /^\+/);
      if ($sID =~ /\.fn$/){
         print SINGLE join("\n", $sID, $sSeq, $sTmp, $sQual), "\n";
      }
      elsif ($sID =~ /\.[rf]$/){
         my $sTrimmedID = $sID;
         $sTrimmedID =~ s/\.[rf]$//;
         # Make sure we're only printing pairs
         if ($sTrimmedID eq $sLastID){
            $nLastCount++;
         }
         else{
            die "Error: incomplete mate pair for '$sLastID' in file '$sFile'\n" unless ($nLastCount==2);
            $nLastCount = 1;
            $sLastID    = $sTrimmedID;
         }
         
         # Print the reads
         if ($sID =~ /\.f$/){
            print READ1 join("\n", $sID, $sSeq, $sTmp, $sQual), "\n";
         }
         else{
            print READ2 join("\n", $sID, $sSeq, $sTmp, $sQual), "\n";
         }
      }
      else{
         warn "Error: read $sID does not end in '.f', '.r', or '.fn' in file '$sFile', skipping\n";
      }
   }
   close INPUT;
   close READ1;
   close READ2;
   close SINGLE;
}

