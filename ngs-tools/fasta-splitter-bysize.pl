#!/usr/bin/env perl

# 08.27.2014 13:35:20 EST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use POSIX;

# GET PARAMETERS
my $sHelp      = 0;
my $sPrefix    = "";
my $nMinLength   = 2000000000;
GetOptions("prefix:s" => \$sPrefix,
           "size:i"   => \$nMinLength);

# PRINT HELP
$sHelp = 1 unless(@ARGV);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [-n -p] <fasta-file>
    
   Split a file containing multiple fasta sequences in a set
   of smaller files, each containing a subset of sequences that
   add up to the approximate specified length.
    
    -s --size <integer>
      Minimum total length of fasta sequences in each split file
      default: $nMinLength
    -p --prefix <string>
      Output file prefix
   
HELP
}


###########
## START ##
###########

# Check arguments
die "Error: minimum length of fasta sequences per file must be greater than 1\n" unless ($nMinLength > 0);

# Figure out how many files we're going to output
my $nFileSize   = -s $ARGV[0];
my $nTotalFiles = POSIX::ceil($nFileSize / $nMinLength);
my $nPadding    = length($nTotalFiles);
die "Error: split would result in $nTotalFiles output files (max allowed is 2000)\n" if ($nTotalFiles > 2000);

# Start splitting the fasta file
my $nLengthCounter = -1;
my $nIncrement  = 1;
open IN, $ARGV[0] or die "Error: can't open '$ARGV[0]': $!\n";
while (<IN>){
   next if /^\s*$/;
   next if /^\s*#/;
   if ($nLengthCounter == -1){
      if (/^\s*>/){
         my $sIncrement = sprintf("%0${nPadding}d", $nIncrement);
         open OUT, ">${sPrefix}${sIncrement}.fa" or die "Error: can't write to '${sPrefix}${nIncrement}.fa': $!\n";
         print OUT $_;
         $nLengthCounter = 0;
         $nIncrement++;
      }
      else{
         die "Error: file '$ARGV[0]' does not start with a fasta header\n";
      }
   }
   else{
      if (/^\s*>/){
         if ($nLengthCounter > $nMinLength){
            close OUT;
            my $sIncrement = sprintf("%0${nPadding}d", $nIncrement);
            open OUT, ">${sPrefix}${sIncrement}.fa" or die "Error: can't write to '${sPrefix}${nIncrement}.fa': $!\n";
            $nLengthCounter = 0;
            $nIncrement++;
         }
         print OUT $_;
      }
      else{
         print OUT $_;
         $nLengthCounter += length($_)-1;
      }
   }
}
close IN;
close OUT;
