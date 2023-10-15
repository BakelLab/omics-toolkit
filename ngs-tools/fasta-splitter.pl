#!/usr/bin/env perl

# 08.11.2009 13:35:20 EST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use POSIX;

# GET PARAMETERS
my $sHelp      = 0;
my $nMaxCount  = 200;
my $sPrefix    = "";
GetOptions("prefix:s" => \$sPrefix,
           "number:i" => \$nMaxCount);

# PRINT HELP
$sHelp = 1 unless(@ARGV);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [-n -p] <fasta-file>
    
   Split a file containing multiple fasta sequences in a set
   of smaller files, each containing a specified number of sequences
    
    -n --number <integer>
      Maximum number of fasta sequences per file
      default: $nMaxCount
    -p --prefix <string>
      Output file prefix
   
HELP
}


###########
## START ##
###########

# Check arguments
die "Error: maximum number of fasta sequences per file must be greater than 1\n" unless ($nMaxCount > 0);

# Figure out how many files we're going to output
my $nTotalCount = count_fasta_headers($ARGV[0]);
my $nTotalFiles = POSIX::ceil($nTotalCount / $nMaxCount);
my $nPadding    = length($nTotalFiles);
die "Error: file '$ARGV[0]' does not appear to be a fasta file\n" unless ($nTotalCount);
die "Error: split would result in $nTotalFiles output files (max allowed is 2000)\n" if ($nTotalFiles > 2000);

# Start splitting the fasta file
my $nFastaCounter = 0;
my $nIncrement  = 1;
open IN, $ARGV[0] or die "Error: can't open '$ARGV[0]': $!\n";
while (<IN>){
   next if /^\s*$/;
   next if /^\s*#/;
   if ($nFastaCounter == 0){
      if (/^\s*>/){
         my $sIncrement = sprintf("%0${nPadding}d", $nIncrement);
         open OUT, ">${sPrefix}${sIncrement}.fa" or die "Error: can't write to '${sPrefix}${nIncrement}.fa': $!\n";
         print OUT $_;
         $nFastaCounter = 1;
         $nIncrement++;
      }
      else{
         die "Error: file '$ARGV[0]' does not start with a fasta header\n";
      }
   }
   else{
      if (/^\s*>/){
         $nFastaCounter++;
         if ($nFastaCounter > $nMaxCount){
            close OUT;
            my $sIncrement = sprintf("%0${nPadding}d", $nIncrement);
            open OUT, ">${sPrefix}${sIncrement}.fa" or die "Error: can't write to '${sPrefix}${nIncrement}.fa': $!\n";
            $nFastaCounter = 1;
            $nIncrement++;
         }
         print OUT $_;
      }
      else{
         print OUT $_;
      }
   }
}
close IN;
close OUT;

#################
## SUBROUTINES ##
#################

# count_fasta_headers
#
# Returns a count of the total number of fasta headers in a multi-fasta file
sub count_fasta_headers{
   my $sFile = shift @_;
   my $nCount = 0;
   open COUNT, $sFile or die "Error: Can't open '$sFile': $!\n";
   while (<COUNT>){
      $nCount++ if (/^\s*>/);
   }
   close COUNT;
   return $nCount;
}

