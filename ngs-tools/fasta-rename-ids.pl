#!/usr/bin/env perl

# 26.08.2011 13:51:59 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
my $sIDfile      = '';
my $sFastaFile   = '';
GetOptions("help!"    => \$sHelp,
           "ids:s"    => \$sIDfile,
           "fasta:s"  => \$sFastaFile);

# PRINT HELP
$sHelp = 1 unless($sIDfile and $sFastaFile);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName -f <fasta file> -i <id file>

   Renames the identifiers in the fasta headers according to a reference
   file that contains a mapping between the current and a new identifier.
   
   Options:
    -i --ids <string>
      Tab-delimited file with identifier mapping: <current-id> <new-id> 
    -f --fasta <string>
      Fasta file. The first element of the fasta header (before the first 
      whitespace character) will be used for ID matching.
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Read IDs from file
my %hIDs;
if ($sIDfile){
   open IN, "<$sIDfile" or die "Error: can't read ID file\n";
   while (<IN>){
      next if (/^\s*$/);
      next if (/^ *#/);
      s/[\n\r]+$//;
      my ($sOldID, $sNewID) = split /\t/;
      $sOldID =~ s/^>+//;
      $sNewID =~ s/^>+//;
      if (exists $hIDs{$sOldID}){
         die "Error: conflicting mapping data provided for sequence '$sOldID' in '$sIDfile'\n" unless ($hIDs{$sOldID} eq $sNewID);
      }
      else{
         $hIDs{$sOldID}   = $sNewID;
      }
   }
   close IN;
}

# Now process the fasta file
my $nTotCount;
my $nRemappedCount;
open FASTA, $sFastaFile or die "Error: can't read fasta file '$sFastaFile': $!\n";
while (<FASTA>){
   if (/^\s*>/){
      $nTotCount++;
      my ($sID, @asRest) = split /\s/;
      $sID =~ s/^\s*>//;
      $sID =~ s/[\n\r]+$//;
      if (exists $hIDs{$sID}){
         print ">$hIDs{$sID}\n";
         $nRemappedCount++;
      }
      else{
         print;
      }
   }
   else{
      print;
   }
}
close FASTA;

# Print stats
print STDERR "$nRemappedCount out of $nTotCount fasta headers were remapped\n";
