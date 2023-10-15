#!/usr/bin/env perl

# 03.08.2010 13:02:34 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
my $sPrimers     = "";
my $sLibrary     = "";
GetOptions("help!"     => \$sHelp,
           "primers:s" => \$sPrimers,
           "library:s" => \$sLibrary);

# PRINT HELP
$sHelp = 1 unless($sPrimers and $sLibrary);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName
    
   Checks a bunch of primers against a library to see if they
   uniquely amplify only a single library sequence.
    
    -primers
      File with primer sequences
    -library
      File with library sequences
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Read in the sequence library
my @asLibrary = read_library($sLibrary);

# Now match the primers
open PRIMER, "<$sPrimers" or die "Error: can't open '$sPrimers': $!\n";
print "primer\tfwd_seq\trev_seq\tmatch_id\tmatch_fragment_length\n";
while (<PRIMER>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   
   my $flMatched = 0;
   my ($sPrimerID, $sFwd, $sRev) = split /\t/;
   $sFwd = uc($sFwd);
   $sRev = reverse_complement(uc($sRev));
   die "Error: illegal forward primer sequence on line $. in file '$sPrimers'\n" unless $sFwd =~ (/^[ACTG]+$/);
   die "Error: illegal reverse primer sequence on line $. in file '$sPrimers'\n" unless $sRev =~ (/^[ACTG]+$/);
   
   # Check each library seq
   foreach my $rLibrarySeq (@asLibrary){
      my ($sLibraryID, $sLibrarySeq_w) = @$rLibrarySeq;
      my $sLibrarySeq_c = reverse_complement($sLibrarySeq_w);
      
      # Match watson strand
      my ($nFwdMatch, $nRevMatch) = match_primers($sFwd, $sRev, $sLibrarySeq_w);
      
      # Match crick strand
      unless ($nFwdMatch and $nRevMatch){
         ($nFwdMatch, $nRevMatch) = match_primers($sFwd, $sRev, $sLibrarySeq_c);
      }
      
      # Check if we have a valid primer pair
      if ($nFwdMatch and $nRevMatch){
         if ($nFwdMatch < $nRevMatch){
            my $nFragmentLength = length($sFwd) + length($sRev) + ($nRevMatch - $nFwdMatch - 1);
            print join("\t", $sPrimerID, $sFwd, reverse_complement($sRev), $sLibraryID, $nFragmentLength), "\n";
            $flMatched = 1;
         }
      }
   }
   
   # Print something if there was no match in the library
   print join("\t", $sPrimerID, $sFwd, reverse_complement($sRev), 'no match', 'no match'), "\n" unless ($flMatched);
}
close PRIMER;


#################
## SUBROUTINES ##
#################

# read_library
#
# Read in the library of primer sequences
sub read_library {
   my ($sFile) = @_;
   my @asLibrary;
   
   open LIB, "<$sFile" or die "Error: can't open '$sFile': $!\n";
   while (<LIB>){
      next if (/^\s*$/);
      next if (/^ *#/);
      s/[\n\r]+$//;
      my ($sID, $sSeq) = split /\t/;
      $sSeq = uc($sSeq);
      die "Error: illegal sequence characters found on line $. in file $sFile\n" unless $sSeq =~ (/^[ACTG]+$/);
      push @asLibrary, [($sID, $sSeq)];
   }
   close LIB;
   return @asLibrary;
}


# reverse_complement
#
# Returns the reverse-complement of the supplied sequence
sub reverse_complement{
   my $seq = shift(@_);
   my $rev = reverse $seq;
   $rev =~ tr/ACGTacgt/TGCAtgca/;
   return $rev;
}


# match_primers
#
# match primers against a sequence
sub match_primers {
   my ($sFwd, $sRev, $sSeq) = @_;
   my ($nFwdMatch, $nRevMatch) = (0,0);
   $nFwdMatch = $sSeq =~ /$sFwd/g ? pos($sSeq) : 0;
   if ($nFwdMatch){
      while($sSeq =~ /$sRev/g){
         $nRevMatch = pos($sSeq);
      }
   }
   $nRevMatch = $nRevMatch - length($sRev)+1 if ($nRevMatch);
   return ($nFwdMatch, $nRevMatch);
}
