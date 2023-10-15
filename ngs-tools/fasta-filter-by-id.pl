#!/usr/bin/env perl

# 07.09.2010 20:00:48 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp         = 0;
my $sIDfile       = '';
my $sFastaFile    = '';
my $flInvertMatch = 0;
GetOptions("help!"           => \$sHelp,
           "match:s"         => \$sIDfile,
           "fasta:s"         => \$sFastaFile,
           "v|invert-match!" => \$flInvertMatch);

# PRINT HELP
$sHelp = 1 unless($sFastaFile and ($sIDfile or @ARGV));
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName -f <fasta file> [-m <id file>] [<id1> ... <idN>]

   Extract a set of sequences from a multifasta file.
   IDs can be read from a file or specified on the command line.
   
   Options:
    -m --match <string>
      Optional file with a column of fasta identifiers and an optional
      column to indicate whether to reverse-complement sequence (1) or
      keep the original orientation (0).
    -f --fasta <string>
      Fasta file. The first element of the fasta header before 
      the first whitespace character will be used for ID matching
    -v --invert-match
      Invert the sense of matching, to select non-matching entries.
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
      my ($sID, $nRevSeq) = split /\t/;
      $sID =~ s/^>+//;
      $nRevSeq    ||= 0;
      if (exists $hIDs{$sID}){
         die "Error: conflicting orientation data provided for sequence '$sID' in '$sIDfile'\n" unless ($hIDs{$sID} == $nRevSeq);
      }
      else{
         $hIDs{$sID}   = $nRevSeq;
      }
   }
   close IN;
}

# Read IDs from the command line
foreach my $sID (@ARGV){
   $sID =~ s/^>+//;
   $hIDs{$sID}   ||= 0;
}

# Make sure we have IDs at this point
warn "Warning: no fasta IDs were provided\n" unless(%hIDs);

# Process the fasta
my $nKept     = 0;
my $sFastaHeader = '';
my $sFastaSeq    = '';
my $flRevSeq     = 0;
open FASTA, "<$sFastaFile" or die "Error: can't read the fasta file\n";
while (<FASTA>){
   next if (/^\s*$/);
   next if (/^ *#/);
   if (/^>(\S*)(\s+|$)/){
      if (eof){
         die "Error: file ends in fasta header without sequence\n" if (/^>/);
      }
   
      # Print the previous match if there is any
      if ($sFastaHeader){
         if ($flRevSeq){
            $sFastaSeq =~ s/[\n\r]+$//;
            $sFastaSeq =  reverse_complement($sFastaSeq);
            $sFastaSeq .= "\n";
         }
         print join("", $sFastaHeader, $sFastaSeq);
      }
      
      # Reset everything for the new sequence
      my $flExtract = exists($hIDs{$1}) ? 1 : 0;
      $flExtract    = 1- $flExtract if ($flInvertMatch);
      if ($flExtract){
         $sFastaHeader = $_;
         $flRevSeq     = $hIDs{$1};
         $nKept++;
      }
      else{
         $sFastaHeader = '';
         $flRevSeq     = 0;
      }
      $sFastaSeq = '';
   }
   else{
      $sFastaSeq .= $_ if ($sFastaHeader);
   }
}
close FASTA;

# Print the last saved sequence if there is any
if ($sFastaHeader){
   if ($flRevSeq){
      $sFastaSeq =~ s/[\n\r]+$//;
      $sFastaSeq =  reverse_complement($sFastaSeq);
      $sFastaSeq .= "\n";
   }
   print join("", $sFastaHeader, $sFastaSeq);
}

# Report the number of sequences kept
my $sMessage = $nKept == 1 ? "Extracted $nKept sequence\n" : "Extracted $nKept sequences\n";
warn ($sMessage);

#################
## SUBROUTINES ##
#################

# rev_comp
#
# Returns the reverse-complement of the supplied sequence
sub reverse_complement{
   my $seq = shift(@_);
   my $rev = reverse $seq;
   $rev =~ tr/ACGTacgt/TGCAtgca/;
   return $rev;
}
